from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, Sequence, Tuple

import numpy as np
from scipy.interpolate import CubicSpline

from parametrization.shared.dependency_setup import (
    ensure_local_dependency_paths,
    load_pygeo_class,
    load_pyspline_curve,
)

from .laws import InterpolationMethod, InterpolationSpec


@dataclass(frozen=True)
class FuselageSectionSpec:
    section_id: str
    eta: float
    width: float
    height: float
    center_y: float = 0.0
    center_z: float = 0.0
    top_shape_exp: float = 2.0
    bottom_shape_exp: float = 2.0
    side_shape_exp: float = 2.0
    rotation_deg: float = 0.0

    def validate(self) -> None:
        if not self.section_id:
            raise ValueError("section_id must be non-empty")
        if not (0.0 <= self.eta <= 1.0):
            raise ValueError(f"fuselage section {self.section_id!r} eta must lie in [0, 1], got {self.eta}")
        if self.width <= 0.0:
            raise ValueError(f"fuselage section {self.section_id!r} width must be positive, got {self.width}")
        if self.height <= 0.0:
            raise ValueError(f"fuselage section {self.section_id!r} height must be positive, got {self.height}")
        if self.top_shape_exp <= 0.0 or self.bottom_shape_exp <= 0.0 or self.side_shape_exp <= 0.0:
            raise ValueError(
                f"fuselage section {self.section_id!r} shape exponents must be positive, "
                f"got top={self.top_shape_exp}, bottom={self.bottom_shape_exp}, side={self.side_shape_exp}"
            )


@dataclass(frozen=True)
class FuselageSpec:
    fuselage_id: str
    length: float
    sections: Tuple[FuselageSectionSpec, ...]
    nose_x: float = 0.0
    section_interpolation: InterpolationSpec = InterpolationSpec()
    spine_interpolation: InterpolationSpec = InterpolationSpec()
    metadata: Dict[str, str] = field(default_factory=dict)

    def validate(self) -> None:
        if not self.fuselage_id:
            raise ValueError("fuselage_id must be non-empty")
        if self.length <= 0.0:
            raise ValueError(f"length must be positive, got {self.length}")
        if len(self.sections) < 2:
            raise ValueError("a fuselage must define at least 2 sections")

        self.section_interpolation.validate()
        self.spine_interpolation.validate()

        seen = set()
        etas = []
        for section in self.sections:
            section.validate()
            if section.section_id in seen:
                raise ValueError(f"duplicate fuselage section_id detected: {section.section_id!r}")
            seen.add(section.section_id)
            etas.append(float(section.eta))

        if any(left >= right for left, right in zip(etas[:-1], etas[1:])):
            raise ValueError(f"fuselage sections must be strictly increasing in eta, got {tuple(etas)}")
        if abs(etas[0]) > 1e-12:
            raise ValueError(f"the first fuselage section must be at eta=0.0, got {etas[0]}")
        if abs(etas[-1] - 1.0) > 1e-12:
            raise ValueError(f"the last fuselage section must be at eta=1.0, got {etas[-1]}")


@dataclass(frozen=True)
class FuselageBuildOptions:
    station_count: int = 21
    station_etas: Tuple[float, ...] = ()
    include_anchor_sections: bool = True
    perimeter_point_count: int = 181
    ku: int = 4
    kv: int = 4
    n_ctlu: int | None = None
    n_ctlv: int | None = None
    rebuild_dependencies: bool = True

    def validate(self) -> None:
        if self.station_count < 2:
            raise ValueError(f"station_count must be >= 2, got {self.station_count}")
        if self.station_etas:
            etas = np.asarray(self.station_etas, dtype=float)
            if np.any(etas < 0.0) or np.any(etas > 1.0):
                raise ValueError(f"station_etas must lie in [0, 1], got {self.station_etas}")
            if np.any(np.diff(np.sort(etas)) <= 0.0):
                raise ValueError(f"station_etas must be strictly increasing, got {self.station_etas}")
        if self.perimeter_point_count < 21:
            raise ValueError(
                f"perimeter_point_count must be >= 21 for a robust closed section, got {self.perimeter_point_count}"
            )
        if self.ku not in {2, 3, 4}:
            raise ValueError(f"ku must be one of 2, 3 or 4, got {self.ku}")
        if self.kv not in {2, 3, 4}:
            raise ValueError(f"kv must be one of 2, 3 or 4, got {self.kv}")
        if self.n_ctlu is not None and self.n_ctlu < self.ku:
            raise ValueError(f"n_ctlu must be >= ku, got n_ctlu={self.n_ctlu}, ku={self.ku}")
        if self.n_ctlv is not None and self.n_ctlv < self.kv:
            raise ValueError(f"n_ctlv must be >= kv, got n_ctlv={self.n_ctlv}, kv={self.kv}")


@dataclass(frozen=True)
class PreparedFuselageStation:
    section_id: str
    eta: float
    x: float
    center_y: float
    center_z: float
    width: float
    height: float
    top_shape_exp: float
    bottom_shape_exp: float
    side_shape_exp: float
    rotation_deg: float
    loop_y: np.ndarray
    loop_z: np.ndarray


@dataclass(frozen=True)
class PreparedFuselage:
    component_id: str
    stations: Tuple[PreparedFuselageStation, ...]
    ku: int
    kv: int
    n_ctlu: int | None
    n_ctlv: int | None

    @property
    def x_grid(self) -> np.ndarray:
        return np.asarray(
            [np.full(station.loop_y.shape, float(station.x), dtype=float) for station in self.stations],
            dtype=float,
        )

    @property
    def y_grid(self) -> np.ndarray:
        return np.asarray([station.loop_y for station in self.stations], dtype=float)

    @property
    def z_grid(self) -> np.ndarray:
        return np.asarray([station.loop_z for station in self.stations], dtype=float)


@dataclass(frozen=True)
class FuselageExportResult:
    prepared: PreparedFuselage
    section_paths: Tuple[Path, ...]
    iges_path: Path


def _sanitize_label(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in {"-", "_"} else "_" for ch in value)


def _build_scalar_interpolant(
    etas: Sequence[float],
    values: Sequence[float],
    interpolation: InterpolationSpec,
) -> Callable[[float], float]:
    interpolation.validate()
    eta_array = np.asarray(etas, dtype=float)
    value_array = np.asarray(values, dtype=float)

    if eta_array.size != value_array.size:
        raise ValueError(
            f"eta/value size mismatch: got {eta_array.size} eta samples and {value_array.size} values"
        )
    if eta_array.size < 2:
        raise ValueError("at least 2 samples are required to build an interpolant")
    if np.any(np.diff(eta_array) <= 0.0):
        raise ValueError(f"etas must be strictly increasing, got {tuple(eta_array.tolist())}")

    if interpolation.method == InterpolationMethod.LINEAR or value_array.size == 2:
        return lambda eta: float(np.interp(float(eta), eta_array, value_array))

    if interpolation.method == InterpolationMethod.CUBIC:
        spline = CubicSpline(eta_array, value_array, bc_type="natural")
        return lambda eta: float(spline(float(eta)))

    ensure_local_dependency_paths()
    Curve = load_pyspline_curve(rebuild_if_needed=True)
    order = min(interpolation.pyspline_degree(), int(value_array.size))
    curve = Curve(x=value_array, s=eta_array, k=order)
    return lambda eta: float(np.asarray(curve(float(eta))).reshape(-1)[0])


def _build_positive_interpolant(
    etas: Sequence[float],
    values: Sequence[float],
    interpolation: InterpolationSpec,
    floor: float = 1e-8,
) -> Callable[[float], float]:
    value_array = np.asarray(values, dtype=float)
    if np.any(value_array <= 0.0):
        raise ValueError(f"positive interpolation received non-positive values: {tuple(value_array.tolist())}")
    if np.allclose(value_array, value_array[0]):
        constant = max(float(value_array[0]), floor)
        return lambda eta: constant

    safe_values = np.maximum(value_array, floor)
    log_interp = _build_scalar_interpolant(etas, np.log(safe_values), interpolation)
    return lambda eta: float(np.exp(log_interp(float(eta))))


def _resolve_station_etas(
    fuselage: FuselageSpec,
    options: FuselageBuildOptions,
) -> np.ndarray:
    if options.station_etas:
        etas = np.asarray(options.station_etas, dtype=float)
    else:
        etas = np.linspace(0.0, 1.0, int(options.station_count))

    if options.include_anchor_sections:
        anchor_etas = np.asarray([section.eta for section in fuselage.sections], dtype=float)
        etas = np.concatenate([etas, anchor_etas])

    return np.asarray(sorted({round(float(value), 12) for value in etas}), dtype=float)


def _sample_superellipse_loop(
    width: float,
    height: float,
    top_shape_exp: float,
    bottom_shape_exp: float,
    side_shape_exp: float,
    rotation_deg: float,
    point_count: int,
) -> tuple[np.ndarray, np.ndarray]:
    theta = np.linspace(0.0, 2.0 * np.pi, int(point_count), endpoint=True)
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)

    half_width = 0.5 * float(width)
    half_height = 0.5 * float(height)

    z = half_width * np.sign(cos_theta) * np.abs(cos_theta) ** (2.0 / float(side_shape_exp))
    y = np.empty_like(z)
    upper_mask = sin_theta >= 0.0
    y[upper_mask] = half_height * np.abs(sin_theta[upper_mask]) ** (2.0 / float(top_shape_exp))
    y[~upper_mask] = -half_height * np.abs(sin_theta[~upper_mask]) ** (2.0 / float(bottom_shape_exp))

    angle = np.deg2rad(float(rotation_deg))
    c = float(np.cos(angle))
    s = float(np.sin(angle))
    y_rot = c * y - s * z
    z_rot = s * y + c * z
    return y_rot, z_rot


def prepare_fuselage(
    fuselage: FuselageSpec,
    options: FuselageBuildOptions = FuselageBuildOptions(),
) -> PreparedFuselage:
    fuselage.validate()
    options.validate()

    section_etas = [section.eta for section in fuselage.sections]
    width_interp = _build_positive_interpolant(
        section_etas,
        [section.width for section in fuselage.sections],
        fuselage.section_interpolation,
    )
    height_interp = _build_positive_interpolant(
        section_etas,
        [section.height for section in fuselage.sections],
        fuselage.section_interpolation,
    )
    center_y_interp = _build_scalar_interpolant(
        section_etas,
        [section.center_y for section in fuselage.sections],
        fuselage.spine_interpolation,
    )
    center_z_interp = _build_scalar_interpolant(
        section_etas,
        [section.center_z for section in fuselage.sections],
        fuselage.spine_interpolation,
    )
    top_exp_interp = _build_positive_interpolant(
        section_etas,
        [section.top_shape_exp for section in fuselage.sections],
        fuselage.section_interpolation,
    )
    bottom_exp_interp = _build_positive_interpolant(
        section_etas,
        [section.bottom_shape_exp for section in fuselage.sections],
        fuselage.section_interpolation,
    )
    side_exp_interp = _build_positive_interpolant(
        section_etas,
        [section.side_shape_exp for section in fuselage.sections],
        fuselage.section_interpolation,
    )
    rotation_interp = _build_scalar_interpolant(
        section_etas,
        [section.rotation_deg for section in fuselage.sections],
        fuselage.section_interpolation,
    )

    stations = []
    for index, eta in enumerate(_resolve_station_etas(fuselage, options)):
        width = float(width_interp(float(eta)))
        height = float(height_interp(float(eta)))
        center_y = float(center_y_interp(float(eta)))
        center_z = float(center_z_interp(float(eta)))
        top_exp = float(top_exp_interp(float(eta)))
        bottom_exp = float(bottom_exp_interp(float(eta)))
        side_exp = float(side_exp_interp(float(eta)))
        rotation_deg = float(rotation_interp(float(eta)))
        loop_y, loop_z = _sample_superellipse_loop(
            width=width,
            height=height,
            top_shape_exp=top_exp,
            bottom_shape_exp=bottom_exp,
            side_shape_exp=side_exp,
            rotation_deg=rotation_deg,
            point_count=options.perimeter_point_count,
        )
        loop_y = loop_y + center_y
        loop_z = loop_z + center_z

        stations.append(
            PreparedFuselageStation(
                section_id=f"{fuselage.fuselage_id}_station_{index:03d}",
                eta=float(eta),
                x=float(fuselage.nose_x) + float(fuselage.length) * float(eta),
                center_y=center_y,
                center_z=center_z,
                width=width,
                height=height,
                top_shape_exp=top_exp,
                bottom_shape_exp=bottom_exp,
                side_shape_exp=side_exp,
                rotation_deg=rotation_deg,
                loop_y=loop_y,
                loop_z=loop_z,
            )
        )

    return PreparedFuselage(
        component_id=fuselage.fuselage_id,
        stations=tuple(stations),
        ku=int(options.ku),
        kv=int(options.kv),
        n_ctlu=options.n_ctlu,
        n_ctlv=options.n_ctlv,
    )


def write_fuselage_sections(
    prepared: PreparedFuselage,
    out_dir: str | Path,
) -> Tuple[Path, ...]:
    output_dir = Path(out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    component_label = _sanitize_label(prepared.component_id)
    section_paths = []
    for index, station in enumerate(prepared.stations):
        path = output_dir / f"{component_label}_section_{index:03d}_eta{station.eta:.3f}.dat"
        data = np.column_stack(
            [
                np.full(station.loop_y.shape, float(station.x), dtype=float),
                station.loop_y,
                station.loop_z,
            ]
        )
        header = "x y z"
        np.savetxt(path, data, header=header, comments="")
        section_paths.append(path)
    return tuple(section_paths)


def build_fuselage_surface(
    prepared: PreparedFuselage,
    rebuild_if_needed: bool = True,
):
    ensure_local_dependency_paths()
    load_pyspline_curve(rebuild_if_needed=rebuild_if_needed)
    from pyspline import Surface

    kwargs = {
        "x": prepared.x_grid,
        "y": prepared.y_grid,
        "z": prepared.z_grid,
        "ku": prepared.ku,
        "kv": prepared.kv,
    }
    if prepared.n_ctlu is not None:
        kwargs["nCtlu"] = prepared.n_ctlu
    if prepared.n_ctlv is not None:
        kwargs["nCtlv"] = prepared.n_ctlv
    return Surface(**kwargs)


def build_fuselage_pygeo(
    prepared: PreparedFuselage,
    rebuild_if_needed: bool = True,
):
    pyGeo = load_pygeo_class(rebuild_if_needed=rebuild_if_needed)
    surface = build_fuselage_surface(prepared, rebuild_if_needed=rebuild_if_needed)
    geometry = pyGeo("create")
    geometry.surfs = [surface]
    geometry.nSurf = 1
    return geometry


def export_fuselage_iges(
    fuselage: FuselageSpec,
    out_dir: str | Path,
    options: FuselageBuildOptions = FuselageBuildOptions(),
    file_name: str | None = None,
) -> FuselageExportResult:
    prepared = prepare_fuselage(fuselage, options=options)
    output_dir = Path(out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    section_paths = write_fuselage_sections(prepared, output_dir / "sections")
    geometry = build_fuselage_pygeo(prepared, rebuild_if_needed=options.rebuild_dependencies)
    iges_name = file_name or f"{_sanitize_label(fuselage.fuselage_id)}.igs"
    iges_path = output_dir / iges_name
    geometry.writeIGES(str(iges_path))

    return FuselageExportResult(
        prepared=prepared,
        section_paths=section_paths,
        iges_path=iges_path,
    )
