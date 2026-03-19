from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Sequence

import numpy as np
from scipy.interpolate import CubicSpline

from parametrization.core.airfoil_io import write_airfoil_dat
from parametrization.core.cst_sharedle import KulfanCSTAirfoil, cosine_spacing
from parametrization.core.dependency_setup import (
    ensure_local_dependency_paths,
    load_pygeo_class,
    load_pyspline_curve,
)

from .components import ComponentKind, LoftedComponentSpec
from .laws import ContinuityOrder, InterpolationMethod, InterpolationSpec, ScalarLawSpec
from .profiles import ProfileCatalog

SUPPORTED_LIFTING_SURFACE_LAWS = (
    "x",
    "y",
    "z",
    "chord",
    "roll_deg",
    "pitch_deg",
    "twist_deg",
    "thickness_scale",
    "te_thickness",
)


@dataclass(frozen=True)
class LiftingSurfaceBuildOptions:
    station_count: int = 9
    station_etas: tuple[float, ...] = ()
    include_anchor_sections: bool = True
    airfoil_sample_count: int = 241
    fit_n_ctl: int | None = None
    k_span: int | None = None
    tip_style: str = "rounded"
    blunt_te: bool = False
    rounded_te: bool = False
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
        if self.airfoil_sample_count < 21:
            raise ValueError(
                f"airfoil_sample_count must be >= 21 for robust CST sampling, got {self.airfoil_sample_count}"
            )
        if self.fit_n_ctl is not None and self.fit_n_ctl < 4:
            raise ValueError(f"fit_n_ctl must be >= 4 when provided, got {self.fit_n_ctl}")
        if self.k_span is not None and self.k_span not in {2, 3, 4}:
            raise ValueError(f"k_span must be one of 2, 3 or 4, got {self.k_span}")
        if self.tip_style not in {"rounded", "pinched"}:
            raise ValueError(f"unsupported tip_style {self.tip_style!r}")
        if self.rounded_te and not self.blunt_te:
            raise ValueError("rounded_te=True requires blunt_te=True")


@dataclass(frozen=True)
class PreparedLiftingSurfaceStation:
    section_id: str
    eta: float
    x: float
    y: float
    z: float
    chord: float
    roll_deg: float
    pitch_deg: float
    twist_deg: float
    thickness_scale: float
    te_height_scaled: float
    x_air: np.ndarray
    yu: np.ndarray
    yl: np.ndarray


@dataclass(frozen=True)
class PreparedLiftingSurface:
    component_id: str
    stations: tuple[PreparedLiftingSurfaceStation, ...]
    offset: np.ndarray
    k_span: int
    fit_n_ctl: int | None
    tip_style: str
    blunt_te: bool
    rounded_te: bool

    @property
    def station_etas(self) -> np.ndarray:
        return np.asarray([station.eta for station in self.stations], dtype=float)

    @property
    def x(self) -> np.ndarray:
        return np.asarray([station.x for station in self.stations], dtype=float)

    @property
    def y(self) -> np.ndarray:
        return np.asarray([station.y for station in self.stations], dtype=float)

    @property
    def z(self) -> np.ndarray:
        return np.asarray([station.z for station in self.stations], dtype=float)

    @property
    def chord(self) -> np.ndarray:
        return np.asarray([station.chord for station in self.stations], dtype=float)

    @property
    def rot_x_deg(self) -> np.ndarray:
        return np.asarray([station.roll_deg for station in self.stations], dtype=float)

    @property
    def rot_y_deg(self) -> np.ndarray:
        return np.asarray([station.pitch_deg for station in self.stations], dtype=float)

    @property
    def rot_z_deg(self) -> np.ndarray:
        return np.asarray([station.twist_deg for station in self.stations], dtype=float)

    @property
    def te_height_scaled(self) -> np.ndarray:
        return np.asarray([station.te_height_scaled for station in self.stations], dtype=float)


@dataclass(frozen=True)
class LiftingSurfaceExportResult:
    prepared: PreparedLiftingSurface
    airfoil_paths: tuple[Path, ...]
    iges_path: Path


def _sanitize_label(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in {"-", "_"} else "_" for ch in value)


def _segmented_blend_weight(s: float, interpolation: InterpolationSpec) -> float:
    s_clamped = min(max(float(s), 0.0), 1.0)
    if interpolation.continuity == ContinuityOrder.C1:
        return s_clamped * s_clamped * (3.0 - 2.0 * s_clamped)
    if interpolation.continuity == ContinuityOrder.C2:
        return s_clamped**3 * (10.0 + s_clamped * (-15.0 + 6.0 * s_clamped))
    raise ValueError(
        "segmented interpolation only supports C1 or C2 continuity, "
        f"got {interpolation.continuity!r}"
    )


def _build_segmented_interpolant(
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
    if eta_array.size == 2:
        return lambda eta: float(np.interp(float(eta), eta_array, value_array))

    span_widths = np.diff(eta_array)
    slopes = np.diff(value_array) / span_widths
    transitions: list[tuple[int, float, float]] = []
    for idx in range(1, eta_array.size - 1):
        half_width = float(interpolation.blend_fraction) * min(
            float(eta_array[idx] - eta_array[idx - 1]),
            float(eta_array[idx + 1] - eta_array[idx]),
        )
        if half_width <= 0.0:
            continue
        transitions.append((idx, float(eta_array[idx] - half_width), float(eta_array[idx] + half_width)))

    def _linear_segment_value(segment_index: int, eta: float) -> float:
        eta0 = float(eta_array[segment_index])
        return float(value_array[segment_index] + slopes[segment_index] * (float(eta) - eta0))

    def _interpolant(eta: float) -> float:
        eta_value = float(eta)
        if eta_value <= float(eta_array[0]):
            return float(value_array[0])
        if eta_value >= float(eta_array[-1]):
            return float(value_array[-1])

        for anchor_index, eta_left, eta_right in transitions:
            if eta_left <= eta_value <= eta_right:
                left_value = _linear_segment_value(anchor_index - 1, eta_value)
                right_value = _linear_segment_value(anchor_index, eta_value)
                local_eta = (eta_value - eta_left) / (eta_right - eta_left)
                weight = _segmented_blend_weight(local_eta, interpolation)
                return float((1.0 - weight) * left_value + weight * right_value)

        segment_index = int(np.searchsorted(eta_array, eta_value, side="right") - 1)
        segment_index = min(max(segment_index, 0), eta_array.size - 2)
        return _linear_segment_value(segment_index, eta_value)

    return _interpolant


def _segmented_support_etas(
    etas: Sequence[float],
    interpolation: InterpolationSpec,
) -> tuple[float, ...]:
    if interpolation.method != InterpolationMethod.SEGMENTED:
        return ()

    eta_array = np.asarray(etas, dtype=float)
    if eta_array.size < 3:
        return ()

    support: list[float] = []
    for idx in range(1, eta_array.size - 1):
        half_width = float(interpolation.blend_fraction) * min(
            float(eta_array[idx] - eta_array[idx - 1]),
            float(eta_array[idx + 1] - eta_array[idx]),
        )
        if half_width <= 0.0:
            continue
        support.extend(
            [
                float(eta_array[idx] - half_width),
                float(eta_array[idx] - 0.5 * half_width),
                float(eta_array[idx]),
                float(eta_array[idx] + 0.5 * half_width),
                float(eta_array[idx] + half_width),
            ]
        )
    return tuple(support)


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

    if interpolation.method == InterpolationMethod.SEGMENTED:
        return _build_segmented_interpolant(eta_array, value_array, interpolation)

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
    if np.allclose(value_array, value_array[0]):
        constant = max(float(value_array[0]), floor)
        return lambda eta: constant

    if interpolation.method == InterpolationMethod.SEGMENTED:
        raw_interp = _build_segmented_interpolant(etas, np.maximum(value_array, floor), interpolation)
        return lambda eta: max(float(raw_interp(float(eta))), floor)

    safe_values = np.maximum(value_array, floor)
    log_interp = _build_scalar_interpolant(etas, np.log(safe_values), interpolation)
    return lambda eta: float(np.exp(log_interp(float(eta))))


def _build_nonnegative_interpolant(
    etas: Sequence[float],
    values: Sequence[float],
    interpolation: InterpolationSpec,
) -> Callable[[float], float]:
    raw_interp = _build_scalar_interpolant(etas, values, interpolation)
    return lambda eta: max(float(raw_interp(float(eta))), 0.0)


def _component_eta_bounds(component: LoftedComponentSpec) -> tuple[float, float]:
    component.validate()
    return float(component.sections[0].eta), float(component.sections[-1].eta)


def _resolve_station_etas(
    component: LoftedComponentSpec,
    options: LiftingSurfaceBuildOptions,
) -> np.ndarray:
    eta_start, eta_end = _component_eta_bounds(component)
    if options.station_etas:
        etas = np.asarray(options.station_etas, dtype=float)
    else:
        etas = np.linspace(eta_start, eta_end, int(options.station_count))

    if options.include_anchor_sections:
        anchor_etas = np.asarray([section.eta for section in component.sections], dtype=float)
        etas = np.concatenate([etas, anchor_etas])

    if component.section_interpolation.method == InterpolationMethod.SEGMENTED:
        section_etas = [section.eta for section in component.sections]
        etas = np.concatenate([etas, np.asarray(_segmented_support_etas(section_etas, component.section_interpolation))])

    if component.spine_interpolation.method == InterpolationMethod.SEGMENTED:
        section_etas = [section.eta for section in component.sections]
        etas = np.concatenate([etas, np.asarray(_segmented_support_etas(section_etas, component.spine_interpolation))])

    for law in component.scalar_laws:
        if law.interpolation.method == InterpolationMethod.SEGMENTED:
            etas = np.concatenate([etas, np.asarray(_segmented_support_etas(law.anchor_etas(), law.interpolation))])

    etas = np.asarray(sorted({round(float(value), 12) for value in etas}), dtype=float)
    if etas[0] < eta_start - 1e-12 or etas[-1] > eta_end + 1e-12:
        raise ValueError(
            f"station eta grid must stay inside component eta range [{eta_start}, {eta_end}], got {tuple(etas)}"
        )
    return etas


def _component_scalar_defaults(component: LoftedComponentSpec, name: str) -> np.ndarray:
    if name not in SUPPORTED_LIFTING_SURFACE_LAWS:
        raise ValueError(f"unsupported lifting-surface scalar {name!r}")

    if name == "te_thickness":
        raise ValueError("te_thickness is not a placement scalar and must be resolved from the profile family")

    return np.asarray([getattr(section.placement, name) for section in component.sections], dtype=float)


def _component_profile_family(
    component: LoftedComponentSpec,
    catalog: ProfileCatalog,
) -> tuple[KulfanCSTAirfoil, np.ndarray, np.ndarray, np.ndarray]:
    profile_specs = [catalog.resolve(section.profile_id) for section in component.sections]
    reference = profile_specs[0]

    for profile in profile_specs[1:]:
        if profile.degree != reference.degree:
            raise ValueError(
                f"component {component.component_id!r} mixes incompatible CST degrees: "
                f"{reference.degree} vs {profile.degree}"
            )
        if abs(profile.n1 - reference.n1) > 1e-12 or abs(profile.n2 - reference.n2) > 1e-12:
            raise ValueError(
                f"component {component.component_id!r} mixes incompatible CST class exponents "
                f"between profiles {reference.profile_id!r} and {profile.profile_id!r}"
            )
        if tuple(profile.x_tc_window) != tuple(reference.x_tc_window):
            raise ValueError(
                f"component {component.component_id!r} mixes incompatible x_tc_window values "
                f"between profiles {reference.profile_id!r} and {profile.profile_id!r}"
            )

    shape = KulfanCSTAirfoil(
        degree=reference.degree,
        n1=reference.n1,
        n2=reference.n2,
        x_tc_window=reference.x_tc_window,
    )
    upper = np.asarray([profile.upper_coeffs for profile in profile_specs], dtype=float)
    lower = np.asarray([profile.lower_coeffs for profile in profile_specs], dtype=float)
    te = np.asarray([profile.te_thickness for profile in profile_specs], dtype=float)
    return shape, upper, lower, te


def _scalar_law_map(component: LoftedComponentSpec) -> dict[str, ScalarLawSpec]:
    law_map = component.scalar_law_map()
    unknown = sorted(set(law_map).difference(SUPPORTED_LIFTING_SURFACE_LAWS))
    if unknown:
        raise ValueError(
            f"component {component.component_id!r} defines unsupported lifting-surface scalar laws: {unknown}"
        )
    return law_map


def _validate_law_coverage(
    component: LoftedComponentSpec,
    law: ScalarLawSpec,
) -> None:
    eta_start, eta_end = _component_eta_bounds(component)
    law_etas = law.anchor_etas()
    if law_etas[0] > eta_start or law_etas[-1] < eta_end:
        raise ValueError(
            f"law {law.name!r} must cover the full component eta range "
            f"[{eta_start}, {eta_end}], got [{law_etas[0]}, {law_etas[-1]}]"
        )


def _resolve_component_scalar(
    component: LoftedComponentSpec,
    law_map: dict[str, ScalarLawSpec],
    name: str,
) -> Callable[[float], float]:
    if name in law_map:
        law = law_map[name]
        _validate_law_coverage(component, law)
        builder = _build_positive_interpolant if name in {"chord", "thickness_scale"} else _build_scalar_interpolant
        return builder(law.anchor_etas(), law.anchor_values(), law.interpolation)

    default_values = _component_scalar_defaults(component, name)
    builder = _build_positive_interpolant if name in {"chord", "thickness_scale"} else _build_scalar_interpolant
    section_etas = [section.eta for section in component.sections]
    return builder(section_etas, default_values, component.spine_interpolation)


def _resolve_te_thickness(
    component: LoftedComponentSpec,
    law_map: dict[str, ScalarLawSpec],
    te_anchor_values: np.ndarray,
) -> Callable[[float], float]:
    if "te_thickness" in law_map:
        law = law_map["te_thickness"]
        _validate_law_coverage(component, law)
        return _build_nonnegative_interpolant(law.anchor_etas(), law.anchor_values(), law.interpolation)

    section_etas = [section.eta for section in component.sections]
    return _build_nonnegative_interpolant(section_etas, te_anchor_values, component.section_interpolation)


def _resolve_profile_interpolants(
    component: LoftedComponentSpec,
    catalog: ProfileCatalog,
) -> tuple[KulfanCSTAirfoil, list[Callable[[float], float]], list[Callable[[float], float]], np.ndarray]:
    shape, upper_sections, lower_sections, te_sections = _component_profile_family(component, catalog)
    section_etas = [section.eta for section in component.sections]

    upper_interpolants = [
        _build_scalar_interpolant(section_etas, upper_sections[:, idx], component.section_interpolation)
        for idx in range(upper_sections.shape[1])
    ]
    lower_interpolants = [
        _build_scalar_interpolant(section_etas, lower_sections[:, idx], component.section_interpolation)
        for idx in range(lower_sections.shape[1])
    ]
    return shape, upper_interpolants, lower_interpolants, te_sections


def _resolved_k_span(
    component: LoftedComponentSpec,
    station_count: int,
    options: LiftingSurfaceBuildOptions,
) -> int:
    def _order(interpolation: InterpolationSpec) -> int:
        if interpolation.method == InterpolationMethod.LINEAR:
            return 2
        if interpolation.method == InterpolationMethod.SEGMENTED:
            return 3 if interpolation.continuity == ContinuityOrder.C1 else 4
        if interpolation.method == InterpolationMethod.CUBIC:
            return 4
        return interpolation.pyspline_degree()

    if station_count <= 2:
        return 2
    if options.k_span is not None:
        return min(int(options.k_span), int(station_count))

    law_orders = [_order(law.interpolation) for law in component.scalar_laws]
    requested_order = max(
        [_order(component.section_interpolation), _order(component.spine_interpolation), *law_orders],
        default=2,
    )
    return min(int(requested_order), int(station_count))


def prepare_lifting_surface(
    component: LoftedComponentSpec,
    profiles: ProfileCatalog,
    options: LiftingSurfaceBuildOptions = LiftingSurfaceBuildOptions(),
) -> PreparedLiftingSurface:
    component.validate()
    profiles.validate()
    options.validate()

    if component.kind not in {ComponentKind.LIFTING_SURFACE, ComponentKind.GENERIC_LOFT}:
        raise ValueError(
            f"component {component.component_id!r} has kind {component.kind.value!r}, "
            "which is not compatible with lifting-surface export"
        )
    if component.mirrored:
        raise NotImplementedError(
            "mirrored=True is not supported yet in the generic lifting-surface builder. "
            "For now define one physical side per component and mirror later at aircraft-assembly level."
        )

    law_map = _scalar_law_map(component)
    station_etas = _resolve_station_etas(component, options)
    x_air = cosine_spacing(int(options.airfoil_sample_count))

    x_interp = _resolve_component_scalar(component, law_map, "x")
    y_interp = _resolve_component_scalar(component, law_map, "y")
    z_interp = _resolve_component_scalar(component, law_map, "z")
    chord_interp = _resolve_component_scalar(component, law_map, "chord")
    roll_interp = _resolve_component_scalar(component, law_map, "roll_deg")
    pitch_interp = _resolve_component_scalar(component, law_map, "pitch_deg")
    twist_interp = _resolve_component_scalar(component, law_map, "twist_deg")
    thickness_interp = _resolve_component_scalar(component, law_map, "thickness_scale")

    shape, upper_interpolants, lower_interpolants, te_sections = _resolve_profile_interpolants(component, profiles)
    te_interp = _resolve_te_thickness(component, law_map, te_sections)

    stations = []
    for index, eta in enumerate(station_etas):
        upper = np.asarray([interp(float(eta)) for interp in upper_interpolants], dtype=float)
        lower = np.asarray([interp(float(eta)) for interp in lower_interpolants], dtype=float)
        coeffs = np.concatenate([upper, lower])
        te_thickness = float(te_interp(float(eta)))
        yu, yl = shape.evaluate(x_air, coeffs, te_thickness=te_thickness)

        thickness_scale = float(thickness_interp(float(eta)))
        yu = thickness_scale * yu
        yl = thickness_scale * yl

        stations.append(
            PreparedLiftingSurfaceStation(
                section_id=f"{component.component_id}_station_{index:03d}",
                eta=float(eta),
                x=float(x_interp(float(eta))),
                y=float(y_interp(float(eta))),
                z=float(z_interp(float(eta))),
                chord=float(chord_interp(float(eta))),
                roll_deg=float(roll_interp(float(eta))),
                pitch_deg=float(pitch_interp(float(eta))),
                twist_deg=float(twist_interp(float(eta))),
                thickness_scale=thickness_scale,
                te_height_scaled=float(yu[-1] - yl[-1]),
                x_air=x_air.copy(),
                yu=yu,
                yl=yl,
            )
        )

    offset = np.zeros((len(stations), 2), dtype=float)
    return PreparedLiftingSurface(
        component_id=component.component_id,
        stations=tuple(stations),
        offset=offset,
        k_span=_resolved_k_span(component, len(stations), options),
        fit_n_ctl=options.fit_n_ctl,
        tip_style=options.tip_style,
        blunt_te=options.blunt_te,
        rounded_te=options.rounded_te,
    )


def write_lifting_surface_airfoils(
    prepared: PreparedLiftingSurface,
    out_dir: str | Path,
) -> tuple[Path, ...]:
    output_dir = Path(out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    component_label = _sanitize_label(prepared.component_id)
    airfoil_paths = []
    for index, station in enumerate(prepared.stations):
        file_name = f"{component_label}_station_{index:03d}_eta{station.eta:.3f}.dat"
        path = output_dir / file_name
        write_airfoil_dat(
            str(path),
            station.x_air,
            station.yu,
            station.yl,
            name=f"{component_label.upper()}_{index:03d}",
        )
        airfoil_paths.append(path)
    return tuple(airfoil_paths)


def build_lifting_surface_pygeo(
    prepared: PreparedLiftingSurface,
    airfoil_paths: Iterable[str | Path],
    rebuild_if_needed: bool = True,
):
    ensure_local_dependency_paths()
    pyGeo = load_pygeo_class(rebuild_if_needed=rebuild_if_needed)

    kwargs = {
        "xsections": [str(Path(path)) for path in airfoil_paths],
        "scale": prepared.chord,
        "offset": prepared.offset,
        "x": prepared.x,
        "y": prepared.y,
        "z": prepared.z,
        "rotX": prepared.rot_x_deg,
        "rotY": prepared.rot_y_deg,
        "rotZ": prepared.rot_z_deg,
        "kSpan": prepared.k_span,
        "tip": prepared.tip_style,
        "bluntTe": prepared.blunt_te,
        "roundedTe": prepared.rounded_te,
    }
    if prepared.fit_n_ctl is not None:
        kwargs["nCtl"] = prepared.fit_n_ctl
    if prepared.blunt_te:
        kwargs["teHeightScaled"] = prepared.te_height_scaled.tolist()

    return pyGeo("liftingSurface", **kwargs)


def export_lifting_surface_iges(
    component: LoftedComponentSpec,
    profiles: ProfileCatalog,
    out_dir: str | Path,
    options: LiftingSurfaceBuildOptions = LiftingSurfaceBuildOptions(),
    file_name: str | None = None,
) -> LiftingSurfaceExportResult:
    prepared = prepare_lifting_surface(component, profiles, options=options)
    output_dir = Path(out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    airfoil_paths = write_lifting_surface_airfoils(prepared, output_dir / "profiles")
    geometry = build_lifting_surface_pygeo(
        prepared,
        airfoil_paths,
        rebuild_if_needed=options.rebuild_dependencies,
    )

    iges_name = file_name or f"{_sanitize_label(component.component_id)}.igs"
    iges_path = output_dir / iges_name
    geometry.writeIGES(str(iges_path))

    return LiftingSurfaceExportResult(
        prepared=prepared,
        airfoil_paths=airfoil_paths,
        iges_path=iges_path,
    )
