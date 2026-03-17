from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

from ..cst_sharedle import KulfanCSTAirfoil, cosine_spacing
from ..dependency_setup import ensure_local_dependency_paths, load_pyspline_curve
from ..design_variables import SectionedBWBDesignVariables
from ..planform import build_sectioned_bwb_planform
from ..sections import apply_camber_delta
from ..sections import build_section_model
from ..spanwise_laws import resolve_spanwise_laws


@dataclass(frozen=True)
class FFDBoxSpec:
    n_streamwise: int = 10
    n_vertical: int = 2
    n_spanwise: int = 17
    x_margin_le_m: float = 0.42
    x_margin_te_m: float = 0.34
    y_margin_m: float = 0.24
    z_margin_root_m: float = 0.24
    z_margin_tip_m: float = 0.16
    contour_fit_n_spanwise: int = 121
    contour_fit_n_chordwise: int = 161


@dataclass(frozen=True)
class FFDBoxSummary:
    xyz_path: Path
    n_streamwise: int
    n_vertical: int
    n_spanwise: int
    full_wing: bool
    x_min: float
    x_max: float
    y_min: float
    y_max: float
    z_min: float
    z_max: float
    dvgeo_load_ok: Optional[bool]


def _write_plot3d_ffd(path: Path, points: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    ni, nj, nk, _ = points.shape
    with path.open("w", encoding="utf-8") as stream:
        stream.write("1\n")
        stream.write(f"{ni} {nj} {nk}\n")
        for axis in range(3):
            for k in range(nk):
                for j in range(nj):
                    for i in range(ni):
                        stream.write(f"{points[i, j, k, axis]:.15f} ")
                    stream.write("\n")


def read_plot3d_ffd(path: Path) -> np.ndarray:
    values = path.read_text(encoding="utf-8").split()
    if len(values) < 4:
        raise ValueError(f"Invalid Plot3D FFD file: {path}")

    nblocks = int(values[0])
    if nblocks != 1:
        raise ValueError(
            f"Only single-block Plot3D FFD files are supported here, got {nblocks} blocks"
        )

    ni = int(values[1])
    nj = int(values[2])
    nk = int(values[3])
    coord_count = 3 * ni * nj * nk
    if len(values) != 4 + coord_count:
        raise ValueError(
            f"Unexpected Plot3D payload size in {path}: expected {coord_count} coordinates, "
            f"got {len(values) - 4}"
        )

    coords = np.asarray([float(value) for value in values[4:]], dtype=float).reshape(3, nk, nj, ni)
    points = np.zeros((ni, nj, nk, 3), dtype=float)
    for axis in range(3):
        points[..., axis] = np.transpose(coords[axis], (2, 1, 0))
    return points


def _build_pyspline_scalar_interpolant(
    span_values: np.ndarray,
    data_values: np.ndarray,
):
    span_values = np.asarray(span_values, dtype=float)
    data_values = np.asarray(data_values, dtype=float)
    if span_values.size == 1:
        return lambda yy: float(data_values[0])
    Curve = load_pyspline_curve(rebuild_if_needed=True)
    order = min(4, int(span_values.size))
    curve = Curve(x=data_values, s=span_values, k=order)
    return lambda yy: float(np.asarray(curve(float(yy))).reshape(-1)[0])


def _cosine_span_stations(span: float, n_spanwise: int, full_wing: bool = False) -> np.ndarray:
    semispan = float(span) * cosine_spacing(int(n_spanwise))
    if not full_wing:
        return semispan
    return np.concatenate([-semispan[:0:-1], semispan])


def _build_reference_section_data(
    design: SectionedBWBDesignVariables,
):
    config = design.to_model_config(profile_generation_mode="cst_only")
    topology = config.topology
    planform = config.planform
    sections = config.sections

    y_sections = topology.y_sections_array
    le_sections = planform.leading_edge_x_sections(topology)
    chord_sections = planform.section_chords()
    vertical_offsets = np.tan(np.deg2rad(config.spanwise.dihedral_deg)) * y_sections

    upper_coeffs = np.asarray([spec.upper_coeffs for spec in sections.section_specs], dtype=float)
    lower_coeffs = np.asarray([spec.lower_coeffs for spec in sections.section_specs], dtype=float)
    te_thickness = np.asarray([spec.te_thickness for spec in sections.section_specs], dtype=float)

    twist_section_indices = np.asarray(config.spanwise.twist_deg.section_indices, dtype=int)
    twist_anchor_y = y_sections[twist_section_indices]
    twist_anchor_values = np.asarray(config.spanwise.twist_deg.values, dtype=float)

    camber_section_indices = np.asarray(config.spanwise.camber_delta.section_indices, dtype=int)
    camber_anchor_y = y_sections[camber_section_indices]
    camber_anchor_values = np.asarray(config.spanwise.camber_delta.values, dtype=float)

    return {
        "config": config,
        "topology": topology,
        "sections": sections,
        "y_sections": y_sections,
        "le_sections": le_sections,
        "chord_sections": chord_sections,
        "vertical_offsets": vertical_offsets,
        "upper_coeffs": upper_coeffs,
        "lower_coeffs": lower_coeffs,
        "te_thickness": te_thickness,
        "twist_anchor_y": twist_anchor_y,
        "twist_anchor_values": twist_anchor_values,
        "camber_anchor_y": camber_anchor_y,
        "camber_anchor_values": camber_anchor_values,
    }


def _interpolate_section_coefficients(
    span_abs: float,
    section_data,
) -> tuple[np.ndarray, float]:
    y_sections = section_data["y_sections"]
    upper_sections = section_data["upper_coeffs"]
    lower_sections = section_data["lower_coeffs"]
    te_sections = section_data["te_thickness"]
    sections = section_data["sections"]

    upper = np.asarray(
        [np.interp(span_abs, y_sections, upper_sections[:, idx]) for idx in range(upper_sections.shape[1])],
        dtype=float,
    )
    lower = np.asarray(
        [np.interp(span_abs, y_sections, lower_sections[:, idx]) for idx in range(lower_sections.shape[1])],
        dtype=float,
    )
    coeffs = np.concatenate([upper, lower])

    camber_delta = float(
        np.interp(span_abs, section_data["camber_anchor_y"], section_data["camber_anchor_values"])
    )
    coeffs = apply_camber_delta(
        coeffs,
        camber_delta=camber_delta,
        degree=sections.cst_degree,
        center=sections.camber_mode_center,
        width=sections.camber_mode_width,
    )
    te_thickness = float(np.interp(span_abs, y_sections, te_sections))
    return coeffs, te_thickness


def _twist_at_span(span_abs: float, section_data) -> float:
    return float(
        np.interp(
            span_abs,
            section_data["twist_anchor_y"],
            section_data["twist_anchor_values"],
        )
    )


def build_reference_surface_grids(
    design: Optional[SectionedBWBDesignVariables] = None,
    n_spanwise: int = 17,
    n_chordwise: int = 81,
    full_wing: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    if design is None:
        design = SectionedBWBDesignVariables.reference_seed()

    config = design.to_model_config(profile_generation_mode="cst_only")
    config.sampling.num_airfoil_points = int(n_chordwise)
    topology = config.topology
    section_data = _build_reference_section_data(design)
    planform = build_sectioned_bwb_planform(topology, config.planform)
    laws = resolve_spanwise_laws(config)
    section_model = build_section_model(config, laws)

    span_stations = _cosine_span_stations(
        topology.span,
        int(n_spanwise),
        full_wing=bool(full_wing),
    )

    upper_grid = np.zeros((span_stations.size, section_model.x_air.size, 3), dtype=float)
    lower_grid = np.zeros((span_stations.size, section_model.x_air.size, 3), dtype=float)

    for idx, span_station in enumerate(span_stations):
        span_abs = abs(float(span_station))
        le_x = float(planform.le_x(span_abs))
        te_x = float(planform.te_x(span_abs))
        chord = float(te_x - le_x)
        y_offset = float(np.tan(np.deg2rad(config.spanwise.dihedral_deg)) * span_abs)
        twist_deg = float(laws.twist_deg(span_abs))
        yu, yl, _ = section_model.coordinates_at_y(span_abs)

        x_local = chord * section_model.x_air
        y_upper_local = chord * yu
        y_lower_local = chord * yl
        theta = np.deg2rad(twist_deg)
        cos_theta = float(np.cos(theta))
        sin_theta = float(np.sin(theta))

        upper_grid[idx, :, 0] = le_x + cos_theta * x_local - sin_theta * y_upper_local
        upper_grid[idx, :, 1] = y_offset + sin_theta * x_local + cos_theta * y_upper_local
        upper_grid[idx, :, 2] = span_station

        lower_grid[idx, :, 0] = le_x + cos_theta * x_local - sin_theta * y_lower_local
        lower_grid[idx, :, 1] = y_offset + sin_theta * x_local + cos_theta * y_lower_local
        lower_grid[idx, :, 2] = span_station

    return upper_grid, lower_grid


def _surface_envelope_curves(
    design: SectionedBWBDesignVariables,
    spec: FFDBoxSpec,
):
    upper_grid, lower_grid = build_reference_surface_grids(
        design=design,
        n_spanwise=int(spec.contour_fit_n_spanwise),
        n_chordwise=int(spec.contour_fit_n_chordwise),
        full_wing=False,
    )
    surface_points = np.concatenate([upper_grid, lower_grid], axis=1)
    span_dense = upper_grid[:, 0, 2].astype(float)
    x_min = np.min(surface_points[:, :, 0], axis=1)
    x_max = np.max(surface_points[:, :, 0], axis=1)
    local_chord = x_max - x_min
    span_max = float(span_dense[-1])
    relief_width = max(0.8, 0.08 * span_max)
    root_tip_relief = np.exp(-((span_dense - span_dense[0]) / relief_width) ** 2)
    root_tip_relief += np.exp(-((span_dense - span_max) / relief_width) ** 2)
    inboard_transition_relief = np.exp(-((span_dense - 0.18 * span_max) / max(0.8, 0.06 * span_max)) ** 2)

    streamwise_eta = cosine_spacing(int(spec.n_streamwise))
    chord_indices = np.clip(
        np.round(streamwise_eta * float(upper_grid.shape[1] - 1)).astype(int),
        0,
        upper_grid.shape[1] - 1,
    )

    x_center_curves = []
    y_lower_curves = []
    y_upper_curves = []
    for i_idx, chord_idx in enumerate(chord_indices):
        eta = float(streamwise_eta[i_idx])
        x_upper = upper_grid[:, chord_idx, 0]
        x_lower = lower_grid[:, chord_idx, 0]
        y_upper_surface = np.maximum(upper_grid[:, chord_idx, 1], lower_grid[:, chord_idx, 1])
        y_lower_surface = np.minimum(upper_grid[:, chord_idx, 1], lower_grid[:, chord_idx, 1])
        local_thickness = y_upper_surface - y_lower_surface

        x_center = 0.5 * (x_upper + x_lower)
        x_center = x_center - float(spec.x_margin_le_m) * (1.0 - eta) ** 2
        x_center = x_center + float(spec.x_margin_te_m) * eta**2
        x_center = x_center - 0.020 * (1.0 - eta) * local_chord
        x_center = x_center + (0.020 + 0.065 * inboard_transition_relief) * eta * local_chord

        y_pad = (
            float(spec.y_margin_m)
            + 0.10 * local_thickness
            + 0.02 * root_tip_relief * local_chord
            + 0.012 * inboard_transition_relief * local_chord
        )
        y_lower = y_lower_surface - y_pad
        y_upper = y_upper_surface + y_pad

        x_center_curves.append(_build_pyspline_scalar_interpolant(span_dense, x_center))
        y_lower_curves.append(_build_pyspline_scalar_interpolant(span_dense, y_lower))
        y_upper_curves.append(_build_pyspline_scalar_interpolant(span_dense, y_upper))

    return {
        "span_dense": span_dense,
        "x_center_curves": x_center_curves,
        "y_lower_curves": y_lower_curves,
        "y_upper_curves": y_upper_curves,
    }


def build_reference_surface_pointset(
    design: Optional[SectionedBWBDesignVariables] = None,
    n_spanwise: int = 17,
    n_chordwise: int = 81,
    full_wing: bool = False,
) -> np.ndarray:
    upper_grid, lower_grid = build_reference_surface_grids(
        design=design,
        n_spanwise=n_spanwise,
        n_chordwise=n_chordwise,
        full_wing=full_wing,
    )
    return np.vstack([upper_grid.reshape(-1, 3), lower_grid.reshape(-1, 3)])


def reshape_surface_pointset(
    pointset: np.ndarray,
    n_spanwise: int,
    n_chordwise: int,
) -> tuple[np.ndarray, np.ndarray]:
    pointset = np.asarray(pointset, dtype=float)
    expected_count = 2 * int(n_spanwise) * int(n_chordwise)
    if pointset.shape != (expected_count, 3):
        raise ValueError(
            f"Expected pointset shape {(expected_count, 3)}, got {pointset.shape}"
        )
    half = expected_count // 2
    upper = pointset[:half].reshape(int(n_spanwise), int(n_chordwise), 3)
    lower = pointset[half:].reshape(int(n_spanwise), int(n_chordwise), 3)
    return upper, lower


def build_semispan_ffd_box_points(
    design: SectionedBWBDesignVariables,
    spec: FFDBoxSpec,
) -> np.ndarray:
    config = design.to_model_config(profile_generation_mode="cst_only")
    topology = config.topology
    envelopes = _surface_envelope_curves(design, spec)

    z_span = _cosine_span_stations(topology.span, int(spec.n_spanwise), full_wing=False)
    points = np.zeros((int(spec.n_streamwise), int(spec.n_vertical), int(spec.n_spanwise), 3), dtype=float)

    for k, zz in enumerate(z_span):
        for i_idx in range(int(spec.n_streamwise)):
            x_center = float(envelopes["x_center_curves"][i_idx](zz))
            y_lower = float(envelopes["y_lower_curves"][i_idx](zz))
            y_upper = float(envelopes["y_upper_curves"][i_idx](zz))
            points[i_idx, 0, k, 0] = x_center
            points[i_idx, 1, k, 0] = x_center
            points[i_idx, 0, k, 1] = y_lower
            points[i_idx, 1, k, 1] = y_upper
        points[:, :, k, 2] = zz

    points[:, :, 0, 2] -= float(spec.z_margin_root_m)
    points[:, :, -1, 2] += float(spec.z_margin_tip_m)
    return points


def build_fullwing_ffd_box_points(
    design: SectionedBWBDesignVariables,
    spec: FFDBoxSpec,
) -> np.ndarray:
    config = design.to_model_config(profile_generation_mode="cst_only")
    topology = config.topology
    envelopes = _surface_envelope_curves(design, spec)

    z_core = _cosine_span_stations(topology.span, int(spec.n_spanwise), full_wing=True)
    z_span = np.asarray(z_core, dtype=float)
    z_span[0] -= float(spec.z_margin_tip_m)
    z_span[-1] += float(spec.z_margin_tip_m)
    points = np.zeros((int(spec.n_streamwise), int(spec.n_vertical), z_span.size, 3), dtype=float)

    for k, zz in enumerate(z_span):
        span_abs = abs(float(zz))
        for i_idx in range(int(spec.n_streamwise)):
            x_center = float(envelopes["x_center_curves"][i_idx](span_abs))
            y_lower = float(envelopes["y_lower_curves"][i_idx](span_abs))
            y_upper = float(envelopes["y_upper_curves"][i_idx](span_abs))
            points[i_idx, 0, k, 0] = x_center
            points[i_idx, 1, k, 0] = x_center
            points[i_idx, 0, k, 1] = y_lower
            points[i_idx, 1, k, 1] = y_upper
        points[:, :, k, 2] = zz

    return points


def load_dvgeometry(
    ffd_path: Path,
    rebuild_if_needed: bool = True,
):
    ensure_local_dependency_paths()
    load_pyspline_curve(rebuild_if_needed=rebuild_if_needed)
    try:
        from pygeo import DVGeometry
    except Exception as exc:
        raise RuntimeError(
            "DVGeometry could not be imported. Make sure the local pySpline/pyGeo stack is available."
        ) from exc
    return DVGeometry(str(ffd_path))


def _try_load_dvgeometry(ffd_path: Path) -> Optional[bool]:
    try:
        load_dvgeometry(ffd_path, rebuild_if_needed=True)
    except RuntimeError:
        return None
    return True


def _summary_from_points(
    out_path: Path,
    points: np.ndarray,
    full_wing: bool,
    dvgeo_ok: Optional[bool],
) -> FFDBoxSummary:
    return FFDBoxSummary(
        xyz_path=out_path,
        n_streamwise=int(points.shape[0]),
        n_vertical=int(points.shape[1]),
        n_spanwise=int(points.shape[2]),
        full_wing=bool(full_wing),
        x_min=float(np.min(points[..., 0])),
        x_max=float(np.max(points[..., 0])),
        y_min=float(np.min(points[..., 1])),
        y_max=float(np.max(points[..., 1])),
        z_min=float(np.min(points[..., 2])),
        z_max=float(np.max(points[..., 2])),
        dvgeo_load_ok=dvgeo_ok,
    )


def build_reference_ffd_box(
    out_path: Path,
    spec: Optional[FFDBoxSpec] = None,
    design: Optional[SectionedBWBDesignVariables] = None,
    full_wing: bool = False,
    check_dvgeo: bool = True,
) -> FFDBoxSummary:
    if spec is None:
        spec = FFDBoxSpec()
    if design is None:
        design = SectionedBWBDesignVariables.reference_seed()

    if full_wing:
        points = build_fullwing_ffd_box_points(design, spec)
    else:
        points = build_semispan_ffd_box_points(design, spec)
    _write_plot3d_ffd(out_path, points)
    dvgeo_ok = _try_load_dvgeometry(out_path) if check_dvgeo else None
    return _summary_from_points(out_path, points, full_wing=full_wing, dvgeo_ok=dvgeo_ok)
