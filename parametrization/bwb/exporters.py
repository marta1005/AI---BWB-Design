from typing import List, Optional

import numpy as np

from parametrization.shared.airfoil_io import write_airfoil_dat
from parametrization.shared.dependency_setup import ensure_local_dependency_paths, load_pygeo_class
from parametrization.shared.pyspline_shim import patch_pyspline_for_pygeo
from .specs import SectionedBWBModelConfig
from .validation import LoftDefinition

def build_te_height_scaled(section_model, loft: LoftDefinition) -> np.ndarray:
    te_height_scaled = np.zeros(loft.span_stations.size, dtype=float)
    for idx, yy in enumerate(loft.span_stations):
        params = section_model.params_at_y(float(yy))
        te_height_scaled[idx] = float(params.te_thickness)
    return te_height_scaled


def write_station_airfoils(
    config: SectionedBWBModelConfig,
    section_model,
    loft: LoftDefinition,
) -> List[Optional[str]]:
    out_dir = config.export.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    if config.sampling.airfoil_distribution_mode == "all":
        airfoil_list: List[Optional[str]] = []
        for idx, yy in enumerate(loft.span_stations):
            yu, yl, _ = section_model.coordinates_at_y(float(yy))
            path = out_dir / f"station_{idx:03d}_y{yy:.3f}.dat"
            write_airfoil_dat(str(path), section_model.x_air, yu, yl, name=f"CST_STATION_{idx:03d}")
            airfoil_list.append(str(path))
        return airfoil_list

    airfoil_list = [None] * loft.span_stations.size
    for idx, yy in enumerate(config.topology.anchor_y_array):
        station_idx = int(np.argmin(np.abs(loft.span_stations - yy)))
        station_distance = abs(float(loft.span_stations[station_idx] - yy))
        if station_distance > 1e-9:
            raise ValueError(
                f"anchor station y={yy:.6f} is not represented in span_stations; "
                f"closest station is {loft.span_stations[station_idx]:.12f}"
            )

        yu, yl, _ = section_model.coordinates_at_y(float(yy))
        path = out_dir / f"anchor_{idx:02d}_y{yy:.3f}.dat"
        write_airfoil_dat(str(path), section_model.x_air, yu, yl, name=f"CST_ANCHOR_{idx:02d}")
        airfoil_list[station_idx] = str(path)
    return airfoil_list


def build_pygeo_surface(
    config: SectionedBWBModelConfig,
    loft: LoftDefinition,
    airfoil_list: List[Optional[str]],
    te_height_scaled: Optional[np.ndarray] = None,
):
    try:
        ensure_local_dependency_paths()
        pyGeo = load_pygeo_class(rebuild_if_needed=True)
    except Exception as exc:
        detail = ""
        exc_text = str(exc)
        if "invalid ELF header" in exc_text or "Mach-O" in exc_text:
            detail = (
                " Detected platform mismatch in local binaries (for example macOS binaries used in Linux)."
            )
        raise RuntimeError(
            "pyGeo export is unavailable because the local pygeo/pyspline stack could not be imported. "
            "The parametric geometry can still be prepared, but IGES export needs a working pygeo installation."
            f"{detail}"
        ) from exc

    xsections = airfoil_list
    scale = np.asarray(loft.chord, dtype=float)
    offset = np.asarray(loft.offset, dtype=float)
    x = np.asarray(loft.leading_edge_x, dtype=float)
    y = np.asarray(loft.vertical_y, dtype=float)
    z = np.asarray(loft.span_z, dtype=float)
    rot_z = np.asarray(loft.twist_deg, dtype=float)
    te_height_local = None if te_height_scaled is None else np.asarray(te_height_scaled, dtype=float)

    if config.export.symmetric:
        mirror_slice = slice(1, None)
        xsections = list(xsections[mirror_slice][::-1]) + list(xsections)
        scale = np.concatenate([scale[mirror_slice][::-1], scale], axis=0)
        offset = np.vstack([offset[mirror_slice][::-1], offset])
        x = np.concatenate([x[mirror_slice][::-1], x], axis=0)
        y = np.concatenate([y[mirror_slice][::-1], y], axis=0)
        z = np.concatenate([-z[mirror_slice][::-1], z], axis=0)
        rot_z = np.concatenate([rot_z[mirror_slice][::-1], rot_z], axis=0)
        if te_height_local is not None:
            te_height_local = np.concatenate([te_height_local[mirror_slice][::-1], te_height_local], axis=0)

    kwargs = {}
    if config.export.blunt_te:
        if te_height_local is None:
            raise ValueError("blunt_te=True requires te_height_scaled for pyGeo export")
        kwargs["teHeightScaled"] = te_height_local.tolist()

    return pyGeo(
        "liftingSurface",
        xsections=xsections,
        scale=scale,
        offset=offset,
        x=x,
        y=y,
        z=z,
        rotZ=rot_z,
        nCtl=config.sampling.section_curve_n_ctl,
        kSpan=config.sampling.k_span,
        bluntTe=config.export.blunt_te,
        tip=config.export.tip_style,
        **kwargs,
    )


def _load_pyspline_surface_class():
    ensure_local_dependency_paths()
    patch_pyspline_for_pygeo()
    from pyspline import Surface

    return Surface


def _build_planar_surface_patch(
    x_outer: np.ndarray,
    y_outer: np.ndarray,
    x_inner: np.ndarray,
    y_inner: np.ndarray,
):
    Surface = _load_pyspline_surface_class()
    x_outer = np.asarray(x_outer, dtype=float)
    y_outer = np.asarray(y_outer, dtype=float)
    x_inner = np.asarray(x_inner, dtype=float)
    y_inner = np.asarray(y_inner, dtype=float)
    if x_outer.shape != x_inner.shape or y_outer.shape != y_inner.shape:
        raise ValueError("planar surface patch boundaries must share the same discretization")

    nu = int(x_outer.size)
    coef = np.zeros((nu, 2, 3), dtype=float)
    coef[:, 0, 0] = x_outer
    coef[:, 0, 1] = y_outer
    coef[:, 1, 0] = x_inner
    coef[:, 1, 1] = y_inner

    interior = np.arange(1, nu - 1, dtype=float) / max(nu - 1, 1)
    tu = np.concatenate([[0.0, 0.0], interior, [1.0, 1.0]]).astype(float)
    tv = np.asarray([0.0, 0.0, 1.0, 1.0], dtype=float)
    return Surface(
        coef=coef,
        ku=2,
        kv=2,
        tu=tu,
        tv=tv,
    )


def build_xy_symmetry_frame_surfaces(
    prepared,
    x_margin_factor: float = 0.08,
    y_margin_factor: float = 0.10,
    min_x_margin_m: float = 0.50,
    min_y_margin_m: float = 0.25,
    side_point_count: int = 121,
    x_min_m: float | None = None,
    x_max_m: float | None = None,
    y_min_m: float | None = None,
    y_max_m: float | None = None,
):
    root_span = 0.0
    x_air = np.asarray(prepared.section_model.x_air, dtype=float)
    upper, lower, _ = prepared.section_model.coordinates_at_y(root_span)
    chord = float(np.interp(root_span, prepared.loft.span_stations, prepared.loft.chord))
    le_x = float(np.interp(root_span, prepared.loft.span_stations, prepared.loft.leading_edge_x))
    vertical_y = float(np.interp(root_span, prepared.loft.span_stations, prepared.loft.vertical_y))
    twist_deg = float(prepared.spanwise_laws.twist_deg(root_span))

    x_local = x_air * chord
    y_upper_local = np.asarray(upper, dtype=float) * chord
    y_lower_local = np.asarray(lower, dtype=float) * chord
    twist_rad = np.deg2rad(twist_deg)
    cos_twist = float(np.cos(twist_rad))
    sin_twist = float(np.sin(twist_rad))

    x_upper = le_x + x_local * cos_twist - y_upper_local * sin_twist
    y_upper = vertical_y + x_local * sin_twist + y_upper_local * cos_twist
    x_lower = le_x + x_local * cos_twist - y_lower_local * sin_twist
    y_lower = vertical_y + x_local * sin_twist + y_lower_local * cos_twist

    x_min_hole = float(min(np.min(x_upper), np.min(x_lower)))
    x_max_hole = float(max(np.max(x_upper), np.max(x_lower)))
    y_min_hole = float(min(np.min(y_upper), np.min(y_lower)))
    y_max_hole = float(max(np.max(y_upper), np.max(y_lower)))
    root_height = max(y_max_hole - y_min_hole, 1.0e-6)

    x_margin = max(float(x_margin_factor) * chord, float(min_x_margin_m))
    y_margin = max(float(y_margin_factor) * root_height, float(min_y_margin_m))
    x_min = x_min_hole - x_margin if x_min_m is None else float(x_min_m)
    x_max = x_max_hole + x_margin if x_max_m is None else float(x_max_m)
    y_min = y_min_hole - y_margin if y_min_m is None else float(y_min_m)
    y_max = y_max_hole + y_margin if y_max_m is None else float(y_max_m)

    outer_top_x = np.linspace(x_min, x_max, x_upper.size, dtype=float)
    outer_top_y = np.full_like(outer_top_x, y_max)
    outer_bottom_x = np.linspace(x_min, x_max, x_lower.size, dtype=float)
    outer_bottom_y = np.full_like(outer_bottom_x, y_min)

    outer_side_y = np.linspace(y_max, y_min, int(side_point_count), dtype=float)
    outer_left_x = np.full_like(outer_side_y, x_min)
    outer_right_x = np.full_like(outer_side_y, x_max)

    inner_left_x = np.linspace(float(x_upper[0]), float(x_lower[0]), int(side_point_count), dtype=float)
    inner_left_y = np.linspace(float(y_upper[0]), float(y_lower[0]), int(side_point_count), dtype=float)
    inner_right_x = np.linspace(float(x_upper[-1]), float(x_lower[-1]), int(side_point_count), dtype=float)
    inner_right_y = np.linspace(float(y_upper[-1]), float(y_lower[-1]), int(side_point_count), dtype=float)

    return [
        _build_planar_surface_patch(outer_top_x, outer_top_y, x_upper, y_upper),
        _build_planar_surface_patch(outer_bottom_x, outer_bottom_y, x_lower, y_lower),
        _build_planar_surface_patch(outer_left_x, outer_side_y, inner_left_x, inner_left_y),
        _build_planar_surface_patch(outer_right_x, outer_side_y, inner_right_x, inner_right_y),
    ]


def append_xy_symmetry_frame_surfaces(
    geometry,
    prepared,
    **kwargs,
):
    frame_surfaces = build_xy_symmetry_frame_surfaces(prepared, **kwargs)
    geometry.surfs.extend(frame_surfaces)
    geometry.nSurf = len(geometry.surfs)
    return geometry
