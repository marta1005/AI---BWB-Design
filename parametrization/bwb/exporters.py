from typing import List, Optional

import numpy as np

from parametrization.shared.airfoil_io import write_airfoil_dat
from parametrization.shared.dependency_setup import ensure_local_dependency_paths, load_pygeo_class
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
            path = out_dir / f"station_{idx:03d}_y{yy:.2f}.dat"
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
        path = out_dir / f"anchor_{idx:02d}_y{yy:.2f}.dat"
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
