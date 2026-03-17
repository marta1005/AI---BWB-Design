from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..design_variables import SectionedBWBDesignVariables
from .ffd_box import (
    FFDBoxSpec,
    FFDBoxSummary,
    build_reference_ffd_box,
    build_reference_surface_grids,
    build_reference_surface_pointset,
    load_dvgeometry,
    read_plot3d_ffd,
    reshape_surface_pointset,
)


@dataclass(frozen=True)
class DVGeoReferenceCase:
    ffd_summary: FFDBoxSummary
    pointset_name: str
    pointset: np.ndarray
    upper_grid: np.ndarray
    lower_grid: np.ndarray


def create_reference_dvgeo_case(
    output_ffd_path: Path,
    design: SectionedBWBDesignVariables | None = None,
    ffd_spec: FFDBoxSpec | None = None,
    full_wing: bool = False,
    surface_n_spanwise: int = 19,
    surface_n_chordwise: int = 81,
    pointset_name: str = "wing_surface",
):
    if design is None:
        design = SectionedBWBDesignVariables.reference_seed()
    if ffd_spec is None:
        ffd_spec = FFDBoxSpec()

    summary = build_reference_ffd_box(
        out_path=output_ffd_path,
        spec=ffd_spec,
        design=design,
        full_wing=full_wing,
        check_dvgeo=False,
    )
    dvgeo = load_dvgeometry(output_ffd_path, rebuild_if_needed=True)

    upper_grid, lower_grid = build_reference_surface_grids(
        design=design,
        n_spanwise=int(surface_n_spanwise),
        n_chordwise=int(surface_n_chordwise),
        full_wing=full_wing,
    )
    pointset = build_reference_surface_pointset(
        design=design,
        n_spanwise=int(surface_n_spanwise),
        n_chordwise=int(surface_n_chordwise),
        full_wing=full_wing,
    )
    dvgeo.addPointSet(pointset, pointset_name)
    return DVGeoReferenceCase(
        ffd_summary=summary,
        pointset_name=pointset_name,
        pointset=pointset,
        upper_grid=upper_grid,
        lower_grid=lower_grid,
    ), dvgeo


def write_dvgeo_baseline_outputs(
    dvgeo,
    pointset_name: str,
    output_dir: Path,
    stem: str,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    dvgeo.writeTecplot(str(output_dir / f"{stem}_ffd.dat"))
    dvgeo.writePointSet(pointset_name, str(output_dir / f"{stem}_surface"))


def extract_ffd_control_lattice(dvgeo, volume_index: int = 0) -> np.ndarray:
    local_index = dvgeo.getLocalIndex(int(volume_index))
    return np.asarray(dvgeo.FFD.coef[local_index], dtype=float)


def build_demo_local_dv_target(
    dvgeo,
    dv_scale: float = 0.68,
) -> dict[str, np.ndarray]:
    current_values = dvgeo.getValues()
    if "shape_up" not in current_values:
        dvgeo.addLocalDV("shape_up", lower=-1.60, upper=1.60, axis="y", scale=1.0)
        current_values = dvgeo.getValues()
    if "shape_aft" not in current_values:
        dvgeo.addLocalDV("shape_aft", lower=-1.60, upper=1.60, axis="x", scale=1.0)
        current_values = dvgeo.getValues()

    dv_values_up = np.zeros_like(current_values["shape_up"])
    dv_values_aft = np.zeros_like(current_values["shape_aft"])
    local_index = dvgeo.getLocalIndex(0)
    control_points = np.asarray(dvgeo.FFD.coef[local_index], dtype=float)
    z_span = np.abs(control_points[0, 0, :, 2])
    z_max = max(float(np.max(z_span)), 1.0e-12)
    z_norm = z_span / z_max

    root_start = 0.03
    span_weight = np.clip((z_norm - root_start) / max(1.0e-12, 1.0 - root_start), 0.0, 1.0)
    span_weight = span_weight ** 1.15
    tip_relax = 1.0 - 0.08 * np.exp(-((1.0 - z_norm) / 0.08) ** 2)
    span_weight *= tip_relax

    for k_idx, weight in enumerate(span_weight):
        if float(weight) < 1.0e-6:
            continue
        up_value = 0.95 * float(dv_scale) * float(weight)
        aft_value = 1.45 * float(dv_scale) * float(weight)
        dv_values_up[local_index[:, :, k_idx]] += up_value
        dv_values_aft[local_index[:, :, k_idx]] += aft_value
    return {"shape_up": dv_values_up, "shape_aft": dv_values_aft}


def evaluate_demo_local_deformation(
    dvgeo,
    pointset_name: str,
    target_shape_values: dict[str, np.ndarray],
    fraction: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    scaled_values = {
        key: float(fraction) * np.asarray(values, dtype=float)
        for key, values in target_shape_values.items()
    }
    dvgeo.setDesignVars(scaled_values)
    deformed_pointset = dvgeo.update(pointset_name)
    return deformed_pointset, extract_ffd_control_lattice(dvgeo)


def apply_demo_local_deformation(
    dvgeo,
    pointset_name: str,
    output_dir: Path,
    stem: str,
    dv_scale: float = 0.68,
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    target_shape_values = build_demo_local_dv_target(dvgeo, dv_scale=dv_scale)

    deformed_pointset, deformed_ffd_points = evaluate_demo_local_deformation(
        dvgeo,
        pointset_name=pointset_name,
        target_shape_values=target_shape_values,
        fraction=1.0,
    )
    deformed_ffd_path = output_dir / f"{stem}_ffd.xyz"
    dvgeo.writePlot3d(str(deformed_ffd_path))
    dvgeo.writeTecplot(str(output_dir / f"{stem}_ffd.dat"))
    dvgeo.writePointSet(pointset_name, str(output_dir / f"{stem}_surface"))
    return deformed_pointset, deformed_ffd_points, target_shape_values


def pointset_to_surface_grids(
    pointset: np.ndarray,
    reference_upper_grid: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    return reshape_surface_pointset(
        pointset,
        n_spanwise=reference_upper_grid.shape[0],
        n_chordwise=reference_upper_grid.shape[1],
    )
