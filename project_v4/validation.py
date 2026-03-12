from dataclasses import dataclass

import numpy as np

from .sections import SectionModel
from .topology import SectionedBWBTopologySpec


@dataclass
class LoftDefinition:
    span_stations: np.ndarray
    leading_edge_x: np.ndarray
    chord: np.ndarray
    vertical_y: np.ndarray
    span_z: np.ndarray
    twist_deg: np.ndarray
    offset: np.ndarray


@dataclass
class ValidationSummary:
    min_inner_tc: float
    min_inner_tc_y: float
    min_inner_tc_xc: float
    num_samples: int


def validate_loft_definition(loft: LoftDefinition) -> None:
    if not np.all(np.isfinite(loft.leading_edge_x)):
        raise ValueError("leading_edge_x contains non-finite values")
    if not np.all(np.isfinite(loft.chord)):
        raise ValueError("chord contains non-finite values")

    invalid = np.where(loft.chord <= 0.0)[0]
    if invalid.size:
        idx = int(invalid[0])
        raise ValueError(
            f"non-positive chord at y={loft.span_stations[idx]:.6f}: chord={loft.chord[idx]:.6f}"
        )


def build_validation_stations(
    topology: SectionedBWBTopologySpec,
    loft: LoftDefinition,
) -> np.ndarray:
    dense_count = max(101, 5 * loft.span_stations.size)
    dense_stations = np.linspace(0.0, topology.span, dense_count)
    stations = np.unique(
        np.concatenate(
            [
                dense_stations,
                topology.y_sections_array,
                topology.anchor_y_array,
                loft.span_stations,
            ]
        )
    )
    return stations.astype(float)


def validate_section_geometry(
    section_model: SectionModel,
    sample_y: np.ndarray,
) -> ValidationSummary:
    global_min_tc = np.inf
    global_min_y = np.nan
    global_min_xc = np.nan
    valid_indices = np.where(section_model.mask_valid_window)[0]

    for yy in sample_y:
        params = section_model.params_at_y(float(yy))
        if not np.all(np.isfinite(params.coeffs)):
            raise ValueError(f"non-finite CST coefficients at y={yy:.6f}")
        if not np.isfinite(params.te_thickness) or params.te_thickness < 0.0:
            raise ValueError(f"invalid trailing-edge thickness at y={yy:.6f}: {params.te_thickness}")
        if not np.isfinite(params.x_tmax):
            raise ValueError(f"invalid x_tmax at y={yy:.6f}: {params.x_tmax}")

        yu, yl = section_model.shape.evaluate(
            section_model.x_air,
            params.coeffs,
            te_thickness=params.te_thickness,
            tc_target=params.tc_max,
            x_tmax=params.x_tmax,
        )
        if not np.all(np.isfinite(yu)) or not np.all(np.isfinite(yl)):
            raise ValueError(f"non-finite airfoil coordinates at y={yy:.6f}")

        thickness = yu - yl
        min_local_tc = float(np.min(thickness[section_model.mask_valid_window]))
        if min_local_tc <= 0.0:
            bad_local_idx = int(np.argmin(thickness[section_model.mask_valid_window]))
            bad_idx = int(valid_indices[bad_local_idx])
            raise ValueError(
                f"thickness collapse at y={yy:.6f}, x/c={section_model.x_air[bad_idx]:.6f}, "
                f"t/c={min_local_tc:.6e}"
            )

        if min_local_tc < global_min_tc:
            bad_local_idx = int(np.argmin(thickness[section_model.mask_valid_window]))
            bad_idx = int(valid_indices[bad_local_idx])
            global_min_tc = min_local_tc
            global_min_y = float(yy)
            global_min_xc = float(section_model.x_air[bad_idx])

    return ValidationSummary(
        min_inner_tc=global_min_tc,
        min_inner_tc_y=global_min_y,
        min_inner_tc_xc=global_min_xc,
        num_samples=int(sample_y.size),
    )
