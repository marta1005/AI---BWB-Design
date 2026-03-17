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


@dataclass
class ConstraintSample:
    y: float
    max_tc: float
    x_tmax: float
    te_thickness: float
    min_inner_tc: float
    min_inner_tc_xc: float
    max_camber: float
    reference_tc_max: float
    reference_x_tmax: float
    reference_te_thickness: float


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

    for yy in sample_y:
        params = section_model.params_at_y(float(yy))
        if not np.all(np.isfinite(params.coeffs)):
            raise ValueError(f"non-finite CST coefficients at y={yy:.6f}")
        if not np.isfinite(params.te_thickness) or params.te_thickness < 0.0:
            raise ValueError(f"invalid trailing-edge thickness at y={yy:.6f}: {params.te_thickness}")
        if not np.isfinite(params.x_tmax):
            raise ValueError(f"invalid x_tmax at y={yy:.6f}: {params.x_tmax}")

        yu, yl, _ = section_model.coordinates_at_y(float(yy))
        if not np.all(np.isfinite(yu)) or not np.all(np.isfinite(yl)):
            raise ValueError(f"non-finite airfoil coordinates at y={yy:.6f}")

        metrics, _ = section_model.geometry_metrics_at_y(float(yy))
        min_local_tc = metrics.min_inner_tc
        if min_local_tc <= 0.0:
            raise ValueError(
                f"thickness collapse at y={yy:.6f}, x/c={metrics.min_inner_tc_xc:.6f}, "
                f"t/c={min_local_tc:.6e}"
            )

        if min_local_tc < global_min_tc:
            global_min_tc = min_local_tc
            global_min_y = float(yy)
            global_min_xc = float(metrics.min_inner_tc_xc)

    return ValidationSummary(
        min_inner_tc=global_min_tc,
        min_inner_tc_y=global_min_y,
        min_inner_tc_xc=global_min_xc,
        num_samples=int(sample_y.size),
    )


def evaluate_section_constraints(
    section_model: SectionModel,
    sample_y: np.ndarray,
) -> list[ConstraintSample]:
    samples: list[ConstraintSample] = []
    for yy in np.asarray(sample_y, dtype=float):
        metrics, params = section_model.geometry_metrics_at_y(float(yy))
        samples.append(
            ConstraintSample(
                y=float(yy),
                max_tc=metrics.max_tc,
                x_tmax=metrics.x_tmax,
                te_thickness=metrics.te_thickness,
                min_inner_tc=metrics.min_inner_tc,
                min_inner_tc_xc=metrics.min_inner_tc_xc,
                max_camber=metrics.max_camber,
                reference_tc_max=float(params.tc_max),
                reference_x_tmax=float(params.x_tmax),
                reference_te_thickness=float(params.te_thickness),
            )
        )
    return samples
