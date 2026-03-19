from dataclasses import dataclass

import numpy as np

from .planform import SectionedPlanform
from .sections import SectionModel
from .specs import SectionedBWBModelConfig


@dataclass
class VolumeConstraintSummary:
    enabled: bool
    satisfied: bool
    enclosed_volume_m3: float
    required_volume_m3: float
    volume_margin_m3: float
    volume_ratio: float
    mean_cross_section_area_m2: float
    max_cross_section_area_m2: float
    max_cross_section_area_y: float
    span_samples: int


def _integrate(values: np.ndarray, coordinates: np.ndarray) -> float:
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(values, coordinates))
    return float(np.trapz(values, coordinates))


def evaluate_volume_constraint(
    config: SectionedBWBModelConfig,
    planform: SectionedPlanform,
    section_model: SectionModel,
) -> VolumeConstraintSummary:
    volume = config.volume
    if not volume.enabled:
        return VolumeConstraintSummary(
            enabled=False,
            satisfied=True,
            enclosed_volume_m3=0.0,
            required_volume_m3=0.0,
            volume_margin_m3=0.0,
            volume_ratio=1.0,
            mean_cross_section_area_m2=0.0,
            max_cross_section_area_m2=0.0,
            max_cross_section_area_y=0.0,
            span_samples=0,
        )

    sample_count = max(int(volume.span_samples), 5 * int(config.sampling.num_base_stations))
    y_samples = np.linspace(0.0, config.topology.span, sample_count, dtype=float)
    local_section_area = np.zeros_like(y_samples)

    for idx, yy in enumerate(y_samples):
        le_x = float(planform.le_x(float(yy)))
        te_x = float(planform.te_x(float(yy)))
        chord = te_x - le_x
        if chord <= 1e-12:
            raise ValueError(
                f"invalid chord while evaluating volume at y={yy:.6f}: chord={chord:.6e}"
            )

        yu, yl, _ = section_model.coordinates_at_y(float(yy))
        thickness = yu - yl
        area_coeff = _integrate(thickness, section_model.x_air)
        local_section_area[idx] = float(area_coeff * chord * chord)

    one_side_volume = _integrate(local_section_area, y_samples)
    enclosed_volume = 2.0 * one_side_volume
    mean_cross_section_area = enclosed_volume / max(2.0 * config.topology.span, 1e-12)
    max_idx = int(np.argmax(local_section_area))
    max_cross_section_area = 2.0 * float(local_section_area[max_idx])
    max_cross_section_area_y = float(y_samples[max_idx])
    volume_margin = float(enclosed_volume - volume.required_volume_m3)
    volume_ratio = float(enclosed_volume / max(volume.required_volume_m3, 1e-12))

    return VolumeConstraintSummary(
        enabled=True,
        satisfied=bool(volume_margin >= 0.0),
        enclosed_volume_m3=float(enclosed_volume),
        required_volume_m3=float(volume.required_volume_m3),
        volume_margin_m3=volume_margin,
        volume_ratio=volume_ratio,
        mean_cross_section_area_m2=float(mean_cross_section_area),
        max_cross_section_area_m2=max_cross_section_area,
        max_cross_section_area_y=max_cross_section_area_y,
        span_samples=int(sample_count),
    )
