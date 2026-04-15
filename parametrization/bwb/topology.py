from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np


@dataclass
class SectionedBWBTopologySpec:
    """Topology driven by three spanwise segments over the semispan."""

    span: float = 39.5
    b1_span_ratio: float = 0.18
    b2_span_ratio: float = 0.12
    b3_span_ratio: float = 0.70
    section_y: Optional[Tuple[float, ...]] = None
    anchor_y: Optional[Tuple[float, ...]] = None
    topology_name: str = "cbs_bwb_v4"

    @property
    def b_ratio_sum(self) -> float:
        return float(self.b1_span_ratio + self.b2_span_ratio + self.b3_span_ratio)

    @property
    def y1(self) -> float:
        if self.section_y is not None:
            if len(self.section_y) < 2:
                raise ValueError("section_y must contain at least two stations to resolve y1")
            return float(self.section_y[1])
        return float(self.span * self.b1_span_ratio)

    @property
    def y2(self) -> float:
        if self.section_y is not None:
            if len(self.section_y) < 3:
                raise ValueError("section_y must contain at least three stations to resolve y2")
            return float(self.section_y[2])
        return float(self.span * (self.b1_span_ratio + self.b2_span_ratio))

    @property
    def y_sections(self) -> Tuple[float, ...]:
        if self.section_y is not None:
            return tuple(float(value) for value in self.section_y)
        return (0.0, self.y1, self.y2, self.span)

    @property
    def y_sections_array(self) -> np.ndarray:
        return np.asarray(self.y_sections, dtype=float)

    @property
    def anchor_y_array(self) -> np.ndarray:
        values = self.y_sections if self.anchor_y is None else self.anchor_y
        return np.asarray(values, dtype=float)

    def validate(self) -> None:
        y_sections = self.y_sections_array
        anchor_y = self.anchor_y_array

        if self.span <= 0.0:
            raise ValueError(f"span must be positive, got {self.span:.6f}")

        if self.section_y is None:
            ratios = (self.b1_span_ratio, self.b2_span_ratio, self.b3_span_ratio)
            if not all(value > 0.0 for value in ratios):
                raise ValueError(f"B ratios must be positive, got {ratios}")
            if abs(self.b_ratio_sum - 1.0) > 1e-9:
                raise ValueError(
                    "b1_span_ratio + b2_span_ratio + b3_span_ratio must equal 1.0 "
                    f"for semispan-based segment ratios, got {self.b_ratio_sum:.12f}"
                )
        elif len(self.section_y) < 2:
            raise ValueError(
                f"section_y must contain at least 2 spanwise stations, got {self.section_y}"
            )

        if not np.all(np.diff(y_sections) > 0.0):
            raise ValueError(f"section stations must be strictly increasing, got {tuple(y_sections)}")
        if not np.isclose(y_sections[0], 0.0):
            raise ValueError(f"first section station must be 0.0, got {y_sections[0]:.6f}")
        if not np.isclose(y_sections[-1], self.span):
            raise ValueError(
                f"last section station must equal span={self.span:.6f}, got {y_sections[-1]:.6f}"
            )

        if anchor_y.size == 0 or not np.all(np.diff(anchor_y) > 0.0):
            raise ValueError("anchor_y must be strictly increasing and non-empty")
        if anchor_y[0] < 0.0 or anchor_y[-1] > self.span:
            raise ValueError("anchor_y must remain inside [0, span]")
        if not np.isclose(anchor_y[0], 0.0):
            raise ValueError(f"anchor_y must start at 0.0, got {anchor_y[0]:.6f}")
        if not np.isclose(anchor_y[-1], self.span):
            raise ValueError(
                f"anchor_y must end at span={self.span:.6f}, got {anchor_y[-1]:.6f}"
            )
