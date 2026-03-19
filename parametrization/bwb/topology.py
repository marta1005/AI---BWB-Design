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
    anchor_y: Optional[Tuple[float, ...]] = None
    topology_name: str = "cbs_bwb_v4"

    @property
    def b_ratio_sum(self) -> float:
        return float(self.b1_span_ratio + self.b2_span_ratio + self.b3_span_ratio)

    @property
    def y1(self) -> float:
        return float(self.span * self.b1_span_ratio)

    @property
    def y2(self) -> float:
        return float(self.span * (self.b1_span_ratio + self.b2_span_ratio))

    @property
    def y_sections(self) -> Tuple[float, ...]:
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

        ratios = (self.b1_span_ratio, self.b2_span_ratio, self.b3_span_ratio)
        if not all(value > 0.0 for value in ratios):
            raise ValueError(f"B ratios must be positive, got {ratios}")
        if abs(self.b_ratio_sum - 1.0) > 1e-9:
            raise ValueError(
                "b1_span_ratio + b2_span_ratio + b3_span_ratio must equal 1.0 "
                f"for semispan-based segment ratios, got {self.b_ratio_sum:.12f}"
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
