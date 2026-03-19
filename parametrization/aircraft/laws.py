from dataclasses import dataclass
from enum import Enum
from typing import Tuple


class InterpolationMethod(str, Enum):
    LINEAR = "linear"
    SEGMENTED = "segmented"
    CUBIC = "cubic"
    PYSPLINE = "pyspline"


class ContinuityOrder(int, Enum):
    C0 = 0
    C1 = 1
    C2 = 2


@dataclass(frozen=True)
class InterpolationSpec:
    method: InterpolationMethod = InterpolationMethod.PYSPLINE
    continuity: ContinuityOrder = ContinuityOrder.C2
    blend_fraction: float = 0.18

    def validate(self) -> None:
        if self.method == InterpolationMethod.SEGMENTED:
            if self.continuity not in {ContinuityOrder.C1, ContinuityOrder.C2}:
                raise ValueError(
                    "segmented interpolation only supports C1 or C2 local blends; "
                    f"got continuity={self.continuity!r}"
                )
            if not (0.0 < float(self.blend_fraction) < 0.5):
                raise ValueError(
                    "blend_fraction must lie in (0, 0.5) for segmented interpolation, "
                    f"got {self.blend_fraction}"
                )
            return

        if self.method != InterpolationMethod.PYSPLINE and self.continuity != ContinuityOrder.C2:
            raise ValueError(
                "continuity is only meaningful for pySpline or segmented interpolation; "
                f"got method={self.method!r}, continuity={self.continuity!r}"
            )

    def pyspline_degree(self) -> int:
        self.validate()
        if self.method != InterpolationMethod.PYSPLINE:
            raise ValueError("pyspline_degree() is only valid for pySpline interpolation")
        return {
            ContinuityOrder.C0: 2,
            ContinuityOrder.C1: 3,
            ContinuityOrder.C2: 4,
        }[self.continuity]


@dataclass(frozen=True)
class ScalarLawAnchor:
    eta: float
    value: float

    def validate(self) -> None:
        if not (0.0 <= self.eta <= 1.0):
            raise ValueError(f"eta must lie in [0, 1], got {self.eta}")


@dataclass(frozen=True)
class ScalarLawSpec:
    name: str
    anchors: Tuple[ScalarLawAnchor, ...]
    interpolation: InterpolationSpec = InterpolationSpec()

    def validate(self) -> None:
        if not self.name:
            raise ValueError("law name must be non-empty")
        if len(self.anchors) < 2:
            raise ValueError(f"law {self.name!r} must define at least 2 anchors")
        self.interpolation.validate()

        etas = []
        for anchor in self.anchors:
            anchor.validate()
            etas.append(float(anchor.eta))

        if any(left >= right for left, right in zip(etas[:-1], etas[1:])):
            raise ValueError(
                f"law {self.name!r} anchors must be strictly increasing in eta, got {tuple(etas)}"
            )

    def anchor_etas(self) -> Tuple[float, ...]:
        self.validate()
        return tuple(float(anchor.eta) for anchor in self.anchors)

    def anchor_values(self) -> Tuple[float, ...]:
        self.validate()
        return tuple(float(anchor.value) for anchor in self.anchors)
