"""Section anchors for generic aircraft components.

Axis convention for 3D placement:

- ``x``: longitudinal / streamwise
- ``y``: vertical
- ``z``: spanwise / lateral

Rotation convention for lifting-surface export:

- ``roll_deg``  -> rotation about ``x``
- ``pitch_deg`` -> rotation about ``y``
- ``twist_deg`` -> rotation about ``z``
"""

from dataclasses import dataclass
from typing import Tuple


@dataclass(frozen=True)
class SectionPlacement:
    x: float
    y: float
    z: float
    chord: float
    twist_deg: float = 0.0
    pitch_deg: float = 0.0
    roll_deg: float = 0.0
    thickness_scale: float = 1.0

    def validate(self) -> None:
        if self.chord <= 0.0:
            raise ValueError(f"section chord must be positive, got {self.chord}")
        if self.thickness_scale <= 0.0:
            raise ValueError(
                f"section thickness_scale must be positive, got {self.thickness_scale}"
            )

    def origin_xyz(self) -> Tuple[float, float, float]:
        self.validate()
        return (float(self.x), float(self.y), float(self.z))

    def rotation_xyz_deg(self) -> Tuple[float, float, float]:
        self.validate()
        return (float(self.roll_deg), float(self.pitch_deg), float(self.twist_deg))


@dataclass(frozen=True)
class SectionAnchorSpec:
    section_id: str
    eta: float
    profile_id: str
    placement: SectionPlacement

    def validate(self) -> None:
        if not self.section_id:
            raise ValueError("section_id must be non-empty")
        if not (0.0 <= self.eta <= 1.0):
            raise ValueError(f"section {self.section_id!r} eta must lie in [0, 1], got {self.eta}")
        if not self.profile_id:
            raise ValueError(f"section {self.section_id!r} profile_id must be non-empty")
        self.placement.validate()
