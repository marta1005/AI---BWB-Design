from __future__ import annotations

from dataclasses import dataclass, field
from math import radians, tan
from typing import Dict, Tuple

from .components import ComponentKind, LoftedComponentSpec
from .laws import ContinuityOrder, InterpolationMethod, InterpolationSpec, ScalarLawAnchor, ScalarLawSpec
from .sections import SectionAnchorSpec, SectionPlacement


@dataclass(frozen=True)
class WingStationSpec:
    station_id: str
    eta: float
    profile_id: str
    chord: float
    twist_deg: float = 0.0
    pitch_deg: float = 0.0
    roll_deg: float = 0.0
    thickness_scale: float = 1.0
    x_le: float | None = None
    vertical_y: float | None = None
    sweep_le_deg: float | None = None
    dihedral_deg: float | None = None

    def validate(self) -> None:
        if not self.station_id:
            raise ValueError("station_id must be non-empty")
        if not (0.0 <= self.eta <= 1.0):
            raise ValueError(f"station {self.station_id!r} eta must lie in [0, 1], got {self.eta}")
        if not self.profile_id:
            raise ValueError(f"station {self.station_id!r} profile_id must be non-empty")
        if self.chord <= 0.0:
            raise ValueError(f"station {self.station_id!r} chord must be positive, got {self.chord}")
        if self.thickness_scale <= 0.0:
            raise ValueError(
                f"station {self.station_id!r} thickness_scale must be positive, got {self.thickness_scale}"
            )
        if self.x_le is not None and self.sweep_le_deg is not None:
            raise ValueError(
                f"station {self.station_id!r} cannot define both x_le and sweep_le_deg; "
                "use one absolute or one segment-based layout input"
            )
        if self.vertical_y is not None and self.dihedral_deg is not None:
            raise ValueError(
                f"station {self.station_id!r} cannot define both vertical_y and dihedral_deg; "
                "use one absolute or one segment-based layout input"
            )


@dataclass(frozen=True)
class WingSpec:
    wing_id: str
    semispan: float
    stations: Tuple[WingStationSpec, ...]
    root_x: float = 0.0
    root_y: float = 0.0
    root_z: float = 0.0
    section_interpolation: InterpolationSpec = InterpolationSpec(
        method=InterpolationMethod.SEGMENTED,
        continuity=ContinuityOrder.C2,
        blend_fraction=0.18,
    )
    spine_interpolation: InterpolationSpec = InterpolationSpec(
        method=InterpolationMethod.SEGMENTED,
        continuity=ContinuityOrder.C2,
        blend_fraction=0.18,
    )
    scalar_laws: Tuple[ScalarLawSpec, ...] = ()
    mirrored: bool = False
    metadata: Dict[str, str] = field(default_factory=dict)

    def validate(self) -> None:
        if not self.wing_id:
            raise ValueError("wing_id must be non-empty")
        if self.semispan <= 0.0:
            raise ValueError(f"semispan must be positive, got {self.semispan}")
        if len(self.stations) < 2:
            raise ValueError("a wing must define at least 2 stations")

        self.section_interpolation.validate()
        self.spine_interpolation.validate()

        seen = set()
        etas = []
        for station in self.stations:
            station.validate()
            if station.station_id in seen:
                raise ValueError(f"duplicate wing station_id detected: {station.station_id!r}")
            seen.add(station.station_id)
            etas.append(float(station.eta))

        if any(left >= right for left, right in zip(etas[:-1], etas[1:])):
            raise ValueError(f"wing stations must be strictly increasing in eta, got {tuple(etas)}")
        if abs(etas[0]) > 1e-12:
            raise ValueError(
                f"the first wing station must be at eta=0.0 so root_x/root_y/root_z are well-defined, got {etas[0]}"
            )
        if abs(etas[-1] - 1.0) > 1e-12:
            raise ValueError(
                f"the last wing station must be at eta=1.0 so semispan reaches the tip, got {etas[-1]}"
            )

        root_station = self.stations[0]
        if root_station.x_le is not None or root_station.sweep_le_deg is not None:
            raise ValueError(
                "the root station layout is defined by root_x; do not set x_le or sweep_le_deg on the first station"
            )
        if root_station.vertical_y is not None or root_station.dihedral_deg is not None:
            raise ValueError(
                "the root station layout is defined by root_y; do not set vertical_y or dihedral_deg on the first station"
            )

    def _resolved_station_placements(self) -> Tuple[SectionPlacement, ...]:
        self.validate()

        placements = []
        prev_x = float(self.root_x)
        prev_y = float(self.root_y)
        prev_z = float(self.root_z)

        for index, station in enumerate(self.stations):
            z = float(self.root_z) + float(self.semispan) * float(station.eta)

            if index == 0:
                x = float(self.root_x)
                y = float(self.root_y)
            else:
                dz = z - prev_z
                if station.x_le is not None:
                    x = float(station.x_le)
                elif station.sweep_le_deg is not None:
                    x = prev_x + dz * tan(radians(float(station.sweep_le_deg)))
                else:
                    x = prev_x

                if station.vertical_y is not None:
                    y = float(station.vertical_y)
                elif station.dihedral_deg is not None:
                    y = prev_y + dz * tan(radians(float(station.dihedral_deg)))
                else:
                    y = prev_y

            placements.append(
                SectionPlacement(
                    x=x,
                    y=y,
                    z=z,
                    chord=float(station.chord),
                    twist_deg=float(station.twist_deg),
                    pitch_deg=float(station.pitch_deg),
                    roll_deg=float(station.roll_deg),
                    thickness_scale=float(station.thickness_scale),
                )
            )
            prev_x = x
            prev_y = y
            prev_z = z

        return tuple(placements)

    def to_component_spec(self) -> LoftedComponentSpec:
        placements = self._resolved_station_placements()
        sections = tuple(
            SectionAnchorSpec(
                section_id=station.station_id,
                eta=float(station.eta),
                profile_id=station.profile_id,
                placement=placement,
            )
            for station, placement in zip(self.stations, placements)
        )
        return LoftedComponentSpec(
            component_id=self.wing_id,
            kind=ComponentKind.LIFTING_SURFACE,
            sections=sections,
            section_interpolation=self.section_interpolation,
            spine_interpolation=self.spine_interpolation,
            scalar_laws=self.scalar_laws,
            mirrored=self.mirrored,
            metadata=dict(self.metadata),
        )

    def _mirrored_scalar_law(self, law: ScalarLawSpec) -> ScalarLawSpec:
        law.validate()
        if law.name == "z":
            raise ValueError(
                "full-span symmetric wing construction cannot mirror an explicit 'z' scalar law; "
                "the spanwise coordinate is generated from symmetry"
            )

        anchors = tuple(law.anchors)
        if abs(float(anchors[0].eta)) > 1e-12 or abs(float(anchors[-1].eta) - 1.0) > 1e-12:
            raise ValueError(
                f"symmetric wing law {law.name!r} must span the full half-wing range [0, 1], "
                f"got [{anchors[0].eta}, {anchors[-1].eta}]"
            )

        left = tuple(
            ScalarLawAnchor(eta=0.5 * (1.0 - float(anchor.eta)), value=float(anchor.value))
            for anchor in reversed(anchors[1:])
        )
        center = (ScalarLawAnchor(eta=0.5, value=float(anchors[0].value)),)
        right = tuple(
            ScalarLawAnchor(eta=0.5 * (1.0 + float(anchor.eta)), value=float(anchor.value))
            for anchor in anchors[1:]
        )
        return ScalarLawSpec(
            name=law.name,
            anchors=left + center + right,
            interpolation=law.interpolation,
        )

    def to_full_span_component_spec(self) -> LoftedComponentSpec:
        placements = self._resolved_station_placements()
        paired = tuple(zip(self.stations, placements))

        left_sections = tuple(
            SectionAnchorSpec(
                section_id=f"left_{station.station_id}",
                eta=0.5 * (1.0 - float(station.eta)),
                profile_id=station.profile_id,
                placement=SectionPlacement(
                    x=float(placement.x),
                    y=float(placement.y),
                    z=2.0 * float(self.root_z) - float(placement.z),
                    chord=float(placement.chord),
                    twist_deg=float(placement.twist_deg),
                    pitch_deg=float(placement.pitch_deg),
                    roll_deg=float(placement.roll_deg),
                    thickness_scale=float(placement.thickness_scale),
                ),
            )
            for station, placement in reversed(paired[1:])
        )

        root_station, root_placement = paired[0]
        center_section = (
            SectionAnchorSpec(
                section_id=f"center_{root_station.station_id}",
                eta=0.5,
                profile_id=root_station.profile_id,
                placement=SectionPlacement(
                    x=float(root_placement.x),
                    y=float(root_placement.y),
                    z=float(root_placement.z),
                    chord=float(root_placement.chord),
                    twist_deg=float(root_placement.twist_deg),
                    pitch_deg=float(root_placement.pitch_deg),
                    roll_deg=float(root_placement.roll_deg),
                    thickness_scale=float(root_placement.thickness_scale),
                ),
            ),
        )

        right_sections = tuple(
            SectionAnchorSpec(
                section_id=f"right_{station.station_id}",
                eta=0.5 * (1.0 + float(station.eta)),
                profile_id=station.profile_id,
                placement=SectionPlacement(
                    x=float(placement.x),
                    y=float(placement.y),
                    z=float(placement.z),
                    chord=float(placement.chord),
                    twist_deg=float(placement.twist_deg),
                    pitch_deg=float(placement.pitch_deg),
                    roll_deg=float(placement.roll_deg),
                    thickness_scale=float(placement.thickness_scale),
                ),
            )
            for station, placement in paired[1:]
        )

        return LoftedComponentSpec(
            component_id=self.wing_id,
            kind=ComponentKind.LIFTING_SURFACE,
            sections=left_sections + center_section + right_sections,
            section_interpolation=self.section_interpolation,
            spine_interpolation=self.spine_interpolation,
            scalar_laws=tuple(self._mirrored_scalar_law(law) for law in self.scalar_laws),
            mirrored=False,
            metadata={**dict(self.metadata), "symmetry": "full_span"},
        )
