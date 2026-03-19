from __future__ import annotations

from dataclasses import dataclass, field
from math import radians, tan
from typing import Dict, Tuple

from .components import ComponentKind, LoftedComponentSpec
from .laws import InterpolationMethod, InterpolationSpec, ScalarLawSpec
from .sections import SectionAnchorSpec, SectionPlacement


@dataclass(frozen=True)
class VerticalTailStationSpec:
    station_id: str
    eta: float
    profile_id: str
    chord: float
    twist_deg: float = 0.0
    pitch_deg: float = 0.0
    roll_deg: float = 0.0
    thickness_scale: float = 1.0
    x_le: float | None = None
    lateral_z: float | None = None
    sweep_le_deg: float | None = None
    cant_deg: float | None = None

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
                f"station {self.station_id!r} cannot define both x_le and sweep_le_deg"
            )
        if self.lateral_z is not None and self.cant_deg is not None:
            raise ValueError(
                f"station {self.station_id!r} cannot define both lateral_z and cant_deg"
            )


@dataclass(frozen=True)
class VerticalTailSpec:
    tail_id: str
    span: float
    stations: Tuple[VerticalTailStationSpec, ...]
    root_x: float = 0.0
    root_y: float = 0.0
    root_z: float = 0.0
    section_interpolation: InterpolationSpec = InterpolationSpec(method=InterpolationMethod.LINEAR)
    spine_interpolation: InterpolationSpec = InterpolationSpec(method=InterpolationMethod.LINEAR)
    scalar_laws: Tuple[ScalarLawSpec, ...] = ()
    base_roll_deg: float = -90.0
    metadata: Dict[str, str] = field(default_factory=dict)

    def validate(self) -> None:
        if not self.tail_id:
            raise ValueError("tail_id must be non-empty")
        if self.span <= 0.0:
            raise ValueError(f"span must be positive, got {self.span}")
        if len(self.stations) < 2:
            raise ValueError("a vertical tail must define at least 2 stations")
        self.section_interpolation.validate()
        self.spine_interpolation.validate()

        seen = set()
        etas = []
        for station in self.stations:
            station.validate()
            if station.station_id in seen:
                raise ValueError(f"duplicate vertical-tail station_id detected: {station.station_id!r}")
            seen.add(station.station_id)
            etas.append(float(station.eta))

        if any(left >= right for left, right in zip(etas[:-1], etas[1:])):
            raise ValueError(f"vertical-tail stations must be strictly increasing in eta, got {tuple(etas)}")
        if abs(etas[0]) > 1e-12:
            raise ValueError(f"the first vertical-tail station must be at eta=0.0, got {etas[0]}")
        if abs(etas[-1] - 1.0) > 1e-12:
            raise ValueError(f"the last vertical-tail station must be at eta=1.0, got {etas[-1]}")

    def _resolved_station_placements(self) -> Tuple[SectionPlacement, ...]:
        self.validate()

        placements = []
        prev_x = float(self.root_x)
        prev_y = float(self.root_y)
        prev_z = float(self.root_z)

        for index, station in enumerate(self.stations):
            y = float(self.root_y) + float(self.span) * float(station.eta)

            if index == 0:
                x = float(self.root_x)
                z = float(self.root_z)
            else:
                dy = y - prev_y
                if station.x_le is not None:
                    x = float(station.x_le)
                elif station.sweep_le_deg is not None:
                    x = prev_x + dy * tan(radians(float(station.sweep_le_deg)))
                else:
                    x = prev_x

                if station.lateral_z is not None:
                    z = float(station.lateral_z)
                elif station.cant_deg is not None:
                    z = prev_z + dy * tan(radians(float(station.cant_deg)))
                else:
                    z = prev_z

            placements.append(
                SectionPlacement(
                    x=x,
                    y=y,
                    z=z,
                    chord=float(station.chord),
                    twist_deg=float(station.twist_deg),
                    pitch_deg=float(station.pitch_deg),
                    roll_deg=float(self.base_roll_deg) + float(station.roll_deg),
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
            component_id=self.tail_id,
            kind=ComponentKind.LIFTING_SURFACE,
            sections=sections,
            section_interpolation=self.section_interpolation,
            spine_interpolation=self.spine_interpolation,
            scalar_laws=self.scalar_laws,
            mirrored=False,
            metadata={**dict(self.metadata), "surface_role": "vtp"},
        )
