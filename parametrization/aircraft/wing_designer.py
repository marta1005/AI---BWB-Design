from __future__ import annotations

from dataclasses import dataclass, field
from math import radians, tan
from typing import Dict, Tuple

from .components import LoftedComponentSpec
from .laws import ContinuityOrder, InterpolationMethod, InterpolationSpec, ScalarLawSpec
from .wing import WingSpec, WingStationSpec


@dataclass(frozen=True)
class WingTransitionSpec:
    method: InterpolationMethod = InterpolationMethod.SEGMENTED
    continuity: ContinuityOrder = ContinuityOrder.C2
    blend_fraction: float = 0.18

    def validate(self) -> None:
        self.to_interpolation_spec().validate()

    def to_interpolation_spec(self) -> InterpolationSpec:
        return InterpolationSpec(
            method=self.method,
            continuity=self.continuity,
            blend_fraction=self.blend_fraction,
        )


@dataclass(frozen=True)
class WingStationDesignerSpec:
    station_id: str
    eta: float
    profile_id: str
    chord_ratio: float
    twist_deg: float = 0.0
    pitch_deg: float = 0.0
    roll_deg: float = 0.0
    thickness_scale: float = 1.0
    sweep_qc_deg: float = 0.0
    dihedral_deg: float = 0.0

    def validate(self) -> None:
        if not self.station_id:
            raise ValueError("station_id must be non-empty")
        if not (0.0 <= float(self.eta) <= 1.0):
            raise ValueError(f"station {self.station_id!r} eta must lie in [0, 1], got {self.eta}")
        if not self.profile_id:
            raise ValueError(f"station {self.station_id!r} profile_id must be non-empty")
        if float(self.chord_ratio) <= 0.0:
            raise ValueError(
                f"station {self.station_id!r} chord_ratio must be positive, got {self.chord_ratio}"
            )
        if float(self.thickness_scale) <= 0.0:
            raise ValueError(
                f"station {self.station_id!r} thickness_scale must be positive, got {self.thickness_scale}"
            )


@dataclass(frozen=True)
class WingDesignerSpec:
    wing_id: str
    span: float
    root_chord: float
    stations: Tuple[WingStationDesignerSpec, ...]
    symmetric: bool = True
    root_le_x: float = 0.0
    root_vertical_y: float = 0.0
    root_z: float = 0.0
    spine_transition: WingTransitionSpec = WingTransitionSpec()
    section_transition: WingTransitionSpec | None = None
    scalar_laws: Tuple[ScalarLawSpec, ...] = ()
    metadata: Dict[str, str] = field(default_factory=dict)

    def validate(self) -> None:
        if not self.wing_id:
            raise ValueError("wing_id must be non-empty")
        if float(self.span) <= 0.0:
            raise ValueError(f"span must be positive, got {self.span}")
        if float(self.root_chord) <= 0.0:
            raise ValueError(f"root_chord must be positive, got {self.root_chord}")
        if len(self.stations) < 2:
            raise ValueError("a designer wing must define at least 2 stations")

        self.spine_transition.validate()
        if self.section_transition is not None:
            self.section_transition.validate()

        seen_ids = set()
        etas = []
        for station in self.stations:
            station.validate()
            if station.station_id in seen_ids:
                raise ValueError(f"duplicate station_id detected: {station.station_id!r}")
            seen_ids.add(station.station_id)
            etas.append(float(station.eta))

        if any(left >= right for left, right in zip(etas[:-1], etas[1:])):
            raise ValueError(
                f"designer wing stations must be strictly increasing in eta, got {tuple(etas)}"
            )
        if abs(etas[0]) > 1e-12:
            raise ValueError(f"the first designer station must be at eta=0.0, got {etas[0]}")
        if abs(etas[-1] - 1.0) > 1e-12:
            raise ValueError(f"the last designer station must be at eta=1.0, got {etas[-1]}")

        root_station = self.stations[0]
        if abs(float(root_station.chord_ratio) - 1.0) > 1e-12:
            raise ValueError(
                "the root designer station must use chord_ratio=1.0 because root_chord is defined separately"
            )
        if abs(float(root_station.sweep_qc_deg)) > 1e-12:
            raise ValueError("the root designer station cannot define sweep_qc_deg")
        if abs(float(root_station.dihedral_deg)) > 1e-12:
            raise ValueError("the root designer station cannot define dihedral_deg")

        for law in self.scalar_laws:
            law.validate()

    def semispan(self) -> float:
        self.validate()
        if self.symmetric:
            return 0.5 * float(self.span)
        return float(self.span)

    def section_interpolation(self) -> InterpolationSpec:
        transition = self.spine_transition if self.section_transition is None else self.section_transition
        return transition.to_interpolation_spec()

    def spine_interpolation(self) -> InterpolationSpec:
        return self.spine_transition.to_interpolation_spec()

    def _build_wing_stations(self) -> Tuple[WingStationSpec, ...]:
        self.validate()

        semispan = self.semispan()
        wing_stations = []
        prev_x_le = float(self.root_le_x)
        prev_y = float(self.root_vertical_y)
        prev_chord = float(self.root_chord)
        prev_eta = 0.0

        for index, station in enumerate(self.stations):
            chord = float(self.root_chord) * float(station.chord_ratio)
            eta = float(station.eta)

            if index == 0:
                wing_stations.append(
                    WingStationSpec(
                        station_id=station.station_id,
                        eta=eta,
                        profile_id=station.profile_id,
                        chord=chord,
                        twist_deg=float(station.twist_deg),
                        pitch_deg=float(station.pitch_deg),
                        roll_deg=float(station.roll_deg),
                        thickness_scale=float(station.thickness_scale),
                    )
                )
            else:
                span_delta = semispan * (eta - prev_eta)
                prev_qc_x = prev_x_le + 0.25 * prev_chord
                next_qc_x = prev_qc_x + span_delta * tan(radians(float(station.sweep_qc_deg)))
                x_le = next_qc_x - 0.25 * chord
                vertical_y = prev_y + span_delta * tan(radians(float(station.dihedral_deg)))

                wing_stations.append(
                    WingStationSpec(
                        station_id=station.station_id,
                        eta=eta,
                        profile_id=station.profile_id,
                        chord=chord,
                        twist_deg=float(station.twist_deg),
                        pitch_deg=float(station.pitch_deg),
                        roll_deg=float(station.roll_deg),
                        thickness_scale=float(station.thickness_scale),
                        x_le=x_le,
                        vertical_y=vertical_y,
                    )
                )
                prev_x_le = x_le
                prev_y = vertical_y
                prev_chord = chord
                prev_eta = eta
                continue

            prev_chord = chord
            prev_eta = eta

        return tuple(wing_stations)

    def to_wing_spec(self) -> WingSpec:
        wing = WingSpec(
            wing_id=self.wing_id,
            semispan=self.semispan(),
            stations=self._build_wing_stations(),
            root_x=float(self.root_le_x),
            root_y=float(self.root_vertical_y),
            root_z=float(self.root_z),
            section_interpolation=self.section_interpolation(),
            spine_interpolation=self.spine_interpolation(),
            scalar_laws=self.scalar_laws,
            mirrored=False,
            metadata={
                **dict(self.metadata),
                "designer_api": "wing_designer",
                "designer_span": str(float(self.span)),
                "designer_symmetric": str(bool(self.symmetric)).lower(),
            },
        )
        wing.validate()
        return wing

    def to_component_spec(self) -> LoftedComponentSpec:
        return self.to_wing_spec().to_component_spec()

    def to_full_span_component_spec(self) -> LoftedComponentSpec:
        return self.to_wing_spec().to_full_span_component_spec()
