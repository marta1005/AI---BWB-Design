from dataclasses import dataclass, field
from typing import Dict, Tuple

from .fuselage import FuselageBuildOptions, FuselageSpec, PreparedFuselage, prepare_fuselage
from .lifting_surface import (
    LiftingSurfaceBuildOptions,
    PreparedLiftingSurface,
    prepare_lifting_surface,
)
from .components import LoftedComponentSpec
from .profiles import ProfileCatalog
from .wing import WingSpec
from .vertical_tail import VerticalTailSpec


@dataclass(frozen=True)
class AircraftGeometrySpec:
    aircraft_id: str
    profiles: ProfileCatalog
    components: Tuple[LoftedComponentSpec, ...]

    def validate(self) -> None:
        if not self.aircraft_id:
            raise ValueError("aircraft_id must be non-empty")
        self.profiles.validate()
        if not self.components:
            raise ValueError("an aircraft geometry must contain at least one component")

        profile_catalog = self.profiles.as_dict()
        component_ids = set()
        for component in self.components:
            component.validate()
            if component.component_id in component_ids:
                raise ValueError(
                    f"duplicate component_id detected in aircraft {self.aircraft_id!r}: "
                    f"{component.component_id!r}"
                )
            component_ids.add(component.component_id)

            for profile_id in component.profile_ids():
                if profile_id not in profile_catalog:
                    raise ValueError(
                        f"component {component.component_id!r} references unknown profile_id "
                        f"{profile_id!r}"
                    )

    def component_ids(self) -> Tuple[str, ...]:
        self.validate()
        return tuple(component.component_id for component in self.components)


@dataclass(frozen=True)
class AircraftWingEntry:
    wing: WingSpec
    options: LiftingSurfaceBuildOptions = LiftingSurfaceBuildOptions()
    mirror: bool = True

    def validate(self) -> None:
        self.wing.validate()
        self.options.validate()


@dataclass(frozen=True)
class AircraftFuselageEntry:
    fuselage: FuselageSpec
    options: FuselageBuildOptions = FuselageBuildOptions()

    def validate(self) -> None:
        self.fuselage.validate()
        self.options.validate()


@dataclass(frozen=True)
class AircraftVerticalTailEntry:
    tail: VerticalTailSpec
    options: LiftingSurfaceBuildOptions = LiftingSurfaceBuildOptions()

    def validate(self) -> None:
        self.tail.validate()
        self.options.validate()


@dataclass(frozen=True)
class AircraftAssemblySpec:
    aircraft_id: str
    profiles: ProfileCatalog
    wings: Tuple[AircraftWingEntry, ...] = ()
    vertical_tails: Tuple[AircraftVerticalTailEntry, ...] = ()
    fuselages: Tuple[AircraftFuselageEntry, ...] = ()
    metadata: Dict[str, str] = field(default_factory=dict)

    def validate(self) -> None:
        if not self.aircraft_id:
            raise ValueError("aircraft_id must be non-empty")
        self.profiles.validate()
        if not self.wings and not self.vertical_tails and not self.fuselages:
            raise ValueError("an aircraft assembly must contain at least one lifting surface or one fuselage")

        profile_catalog = self.profiles.as_dict()
        component_ids = set()
        for entry in self.wings:
            entry.validate()
            wing_id = entry.wing.wing_id
            if wing_id in component_ids:
                raise ValueError(f"duplicate aircraft component id detected: {wing_id!r}")
            component_ids.add(wing_id)
            for profile_id in (station.profile_id for station in entry.wing.stations):
                if profile_id not in profile_catalog:
                    raise ValueError(f"wing {wing_id!r} references unknown profile_id {profile_id!r}")

        for entry in self.vertical_tails:
            entry.validate()
            tail_id = entry.tail.tail_id
            if tail_id in component_ids:
                raise ValueError(f"duplicate aircraft component id detected: {tail_id!r}")
            component_ids.add(tail_id)
            for profile_id in (station.profile_id for station in entry.tail.stations):
                if profile_id not in profile_catalog:
                    raise ValueError(f"vertical tail {tail_id!r} references unknown profile_id {profile_id!r}")

        for entry in self.fuselages:
            entry.validate()
            fuselage_id = entry.fuselage.fuselage_id
            if fuselage_id in component_ids:
                raise ValueError(f"duplicate aircraft component id detected: {fuselage_id!r}")
            component_ids.add(fuselage_id)


@dataclass(frozen=True)
class PreparedAircraftWing:
    wing_id: str
    prepared: PreparedLiftingSurface


@dataclass(frozen=True)
class PreparedAircraftVerticalTail:
    tail_id: str
    prepared: PreparedLiftingSurface


@dataclass(frozen=True)
class PreparedAircraftFuselage:
    fuselage_id: str
    prepared: PreparedFuselage


@dataclass(frozen=True)
class PreparedAircraftGeometry:
    aircraft_id: str
    wings: Tuple[PreparedAircraftWing, ...]
    vertical_tails: Tuple[PreparedAircraftVerticalTail, ...]
    fuselages: Tuple[PreparedAircraftFuselage, ...]


def _full_span_build_options(options: LiftingSurfaceBuildOptions) -> LiftingSurfaceBuildOptions:
    options.validate()
    if options.station_etas:
        half_etas = tuple(sorted({float(value) for value in options.station_etas}))
        positive = tuple(value for value in half_etas if value > 1e-12)
        full_station_etas = tuple(
            sorted(
                {
                    0.5,
                    *[0.5 * (1.0 - value) for value in positive],
                    *[0.5 * (1.0 + value) for value in positive],
                }
            )
        )
    else:
        full_station_etas = ()

    return LiftingSurfaceBuildOptions(
        station_count=2 * int(options.station_count) - 1,
        station_etas=full_station_etas,
        include_anchor_sections=options.include_anchor_sections,
        airfoil_sample_count=options.airfoil_sample_count,
        fit_n_ctl=options.fit_n_ctl,
        k_span=options.k_span,
        tip_style=options.tip_style,
        blunt_te=options.blunt_te,
        rounded_te=options.rounded_te,
        rebuild_dependencies=options.rebuild_dependencies,
    )


def prepare_aircraft_geometry(aircraft: AircraftAssemblySpec) -> PreparedAircraftGeometry:
    aircraft.validate()

    prepared_wings = []
    for entry in aircraft.wings:
        if entry.mirror:
            component = entry.wing.to_full_span_component_spec()
            options = _full_span_build_options(entry.options)
        else:
            component = entry.wing.to_component_spec()
            options = entry.options
        prepared_wings.append(
            PreparedAircraftWing(
                wing_id=entry.wing.wing_id,
                prepared=prepare_lifting_surface(component, aircraft.profiles, options=options),
            )
        )

    prepared_vertical_tails = []
    for entry in aircraft.vertical_tails:
        component = entry.tail.to_component_spec()
        prepared_vertical_tails.append(
            PreparedAircraftVerticalTail(
                tail_id=entry.tail.tail_id,
                prepared=prepare_lifting_surface(component, aircraft.profiles, options=entry.options),
            )
        )

    prepared_fuselages = []
    for entry in aircraft.fuselages:
        prepared_fuselages.append(
            PreparedAircraftFuselage(
                fuselage_id=entry.fuselage.fuselage_id,
                prepared=prepare_fuselage(entry.fuselage, options=entry.options),
            )
        )

    return PreparedAircraftGeometry(
        aircraft_id=aircraft.aircraft_id,
        wings=tuple(prepared_wings),
        vertical_tails=tuple(prepared_vertical_tails),
        fuselages=tuple(prepared_fuselages),
    )
