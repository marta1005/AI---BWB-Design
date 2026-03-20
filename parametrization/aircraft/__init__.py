"""Generic aircraft parametrization package."""

from .aircraft import (
    AircraftAssemblySpec,
    AircraftFuselageEntry,
    AircraftGeometrySpec,
    AircraftVerticalTailEntry,
    AircraftWingEntry,
    PreparedAircraftGeometry,
    PreparedAircraftFuselage,
    PreparedAircraftVerticalTail,
    PreparedAircraftWing,
    prepare_aircraft_geometry,
)
from .components import ComponentKind, LoftedComponentSpec
from .fuselage import (
    FuselageBuildOptions,
    FuselageExportResult,
    FuselageSectionSpec,
    FuselageSpec,
    PreparedFuselage,
    PreparedFuselageStation,
    build_fuselage_pygeo,
    build_fuselage_surface,
    export_fuselage_iges,
    prepare_fuselage,
    write_fuselage_sections,
)
from .lifting_surface import (
    LiftingSurfaceBuildOptions,
    LiftingSurfaceExportResult,
    PreparedLiftingSurface,
    PreparedLiftingSurfaceStation,
    SUPPORTED_LIFTING_SURFACE_LAWS,
    build_lifting_surface_pygeo,
    export_lifting_surface_iges,
    prepare_lifting_surface,
    write_lifting_surface_airfoils,
)
from .laws import ContinuityOrder, InterpolationMethod, InterpolationSpec, ScalarLawAnchor, ScalarLawSpec
from .profiles import (
    CSTAirfoilProfileSpec,
    ICSTAirfoilProfileSpec,
    ICSTConstraint,
    IntuitiveAirfoilProfileSpec,
    ProfileCatalog,
    ProfileSample,
)
from .sections import SectionAnchorSpec, SectionPlacement
from .wing import WingSpec, WingSpanStationSpec, WingStationSpec
from .wing_designer import WingDesignerSpec, WingStationDesignerSpec, WingTransitionSpec
from .vertical_tail import VerticalTailSpec, VerticalTailStationSpec


def build_lifting_surface_mesh(*args, **kwargs):
    from .plotting import build_lifting_surface_mesh as _build_lifting_surface_mesh

    return _build_lifting_surface_mesh(*args, **kwargs)


def create_lifting_surface_3d_figure(*args, **kwargs):
    from .plotting import create_lifting_surface_3d_figure as _create_lifting_surface_3d_figure

    return _create_lifting_surface_3d_figure(*args, **kwargs)


def create_lifting_surface_overview_figure(*args, **kwargs):
    from .plotting import create_lifting_surface_overview_figure as _create_lifting_surface_overview_figure

    return _create_lifting_surface_overview_figure(*args, **kwargs)


def save_lifting_surface_3d(*args, **kwargs):
    from .plotting import save_lifting_surface_3d as _save_lifting_surface_3d

    return _save_lifting_surface_3d(*args, **kwargs)


def save_lifting_surface_overview(*args, **kwargs):
    from .plotting import save_lifting_surface_overview as _save_lifting_surface_overview

    return _save_lifting_surface_overview(*args, **kwargs)


def build_fuselage_mesh(*args, **kwargs):
    from .plotting import build_fuselage_mesh as _build_fuselage_mesh

    return _build_fuselage_mesh(*args, **kwargs)


def create_fuselage_3d_figure(*args, **kwargs):
    from .plotting import create_fuselage_3d_figure as _create_fuselage_3d_figure

    return _create_fuselage_3d_figure(*args, **kwargs)


def create_fuselage_overview_figure(*args, **kwargs):
    from .plotting import create_fuselage_overview_figure as _create_fuselage_overview_figure

    return _create_fuselage_overview_figure(*args, **kwargs)


def save_fuselage_3d(*args, **kwargs):
    from .plotting import save_fuselage_3d as _save_fuselage_3d

    return _save_fuselage_3d(*args, **kwargs)


def save_fuselage_overview(*args, **kwargs):
    from .plotting import save_fuselage_overview as _save_fuselage_overview

    return _save_fuselage_overview(*args, **kwargs)


def create_aircraft_3d_figure(*args, **kwargs):
    from .plotting import create_aircraft_3d_figure as _create_aircraft_3d_figure

    return _create_aircraft_3d_figure(*args, **kwargs)


def create_aircraft_overview_figure(*args, **kwargs):
    from .plotting import create_aircraft_overview_figure as _create_aircraft_overview_figure

    return _create_aircraft_overview_figure(*args, **kwargs)


def save_aircraft_3d(*args, **kwargs):
    from .plotting import save_aircraft_3d as _save_aircraft_3d

    return _save_aircraft_3d(*args, **kwargs)


def save_aircraft_overview(*args, **kwargs):
    from .plotting import save_aircraft_overview as _save_aircraft_overview

    return _save_aircraft_overview(*args, **kwargs)

__all__ = [
    "AircraftGeometrySpec",
    "AircraftAssemblySpec",
    "AircraftFuselageEntry",
    "AircraftVerticalTailEntry",
    "AircraftWingEntry",
    "CSTAirfoilProfileSpec",
    "ComponentKind",
    "ContinuityOrder",
    "FuselageBuildOptions",
    "FuselageExportResult",
    "FuselageSectionSpec",
    "FuselageSpec",
    "ICSTAirfoilProfileSpec",
    "ICSTConstraint",
    "IntuitiveAirfoilProfileSpec",
    "InterpolationMethod",
    "InterpolationSpec",
    "LiftingSurfaceBuildOptions",
    "LiftingSurfaceExportResult",
    "LoftedComponentSpec",
    "PreparedFuselage",
    "PreparedFuselageStation",
    "PreparedAircraftGeometry",
    "PreparedAircraftFuselage",
    "PreparedAircraftVerticalTail",
    "PreparedAircraftWing",
    "PreparedLiftingSurface",
    "PreparedLiftingSurfaceStation",
    "ProfileCatalog",
    "ProfileSample",
    "ScalarLawAnchor",
    "ScalarLawSpec",
    "SectionAnchorSpec",
    "SectionPlacement",
    "SUPPORTED_LIFTING_SURFACE_LAWS",
    "VerticalTailSpec",
    "VerticalTailStationSpec",
    "WingSpec",
    "WingDesignerSpec",
    "WingSpanStationSpec",
    "WingStationDesignerSpec",
    "WingStationSpec",
    "WingTransitionSpec",
    "build_fuselage_mesh",
    "build_fuselage_pygeo",
    "build_fuselage_surface",
    "build_lifting_surface_mesh",
    "build_lifting_surface_pygeo",
    "create_fuselage_3d_figure",
    "create_fuselage_overview_figure",
    "create_aircraft_3d_figure",
    "create_aircraft_overview_figure",
    "create_lifting_surface_3d_figure",
    "create_lifting_surface_overview_figure",
    "export_fuselage_iges",
    "export_lifting_surface_iges",
    "prepare_aircraft_geometry",
    "prepare_fuselage",
    "prepare_lifting_surface",
    "save_aircraft_3d",
    "save_aircraft_overview",
    "save_fuselage_3d",
    "save_fuselage_overview",
    "save_lifting_surface_3d",
    "save_lifting_surface_overview",
    "write_fuselage_sections",
    "write_lifting_surface_airfoils",
]
