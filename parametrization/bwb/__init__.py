"""Generic sectioned BWB parametrization tools."""

from .builder import PreparedGeometry, build_loft_definition, build_span_stations, build_surface, export_iges, prepare_geometry
from .design_space import DesignSpace, flatten_design, parameter_groups, parameter_info
from .design_variables import SectionedBWBDesignVariables
from .planform import SectionedPlanform, build_sectioned_bwb_planform
from .specs import SectionedBWBModelConfig
from .topology import SectionedBWBTopologySpec

__all__ = [
    "DesignSpace",
    "PreparedGeometry",
    "SectionedBWBDesignVariables",
    "SectionedBWBModelConfig",
    "SectionedBWBTopologySpec",
    "SectionedPlanform",
    "build_loft_definition",
    "build_sectioned_bwb_planform",
    "build_span_stations",
    "build_surface",
    "export_iges",
    "flatten_design",
    "parameter_groups",
    "parameter_info",
    "prepare_geometry",
]
