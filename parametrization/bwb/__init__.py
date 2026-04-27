"""Generic sectioned BWB parametrization tools."""

from .builder import PreparedGeometry, build_loft_definition, build_span_stations, build_surface, export_iges, prepare_geometry
from .design_space import DesignSpace, flatten_design, parameter_groups, parameter_info
from .design_variables import SectionedBWBDesignVariables
from .internal_volume_constraints import (
    CadReferenceFrame,
    GeometryEnvelopeEvaluator,
    IndicatorSurfaceResult,
    IndicatorSurfaceSpec,
    InternalVolumeConstraintResult,
    InternalVolumeConstraintSet,
    evaluate_indicator_surface,
    evaluate_internal_volume_constraints,
    evaluate_internal_volume_constraints_from_config,
    load_internal_volume_constraint_set,
)
from .planform import SectionedPlanform, build_sectioned_bwb_planform
from .specs import SectionedBWBModelConfig
from .topology import SectionedBWBTopologySpec

__all__ = [
    "CadReferenceFrame",
    "DesignSpace",
    "GeometryEnvelopeEvaluator",
    "IndicatorSurfaceResult",
    "IndicatorSurfaceSpec",
    "InternalVolumeConstraintResult",
    "InternalVolumeConstraintSet",
    "PreparedGeometry",
    "SectionedBWBDesignVariables",
    "SectionedBWBModelConfig",
    "SectionedBWBTopologySpec",
    "SectionedPlanform",
    "build_loft_definition",
    "build_sectioned_bwb_planform",
    "build_span_stations",
    "build_surface",
    "evaluate_indicator_surface",
    "evaluate_internal_volume_constraints",
    "evaluate_internal_volume_constraints_from_config",
    "export_iges",
    "flatten_design",
    "load_internal_volume_constraint_set",
    "parameter_groups",
    "parameter_info",
    "prepare_geometry",
]
