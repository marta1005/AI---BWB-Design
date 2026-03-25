"""CTA-specific BWB parametrization package."""

from .design_space import (
    CTA_ACTIVE_VARIABLES,
    CTA_FIXED_PARAMETERS,
    build_cta_design_space,
    cta_fixed_parameters,
    cta_parameter_metadata,
    sample_cta_designs,
)
from .reference import (
    FIXED_TE_EXACT_SEGMENTS,
    SWEEP_NAME_TO_VARIABLE,
    VARIABLE_TO_SWEEP_NAME,
    apply_cta_fixed_parameters,
    build_reference_design,
    to_cta_model_config,
)
from .gemseo_space import (
    CTAGemseoDesignSpaceAdapter,
    CTAGemseoVariableSpec,
    CTASampleGeometryEvaluation,
    available_cta_gemseo_doe_algorithms,
    build_cta_gemseo_design_space,
    build_cta_gemseo_design_space_definition,
    evaluate_cta_gemseo_sample_geometry,
    sample_cta_gemseo_doe,
)

__all__ = [
    "CTA_ACTIVE_VARIABLES",
    "CTA_FIXED_PARAMETERS",
    "FIXED_TE_EXACT_SEGMENTS",
    "SWEEP_NAME_TO_VARIABLE",
    "VARIABLE_TO_SWEEP_NAME",
    "apply_cta_fixed_parameters",
    "build_cta_design_space",
    "build_reference_design",
    "CTAGemseoDesignSpaceAdapter",
    "CTAGemseoVariableSpec",
    "CTASampleGeometryEvaluation",
    "available_cta_gemseo_doe_algorithms",
    "build_cta_gemseo_design_space",
    "build_cta_gemseo_design_space_definition",
    "cta_fixed_parameters",
    "cta_parameter_metadata",
    "evaluate_cta_gemseo_sample_geometry",
    "sample_cta_designs",
    "sample_cta_gemseo_doe",
    "to_cta_model_config",
]
