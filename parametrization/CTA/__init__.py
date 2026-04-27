"""CTA-specific BWB parametrization package."""

from .design_space import (
    CTA_ACTIVE_VARIABLES,
    CTA_FIXED_PARAMETERS,
    build_cta_design_space,
    cta_fixed_parameters,
    cta_parameter_metadata,
    sample_cta_designs,
)
from .case import (
    FIXED_TE_EXACT_SEGMENTS,
    SWEEP_NAME_TO_VARIABLE,
    VARIABLE_TO_SWEEP_NAME,
    apply_cta_fixed_parameters,
    build_cta_design,
    to_cta_model_config,
)
from .internal_volume_constraints import (
    CTA_CAD_REFERENCE_FRAME,
    CTA_INTERNAL_VOLUME_CONSTRAINTS_PATH,
    evaluate_cta_internal_volume_constraints,
    load_cta_internal_volume_constraint_set,
)

__all__ = [
    "CTA_ACTIVE_VARIABLES",
    "CTA_CAD_REFERENCE_FRAME",
    "CTA_FIXED_PARAMETERS",
    "CTA_INTERNAL_VOLUME_CONSTRAINTS_PATH",
    "FIXED_TE_EXACT_SEGMENTS",
    "SWEEP_NAME_TO_VARIABLE",
    "VARIABLE_TO_SWEEP_NAME",
    "apply_cta_fixed_parameters",
    "build_cta_design_space",
    "build_cta_design",
    "cta_fixed_parameters",
    "cta_parameter_metadata",
    "evaluate_cta_internal_volume_constraints",
    "load_cta_internal_volume_constraint_set",
    "sample_cta_designs",
    "to_cta_model_config",
]
