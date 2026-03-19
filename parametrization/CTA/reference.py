"""Reference CTA parametrization built on top of the local BWB package."""

from dataclasses import replace
from typing import Dict, Optional, Tuple

from parametrization.bwb.design_variables import SectionedBWBDesignVariables
from parametrization.bwb.specs import SectionedBWBModelConfig

# Requested sweep renaming for the CTA package:
# - legacy S1 -> S (fixed)
# - legacy S2 -> S1 (variable)
# - legacy S3 -> S2 (variable)
SWEEP_NAME_TO_VARIABLE: Dict[str, str] = {
    "S": "s1_deg",
    "S1": "s2_deg",
    "S2": "s3_deg",
}
VARIABLE_TO_SWEEP_NAME: Dict[str, str] = {value: key for key, value in SWEEP_NAME_TO_VARIABLE.items()}

# CTA planform straight TE segments requested by design:
# - C0 -> C1
# - C3 -> C4
FIXED_TE_EXACT_SEGMENTS: Tuple[int, int] = (0, 3)


def build_reference_design() -> SectionedBWBDesignVariables:
    """Return the CTA reference design in the local parametrization core."""
    return SectionedBWBDesignVariables.reference_seed()


def cta_fixed_values(
    reference_design: Optional[SectionedBWBDesignVariables] = None,
) -> Dict[str, object]:
    """Fixed CTA values used by the AI design-space definition."""
    reference = build_reference_design() if reference_design is None else reference_design
    return {
        "s_deg": float(reference.s1_deg),
        "c1_root_chord": float(reference.c1_root_chord),
        "c2_c1_ratio": float(reference.c2_c1_ratio),
        "te_exact_segments": FIXED_TE_EXACT_SEGMENTS,
    }


def apply_cta_fixed_parameters(
    design: SectionedBWBDesignVariables,
    reference_design: Optional[SectionedBWBDesignVariables] = None,
) -> SectionedBWBDesignVariables:
    """Force CTA fixed parameters on a design vector."""
    fixed = cta_fixed_values(reference_design=reference_design)
    return replace(
        design,
        s1_deg=float(fixed["s_deg"]),
        c1_root_chord=float(fixed["c1_root_chord"]),
        c2_c1_ratio=float(fixed["c2_c1_ratio"]),
    )


def to_cta_model_config(
    design: SectionedBWBDesignVariables,
    reference_design: Optional[SectionedBWBDesignVariables] = None,
) -> SectionedBWBModelConfig:
    """
    Build a model config enforcing CTA fixed choices:
    fixed S, fixed C1/C2 and fixed straight TE segments.
    """
    fixed_design = apply_cta_fixed_parameters(design, reference_design=reference_design)
    config = fixed_design.to_model_config()
    config.planform.te_exact_segments = FIXED_TE_EXACT_SEGMENTS
    return config
