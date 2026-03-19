"""CTA design-space definition for AI sampling and optimization."""

from dataclasses import asdict
from typing import Dict, List, Optional, Tuple

from parametrization.bwb.design_space import DesignSpace, parameter_info
from parametrization.bwb.design_variables import SectionedBWBDesignVariables

from .reference import (
    SWEEP_NAME_TO_VARIABLE,
    VARIABLE_TO_SWEEP_NAME,
    apply_cta_fixed_parameters,
    build_reference_design,
    cta_fixed_values,
)

CTA_FIXED_PARAMETERS: Tuple[str, ...] = (
    "s1_deg",
    "c1_root_chord",
    "c2_c1_ratio",
)


def _cst_names(section_index: int) -> Tuple[str, ...]:
    return tuple(
        [f"c{section_index}_upper_cst_{idx}" for idx in range(6)]
        + [f"c{section_index}_lower_cst_{idx}" for idx in range(6)]
    )


CTA_ACTIVE_VARIABLES: Tuple[str, ...] = (
    "span",
    "b1_span_ratio",
    "b2_span_ratio",
    "b3_span_ratio",
    "s2_deg",  # CTA label S1
    "s3_deg",  # CTA label S2
    "c3_c1_ratio",
    "c4_c1_ratio",
    "twist_c1_deg",
    "twist_c2_deg",
    "twist_c3_deg",
    "twist_c4_deg",
    *_cst_names(1),
    *_cst_names(2),
    *_cst_names(3),
    *_cst_names(4),
)


def build_cta_design_space(
    reference_design: Optional[SectionedBWBDesignVariables] = None,
) -> DesignSpace:
    """
    Build CTA AI design space with requested fixed/variable split:
    - fixed: S, C1, C2
    - variable: semi-span, B ratios, following sweeps, twists, C3/C4/C5 chord family, CST.
    """
    reference = build_reference_design() if reference_design is None else reference_design
    reference = apply_cta_fixed_parameters(reference, reference_design=reference)
    bounds = SectionedBWBDesignVariables.default_bounds()
    return DesignSpace(
        preset_name="cta_reference_ai_core",
        active_groups=("cta_core",),
        active_variables=CTA_ACTIVE_VARIABLES,
        reference_design=reference,
        bounds=bounds,
    )


def cta_fixed_parameters(
    reference_design: Optional[SectionedBWBDesignVariables] = None,
) -> Dict[str, object]:
    """Fixed CTA values plus naming aliases for presentation."""
    reference = build_reference_design() if reference_design is None else reference_design
    fixed = cta_fixed_values(reference_design=reference)
    fixed.update(
        {
            "sweep_naming": dict(SWEEP_NAME_TO_VARIABLE),
            "legacy_reference": {
                "S(legacy S1)": fixed["s_deg"],
                "S1(legacy S2)": float(reference.s2_deg),
                "S2(legacy S3)": float(reference.s3_deg),
            },
            "notes": (
                "C5 chord in CTA plot notation is tied to the current tip-chord control "
                "(core c4_c1_ratio)."
            ),
        }
    )
    return fixed


def cta_parameter_metadata(
    reference_design: Optional[SectionedBWBDesignVariables] = None,
) -> List[Dict[str, object]]:
    """Return active CTA variables with metadata and renamed sweep labels."""
    design_space = build_cta_design_space(reference_design=reference_design)
    reference_flat = design_space.reference_flat()
    rows: List[Dict[str, object]] = []
    for name in design_space.active_variables:
        info = dict(parameter_info(name))
        if name in VARIABLE_TO_SWEEP_NAME:
            cta_name = VARIABLE_TO_SWEEP_NAME[name]
            info["display_name"] = f"Sweep {cta_name}"
            info["symbol"] = cta_name
            info["description"] = f"CTA sweep {cta_name} (mapped from {name})."
        if name == "c3_c1_ratio":
            info["description"] = "Chord ratio for CTA chord family C3/C1."
        if name == "c4_c1_ratio":
            info["description"] = "Chord ratio for CTA outboard chord family C4/C5 relative to C1."
        lower, upper = design_space.bounds[name]
        rows.append(
            {
                "parameter": name,
                "display_name": info["display_name"],
                "symbol": info["symbol"],
                "units": info["units"],
                "normalization": info["normalization"],
                "description": info["description"],
                "reference": float(reference_flat[name]),
                "lower_bound": float(lower),
                "upper_bound": float(upper),
            }
        )
    return rows


def sample_cta_designs(
    count: int,
    seed: int = 7,
    variation_scale: float = 0.25,
    reference_design: Optional[SectionedBWBDesignVariables] = None,
) -> List[SectionedBWBDesignVariables]:
    """
    Sample CTA designs from the active AI design-space variables only.
    Fixed parameters remain locked by construction.
    """
    design_space = build_cta_design_space(reference_design=reference_design)
    return design_space.sample_designs(count=count, seed=seed, variation_scale=variation_scale)


def cta_design_space_summary(
    reference_design: Optional[SectionedBWBDesignVariables] = None,
) -> Dict[str, object]:
    """Compact serializable summary used by docs/scripts."""
    ds = build_cta_design_space(reference_design=reference_design)
    return {
        "preset_name": ds.preset_name,
        "active_variable_count": len(ds.active_variables),
        "active_variables": list(ds.active_variables),
        "fixed_parameters": cta_fixed_parameters(reference_design=ds.reference_design),
        "reference_design": asdict(ds.reference_design),
    }
