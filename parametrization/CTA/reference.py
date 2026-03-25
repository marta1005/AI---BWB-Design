"""Reference CTA parametrization built on top of the local BWB package."""

from dataclasses import replace
from math import atan2, degrees, radians, tan
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

# CTA planform trailing-edge definition requested by design:
# - C0 -> C1 straight exact segment
# - C1 -> C3 spline blend
# - C3 -> C4 smooth transition
# - C4 -> C5 straight exact segment
FIXED_TE_EXACT_SEGMENTS: Tuple[int, ...] = (0, 4)

CTA_REFERENCE_SPAN_M = 39.5
CTA_REFERENCE_B1_M = 8.041
CTA_REFERENCE_B2_M = 4.468
CTA_REFERENCE_B3_M = 26.991

CTA_REFERENCE_C0_BODY_CHORD_M = 41.203
CTA_REFERENCE_C3_TRANSITION_CHORD_M = 13.927
CTA_REFERENCE_C4_WING_CHORD_M = 7.768
CTA_REFERENCE_C5_WING_TIP_M = 0.8
CTA_REFERENCE_NOSE_HELPER_1_Y_M = 1.90
CTA_REFERENCE_NOSE_HELPER_1_X_M = 1.905
CTA_REFERENCE_C1_Y_M = 5.694
CTA_REFERENCE_C1_LE_HELPER_X_M = 10.787
CTA_REFERENCE_TE_INBOARD_BLEND_DX_M = 6.7
CTA_REFERENCE_TE_INBOARD_BLEND_Y_M = 7.50
CTA_REFERENCE_TE_INBOARD_RADIUS_FACTOR = 1.5

CTA_REFERENCE_S_DEG = 66.87
CTA_REFERENCE_C1_SWEEP_DEG = 64.85
CTA_REFERENCE_S1_DEG = 54.059
CTA_REFERENCE_S2_DEG = 27.71

CTA_REFERENCE_PROFILE_GENERATION_MODE = "enforce_targets"

# CTA reference airfoil family re-tuned toward a more subsonic-laminar look:
# - lower camber than the original family
# - aft-shifted thickness peak
# - less razor-thin trailing edge
CTA_REFERENCE_C0_UPPER_CST: Tuple[float, ...] = (0.233948, 0.313004, 0.124486, 0.451290, 0.075379, 0.359262)
CTA_REFERENCE_C0_LOWER_CST: Tuple[float, ...] = (0.221796, 0.164675, 0.060559, 0.310561, -0.039576, 0.224411)
CTA_REFERENCE_C3_UPPER_CST: Tuple[float, ...] = (0.205464, 0.283149, 0.112921, 0.403674, 0.073142, 0.322782)
CTA_REFERENCE_C3_LOWER_CST: Tuple[float, ...] = (0.193312, 0.134820, 0.048994, 0.262945, -0.041814, 0.187931)
CTA_REFERENCE_C4_UPPER_CST: Tuple[float, ...] = (0.176980, 0.253294, 0.101356, 0.356058, 0.070904, 0.286303)
CTA_REFERENCE_C4_LOWER_CST: Tuple[float, ...] = (0.164828, 0.104965, 0.037429, 0.215330, -0.044051, 0.151452)
CTA_REFERENCE_C5_UPPER_CST: Tuple[float, ...] = (0.176980, 0.253294, 0.101356, 0.356058, 0.070904, 0.286303)
CTA_REFERENCE_C5_LOWER_CST: Tuple[float, ...] = (0.164828, 0.104965, 0.037429, 0.215330, -0.044051, 0.151452)

CTA_REFERENCE_C0_TC_MAX = 0.18
CTA_REFERENCE_C3_TC_MAX = 0.16
CTA_REFERENCE_C4_TC_MAX = 0.14
CTA_REFERENCE_C5_TC_MAX = 0.115

CTA_REFERENCE_C0_X_TMAX = 0.30
CTA_REFERENCE_C3_X_TMAX = 0.30
CTA_REFERENCE_C4_X_TMAX = 0.31
CTA_REFERENCE_C5_X_TMAX = 0.31

CTA_REFERENCE_C0_TE_THICKNESS = 0.0032
CTA_REFERENCE_C3_TE_THICKNESS = 0.0030
CTA_REFERENCE_C4_TE_THICKNESS = 0.0028
CTA_REFERENCE_C5_TE_THICKNESS = 0.0025

CTA_REFERENCE_C0_CAMBER_DELTA = 0.0
CTA_REFERENCE_C3_CAMBER_DELTA = 0.0
CTA_REFERENCE_C4_CAMBER_DELTA = 0.0
CTA_REFERENCE_C5_CAMBER_DELTA = 0.0


def _cta_reference_c3_transition_chord_m() -> float:
    """Return the public CTA C3 chord from the fixed Airbus-like reference."""
    return float(CTA_REFERENCE_C3_TRANSITION_CHORD_M)


def _cta_reference_le_helper_1_x_m() -> float:
    return float(CTA_REFERENCE_NOSE_HELPER_1_X_M)


def _cta_reference_le_helper_2_x_m() -> float:
    return float(CTA_REFERENCE_C1_LE_HELPER_X_M)


def _cta_reference_le_c3_x_m() -> float:
    dy = float(CTA_REFERENCE_B1_M - CTA_REFERENCE_C1_Y_M)
    return float(_cta_reference_le_helper_2_x_m() + tan(radians(CTA_REFERENCE_C1_SWEEP_DEG)) * dy)


def _cta_reference_internal_secant_s_deg() -> float:
    return float(degrees(atan2(_cta_reference_le_c3_x_m(), CTA_REFERENCE_B1_M)))


def build_reference_design() -> SectionedBWBDesignVariables:
    """Return the CTA reference design in the local parametrization core."""
    base = SectionedBWBDesignVariables.reference_seed()
    return replace(
        base,
        span=CTA_REFERENCE_SPAN_M,
        b1_span_ratio=CTA_REFERENCE_B1_M / CTA_REFERENCE_SPAN_M,
        b2_span_ratio=CTA_REFERENCE_B2_M / CTA_REFERENCE_SPAN_M,
        b3_span_ratio=CTA_REFERENCE_B3_M / CTA_REFERENCE_SPAN_M,
        c1_root_chord=CTA_REFERENCE_C0_BODY_CHORD_M,
        c2_c1_ratio=_cta_reference_c3_transition_chord_m() / CTA_REFERENCE_C0_BODY_CHORD_M,
        c3_c1_ratio=CTA_REFERENCE_C4_WING_CHORD_M / CTA_REFERENCE_C0_BODY_CHORD_M,
        c4_c1_ratio=CTA_REFERENCE_C5_WING_TIP_M / CTA_REFERENCE_C0_BODY_CHORD_M,
        c1_tc_max=CTA_REFERENCE_C0_TC_MAX,
        c2_tc_max=CTA_REFERENCE_C3_TC_MAX,
        c3_tc_max=CTA_REFERENCE_C4_TC_MAX,
        c4_tc_max=CTA_REFERENCE_C5_TC_MAX,
        c1_x_tmax=CTA_REFERENCE_C0_X_TMAX,
        c2_x_tmax=CTA_REFERENCE_C3_X_TMAX,
        c3_x_tmax=CTA_REFERENCE_C4_X_TMAX,
        c4_x_tmax=CTA_REFERENCE_C5_X_TMAX,
        c1_te_thickness=CTA_REFERENCE_C0_TE_THICKNESS,
        c2_te_thickness=CTA_REFERENCE_C3_TE_THICKNESS,
        c3_te_thickness=CTA_REFERENCE_C4_TE_THICKNESS,
        c4_te_thickness=CTA_REFERENCE_C5_TE_THICKNESS,
        camber_c1=CTA_REFERENCE_C0_CAMBER_DELTA,
        camber_c2=CTA_REFERENCE_C3_CAMBER_DELTA,
        camber_c3=CTA_REFERENCE_C4_CAMBER_DELTA,
        camber_c4=CTA_REFERENCE_C5_CAMBER_DELTA,
        c1_upper_cst=CTA_REFERENCE_C0_UPPER_CST,
        c1_lower_cst=CTA_REFERENCE_C0_LOWER_CST,
        c2_upper_cst=CTA_REFERENCE_C3_UPPER_CST,
        c2_lower_cst=CTA_REFERENCE_C3_LOWER_CST,
        c3_upper_cst=CTA_REFERENCE_C4_UPPER_CST,
        c3_lower_cst=CTA_REFERENCE_C4_LOWER_CST,
        c4_upper_cst=CTA_REFERENCE_C5_UPPER_CST,
        c4_lower_cst=CTA_REFERENCE_C5_LOWER_CST,
        twist_c2_deg=base.twist_c1_deg,
        s1_deg=_cta_reference_internal_secant_s_deg(),
        s2_deg=CTA_REFERENCE_S1_DEG,
        s3_deg=CTA_REFERENCE_S2_DEG,
    )


def _cta_fixed_b1_length(
    reference_design: Optional[SectionedBWBDesignVariables] = None,
) -> float:
    reference = build_reference_design() if reference_design is None else reference_design
    return float(reference.span * reference.b1_span_ratio)


def _resolve_cta_span_partition(
    design: SectionedBWBDesignVariables,
    fixed_b1_length: float,
) -> Tuple[float, float, float]:
    if design.span <= fixed_b1_length:
        raise ValueError(
            "CTA fixed B1 length must remain smaller than the semi-span, "
            f"got B1={fixed_b1_length:.6f} m and span={design.span:.6f} m"
        )

    b1_span_ratio = float(fixed_b1_length / design.span)
    b2_span_ratio = float(design.b2_span_ratio)
    b3_span_ratio = float(1.0 - b1_span_ratio - b2_span_ratio)

    if b2_span_ratio <= 0.0:
        raise ValueError(f"CTA b2_span_ratio must stay positive, got {b2_span_ratio:.6f}")
    if b3_span_ratio <= 0.0:
        raise ValueError(
            "CTA span partition leaves no room for B3 after enforcing the fixed B1 length, "
            f"got span={design.span:.6f} m, B1={fixed_b1_length:.6f} m, "
            f"b2_span_ratio={b2_span_ratio:.6f}, resulting b3_span_ratio={b3_span_ratio:.6f}"
        )
    return b1_span_ratio, b2_span_ratio, b3_span_ratio


def cta_fixed_values(
    reference_design: Optional[SectionedBWBDesignVariables] = None,
) -> Dict[str, object]:
    """Fixed CTA values used by the AI design-space definition."""
    reference = build_reference_design() if reference_design is None else reference_design
    b1_fixed_m = _cta_fixed_b1_length(reference)
    c0_reference_body_chord_m = float(reference.c1_root_chord)
    c3_reference_m = float(reference.c2_c1_ratio * reference.c1_root_chord)
    c4_wing_chord_m = float(reference.c3_c1_ratio * reference.c1_root_chord)
    c5_wing_tip_m = float(reference.c4_c1_ratio * reference.c1_root_chord)
    wing_span_m = float(reference.span - b1_fixed_m)
    b2_wing_span_ratio = float((reference.span * reference.b2_span_ratio) / max(wing_span_m, 1e-12))
    return {
        "s_deg": float(CTA_REFERENCE_S_DEG),
        "c1_sweep_deg": float(CTA_REFERENCE_C1_SWEEP_DEG),
        "s_internal_deg": float(reference.s1_deg),
        "b1_fixed_m": b1_fixed_m,
        "c0_body_chord_m": c0_reference_body_chord_m,
        "c0_reference_body_chord_m": c0_reference_body_chord_m,
        "c3_reference_m": c3_reference_m,
        "c4_wing_chord_m": c4_wing_chord_m,
        "c5_wing_tip_m": c5_wing_tip_m,
        "wing_span_m": wing_span_m,
        "b2_wing_span_ratio": b2_wing_span_ratio,
        "te_exact_segments": FIXED_TE_EXACT_SEGMENTS,
    }


def apply_cta_fixed_parameters(
    design: SectionedBWBDesignVariables,
    reference_design: Optional[SectionedBWBDesignVariables] = None,
) -> SectionedBWBDesignVariables:
    """Force CTA fixed parameters on a design vector."""
    fixed = cta_fixed_values(reference_design=reference_design)
    current_c0_body_chord_m = float(design.c1_root_chord)
    if current_c0_body_chord_m <= 0.0:
        raise ValueError(f"CTA C0/body chord must stay positive, got {current_c0_body_chord_m:.6f} m")
    current_c3_transition_chord_m = float(design.c2_c1_ratio * current_c0_body_chord_m)
    if current_c3_transition_chord_m <= 0.0:
        raise ValueError(
            f"CTA C3/transition-wing chord must stay positive, got {current_c3_transition_chord_m:.6f} m"
        )
    current_c5_wing_tip_m = float(design.c4_c1_ratio * current_c0_body_chord_m)

    b1_span_ratio, b2_span_ratio, b3_span_ratio = _resolve_cta_span_partition(
        design=design,
        fixed_b1_length=float(fixed["b1_fixed_m"]),
    )
    current_b2_m = float(design.span * b2_span_ratio)
    current_c4_wing_chord_m = float(design.c3_c1_ratio * current_c0_body_chord_m)
    if current_c4_wing_chord_m <= 0.0:
        raise ValueError(
            f"CTA C4/outer-wing chord must stay positive, got {current_c4_wing_chord_m:.6f} m"
        )
    return replace(
        design,
        b1_span_ratio=b1_span_ratio,
        b2_span_ratio=b2_span_ratio,
        b3_span_ratio=b3_span_ratio,
        s1_deg=float(fixed["s_internal_deg"]),
        c2_c1_ratio=float(current_c3_transition_chord_m / current_c0_body_chord_m),
        c3_c1_ratio=float(current_c4_wing_chord_m / current_c0_body_chord_m),
        c4_c1_ratio=float(current_c5_wing_tip_m / current_c0_body_chord_m),
        twist_c2_deg=float(design.twist_c1_deg),
    )


def to_cta_model_config(
    design: SectionedBWBDesignVariables,
    reference_design: Optional[SectionedBWBDesignVariables] = None,
) -> SectionedBWBModelConfig:
    """
    Build a model config enforcing CTA fixed choices:
    fixed sweep S, fixed CTA chord family, fixed B1 in meters,
    and fixed straight TE segments.
    """
    fixed_design = apply_cta_fixed_parameters(design, reference_design=reference_design)
    config = fixed_design.to_model_config(profile_generation_mode=CTA_REFERENCE_PROFILE_GENERATION_MODE)
    config.sections.shared_leading_edge = True
    config.planform.te_exact_segments = FIXED_TE_EXACT_SEGMENTS
    config.planform.te_c1_span_fraction = CTA_REFERENCE_C1_Y_M / CTA_REFERENCE_B1_M
    config.planform.te_inboard_blend_fraction = CTA_REFERENCE_TE_INBOARD_BLEND_Y_M / CTA_REFERENCE_B1_M
    config.planform.te_inboard_blend_dx = CTA_REFERENCE_TE_INBOARD_BLEND_DX_M
    config.planform.te_inboard_radius_factor = CTA_REFERENCE_TE_INBOARD_RADIUS_FACTOR
    config.planform.te_outer_blend_fraction = 0.0
    # Keep the public CTA geometry but make the C1->C3 TE transition read as
    # a genuine spline instead of a mostly linear break between hidden helpers.
    config.planform.te_blend_fraction = 0.44
    config.planform.te_min_linear_core_fraction = 0.08
    config.planform.te_spline_bridge = (1, 3)
    # Keep the LE as a segmented curve with local smooth blends at the helper
    # and section junctions, while preserving straight cores between them.
    config.planform.blend_fraction = 0.24
    config.planform.min_linear_core_fraction = 0.58
    config.planform.body_le_fixed_points = (
        (_cta_reference_le_helper_1_x_m(), CTA_REFERENCE_NOSE_HELPER_1_Y_M),
        (_cta_reference_le_helper_2_x_m(), CTA_REFERENCE_C1_Y_M),
    )
    config.planform.symmetry_blend_y = CTA_REFERENCE_NOSE_HELPER_1_Y_M
    # pyGeo export is more robust for this CTA when the lifting surface is
    # reconstructed from the anchor sections instead of refitting every
    # already-interpolated station airfoil independently.
    config.sampling.airfoil_distribution_mode = "anchors"
    # Keep sections 0-3 with no twist variation: the first actual twist change
    # is allowed only after C3.
    config.spanwise.twist_deg.interpolation = "pchip"
    config.spanwise.twist_deg.values = (
        float(fixed_design.twist_c1_deg),
        float(fixed_design.twist_c1_deg),
        float(fixed_design.twist_c3_deg),
        float(fixed_design.twist_c4_deg),
    )
    return config
