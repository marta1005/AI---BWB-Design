"""Reference CTA parametrization built on top of the local BWB package."""

from dataclasses import replace
from math import atan2, degrees, radians, tan
from typing import Dict, Optional, Tuple

from parametrization.bwb.design_variables import SectionedBWBDesignVariables
from parametrization.bwb.specs import (
    AnchoredSpanwiseLaw,
    SectionCSTSpec,
    SectionFamilySpec,
    SectionProfileRelationSpec,
    SectionedBWBModelConfig,
)

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

# CTA reference airfoils fitted from the glider geometry for the stations that
# physically exist in the reference half-wing. The two extra outer-wing glider
# sections at y=22.5 and y=32.5 are intentionally not used.
CTA_REFERENCE_SECTION0_UPPER_CST: Tuple[float, ...] = (
    0.3214813601757442,
    0.13541796839093923,
    0.200324122975498,
    0.1038325369718894,
    0.21630141178783335,
    0.16283323861908938,
)
CTA_REFERENCE_SECTION0_LOWER_CST: Tuple[float, ...] = (
    0.11795970916265147,
    0.2382356713547635,
    0.03847969918611987,
    0.29850032826444445,
    0.29406076166617084,
    0.08245081524365487,
)
CTA_REFERENCE_SECTION1_UPPER_CST: Tuple[float, ...] = (
    0.1694179199240756,
    0.2602461265303534,
    0.036871399798569314,
    0.1432700871743449,
    0.1885179368691094,
    0.143355404059814,
)
CTA_REFERENCE_SECTION1_LOWER_CST: Tuple[float, ...] = (
    0.3229822425397869,
    0.10215981482696812,
    0.2412602101433208,
    0.24903666719690892,
    0.2920391779108863,
    0.15126112540981526,
)
CTA_REFERENCE_SECTION2_UPPER_CST: Tuple[float, ...] = (
    0.28475162829971634,
    0.12091204202425475,
    0.2055577386527395,
    0.13780299593578327,
    0.2226651723027435,
    0.16907106535130556,
)
CTA_REFERENCE_SECTION2_LOWER_CST: Tuple[float, ...] = (
    0.11356435505298283,
    0.2083557898981848,
    0.06108847434064107,
    0.09710808414935262,
    0.33228296593192536,
    0.15060629926826177,
)
CTA_REFERENCE_SECTION3_UPPER_CST: Tuple[float, ...] = (
    0.26599421033322984,
    0.2639179005615584,
    0.2826156801348145,
    0.16904898778928093,
    0.17183454718969435,
    0.06761012915760517,
)
CTA_REFERENCE_SECTION3_LOWER_CST: Tuple[float, ...] = (
    0.1654346986913471,
    0.21682216708613417,
    0.09748225097313096,
    0.2483137913575583,
    0.1208849451948237,
    0.09993190335101831,
)
CTA_REFERENCE_SECTION4_UPPER_CST: Tuple[float, ...] = (
    0.192497208133799,
    0.05467717490105108,
    0.24610723721114705,
    0.14203567893957236,
    0.15639758085383196,
    0.05864040954593154,
)
CTA_REFERENCE_SECTION4_LOWER_CST: Tuple[float, ...] = (
    0.07077840726077685,
    0.13093012705982235,
    0.0749666051444539,
    0.1430579993393733,
    0.13613778938194562,
    0.023319557557844452,
)
CTA_REFERENCE_SECTION5_UPPER_CST: Tuple[float, ...] = (
    0.052687919761369426,
    0.1286357146619804,
    0.08772041301281408,
    0.14811826206192863,
    0.1279116086971417,
    0.02496608975723723,
)
CTA_REFERENCE_SECTION5_LOWER_CST: Tuple[float, ...] = (
    0.19542915257084556,
    0.03906111986960669,
    0.2129148491134933,
    0.11474070108949,
    0.14528966197307466,
    0.02952858971421116,
)

CTA_REFERENCE_C0_UPPER_CST = CTA_REFERENCE_SECTION0_UPPER_CST
CTA_REFERENCE_C0_LOWER_CST = CTA_REFERENCE_SECTION0_LOWER_CST
CTA_REFERENCE_C3_UPPER_CST = CTA_REFERENCE_SECTION3_UPPER_CST
CTA_REFERENCE_C3_LOWER_CST = CTA_REFERENCE_SECTION3_LOWER_CST
CTA_REFERENCE_C4_UPPER_CST = CTA_REFERENCE_SECTION4_UPPER_CST
CTA_REFERENCE_C4_LOWER_CST = CTA_REFERENCE_SECTION4_LOWER_CST
CTA_REFERENCE_C5_UPPER_CST = CTA_REFERENCE_SECTION5_UPPER_CST
CTA_REFERENCE_C5_LOWER_CST = CTA_REFERENCE_SECTION5_LOWER_CST

CTA_REFERENCE_SECTION0_TC_MAX = 0.13610295872398795
CTA_REFERENCE_SECTION1_TC_MAX = 0.1424495843805421
CTA_REFERENCE_SECTION2_TC_MAX = 0.1224459047766308
CTA_REFERENCE_SECTION3_TC_MAX = 0.16408005340327847
CTA_REFERENCE_SECTION4_TC_MAX = 0.10318708802443649
CTA_REFERENCE_SECTION5_TC_MAX = 0.09558427385019233

CTA_REFERENCE_SECTION0_X_TMAX = 0.2776824104075363
CTA_REFERENCE_SECTION1_X_TMAX = 0.23208660251050162
CTA_REFERENCE_SECTION2_X_TMAX = 0.2776824104075363
CTA_REFERENCE_SECTION3_X_TMAX = 0.2918596038697993
CTA_REFERENCE_SECTION4_X_TMAX = 0.38582456494467204
CTA_REFERENCE_SECTION5_X_TMAX = 0.38582456494467204

CTA_REFERENCE_SECTION0_TE_THICKNESS = 0.0017868134511835902
CTA_REFERENCE_SECTION1_TE_THICKNESS = 0.0018257217427832245
CTA_REFERENCE_SECTION2_TE_THICKNESS = 0.001853538450733595
CTA_REFERENCE_SECTION3_TE_THICKNESS = 0.003047758054176927
CTA_REFERENCE_SECTION4_TE_THICKNESS = 0.0023630630381066327
CTA_REFERENCE_SECTION5_TE_THICKNESS = 0.0021858619068493594

CTA_REFERENCE_C0_TC_MAX = CTA_REFERENCE_SECTION0_TC_MAX
CTA_REFERENCE_C3_TC_MAX = CTA_REFERENCE_SECTION3_TC_MAX
CTA_REFERENCE_C4_TC_MAX = CTA_REFERENCE_SECTION4_TC_MAX
CTA_REFERENCE_C5_TC_MAX = CTA_REFERENCE_SECTION5_TC_MAX

CTA_REFERENCE_C0_X_TMAX = CTA_REFERENCE_SECTION0_X_TMAX
CTA_REFERENCE_C3_X_TMAX = CTA_REFERENCE_SECTION3_X_TMAX
CTA_REFERENCE_C4_X_TMAX = CTA_REFERENCE_SECTION4_X_TMAX
CTA_REFERENCE_C5_X_TMAX = CTA_REFERENCE_SECTION5_X_TMAX

CTA_REFERENCE_C0_TE_THICKNESS = CTA_REFERENCE_SECTION0_TE_THICKNESS
CTA_REFERENCE_C3_TE_THICKNESS = CTA_REFERENCE_SECTION3_TE_THICKNESS
CTA_REFERENCE_C4_TE_THICKNESS = CTA_REFERENCE_SECTION4_TE_THICKNESS
CTA_REFERENCE_C5_TE_THICKNESS = CTA_REFERENCE_SECTION5_TE_THICKNESS

CTA_REFERENCE_C0_CAMBER_DELTA = 0.0
CTA_REFERENCE_C3_CAMBER_DELTA = 0.0
CTA_REFERENCE_C4_CAMBER_DELTA = 0.0
CTA_REFERENCE_C5_CAMBER_DELTA = 0.0

CTA_REFERENCE_PROFILE_ANCHOR_Y: Tuple[float, ...] = (
    0.0,
    CTA_REFERENCE_NOSE_HELPER_1_Y_M,
    CTA_REFERENCE_C1_Y_M,
    CTA_REFERENCE_B1_M,
    CTA_REFERENCE_B1_M + CTA_REFERENCE_B2_M,
    CTA_REFERENCE_SPAN_M,
)
CTA_REFERENCE_PROFILE_ANCHOR_LE_Z_M: Tuple[float, ...] = (
    0.25865,
    1.03474,
    0.48145,
    0.67693,
    0.72992,
    1.97538,
)
CTA_REFERENCE_PROFILE_ANCHOR_TWIST_DEG: Tuple[float, ...] = (
    0.778,
    -0.342,
    0.371,
    -0.249,
    0.483,
    3.177,
)


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


def _cta_reference_anchor_section_specs() -> Tuple[SectionCSTSpec, ...]:
    return (
        SectionCSTSpec(
            upper_coeffs=CTA_REFERENCE_SECTION0_UPPER_CST,
            lower_coeffs=CTA_REFERENCE_SECTION0_LOWER_CST,
            tc_max=CTA_REFERENCE_SECTION0_TC_MAX,
            x_tmax=CTA_REFERENCE_SECTION0_X_TMAX,
            te_thickness=CTA_REFERENCE_SECTION0_TE_THICKNESS,
        ),
        SectionCSTSpec(
            upper_coeffs=CTA_REFERENCE_SECTION1_UPPER_CST,
            lower_coeffs=CTA_REFERENCE_SECTION1_LOWER_CST,
            tc_max=CTA_REFERENCE_SECTION1_TC_MAX,
            x_tmax=CTA_REFERENCE_SECTION1_X_TMAX,
            te_thickness=CTA_REFERENCE_SECTION1_TE_THICKNESS,
        ),
        SectionCSTSpec(
            upper_coeffs=CTA_REFERENCE_SECTION2_UPPER_CST,
            lower_coeffs=CTA_REFERENCE_SECTION2_LOWER_CST,
            tc_max=CTA_REFERENCE_SECTION2_TC_MAX,
            x_tmax=CTA_REFERENCE_SECTION2_X_TMAX,
            te_thickness=CTA_REFERENCE_SECTION2_TE_THICKNESS,
        ),
        SectionCSTSpec(
            upper_coeffs=CTA_REFERENCE_SECTION3_UPPER_CST,
            lower_coeffs=CTA_REFERENCE_SECTION3_LOWER_CST,
            tc_max=CTA_REFERENCE_SECTION3_TC_MAX,
            x_tmax=CTA_REFERENCE_SECTION3_X_TMAX,
            te_thickness=CTA_REFERENCE_SECTION3_TE_THICKNESS,
        ),
        SectionCSTSpec(
            upper_coeffs=CTA_REFERENCE_SECTION4_UPPER_CST,
            lower_coeffs=CTA_REFERENCE_SECTION4_LOWER_CST,
            tc_max=CTA_REFERENCE_SECTION4_TC_MAX,
            x_tmax=CTA_REFERENCE_SECTION4_X_TMAX,
            te_thickness=CTA_REFERENCE_SECTION4_TE_THICKNESS,
        ),
        SectionCSTSpec(
            upper_coeffs=CTA_REFERENCE_SECTION5_UPPER_CST,
            lower_coeffs=CTA_REFERENCE_SECTION5_LOWER_CST,
            tc_max=CTA_REFERENCE_SECTION5_TC_MAX,
            x_tmax=CTA_REFERENCE_SECTION5_X_TMAX,
            te_thickness=CTA_REFERENCE_SECTION5_TE_THICKNESS,
        ),
    )


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
    use_reference_anchor_twist: bool = False,
) -> SectionedBWBModelConfig:
    """
    Build a model config enforcing CTA fixed choices:
    fixed sweep S, fixed CTA chord family, fixed B1 in meters,
    and fixed straight TE segments.
    """
    fixed_design = apply_cta_fixed_parameters(design, reference_design=reference_design)
    config = fixed_design.to_model_config(profile_generation_mode=CTA_REFERENCE_PROFILE_GENERATION_MODE)
    config.topology.anchor_y = CTA_REFERENCE_PROFILE_ANCHOR_Y
    config.sections = replace(
        config.sections,
        shared_leading_edge=False,
        section_specs_override=_cta_reference_anchor_section_specs(),
        profile_relations=tuple(
            SectionProfileRelationSpec() for _ in range(len(CTA_REFERENCE_PROFILE_ANCHOR_Y))
        ),
    )
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
    # Use linear profile interpolation between anchor sections for the CTA.
    # With the fitted glider sections, higher-order interpolation was creating
    # visible intermediate bumps that do not belong to the reference geometry.
    config.sampling.section_interpolation = "linear"
    # Keep sections 0-3 with no twist variation: the first actual twist change
    # is allowed only after C3.
    if use_reference_anchor_twist:
        # Reference-only mode: use the real chord-line inclinations extracted
        # from the glider sections so the z-y evolution is transferred to the
        # CTA planform consistently in the front view, 3D loft, and IGES.
        twist_values = tuple(float(value) for value in CTA_REFERENCE_PROFILE_ANCHOR_TWIST_DEG)
    else:
        # Design-space mode: keep the original grouped twist controls where
        # sections 0-3 move together and the first independent twist change
        # starts after C3.
        twist_values = (
            float(fixed_design.twist_c1_deg),
            float(fixed_design.twist_c1_deg),
            float(fixed_design.twist_c1_deg),
            float(fixed_design.twist_c1_deg),
            float(fixed_design.twist_c3_deg),
            float(fixed_design.twist_c4_deg),
        )
    config.spanwise.twist_deg = AnchoredSpanwiseLaw(
        section_indices=tuple(range(len(CTA_REFERENCE_PROFILE_ANCHOR_Y))),
        values=twist_values,
        interpolation="pchip",
    )
    config.spanwise.camber_delta = AnchoredSpanwiseLaw(
        section_indices=tuple(range(len(CTA_REFERENCE_PROFILE_ANCHOR_Y))),
        values=(
            float(fixed_design.camber_c1),
            float(fixed_design.camber_c1),
            float(fixed_design.camber_c1),
            float(fixed_design.camber_c1),
            float(fixed_design.camber_c3),
            float(fixed_design.camber_c4),
        ),
        interpolation="pyspline",
    )
    config.spanwise.vertical_offset_z = AnchoredSpanwiseLaw(
        section_indices=tuple(range(len(CTA_REFERENCE_PROFILE_ANCHOR_Y))),
        values=tuple(float(value) for value in CTA_REFERENCE_PROFILE_ANCHOR_LE_Z_M),
        interpolation="pchip",
    )
    return config
