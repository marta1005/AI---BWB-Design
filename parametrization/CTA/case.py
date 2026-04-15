from __future__ import annotations

from dataclasses import dataclass, replace
import json
from math import atan2, degrees, radians, tan
from pathlib import Path

from parametrization.bwb.case_definition import (
    CaseTemplate,
    CanonicalSectionDeclaration,
    ExplicitSectionSpec,
    Relation,
    build_case_config_from_explicit_sections,
    build_case_config_from_resolved,
    leading_edge_sweep_from_chord_sweep,
    chord_sweep_from_leading_edge_sweep,
    solve_relations,
)
from parametrization.bwb.design_variables import SectionedBWBDesignVariables
from parametrization.bwb.specs import (
    AnchoredSpanwiseLaw,
    SectionCSTSpec,
    SectionFamilySpec,
    SectionProfileRelationSpec,
    SectionedBWBModelConfig,
)


CTA_CASE_INPUT_PATH = Path(__file__).resolve().parent / "data" / "cta_case_inputs.json"

SWEEP_NAME_TO_VARIABLE: dict[str, str] = {
    "S": "s1_deg",
    "S1": "s2_deg",
    "S2": "s3_deg",
}
VARIABLE_TO_SWEEP_NAME: dict[str, str] = {value: key for key, value in SWEEP_NAME_TO_VARIABLE.items()}


@dataclass(frozen=True)
class CTAPublicDriverDeclaration:
    wing_span_m: float
    b2_over_wing_span: float
    c0_m: float
    c3_m: float
    wing_med_3_tr: float
    c5_m: float
    s1_50_deg: float
    s2_25_deg: float
    med_3_te_sweep_deg: float


@dataclass(frozen=True)
class CTAAirfoilDefinition:
    upper_cst: tuple[float, ...]
    lower_cst: tuple[float, ...]
    tc_max: float
    x_tmax: float
    te_thickness: float
    camber_delta: float


def load_cta_case_payload(path: Path | None = None) -> dict[str, object]:
    input_path = CTA_CASE_INPUT_PATH if path is None else Path(path)
    return json.loads(input_path.read_text(encoding="utf-8"))


def _section_spec_from_mapping(mapping: dict[str, object]) -> SectionCSTSpec:
    return SectionCSTSpec(
        upper_coeffs=tuple(float(value) for value in mapping["upper_coeffs"]),
        lower_coeffs=tuple(float(value) for value in mapping["lower_coeffs"]),
        tc_max=float(mapping["tc_max"]),
        x_tmax=float(mapping["x_tmax"]),
        te_thickness=float(mapping["te_thickness"]),
    )


def load_cta_canonical_declaration(payload: dict[str, object] | None = None) -> CanonicalSectionDeclaration:
    payload = load_cta_case_payload() if payload is None else payload
    return CanonicalSectionDeclaration(
        **{key: float(value) for key, value in payload["canonical_sections"].items()}
    )


def load_cta_public_declaration(payload: dict[str, object] | None = None) -> CTAPublicDriverDeclaration:
    payload = load_cta_case_payload() if payload is None else payload
    return CTAPublicDriverDeclaration(
        **{key: float(value) for key, value in payload["public_drivers"].items()}
    )


def load_cta_case_template(
    payload: dict[str, object] | None = None,
    canonical: CanonicalSectionDeclaration | None = None,
) -> CaseTemplate:
    payload = load_cta_case_payload() if payload is None else payload
    canonical = load_cta_canonical_declaration(payload) if canonical is None else canonical
    template_payload = payload["template"]
    return CaseTemplate(
        anchor_y_m=tuple(float(value) for value in template_payload["anchor_y_m"]),
        anchor_le_z_m=tuple(float(value) for value in template_payload["anchor_le_z_m"]),
        anchor_twist_deg=tuple(float(value) for value in template_payload["anchor_twist_deg"]),
        anchor_camber_delta=tuple(float(value) for value in template_payload["anchor_camber_delta"]),
        anchor_sections=tuple(
            _section_spec_from_mapping(section_payload) for section_payload in template_payload["anchor_sections"]
        ),
        cst_n1=float(template_payload["cst_n1"]),
        cst_n2=float(template_payload["cst_n2"]),
        shared_leading_edge=bool(template_payload["shared_leading_edge"]),
        profile_generation_mode=str(template_payload["profile_generation_mode"]),
        body_le_fixed_points=tuple(
            (float(point[0]), float(point[1])) for point in template_payload["body_le_fixed_points"]
        ),
        te_c1_y_m=float(template_payload["te_c1_span_fraction"]) * float(canonical.y_1_m),
        te_inboard_blend_y_m=float(template_payload["te_inboard_blend_fraction"]) * float(canonical.y_1_m),
        te_exact_segments=tuple(int(value) for value in template_payload["te_exact_segments"]),
        te_inboard_blend_dx_m=float(template_payload["te_inboard_blend_dx"]),
        te_inboard_radius_factor=float(template_payload["te_inboard_radius_factor"]),
        med_3_te_helper_fraction=float(template_payload["med_3_te_helper_fraction"]),
        te_outer_blend_fraction=float(template_payload["te_outer_blend_fraction"]),
        blend_fraction=float(template_payload["blend_fraction"]),
        min_linear_core_fraction=float(template_payload["min_linear_core_fraction"]),
        te_blend_fraction=(
            None if template_payload["te_blend_fraction"] is None else float(template_payload["te_blend_fraction"])
        ),
        te_min_linear_core_fraction=(
            None
            if template_payload["te_min_linear_core_fraction"] is None
            else float(template_payload["te_min_linear_core_fraction"])
        ),
        te_spline_bridge=(
            None
            if template_payload["te_spline_bridge"] is None
            else tuple(int(value) for value in template_payload["te_spline_bridge"])
        ),
        symmetry_blend_y_m=float(template_payload["symmetry_blend_y"]),
        section_interpolation=str(template_payload["section_interpolation"]),
        airfoil_distribution_mode=str(template_payload["airfoil_distribution_mode"]),
        num_airfoil_points=int(template_payload["num_airfoil_points"]),
        num_base_stations=int(template_payload["num_base_stations"]),
        section_curve_n_ctl=int(template_payload["section_curve_n_ctl"]),
        k_span=int(template_payload["k_span"]),
    )


def make_cta_anchor_y(template: CaseTemplate, resolved: dict[str, float]) -> tuple[float, ...]:
    return (
        0.0,
        template.anchor_y_m[1],
        template.anchor_y_m[2],
        float(resolved["y_1_m"]),
        float(resolved["y_2_m"]),
        float(resolved["y_3_m"]),
    )


def resolve_cta_from_sections(
    canonical: CanonicalSectionDeclaration,
    template: CaseTemplate,
) -> dict[str, float]:
    initial = canonical.__dict__.copy()
    relations = [
        Relation("span_m", ("y_3_m",), lambda y_3_m: y_3_m, "Semispan equals tip station."),
        Relation("b1_span_ratio", ("y_1_m", "span_m"), lambda y_1_m, span_m: y_1_m / span_m, "B1 ratio."),
        Relation(
            "b2_span_ratio",
            ("y_1_m", "y_2_m", "span_m"),
            lambda y_1_m, y_2_m, span_m: (y_2_m - y_1_m) / span_m,
            "B2 ratio.",
        ),
        Relation(
            "b3_span_ratio",
            ("y_2_m", "y_3_m", "span_m"),
            lambda y_2_m, y_3_m, span_m: (y_3_m - y_2_m) / span_m,
            "B3 ratio.",
        ),
        Relation("c2_c1_ratio", ("chord_0_m", "chord_1_m"), lambda c0, c1: c1 / c0, "Chord ratio root->sec1."),
        Relation("c3_c1_ratio", ("chord_0_m", "chord_2_m"), lambda c0, c2: c2 / c0, "Chord ratio root->sec2."),
        Relation("c4_c1_ratio", ("chord_0_m", "chord_3_m"), lambda c0, c3: c3 / c0, "Chord ratio root->tip."),
        Relation(
            "s1_deg",
            ("le_x_0_m", "le_x_1_m", "y_0_m", "y_1_m"),
            lambda x0, x1, y0, y1: degrees(atan2(x1 - x0, y1 - y0)),
            "Internal LE sweep of segment 0->1.",
        ),
        Relation(
            "s2_deg",
            ("le_x_1_m", "le_x_2_m", "y_1_m", "y_2_m"),
            lambda x1, x2, y1, y2: degrees(atan2(x2 - x1, y2 - y1)),
            "Internal LE sweep of segment 1->2.",
        ),
        Relation(
            "s3_deg",
            ("le_x_2_m", "le_x_3_m", "y_2_m", "y_3_m"),
            lambda x2, x3, y2, y3: degrees(atan2(x3 - x2, y3 - y2)),
            "Internal LE sweep of segment 2->3.",
        ),
        Relation(
            "te_c1_span_fraction",
            ("y_1_m",),
            lambda y_1_m: template.te_c1_y_m / y_1_m,
            "Keep C1 helper at fixed absolute spanwise location.",
        ),
        Relation(
            "te_inboard_blend_fraction",
            ("y_1_m",),
            lambda y_1_m: template.te_inboard_blend_y_m / y_1_m,
            "Keep inboard TE blend helper at fixed absolute spanwise location.",
        ),
    ]
    return solve_relations(initial, relations)


def resolve_cta_from_public(
    public: CTAPublicDriverDeclaration,
    canonical: CanonicalSectionDeclaration,
    template: CaseTemplate,
    fixed_internal_s1_deg: float,
) -> dict[str, float]:
    initial = public.__dict__.copy()
    initial["b1_m"] = float(canonical.y_1_m)
    initial["le_x_0_m"] = float(canonical.le_x_0_m)
    initial["fixed_internal_s1_deg"] = float(fixed_internal_s1_deg)
    relations = [
        Relation("y_0_m", (), lambda: 0.0, "Root section fixed at y=0."),
        Relation("y_1_m", ("b1_m",), lambda b1_m: b1_m, "Fixed B1."),
        Relation("y_3_m", ("b1_m", "wing_span_m"), lambda b1_m, wing_span_m: b1_m + wing_span_m, "Semispan."),
        Relation(
            "y_2_m",
            ("b1_m", "wing_span_m", "b2_over_wing_span"),
            lambda b1_m, wing_span_m, ratio: b1_m + ratio * wing_span_m,
            "Outer-body to wing break.",
        ),
        Relation("span_m", ("y_3_m",), lambda y_3_m: y_3_m, "Semispan equals tip station."),
        Relation("b1_span_ratio", ("y_1_m", "span_m"), lambda y_1_m, span_m: y_1_m / span_m, "B1 ratio."),
        Relation(
            "b2_span_ratio",
            ("y_1_m", "y_2_m", "span_m"),
            lambda y_1_m, y_2_m, span_m: (y_2_m - y_1_m) / span_m,
            "B2 ratio.",
        ),
        Relation(
            "b3_span_ratio",
            ("y_2_m", "y_3_m", "span_m"),
            lambda y_2_m, y_3_m, span_m: (y_3_m - y_2_m) / span_m,
            "B3 ratio.",
        ),
        Relation("chord_0_m", ("c0_m",), lambda c0_m: c0_m, "Root chord."),
        Relation("chord_1_m", ("c3_m",), lambda c3_m: c3_m, "Transition chord at section 1."),
        Relation(
            "chord_2_m",
            ("c3_m", "wing_med_3_tr"),
            lambda c3_m, wing_med_3_tr: c3_m * wing_med_3_tr,
            "Outer-wing root chord from taper ratio.",
        ),
        Relation("chord_3_m", ("c5_m",), lambda c5_m: c5_m, "Tip chord."),
        Relation("c2_c1_ratio", ("chord_0_m", "chord_1_m"), lambda c0, c1: c1 / c0, "Chord ratio root->sec1."),
        Relation("c3_c1_ratio", ("chord_0_m", "chord_2_m"), lambda c0, c2: c2 / c0, "Chord ratio root->sec2."),
        Relation("c4_c1_ratio", ("chord_0_m", "chord_3_m"), lambda c0, c3: c3 / c0, "Chord ratio root->tip."),
        Relation(
            "s1_deg",
            ("fixed_internal_s1_deg",),
            lambda fixed_internal_s1_deg: fixed_internal_s1_deg,
            "Keep the inboard body sweep fixed.",
        ),
        Relation(
            "s2_deg",
            ("s1_50_deg", "y_1_m", "y_2_m", "chord_1_m", "chord_2_m"),
            lambda s1_50_deg, y_1_m, y_2_m, c1, c2: leading_edge_sweep_from_chord_sweep(
                s1_50_deg, 0.50, y_2_m - y_1_m, c1, c2
            ),
            "Transition-wing LE sweep from 50% chord sweep.",
        ),
        Relation(
            "s3_deg",
            ("s2_25_deg", "y_2_m", "y_3_m", "chord_2_m", "chord_3_m"),
            lambda s2_25_deg, y_2_m, y_3_m, c2, c3: leading_edge_sweep_from_chord_sweep(
                s2_25_deg, 0.25, y_3_m - y_2_m, c2, c3
            ),
            "Outer-wing LE sweep from 25% chord sweep.",
        ),
        Relation(
            "le_x_1_m",
            ("le_x_0_m", "s1_deg", "y_0_m", "y_1_m"),
            lambda x0, s1_deg, y0, y1: x0 + tan(radians(s1_deg)) * (y1 - y0),
            "Section 1 LE x.",
        ),
        Relation(
            "le_x_2_m",
            ("le_x_1_m", "s2_deg", "y_1_m", "y_2_m"),
            lambda x1, s2_deg, y1, y2: x1 + tan(radians(s2_deg)) * (y2 - y1),
            "Section 2 LE x.",
        ),
        Relation(
            "le_x_3_m",
            ("le_x_2_m", "s3_deg", "y_2_m", "y_3_m"),
            lambda x2, s3_deg, y2, y3: x2 + tan(radians(s3_deg)) * (y3 - y2),
            "Section 3 LE x.",
        ),
        Relation(
            "te_c1_span_fraction",
            ("y_1_m",),
            lambda y_1_m: template.te_c1_y_m / y_1_m,
            "Keep C1 helper at fixed absolute spanwise location.",
        ),
        Relation(
            "te_inboard_blend_fraction",
            ("y_1_m",),
            lambda y_1_m: template.te_inboard_blend_y_m / y_1_m,
            "Keep inboard TE blend helper at fixed absolute spanwise location.",
        ),
    ]
    return solve_relations(initial, relations)


def build_cta_case_config_from_resolved(
    resolved: dict[str, float],
    template: CaseTemplate,
    topology_name: str = "cta_bwb_case",
) -> SectionedBWBModelConfig:
    return build_case_config_from_resolved(
        resolved=resolved,
        template=template,
        anchor_y=make_cta_anchor_y(template, resolved),
        topology_name=topology_name,
    )


def build_cta_explicit_sections_from_resolved(
    resolved: dict[str, float],
) -> tuple[ExplicitSectionSpec, ...]:
    return (
        ExplicitSectionSpec("C0", float(resolved["y_0_m"]), float(resolved["le_x_0_m"]), float(resolved["chord_0_m"])),
        ExplicitSectionSpec("C3", float(resolved["y_1_m"]), float(resolved["le_x_1_m"]), float(resolved["chord_1_m"])),
        ExplicitSectionSpec("C4", float(resolved["y_2_m"]), float(resolved["le_x_2_m"]), float(resolved["chord_2_m"])),
        ExplicitSectionSpec("C5", float(resolved["y_3_m"]), float(resolved["le_x_3_m"]), float(resolved["chord_3_m"])),
    )


def build_cta_leading_edge_control_points(
    resolved: dict[str, float],
    template: CaseTemplate,
) -> tuple[tuple[float, float], ...]:
    points = [(float(resolved["le_x_0_m"]), float(resolved["y_0_m"]))]
    points.extend((float(x_value), float(y_value)) for x_value, y_value in template.body_le_fixed_points)
    points.extend(
        (
            (float(resolved["le_x_1_m"]), float(resolved["y_1_m"])),
            (float(resolved["le_x_2_m"]), float(resolved["y_2_m"])),
            (float(resolved["le_x_3_m"]), float(resolved["y_3_m"])),
        )
    )
    return tuple(points)


def build_cta_trailing_edge_control_points(
    resolved: dict[str, float],
    template: CaseTemplate,
) -> tuple[tuple[float, float], ...]:
    y_sections = (
        float(resolved["y_0_m"]),
        float(resolved["y_1_m"]),
        float(resolved["y_2_m"]),
        float(resolved["y_3_m"]),
    )
    le_sections = (
        float(resolved["le_x_0_m"]),
        float(resolved["le_x_1_m"]),
        float(resolved["le_x_2_m"]),
        float(resolved["le_x_3_m"]),
    )
    chords = (
        float(resolved["chord_0_m"]),
        float(resolved["chord_1_m"]),
        float(resolved["chord_2_m"]),
        float(resolved["chord_3_m"]),
    )
    te_sections = tuple(le_x + chord for le_x, chord in zip(le_sections, chords))

    y_c3 = float(y_sections[1])
    y_c4 = float(y_sections[2])
    y_tip = float(y_sections[3])
    y_c1 = float(template.te_c1_y_m)
    y_inboard_blend = float(template.te_inboard_blend_y_m)

    te_root = float(te_sections[0])
    te_c3 = float(te_sections[1])
    te_c4 = float(te_sections[2])
    te_tip = float(te_sections[3])

    te_c1 = te_root
    if abs(float(template.te_inboard_blend_dx_m)) > 1e-12:
        te_inboard_blend = te_c3 + float(template.te_inboard_blend_dx_m)
        te_inboard_blend = float(min(max(te_c3, te_inboard_blend), te_c1))
    else:
        blend_ratio = (y_inboard_blend - y_c1) / max(y_c3 - y_c1, 1e-12)
        te_inboard_blend = te_c1 + blend_ratio * (te_c3 - te_c1)

    points = [
        (te_root, 0.0),
        (te_c1, y_c1),
        (float(te_inboard_blend), y_inboard_blend),
        (te_c3, y_c3),
    ]

    med_3_te_sweep_deg = float(resolved.get("med_3_te_sweep_deg", 0.0))
    if abs(med_3_te_sweep_deg) > 1e-12:
        y_med3 = y_c3 + float(template.med_3_te_helper_fraction) * (y_c4 - y_c3)
        te_med3 = te_c3 + tan(radians(med_3_te_sweep_deg)) * (y_med3 - y_c3)
        points.append((float(te_med3), float(y_med3)))

    points.append((te_c4, y_c4))

    if float(template.te_outer_blend_fraction) > 1e-12:
        y_outer_blend = y_c4 + float(template.te_outer_blend_fraction) * (y_tip - y_c4)
        outer_ratio = (y_outer_blend - y_c4) / max(y_tip - y_c4, 1e-12)
        te_outer_blend = te_c4 + outer_ratio * (te_tip - te_c4)
        points.append((float(te_outer_blend), float(y_outer_blend)))

    points.append((te_tip, y_tip))
    return tuple(points)


def build_cta_case_target_config() -> SectionedBWBModelConfig:
    payload = load_cta_case_payload()
    canonical = load_cta_canonical_declaration(payload)
    template = load_cta_case_template(payload, canonical)
    resolved = resolve_cta_from_sections(canonical, template)
    return build_case_config_from_explicit_sections(
        sections_definition=build_cta_explicit_sections_from_resolved(resolved),
        template=template,
        topology_name="cta_bwb_case",
        leading_edge_control_points=build_cta_leading_edge_control_points(resolved, template),
        trailing_edge_control_points=build_cta_trailing_edge_control_points(resolved, template),
        te_exact_segments=tuple(int(value) for value in template.te_exact_segments),
        te_spline_bridge=template.te_spline_bridge,
    )


def derive_cta_public_from_canonical(canonical: CanonicalSectionDeclaration) -> CTAPublicDriverDeclaration:
    wing_span_m = float(canonical.y_3_m - canonical.y_1_m)
    b2_over_wing_span = float((canonical.y_2_m - canonical.y_1_m) / max(wing_span_m, 1e-12))
    s1_50_deg = chord_sweep_from_leading_edge_sweep(
        le_sweep_deg=degrees(atan2(canonical.le_x_2_m - canonical.le_x_1_m, canonical.y_2_m - canonical.y_1_m)),
        chord_fraction=0.50,
        dy_m=canonical.y_2_m - canonical.y_1_m,
        chord_in_m=canonical.chord_1_m,
        chord_out_m=canonical.chord_2_m,
    )
    s2_25_deg = chord_sweep_from_leading_edge_sweep(
        le_sweep_deg=degrees(atan2(canonical.le_x_3_m - canonical.le_x_2_m, canonical.y_3_m - canonical.y_2_m)),
        chord_fraction=0.25,
        dy_m=canonical.y_3_m - canonical.y_2_m,
        chord_in_m=canonical.chord_2_m,
        chord_out_m=canonical.chord_3_m,
    )
    return CTAPublicDriverDeclaration(
        wing_span_m=wing_span_m,
        b2_over_wing_span=b2_over_wing_span,
        c0_m=float(canonical.chord_0_m),
        c3_m=float(canonical.chord_1_m),
        wing_med_3_tr=float(canonical.chord_2_m / max(canonical.chord_1_m, 1e-12)),
        c5_m=float(canonical.chord_3_m),
        s1_50_deg=s1_50_deg,
        s2_25_deg=s2_25_deg,
        med_3_te_sweep_deg=0.0,
    )


_CTA_CASE_PAYLOAD = load_cta_case_payload()
_CTA_CANONICAL = load_cta_canonical_declaration(_CTA_CASE_PAYLOAD)
_CTA_PUBLIC = load_cta_public_declaration(_CTA_CASE_PAYLOAD)
_CTA_TEMPLATE = load_cta_case_template(_CTA_CASE_PAYLOAD, _CTA_CANONICAL)
_CTA_RESOLVED = resolve_cta_from_sections(_CTA_CANONICAL, _CTA_TEMPLATE)

FIXED_TE_EXACT_SEGMENTS: tuple[int, ...] = tuple(int(value) for value in _CTA_TEMPLATE.te_exact_segments)

CTA_SPAN_M = float(_CTA_CANONICAL.y_3_m)
CTA_B1_M = float(_CTA_CANONICAL.y_1_m)
CTA_B2_M = float(_CTA_CANONICAL.y_2_m - _CTA_CANONICAL.y_1_m)
CTA_B3_M = float(_CTA_CANONICAL.y_3_m - _CTA_CANONICAL.y_2_m)

CTA_C0_BODY_CHORD_M = float(_CTA_CANONICAL.chord_0_m)
CTA_C3_TRANSITION_CHORD_M = float(_CTA_CANONICAL.chord_1_m)
CTA_C4_WING_CHORD_M = float(_CTA_CANONICAL.chord_2_m)
CTA_C5_WING_TIP_M = float(_CTA_CANONICAL.chord_3_m)
CTA_NOSE_HELPER_1_Y_M = float(_CTA_TEMPLATE.anchor_y_m[1])
CTA_NOSE_HELPER_1_X_M = float(_CTA_TEMPLATE.body_le_fixed_points[0][0])
CTA_C1_Y_M = float(_CTA_TEMPLATE.anchor_y_m[2])
CTA_C1_LE_HELPER_X_M = float(_CTA_TEMPLATE.body_le_fixed_points[1][0])
CTA_TE_INBOARD_BLEND_DX_M = float(_CTA_TEMPLATE.te_inboard_blend_dx_m)
CTA_TE_INBOARD_BLEND_Y_M = float(_CTA_TEMPLATE.te_inboard_blend_y_m)
CTA_TE_INBOARD_RADIUS_FACTOR = float(_CTA_TEMPLATE.te_inboard_radius_factor)

CTA_S_DEG = float(degrees(atan2(CTA_NOSE_HELPER_1_X_M, CTA_NOSE_HELPER_1_Y_M)))
CTA_C1_SWEEP_DEG = float(
    degrees(
        atan2(
            _CTA_CANONICAL.le_x_1_m - CTA_C1_LE_HELPER_X_M,
            _CTA_CANONICAL.y_1_m - CTA_C1_Y_M,
        )
    )
)
CTA_S1_DEG = float(_CTA_PUBLIC.s1_50_deg)
CTA_S2_DEG = float(_CTA_PUBLIC.s2_25_deg)
CTA_MED_3_TESWP_DEG = float(_CTA_PUBLIC.med_3_te_sweep_deg)
CTA_MED_3_TE_HELPER_FRACTION = float(_CTA_TEMPLATE.med_3_te_helper_fraction)

CTA_PROFILE_GENERATION_MODE = str(_CTA_TEMPLATE.profile_generation_mode)
CTA_CST_N1 = float(_CTA_TEMPLATE.cst_n1)
CTA_CST_N2 = float(_CTA_TEMPLATE.cst_n2)

def _cta_airfoil_definition(anchor_index: int) -> CTAAirfoilDefinition:
    spec = _CTA_TEMPLATE.anchor_sections[int(anchor_index)]
    return CTAAirfoilDefinition(
        upper_cst=tuple(float(value) for value in spec.upper_coeffs),
        lower_cst=tuple(float(value) for value in spec.lower_coeffs),
        tc_max=float(spec.tc_max),
        x_tmax=float(spec.x_tmax),
        te_thickness=float(spec.te_thickness),
        camber_delta=float(_CTA_TEMPLATE.anchor_camber_delta[int(anchor_index)]),
    )


_CTA_ACTIVE_AIRFOIL_ANCHOR_INDICES = (0, 3, 4, 5)
_CTA_C0_AIRFOIL, _CTA_C3_AIRFOIL, _CTA_C4_AIRFOIL, _CTA_C5_AIRFOIL = tuple(
    _cta_airfoil_definition(idx) for idx in _CTA_ACTIVE_AIRFOIL_ANCHOR_INDICES
)

CTA_PROFILE_ANCHOR_Y: tuple[float, ...] = tuple(float(value) for value in _CTA_TEMPLATE.anchor_y_m)
CTA_PROFILE_ANCHOR_LE_Z_M: tuple[float, ...] = tuple(float(value) for value in _CTA_TEMPLATE.anchor_le_z_m)
CTA_PROFILE_ANCHOR_TWIST_DEG: tuple[float, ...] = tuple(float(value) for value in _CTA_TEMPLATE.anchor_twist_deg)

_leading_edge_sweep_from_chord_sweep = leading_edge_sweep_from_chord_sweep
_chord_sweep_from_leading_edge_sweep = chord_sweep_from_leading_edge_sweep


def _cta_c3_transition_chord_m() -> float:
    return float(CTA_C3_TRANSITION_CHORD_M)


def _cta_le_helper_1_x_m() -> float:
    return float(CTA_NOSE_HELPER_1_X_M)


def _cta_le_helper_2_x_m() -> float:
    return float(CTA_C1_LE_HELPER_X_M)


def _cta_le_c3_x_m() -> float:
    dy = float(CTA_B1_M - CTA_C1_Y_M)
    return float(_cta_le_helper_2_x_m() + tan(radians(CTA_C1_SWEEP_DEG)) * dy)


def _cta_internal_secant_s_deg() -> float:
    return float(degrees(atan2(_cta_le_c3_x_m(), CTA_B1_M)))


def _cta_anchor_section_specs() -> tuple[SectionCSTSpec, ...]:
    return tuple(
        SectionCSTSpec(
            upper_coeffs=tuple(float(value) for value in spec.upper_coeffs),
            lower_coeffs=tuple(float(value) for value in spec.lower_coeffs),
            tc_max=float(spec.tc_max),
            x_tmax=float(spec.x_tmax),
            te_thickness=float(spec.te_thickness),
        )
        for spec in _CTA_TEMPLATE.anchor_sections
    )


def _cta_design_airfoil_kwargs() -> dict[str, object]:
    airfoils = (_CTA_C0_AIRFOIL, _CTA_C3_AIRFOIL, _CTA_C4_AIRFOIL, _CTA_C5_AIRFOIL)
    kwargs: dict[str, object] = {}
    for idx, airfoil in enumerate(airfoils, start=1):
        kwargs[f"c{idx}_tc_max"] = airfoil.tc_max
        kwargs[f"c{idx}_x_tmax"] = airfoil.x_tmax
        kwargs[f"c{idx}_te_thickness"] = airfoil.te_thickness
        kwargs[f"camber_c{idx}"] = airfoil.camber_delta
        kwargs[f"c{idx}_upper_cst"] = airfoil.upper_cst
        kwargs[f"c{idx}_lower_cst"] = airfoil.lower_cst
    return kwargs


def build_cta_design() -> SectionedBWBDesignVariables:
    base = SectionedBWBDesignVariables.seed()
    return replace(
        base,
        span=CTA_SPAN_M,
        cst_n1=CTA_CST_N1,
        cst_n2=CTA_CST_N2,
        b1_span_ratio=CTA_B1_M / CTA_SPAN_M,
        b2_span_ratio=CTA_B2_M / CTA_SPAN_M,
        b3_span_ratio=CTA_B3_M / CTA_SPAN_M,
        c1_root_chord=CTA_C0_BODY_CHORD_M,
        c2_c1_ratio=_cta_c3_transition_chord_m() / CTA_C0_BODY_CHORD_M,
        c3_c1_ratio=CTA_C4_WING_CHORD_M / CTA_C0_BODY_CHORD_M,
        c4_c1_ratio=CTA_C5_WING_TIP_M / CTA_C0_BODY_CHORD_M,
        twist_c2_deg=base.twist_c1_deg,
        s1_deg=_cta_internal_secant_s_deg(),
        s2_deg=CTA_S1_DEG,
        s3_deg=CTA_S2_DEG,
        med_3_te_sweep_deg=CTA_MED_3_TESWP_DEG,
        **_cta_design_airfoil_kwargs(),
    )


def _cta_fixed_b1_length(
    cta_design: SectionedBWBDesignVariables | None = None,
) -> float:
    design = build_cta_design() if cta_design is None else cta_design
    return float(design.span * design.b1_span_ratio)


def _resolve_cta_span_partition(
    design: SectionedBWBDesignVariables,
    fixed_b1_length: float,
) -> tuple[float, float, float]:
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
    cta_design: SectionedBWBDesignVariables | None = None,
) -> dict[str, object]:
    design = build_cta_design() if cta_design is None else cta_design
    b1_fixed_m = _cta_fixed_b1_length(design)
    c0_body_chord_m = float(design.c1_root_chord)
    c3_transition_chord_m = float(design.c2_c1_ratio * design.c1_root_chord)
    c4_wing_chord_m = float(design.c3_c1_ratio * design.c1_root_chord)
    c5_wing_tip_m = float(design.c4_c1_ratio * design.c1_root_chord)
    wing_span_m = float(design.span - b1_fixed_m)
    b2_wing_span_ratio = float((design.span * design.b2_span_ratio) / max(wing_span_m, 1e-12))
    return {
        "s_deg": float(CTA_S_DEG),
        "c1_sweep_deg": float(CTA_C1_SWEEP_DEG),
        "s_internal_deg": float(design.s1_deg),
        "b1_fixed_m": b1_fixed_m,
        "c0_body_chord_m": c0_body_chord_m,
        "c3_transition_chord_m": c3_transition_chord_m,
        "c4_wing_chord_m": c4_wing_chord_m,
        "c5_wing_tip_m": c5_wing_tip_m,
        "wing_span_m": wing_span_m,
        "b2_wing_span_ratio": b2_wing_span_ratio,
        "te_exact_segments": FIXED_TE_EXACT_SEGMENTS,
    }


def apply_cta_fixed_parameters(
    design: SectionedBWBDesignVariables,
    cta_design: SectionedBWBDesignVariables | None = None,
) -> SectionedBWBDesignVariables:
    fixed = cta_fixed_values(cta_design=cta_design)
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
    current_b3_m = float(design.span * b3_span_ratio)
    current_c4_wing_chord_m = float(design.c3_c1_ratio * current_c0_body_chord_m)
    if current_c4_wing_chord_m <= 0.0:
        raise ValueError(
            f"CTA C4/outer-wing chord must stay positive, got {current_c4_wing_chord_m:.6f} m"
        )
    internal_s1_deg = leading_edge_sweep_from_chord_sweep(
        chord_sweep_deg=float(design.s2_deg),
        chord_fraction=0.50,
        dy_m=current_b2_m,
        chord_in_m=current_c3_transition_chord_m,
        chord_out_m=current_c4_wing_chord_m,
    )
    internal_s2_deg = leading_edge_sweep_from_chord_sweep(
        chord_sweep_deg=float(design.s3_deg),
        chord_fraction=0.25,
        dy_m=current_b3_m,
        chord_in_m=current_c4_wing_chord_m,
        chord_out_m=current_c5_wing_tip_m,
    )
    return replace(
        design,
        b1_span_ratio=b1_span_ratio,
        b2_span_ratio=b2_span_ratio,
        b3_span_ratio=b3_span_ratio,
        s1_deg=float(fixed["s_internal_deg"]),
        s2_deg=internal_s1_deg,
        s3_deg=internal_s2_deg,
        c2_c1_ratio=float(current_c3_transition_chord_m / current_c0_body_chord_m),
        c3_c1_ratio=float(current_c4_wing_chord_m / current_c0_body_chord_m),
        c4_c1_ratio=float(current_c5_wing_tip_m / current_c0_body_chord_m),
        twist_c2_deg=float(design.twist_c1_deg),
    )


def to_cta_model_config(
    design: SectionedBWBDesignVariables,
    cta_design: SectionedBWBDesignVariables | None = None,
    use_cta_anchor_twist: bool = False,
) -> SectionedBWBModelConfig:
    fixed_design = apply_cta_fixed_parameters(design, cta_design=cta_design)
    config = fixed_design.to_model_config(profile_generation_mode=CTA_PROFILE_GENERATION_MODE)
    config.topology.anchor_y = CTA_PROFILE_ANCHOR_Y
    config.sections = replace(
        config.sections,
        shared_leading_edge=False,
        section_specs_override=_cta_anchor_section_specs(),
        profile_relations=tuple(
            SectionProfileRelationSpec() for _ in range(len(CTA_PROFILE_ANCHOR_Y))
        ),
    )
    config.planform.te_exact_segments = FIXED_TE_EXACT_SEGMENTS
    if abs(float(fixed_design.med_3_te_sweep_deg)) > 1e-12:
        config.planform.te_exact_segments = (0, 5)
    config.planform.te_c1_span_fraction = CTA_C1_Y_M / CTA_B1_M
    config.planform.te_inboard_blend_fraction = CTA_TE_INBOARD_BLEND_Y_M / CTA_B1_M
    config.planform.te_inboard_blend_dx = CTA_TE_INBOARD_BLEND_DX_M
    config.planform.te_inboard_radius_factor = CTA_TE_INBOARD_RADIUS_FACTOR
    config.planform.med_3_te_sweep_deg = float(fixed_design.med_3_te_sweep_deg)
    config.planform.med_3_te_helper_fraction = CTA_MED_3_TE_HELPER_FRACTION
    config.planform.te_outer_blend_fraction = 0.0
    config.planform.te_blend_fraction = 0.44
    config.planform.te_min_linear_core_fraction = 0.08
    config.planform.te_spline_bridge = (1, 3)
    config.planform.blend_fraction = 0.24
    config.planform.min_linear_core_fraction = 0.58
    config.planform.body_le_fixed_points = (
        (_cta_le_helper_1_x_m(), CTA_NOSE_HELPER_1_Y_M),
        (_cta_le_helper_2_x_m(), CTA_C1_Y_M),
    )
    config.planform.symmetry_blend_y = CTA_NOSE_HELPER_1_Y_M
    config.sampling.airfoil_distribution_mode = "anchors"
    config.sampling.section_interpolation = "linear"
    if use_cta_anchor_twist:
        twist_values = tuple(float(value) for value in CTA_PROFILE_ANCHOR_TWIST_DEG)
    else:
        twist_values = (
            float(fixed_design.twist_c1_deg),
            float(fixed_design.twist_c1_deg),
            float(fixed_design.twist_c1_deg),
            float(fixed_design.twist_c1_deg),
            float(fixed_design.twist_c3_deg),
            float(fixed_design.twist_c4_deg),
        )
    config.spanwise.twist_deg = AnchoredSpanwiseLaw(
        section_indices=tuple(range(len(CTA_PROFILE_ANCHOR_Y))),
        values=twist_values,
        interpolation="pchip",
    )
    config.spanwise.camber_delta = AnchoredSpanwiseLaw(
        section_indices=tuple(range(len(CTA_PROFILE_ANCHOR_Y))),
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
        section_indices=tuple(range(len(CTA_PROFILE_ANCHOR_Y))),
        values=tuple(float(value) for value in CTA_PROFILE_ANCHOR_LE_Z_M),
        interpolation="pchip",
    )
    return config
