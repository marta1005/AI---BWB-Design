from __future__ import annotations

from dataclasses import dataclass, replace
from functools import lru_cache
import json
from math import atan2, degrees, radians, tan
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

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
from parametrization.bwb.internal_volume_constraints import (
    CadReferenceFrame,
    load_internal_volume_constraint_set,
    required_constraint_bounds_at_points,
)
from parametrization.bwb.planform import build_sectioned_bwb_planform
from parametrization.bwb.sections import build_section_model
from parametrization.bwb.specs import (
    AnchoredSpanwiseLaw,
    ExactSectionProfileOverrideSpec,
    SectionCSTSpec,
    SectionProfileRelationSpec,
    SectionedBWBModelConfig,
)
from parametrization.bwb.spanwise_laws import (
    build_anchored_interpolant,
    build_scalar_interpolant,
    quintic_c2_transition,
    resolve_spanwise_laws,
    vertical_offsets,
)
from parametrization.shared.airfoil_fit import load_xyz_sections_by_span
from parametrization.shared.cst import KulfanCSTAirfoil, cosine_spacing
CTA_CASE_INPUT_PATH = Path(__file__).resolve().parent / "data" / "cta_case_inputs.json"
CTA_GLIDER_GEO_PATH = Path(__file__).resolve().parent / "airfoils" / "bwb_glider.geo"
CTA_INTERNAL_VOLUME_CONSTRAINTS_PATH = (
    Path(__file__).resolve().parent.parent.parent
    / "BWB B359-V0 - Internal Volume constraint - Set1 - BWB B359-V0 - Internal Volume constraint - Set1.csv"
)
CTA_CAD_REFERENCE_FRAME = CadReferenceFrame(
    name="CTA CAD reference",
    offset_x_m=0.077261,
    offset_y_m=0.0,
    offset_z_m=0.0,
    mirror_about_symmetry_plane=True,
)

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
    fixed_prefix = tuple(float(value) for value in template.anchor_y_m[:-3])
    return fixed_prefix + (
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
CTA_LE_ROOT_X_M = float(_CTA_CANONICAL.le_x_0_m)
CTA_NOSE_HELPER_1_Y_M = float(_CTA_TEMPLATE.body_le_fixed_points[0][1])
CTA_NOSE_HELPER_1_X_M = float(_CTA_TEMPLATE.body_le_fixed_points[0][0])
CTA_C1_Y_M = float(_CTA_TEMPLATE.body_le_fixed_points[1][1])
CTA_C1_LE_HELPER_X_M = float(_CTA_TEMPLATE.body_le_fixed_points[1][0])
CTA_TE_INBOARD_BLEND_DX_M = float(_CTA_TEMPLATE.te_inboard_blend_dx_m)
CTA_TE_INBOARD_BLEND_Y_M = float(_CTA_TEMPLATE.te_inboard_blend_y_m)
CTA_TE_INBOARD_RADIUS_FACTOR = float(_CTA_TEMPLATE.te_inboard_radius_factor)
CTA_ROOT_FRONT_BLEND_FRACTION = 0.74
CTA_ROOT_PROFILE_BLEND_FRACTION = 0.68
CTA_ROOT_FRONT_SOLVE_TWIST_WINDOW_DEG = 6.0
CTA_ROOT_FRONT_SOLVE_GRID_SIZE = 121
CTA_ROOT_FRONT_SLOPE_H_M = 1.0e-5
CTA_ROOT_TE_FLATTEN_X_START = 0.92

CTA_S_DEG = float(degrees(atan2(CTA_NOSE_HELPER_1_X_M - CTA_LE_ROOT_X_M, CTA_NOSE_HELPER_1_Y_M)))
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


_CTA_ACTIVE_AIRFOIL_ANCHOR_INDICES = (
    0,
    len(_CTA_TEMPLATE.anchor_sections) - 3,
    len(_CTA_TEMPLATE.anchor_sections) - 2,
    len(_CTA_TEMPLATE.anchor_sections) - 1,
)
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
    return float(degrees(atan2(_cta_le_c3_x_m() - CTA_LE_ROOT_X_M, CTA_B1_M)))


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


def _cta_anchor_section_specs_from_active(
    active_specs: tuple[SectionCSTSpec, ...],
) -> tuple[SectionCSTSpec, ...]:
    anchor_specs = list(_cta_anchor_section_specs())
    for active_idx, anchor_idx in enumerate(_CTA_ACTIVE_AIRFOIL_ANCHOR_INDICES):
        source = active_specs[active_idx]
        anchor_specs[anchor_idx] = SectionCSTSpec(
            upper_coeffs=tuple(float(value) for value in source.upper_coeffs),
            lower_coeffs=tuple(float(value) for value in source.lower_coeffs),
            tc_max=float(source.tc_max),
            x_tmax=float(source.x_tmax),
            te_thickness=float(source.te_thickness),
        )
    return tuple(anchor_specs)


def _cta_front_local_envelope(
    y: float,
    planform,
    section_model,
    twist_deg: float,
) -> tuple[float, float]:
    yu, yl, _ = section_model.coordinates_at_y(float(y))
    chord = float(planform.te_x(float(y)) - planform.le_x(float(y)))
    tr = np.deg2rad(float(twist_deg))
    x_term = section_model.x_air * chord * np.sin(tr)
    z_upper = x_term + yu * chord * np.cos(tr)
    z_lower = x_term + yl * chord * np.cos(tr)
    return float(np.max(z_upper)), float(np.min(z_lower))


def _cta_front_envelope(
    y: float,
    planform,
    section_model,
    twist_deg: float,
    vertical_z: float,
) -> tuple[float, float]:
    local_upper, local_lower = _cta_front_local_envelope(
        y=float(y),
        planform=planform,
        section_model=section_model,
        twist_deg=float(twist_deg),
    )
    z_shift = float(vertical_z)
    return local_upper + z_shift, local_lower + z_shift


def _cta_te_local_pair(
    y: float,
    planform,
    section_model,
    twist_deg: float,
) -> tuple[float, float]:
    yu, yl, _ = section_model.coordinates_at_y(float(y))
    chord = float(planform.te_x(float(y)) - planform.le_x(float(y)))
    tr = np.deg2rad(float(twist_deg))
    x_te = float(section_model.x_air[-1]) * chord
    z_upper_te = x_te * np.sin(tr) + float(yu[-1]) * chord * np.cos(tr)
    z_lower_te = x_te * np.sin(tr) + float(yl[-1]) * chord * np.cos(tr)
    return float(z_upper_te), float(z_lower_te)


def _cta_central_slope(
    fun,
    y: float,
    h: float = CTA_ROOT_FRONT_SLOPE_H_M,
) -> float:
    yy = float(y)
    hh = float(max(h, 1.0e-8))
    return float((fun(yy + hh) - fun(yy - hh)) / (2.0 * hh))


def _cta_solve_twist_for_target_front_thickness(
    y: float,
    target_thickness_m: float,
    current_twist_deg: float,
    planform,
    section_model,
) -> float:
    yy = float(y)
    current_twist = float(current_twist_deg)

    def residual(twist_deg: float) -> float:
        local_upper, local_lower = _cta_front_local_envelope(
            y=yy,
            planform=planform,
            section_model=section_model,
            twist_deg=float(twist_deg),
        )
        return float((local_upper - local_lower) - target_thickness_m)

    grid = np.linspace(
        current_twist - CTA_ROOT_FRONT_SOLVE_TWIST_WINDOW_DEG,
        current_twist + CTA_ROOT_FRONT_SOLVE_TWIST_WINDOW_DEG,
        CTA_ROOT_FRONT_SOLVE_GRID_SIZE,
    )
    residuals = np.asarray([residual(float(value)) for value in grid], dtype=float)

    candidate_intervals: list[tuple[float, float, float]] = []
    for left, right, res_left, res_right in zip(
        grid[:-1],
        grid[1:],
        residuals[:-1],
        residuals[1:],
    ):
        if res_left == 0.0:
            return float(left)
        if res_left * res_right <= 0.0:
            center = 0.5 * (float(left) + float(right))
            candidate_intervals.append((abs(center - current_twist), float(left), float(right)))

    if candidate_intervals:
        _, left, right = min(candidate_intervals, key=lambda item: item[0])
        return float(brentq(residual, left, right))

    best_idx = int(np.argmin(np.abs(residuals)))
    return float(grid[best_idx])


def _cta_root_front_matched_callables(
    config: SectionedBWBModelConfig,
):
    planform = build_sectioned_bwb_planform(config.topology, config.planform)
    base_laws = resolve_spanwise_laws(config)
    section_model = build_section_model(config, base_laws)
    root_y = 0.0
    join_y = float(config.topology.anchor_y_array[1])
    if config.spanwise.vertical_offset_z is not None:
        base_vertical_interpolant = build_anchored_interpolant(config.topology, config.spanwise.vertical_offset_z)
    else:
        dihedral_tan = float(np.tan(np.deg2rad(config.spanwise.dihedral_deg)))
        vertical_offset_m = float(config.spanwise.vertical_offset_m)

        def base_vertical_interpolant(yy: float) -> float:
            return float(vertical_offset_m + dihedral_tan * float(yy))

    def base_vertical(yy: float) -> float:
        return float(base_vertical_interpolant(float(yy)))

    def base_front_upper(yy: float) -> float:
        upper, _ = _cta_front_envelope(
            y=float(yy),
            planform=planform,
            section_model=section_model,
            twist_deg=float(base_laws.twist_deg(float(yy))),
            vertical_z=base_vertical(float(yy)),
        )
        return float(upper)

    def base_front_lower(yy: float) -> float:
        _, lower = _cta_front_envelope(
            y=float(yy),
            planform=planform,
            section_model=section_model,
            twist_deg=float(base_laws.twist_deg(float(yy))),
            vertical_z=base_vertical(float(yy)),
        )
        return float(lower)

    root_upper = base_front_upper(root_y)
    root_lower = base_front_lower(root_y)
    join_upper = base_front_upper(join_y)
    join_lower = base_front_lower(join_y)
    join_upper_slope = _cta_central_slope(base_front_upper, join_y)
    join_lower_slope = _cta_central_slope(base_front_lower, join_y)

    @lru_cache(maxsize=512)
    def solved_root_state(y_key: float) -> tuple[float, float]:
        yy = float(np.clip(float(y_key), root_y, join_y))
        target_upper = quintic_c2_transition(
            yy,
            root_y,
            join_y,
            root_upper,
            join_upper,
            0.0,
            join_upper_slope,
        )
        target_lower = quintic_c2_transition(
            yy,
            root_y,
            join_y,
            root_lower,
            join_lower,
            0.0,
            join_lower_slope,
        )
        target_thickness = float(target_upper - target_lower)
        current_twist = float(base_laws.twist_deg(yy))
        solved_twist = _cta_solve_twist_for_target_front_thickness(
            y=yy,
            target_thickness_m=target_thickness,
            current_twist_deg=current_twist,
            planform=planform,
            section_model=section_model,
        )
        local_upper, local_lower = _cta_front_local_envelope(
            y=yy,
            planform=planform,
            section_model=section_model,
            twist_deg=solved_twist,
        )
        solved_vertical = 0.5 * (
            (float(target_upper) + float(target_lower))
            - (float(local_upper) + float(local_lower))
        )
        return float(solved_twist), float(solved_vertical)

    def twist_callable(yy: float) -> float:
        y_val = float(yy)
        if y_val >= join_y - 1.0e-12:
            return float(base_laws.twist_deg(y_val))
        return float(solved_root_state(round(max(y_val, root_y), 10))[0])

    def vertical_callable(yy: float) -> float:
        y_val = float(yy)
        if y_val >= join_y - 1.0e-12:
            return float(base_vertical(y_val))
        return float(solved_root_state(round(max(y_val, root_y), 10))[1])

    return twist_callable, vertical_callable


def _cta_root_te_flatten_profile_transform(
    config: SectionedBWBModelConfig,
):
    base_sections = replace(config.sections, profile_post_transform=None)
    base_config = replace(config, sections=base_sections)
    planform = build_sectioned_bwb_planform(base_config.topology, base_config.planform)
    laws = resolve_spanwise_laws(base_config)
    section_model = build_section_model(base_config, laws)
    root_y = 0.0
    join_y = float(base_config.topology.anchor_y_array[1])
    if base_config.spanwise.vertical_offset_callable is not None:
        base_vertical = base_config.spanwise.vertical_offset_callable
    elif base_config.spanwise.vertical_offset_z is not None:
        interpolant = build_anchored_interpolant(base_config.topology, base_config.spanwise.vertical_offset_z)

        def base_vertical(yy: float) -> float:
            return float(interpolant(float(yy)))

    else:
        dihedral_tan = float(np.tan(np.deg2rad(base_config.spanwise.dihedral_deg)))
        vertical_offset_m = float(base_config.spanwise.vertical_offset_m)

        def base_vertical(yy: float) -> float:
            return float(vertical_offset_m + dihedral_tan * float(yy))

    def base_te_center(yy: float) -> float:
        local_upper_te, local_lower_te = _cta_te_local_pair(
            y=float(yy),
            planform=planform,
            section_model=section_model,
            twist_deg=float(laws.twist_deg(float(yy))),
        )
        return 0.5 * (float(local_upper_te) + float(local_lower_te)) + float(base_vertical(float(yy)))

    root_te_center = base_te_center(root_y)
    join_te_center = base_te_center(join_y)
    join_te_center_slope = _cta_central_slope(base_te_center, join_y)
    x_start = float(CTA_ROOT_TE_FLATTEN_X_START)

    def transform(
        y: float,
        x_air: np.ndarray,
        upper: np.ndarray,
        lower: np.ndarray,
        _params,
    ) -> tuple[np.ndarray, np.ndarray]:
        yy = float(y)
        if yy <= root_y + 1.0e-12 or yy >= join_y - 1.0e-12:
            return upper, lower
        twist_deg = float(laws.twist_deg(yy))
        chord = float(planform.te_x(yy) - planform.le_x(yy))
        tr = np.deg2rad(twist_deg)
        cos_tr = float(np.cos(tr))
        if abs(chord * cos_tr) <= 1.0e-12:
            return upper, lower
        x_te = float(x_air[-1]) * chord
        current_te_center = (
            x_te * np.sin(tr)
            + 0.5 * (float(upper[-1]) + float(lower[-1])) * chord * cos_tr
            + float(base_vertical(yy))
        )
        target_te_center = quintic_c2_transition(
            yy,
            root_y,
            join_y,
            root_te_center,
            join_te_center,
            0.0,
            join_te_center_slope,
        )
        delta_world = float(target_te_center - current_te_center)
        if abs(delta_world) <= 1.0e-12:
            return upper, lower
        local_delta = delta_world / (chord * cos_tr)
        t = np.clip((np.asarray(x_air, dtype=float) - x_start) / max(1.0 - x_start, 1.0e-12), 0.0, 1.0)
        aft_weight = t * t * (3.0 - 2.0 * t)
        return upper + local_delta * aft_weight, lower + local_delta * aft_weight

    return transform


@lru_cache(maxsize=1)
def _cta_internal_volume_constraint_set():
    return load_internal_volume_constraint_set(
        CTA_INTERNAL_VOLUME_CONSTRAINTS_PATH,
        reference_frame=CTA_CAD_REFERENCE_FRAME,
        name="CTA B359-V0 internal volume constraints",
    )


@lru_cache(maxsize=1)
def _cta_glider_lower_envelope():
    sections = load_xyz_sections_by_span(CTA_GLIDER_GEO_PATH)
    y_sections = np.asarray(sorted(sections.keys()), dtype=float)
    lower = np.asarray([float(np.min(np.asarray(sections[yy], dtype=float)[:, 1])) for yy in y_sections], dtype=float)
    return y_sections, lower


def _cta_internal_volume_profile_transform(
    config: SectionedBWBModelConfig,
):
    constraint_set = _cta_internal_volume_constraint_set()
    planform = build_sectioned_bwb_planform(config.topology, config.planform)
    laws = resolve_spanwise_laws(config)
    smoothing_kernel_x = np.asarray([1.0, 4.0, 6.0, 4.0, 1.0], dtype=float)
    smoothing_kernel_y = np.asarray([1.0, 3.0, 5.0, 3.0, 1.0], dtype=float)
    smoothing_iterations = 24
    fitting_iterations = 6
    fitting_gain = 1.8

    def smooth_along_axis(field: np.ndarray, kernel: np.ndarray, axis: int) -> np.ndarray:
        pad = kernel.size // 2
        padded = np.pad(field, [(pad, pad) if idx == axis else (0, 0) for idx in range(field.ndim)], mode="edge")
        out = np.zeros_like(field, dtype=float)
        for idx, weight in enumerate(kernel):
            slicer = [slice(None)] * field.ndim
            slicer[axis] = slice(idx, idx + field.shape[axis])
            out += float(weight) * padded[tuple(slicer)]
        return out / float(np.sum(kernel))

    def smooth_with_projection_upper(raw: np.ndarray) -> np.ndarray:
        field = np.asarray(raw, dtype=float).copy()
        for _ in range(smoothing_iterations):
            field = smooth_along_axis(field, smoothing_kernel_x, axis=1)
            field = smooth_along_axis(field, smoothing_kernel_y, axis=0)
            field = np.maximum(field, raw)
        return field

    def smooth_with_projection_lower(raw: np.ndarray) -> np.ndarray:
        field = np.asarray(raw, dtype=float).copy()
        for _ in range(smoothing_iterations):
            field = smooth_along_axis(field, smoothing_kernel_x, axis=1)
            field = smooth_along_axis(field, smoothing_kernel_y, axis=0)
            field = np.minimum(field, raw)
        return field

    def compress_dense_to_fit(field_dense: np.ndarray, x_dense: np.ndarray, x_fit: np.ndarray, mode: str) -> np.ndarray:
        compressed = np.zeros((field_dense.shape[0], x_fit.size), dtype=float)
        for idx, x_here in enumerate(x_fit):
            left = 0.0 if idx == 0 else 0.5 * (float(x_fit[idx - 1]) + float(x_here))
            right = 1.0 if idx == x_fit.size - 1 else 0.5 * (float(x_here) + float(x_fit[idx + 1]))
            mask = (x_dense >= left - 1.0e-12) & (x_dense <= right + 1.0e-12)
            if mode == "upper":
                compressed[:, idx] = np.max(field_dense[:, mask], axis=1)
            else:
                compressed[:, idx] = np.min(field_dense[:, mask], axis=1)
        return compressed

    def constraint_surface_x_probes(y_cad: float) -> np.ndarray:
        x_values: list[float] = []
        tol = 1.0e-9
        for surface in constraint_set.surfaces:
            vertices = np.asarray(surface.vertices_xyz_m, dtype=float)
            n_vertices = vertices.shape[0]
            for idx in range(n_vertices):
                p0 = vertices[idx]
                p1 = vertices[(idx + 1) % n_vertices]
                x0, y0 = float(p0[0]), float(p0[1])
                x1, y1 = float(p1[0]), float(p1[1])
                if abs(y_cad - y0) <= tol:
                    x_values.append(x0)
                if abs(y_cad - y1) <= tol:
                    x_values.append(x1)
                if abs(y1 - y0) <= tol:
                    if abs(y_cad - y0) <= tol:
                        x_values.extend((x0, x1))
                    continue
                if (y_cad - y0) * (y_cad - y1) < 0.0:
                    t_edge = (y_cad - y0) / (y1 - y0)
                    x_values.append(x0 + t_edge * (x1 - x0))
        if not x_values:
            return np.empty(0, dtype=float)
        return np.unique(np.asarray(x_values, dtype=float))

    if config.spanwise.vertical_offset_callable is not None:
        base_vertical = config.spanwise.vertical_offset_callable
    elif config.spanwise.vertical_offset_z is not None:
        interpolant = build_anchored_interpolant(config.topology, config.spanwise.vertical_offset_z)

        def base_vertical(yy: float) -> float:
            return float(interpolant(float(yy)))

    else:
        dihedral_tan = float(np.tan(np.deg2rad(config.spanwise.dihedral_deg)))
        vertical_offset_m = float(config.spanwise.vertical_offset_m)

        def base_vertical(yy: float) -> float:
            return float(vertical_offset_m + dihedral_tan * float(yy))

    max_constraint_y = max(float(np.max(surface.vertices_xyz_m[:, 1])) for surface in constraint_set.surfaces)
    min_constraint_y = min(float(np.min(surface.vertices_xyz_m[:, 1])) for surface in constraint_set.surfaces)
    section_model = build_section_model(config, laws)
    x_air_native = np.asarray(section_model.x_air, dtype=float)
    y_grid_model = np.unique(
        np.concatenate(
            [
                np.linspace(0.0, float(config.topology.span), 241, dtype=float),
                np.asarray(
                    [
                        float(vertex[1]) - float(CTA_CAD_REFERENCE_FRAME.offset_y_m)
                        for surface in constraint_set.surfaces
                        for vertex in np.asarray(surface.vertices_xyz_m, dtype=float)
                    ],
                    dtype=float,
                ),
            ]
        )
    )
    y_grid_model = y_grid_model[(y_grid_model >= 0.0) & (y_grid_model <= float(config.topology.span))]
    x_air_probe_list = [np.linspace(0.0, 1.0, 401, dtype=float)]
    for yy in y_grid_model:
        y_cad = float(yy) + float(CTA_CAD_REFERENCE_FRAME.offset_y_m)
        chord = float(planform.te_x(float(yy)) - planform.le_x(float(yy)))
        if chord <= 1.0e-12:
            continue
        le_x = float(planform.le_x(float(yy)))
        twist_deg = float(laws.twist_deg(float(yy)))
        cos_tr = float(np.cos(np.deg2rad(twist_deg)))
        if abs(chord * cos_tr) <= 1.0e-12:
            continue
        x_probe_cad = constraint_surface_x_probes(y_cad)
        if x_probe_cad.size == 0:
            continue
        approx_x_air = np.clip(
            (x_probe_cad - float(CTA_CAD_REFERENCE_FRAME.offset_x_m) - le_x) / (chord * cos_tr),
            0.0,
            1.0,
        )
        x_air_probe_list.append(np.asarray(approx_x_air, dtype=float))
    dense_x_air_grid = np.unique(np.concatenate(x_air_probe_list))
    fit_x_air_grid = np.unique(cosine_spacing(41))
    anchor_y = np.asarray(config.topology.anchor_y_array, dtype=float)
    dense_upper_base = np.zeros((y_grid_model.size, dense_x_air_grid.size), dtype=float)
    dense_lower_base = np.zeros((y_grid_model.size, dense_x_air_grid.size), dtype=float)
    chord_grid = np.zeros(y_grid_model.size, dtype=float)
    le_x_grid = np.zeros(y_grid_model.size, dtype=float)
    vertical_grid = np.zeros(y_grid_model.size, dtype=float)
    cos_grid = np.ones(y_grid_model.size, dtype=float)
    sin_grid = np.zeros(y_grid_model.size, dtype=float)
    active_row = np.zeros(y_grid_model.size, dtype=bool)

    for iy, yy in enumerate(y_grid_model):
        y_cad = float(yy) + float(CTA_CAD_REFERENCE_FRAME.offset_y_m)
        if y_cad < min_constraint_y - 1.0e-9 or y_cad > max_constraint_y + 1.0e-9:
            continue
        chord = float(planform.te_x(float(yy)) - planform.le_x(float(yy)))
        if chord <= 1.0e-12:
            continue
        le_x = float(planform.le_x(float(yy)))
        twist_deg = float(laws.twist_deg(float(yy)))
        vertical = float(base_vertical(float(yy)))
        tr = np.deg2rad(twist_deg)
        cos_tr = float(np.cos(tr))
        sin_tr = float(np.sin(tr))
        if abs(chord * cos_tr) <= 1.0e-12:
            continue
        upper_native, lower_native, _ = section_model.coordinates_at_y(float(yy))
        dense_upper_base[iy] = np.interp(dense_x_air_grid, x_air_native, np.asarray(upper_native, dtype=float))
        dense_lower_base[iy] = np.interp(dense_x_air_grid, x_air_native, np.asarray(lower_native, dtype=float))
        chord_grid[iy] = chord
        le_x_grid[iy] = le_x
        vertical_grid[iy] = vertical
        cos_grid[iy] = cos_tr
        sin_grid[iy] = sin_tr
        active_row[iy] = True

    upper_delta_local_table = np.zeros((y_grid_model.size, fit_x_air_grid.size), dtype=float)
    lower_delta_local_table = np.zeros((y_grid_model.size, fit_x_air_grid.size), dtype=float)

    for _ in range(fitting_iterations):
        upper_delta_local_target_dense = np.zeros((y_grid_model.size, dense_x_air_grid.size), dtype=float)
        lower_delta_local_target_dense = np.zeros((y_grid_model.size, dense_x_air_grid.size), dtype=float)
        for iy, yy in enumerate(y_grid_model):
            if not active_row[iy]:
                continue
            y_cad = float(yy) + float(CTA_CAD_REFERENCE_FRAME.offset_y_m)
            chord = chord_grid[iy]
            le_x = le_x_grid[iy]
            vertical = vertical_grid[iy]
            cos_tr = cos_grid[iy]
            sin_tr = sin_grid[iy]

            upper_delta_dense = np.interp(dense_x_air_grid, fit_x_air_grid, upper_delta_local_table[iy])
            lower_delta_dense = np.interp(dense_x_air_grid, fit_x_air_grid, lower_delta_local_table[iy])
            upper_delta_local_target_dense[iy] = upper_delta_dense
            lower_delta_local_target_dense[iy] = lower_delta_dense
            dense_upper = dense_upper_base[iy] + upper_delta_dense
            dense_lower = dense_lower_base[iy] + lower_delta_dense
            x_local = dense_x_air_grid * chord
            upper_local = dense_upper * chord
            lower_local = dense_lower * chord
            x_upper_cad = le_x + x_local * cos_tr - upper_local * sin_tr + float(CTA_CAD_REFERENCE_FRAME.offset_x_m)
            x_lower_cad = le_x + x_local * cos_tr - lower_local * sin_tr + float(CTA_CAD_REFERENCE_FRAME.offset_x_m)
            y_cad_arr = np.full_like(x_local, y_cad, dtype=float)
            z_upper_cad = vertical + x_local * sin_tr + upper_local * cos_tr + float(CTA_CAD_REFERENCE_FRAME.offset_z_m)
            z_lower_cad = vertical + x_local * sin_tr + lower_local * cos_tr + float(CTA_CAD_REFERENCE_FRAME.offset_z_m)

            upper_required, _ = required_constraint_bounds_at_points(
                constraint_set.surfaces,
                x_upper_cad,
                y_cad_arr,
            )
            _, lower_required = required_constraint_bounds_at_points(
                constraint_set.surfaces,
                x_lower_cad,
                y_cad_arr,
            )

            valid_upper = np.isfinite(upper_required)
            if np.any(valid_upper):
                delta = (upper_required[valid_upper] - z_upper_cad[valid_upper]) / (chord * cos_tr)
                upper_delta_local_target_dense[iy, valid_upper] = np.maximum(
                    upper_delta_local_target_dense[iy, valid_upper],
                    upper_delta_dense[valid_upper] + fitting_gain * np.maximum(delta, 0.0),
                )

            valid_lower = np.isfinite(lower_required)
            if np.any(valid_lower):
                delta = (lower_required[valid_lower] - z_lower_cad[valid_lower]) / (chord * cos_tr)
                lower_delta_local_target_dense[iy, valid_lower] = np.minimum(
                    lower_delta_local_target_dense[iy, valid_lower],
                    lower_delta_dense[valid_lower] + fitting_gain * np.minimum(delta, 0.0),
                )

        upper_delta_local_target = compress_dense_to_fit(
            upper_delta_local_target_dense,
            dense_x_air_grid,
            fit_x_air_grid,
            mode="upper",
        )
        lower_delta_local_target = compress_dense_to_fit(
            lower_delta_local_target_dense,
            dense_x_air_grid,
            fit_x_air_grid,
            mode="lower",
        )
        upper_delta_local_table = smooth_with_projection_upper(upper_delta_local_target)
        lower_delta_local_table = smooth_with_projection_lower(lower_delta_local_target)

    if anchor_y.size >= 3:
        section1_y = float(anchor_y[1])
        section2_y = float(anchor_y[2])
        if section2_y > section1_y + 1.0e-12:
            lower_delta_before_segment_constraint = lower_delta_local_table.copy()
            lower_profile_section1 = np.asarray(
                [np.interp(section1_y, y_grid_model, lower_delta_local_table[:, idx]) for idx in range(fit_x_air_grid.size)],
                dtype=float,
            )
            lower_profile_section2 = np.asarray(
                [np.interp(section2_y, y_grid_model, lower_delta_local_table[:, idx]) for idx in range(fit_x_air_grid.size)],
                dtype=float,
            )
            segment_mask = (y_grid_model >= section1_y) & (y_grid_model < section2_y)
            if np.any(segment_mask):
                lower_profile_mid_target = np.min(lower_delta_before_segment_constraint[segment_mask], axis=0)
                lower_profile_control = 2.0 * lower_profile_mid_target - 0.5 * (
                    lower_profile_section1 + lower_profile_section2
                )
                t = ((y_grid_model[segment_mask] - section1_y) / (section2_y - section1_y)).reshape(-1, 1)
                lower_delta_local_table[segment_mask] = (
                    ((1.0 - t) ** 2) * lower_profile_section1.reshape(1, -1)
                    + (2.0 * (1.0 - t) * t) * lower_profile_control.reshape(1, -1)
                    + (t**2) * lower_profile_section2.reshape(1, -1)
                )

    lower_target_xc = 0.15
    lower_target_indices = tuple(range(min(4, anchor_y.size)))
    lower_target_y = np.asarray([float(anchor_y[idx]) for idx in lower_target_indices], dtype=float)
    lower_target_z = np.zeros(lower_target_y.size, dtype=float)
    for idx, yy in enumerate(lower_target_y):
        chord = float(planform.te_x(float(yy)) - planform.le_x(float(yy)))
        twist_deg = float(laws.twist_deg(float(yy)))
        vertical = float(base_vertical(float(yy)))
        tr = np.deg2rad(twist_deg)
        cos_tr = float(np.cos(tr))
        sin_tr = float(np.sin(tr))
        _, lower_anchor_base, _ = section_model.coordinates_at_y(float(yy))
        lower_z_local = float(np.interp(lower_target_xc, x_air_native, np.asarray(lower_anchor_base, dtype=float))) * chord
        x_local = lower_target_xc * chord
        lower_target_z[idx] = float(
            vertical
            + x_local * sin_tr
            + lower_z_local * cos_tr
            + float(CTA_CAD_REFERENCE_FRAME.offset_z_m)
        )
    lower_front_target = build_scalar_interpolant(lower_target_y, lower_target_z, "pchip")
    lower_front_target_y_end = float(lower_target_y[-1])

    lower_raise_raw_m = np.zeros(y_grid_model.size, dtype=float)
    for iy, yy in enumerate(y_grid_model):
        if not active_row[iy]:
            continue
        y_cad = float(yy) + float(CTA_CAD_REFERENCE_FRAME.offset_y_m)
        chord = chord_grid[iy]
        le_x = le_x_grid[iy]
        vertical = vertical_grid[iy]
        cos_tr = cos_grid[iy]
        sin_tr = sin_grid[iy]
        lower_delta_dense = np.interp(dense_x_air_grid, fit_x_air_grid, lower_delta_local_table[iy])
        dense_lower = dense_lower_base[iy] + lower_delta_dense
        x_local = dense_x_air_grid * chord
        lower_local = dense_lower * chord
        x_lower_cad = le_x + x_local * cos_tr - lower_local * sin_tr + float(CTA_CAD_REFERENCE_FRAME.offset_x_m)
        z_lower_cad = vertical + x_local * sin_tr + lower_local * cos_tr + float(CTA_CAD_REFERENCE_FRAME.offset_z_m)
        y_cad_arr = np.full_like(x_local, y_cad, dtype=float)
        _, lower_required = required_constraint_bounds_at_points(
            constraint_set.surfaces,
            x_lower_cad,
            y_cad_arr,
        )
        valid_lower = np.isfinite(lower_required)
        allowable_raise_m = np.inf
        if np.any(valid_lower):
            allowable_raise_m = max(0.0, float(np.min(z_lower_cad[valid_lower] - lower_required[valid_lower])))
        if float(yy) <= lower_front_target_y_end + 1.0e-12:
            target_lower_m = float(lower_front_target(float(yy)))
            current_lower_env_m = float(np.min(z_lower_cad))
            desired_raise_m = max(0.0, target_lower_m - current_lower_env_m)
        else:
            desired_raise_m = 0.0
        lower_raise_raw_m[iy] = min(desired_raise_m, allowable_raise_m)

    lower_raise_smooth_m = smooth_along_axis(lower_raise_raw_m[:, None], smoothing_kernel_y, axis=0).ravel()
    lower_raise_smooth_m = np.maximum(lower_raise_smooth_m, 0.94 * lower_raise_raw_m)
    lower_raise_smooth_m = np.minimum(lower_raise_smooth_m, lower_raise_raw_m)
    lower_raise_smooth_m = np.maximum(lower_raise_smooth_m, 0.0)

    for iy in range(y_grid_model.size):
        if not active_row[iy]:
            continue
        local_raise = lower_raise_smooth_m[iy] / max(chord_grid[iy] * cos_grid[iy], 1.0e-12)
        lower_delta_local_table[iy] += local_raise

    def transform(
        y: float,
        x_air: np.ndarray,
        upper: np.ndarray,
        lower: np.ndarray,
        _params,
    ) -> tuple[np.ndarray, np.ndarray]:
        yy = float(y)
        y_cad = yy + float(CTA_CAD_REFERENCE_FRAME.offset_y_m)
        if y_cad < min_constraint_y - 1.0e-9 or y_cad > max_constraint_y + 1.0e-9:
            return upper, lower
        x_air_arr = np.asarray(x_air, dtype=float)
        upper_delta_vs_x = np.asarray(
            [np.interp(yy, y_grid_model, upper_delta_local_table[:, idx]) for idx in range(fit_x_air_grid.size)],
            dtype=float,
        )
        lower_delta_vs_x = np.asarray(
            [np.interp(yy, y_grid_model, lower_delta_local_table[:, idx]) for idx in range(fit_x_air_grid.size)],
            dtype=float,
        )
        upper_delta = np.interp(x_air_arr, fit_x_air_grid, upper_delta_vs_x)
        lower_delta = np.interp(x_air_arr, fit_x_air_grid, lower_delta_vs_x)
        return (
            np.asarray(upper, dtype=float) + upper_delta,
            np.asarray(lower, dtype=float) + lower_delta,
        )

    return transform


def _compose_profile_post_transforms(*transforms):
    active = tuple(transform for transform in transforms if transform is not None)
    if not active:
        return None

    def composed(
        y: float,
        x_air: np.ndarray,
        upper: np.ndarray,
        lower: np.ndarray,
        params,
    ) -> tuple[np.ndarray, np.ndarray]:
        out_upper = np.asarray(upper, dtype=float).copy()
        out_lower = np.asarray(lower, dtype=float).copy()
        for transform in active:
            out_upper, out_lower = transform(y, x_air, out_upper, out_lower, params)
        return out_upper, out_lower

    return composed


def _cta_model_anchor_y_from_topology(config: SectionedBWBModelConfig) -> tuple[float, ...]:
    y_sections = np.asarray(config.topology.y_sections_array, dtype=float)
    return (
        float(CTA_PROFILE_ANCHOR_Y[0]),
        float(CTA_PROFILE_ANCHOR_Y[1]),
        float(CTA_PROFILE_ANCHOR_Y[2]),
        float(y_sections[1]),
        float(y_sections[2]),
        float(y_sections[3]),
    )


def _cta_joined_section3_exact_override(
    config: SectionedBWBModelConfig,
) -> ExactSectionProfileOverrideSpec:
    anchor_specs = _cta_anchor_section_specs_from_active(tuple(config.sections.section_specs))
    source_spec = anchor_specs[3]
    x_air = cosine_spacing(int(_CTA_TEMPLATE.num_airfoil_points))
    shape = KulfanCSTAirfoil(
        degree=len(source_spec.upper_coeffs) - 1,
        n1=float(_CTA_TEMPLATE.cst_n1),
        n2=float(_CTA_TEMPLATE.cst_n2),
        x_tc_window=(0.15, 0.65),
        shared_leading_edge=bool(_CTA_TEMPLATE.shared_leading_edge),
    )
    source_coeffs = np.concatenate([source_spec.upper_coeffs, source_spec.lower_coeffs], dtype=float)
    joined_upper, joined_lower = shape.evaluate(
        x_air,
        source_coeffs,
        te_thickness=float(source_spec.te_thickness),
    )
    thickness = joined_upper - joined_lower

    return ExactSectionProfileOverrideSpec(
        section_index=3,
        x=tuple(float(value) for value in x_air),
        upper=tuple(float(value) for value in joined_upper),
        lower=tuple(float(value) for value in joined_lower),
        tc_max=float(np.max(thickness)),
        x_tmax=float(x_air[int(np.argmax(thickness))]),
        te_thickness=float(thickness[-1]),
    )


def _cta_outward_section1_blend_override(
    config: SectionedBWBModelConfig,
) -> ExactSectionProfileOverrideSpec:
    anchor_specs = _cta_anchor_section_specs_from_active(tuple(config.sections.section_specs))
    source_spec = anchor_specs[1]
    x_air = cosine_spacing(int(_CTA_TEMPLATE.num_airfoil_points))
    shape = KulfanCSTAirfoil(
        degree=len(source_spec.upper_coeffs) - 1,
        n1=float(_CTA_TEMPLATE.cst_n1),
        n2=float(_CTA_TEMPLATE.cst_n2),
        x_tc_window=(0.15, 0.65),
        shared_leading_edge=bool(_CTA_TEMPLATE.shared_leading_edge),
    )
    source_coeffs = np.concatenate([source_spec.upper_coeffs, source_spec.lower_coeffs], dtype=float)
    source_upper, source_lower = shape.evaluate(
        x_air,
        source_coeffs,
        te_thickness=float(source_spec.te_thickness),
    )
    override_upper = source_upper
    override_lower = source_lower
    thickness = override_upper - override_lower
    return ExactSectionProfileOverrideSpec(
        section_index=1,
        x=tuple(float(value) for value in x_air),
        upper=tuple(float(value) for value in override_upper),
        lower=tuple(float(value) for value in override_lower),
        tc_max=float(np.max(thickness)),
        x_tmax=float(x_air[int(np.argmax(thickness))]),
        te_thickness=float(thickness[-1]),
        y_end=float(config.topology.anchor_y_array[2]),
        blend_to_base=True,
    )


def _cta_root_to_section1_blend_override(
    config: SectionedBWBModelConfig,
) -> ExactSectionProfileOverrideSpec:
    anchor_specs = _cta_anchor_section_specs_from_active(tuple(config.sections.section_specs))
    source_spec = anchor_specs[0]
    x_air = cosine_spacing(int(_CTA_TEMPLATE.num_airfoil_points))
    shape = KulfanCSTAirfoil(
        degree=len(source_spec.upper_coeffs) - 1,
        n1=float(_CTA_TEMPLATE.cst_n1),
        n2=float(_CTA_TEMPLATE.cst_n2),
        x_tc_window=(0.15, 0.65),
        shared_leading_edge=bool(_CTA_TEMPLATE.shared_leading_edge),
    )
    source_coeffs = np.concatenate([source_spec.upper_coeffs, source_spec.lower_coeffs], dtype=float)
    source_upper, source_lower = shape.evaluate(
        x_air,
        source_coeffs,
        te_thickness=float(source_spec.te_thickness),
    )
    thickness = source_upper - source_lower
    return ExactSectionProfileOverrideSpec(
        section_index=0,
        x=tuple(float(value) for value in x_air),
        upper=tuple(float(value) for value in source_upper),
        lower=tuple(float(value) for value in source_lower),
        tc_max=float(np.max(thickness)),
        x_tmax=float(x_air[int(np.argmax(thickness))]),
        te_thickness=float(thickness[-1]),
        y_end=float(config.topology.anchor_y_array[1]),
        blend_to_base=True,
        blend_curve="ellipse",
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
        le_root_x=CTA_LE_ROOT_X_M,
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
    model_anchor_y = _cta_model_anchor_y_from_topology(config)
    model_anchor_le_z_m = CTA_PROFILE_ANCHOR_LE_Z_M
    model_anchor_twist_deg = CTA_PROFILE_ANCHOR_TWIST_DEG
    active_section_specs = tuple(config.sections.section_specs)
    config.topology.anchor_y = model_anchor_y
    if use_cta_anchor_twist:
        twist_values = model_anchor_twist_deg
    else:
        twist_values = (float(fixed_design.twist_c1_deg),) * (len(model_anchor_y) - 2) + (
            float(fixed_design.twist_c3_deg),
            float(fixed_design.twist_c4_deg),
        )
    config.sections = replace(
        config.sections,
        shared_leading_edge=False,
        section_specs_override=_cta_anchor_section_specs_from_active(active_section_specs),
        profile_relations=tuple(SectionProfileRelationSpec() for _ in range(len(model_anchor_y))),
        shape_hold_segments=(2,),
        shape_hold_tc_ramp_segments=(),
        exact_profile_overrides=(),
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
    config.planform.le_linear_start_index = 4
    config.planform.te_linear_start_index = 4
    config.planform.body_le_fixed_points = (
        (_cta_le_helper_1_x_m(), CTA_NOSE_HELPER_1_Y_M),
        (_cta_le_helper_2_x_m(), CTA_C1_Y_M),
    )
    config.planform.symmetry_blend_y = CTA_NOSE_HELPER_1_Y_M
    config.sampling.airfoil_distribution_mode = "anchors"
    config.sampling.section_interpolation = "pchip"
    config.sampling.section_linear_start_index = 4
    config.spanwise.twist_deg = AnchoredSpanwiseLaw(
        section_indices=tuple(range(len(model_anchor_y))),
        values=twist_values,
        interpolation="pchip",
        linear_start_index=4,
        root_blend_y=float(CTA_ROOT_FRONT_BLEND_FRACTION * model_anchor_y[1]),
    )
    config.spanwise.camber_delta = AnchoredSpanwiseLaw(
        section_indices=tuple(range(len(model_anchor_y))),
        values=(float(fixed_design.camber_c1),) * (len(model_anchor_y) - 2)
        + (float(fixed_design.camber_c3), float(fixed_design.camber_c4)),
        interpolation="pchip",
        linear_start_index=4,
    )
    config.spanwise.vertical_offset_z = AnchoredSpanwiseLaw(
        section_indices=tuple(range(len(model_anchor_y))),
        values=model_anchor_le_z_m,
        interpolation="pchip",
        linear_start_index=4,
        root_blend_y=float(CTA_ROOT_FRONT_BLEND_FRACTION * model_anchor_y[1]),
    )
    (
        config.spanwise.twist_callable,
        config.spanwise.vertical_offset_callable,
    ) = _cta_root_front_matched_callables(config)
    config.sections = replace(
        config.sections,
        profile_post_transform=None,
    )
    return config
