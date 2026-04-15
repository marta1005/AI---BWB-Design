from __future__ import annotations

from dataclasses import dataclass, replace
from math import atan2, degrees, radians, tan
from typing import Callable

from .specs import (
    AnchoredSpanwiseLaw,
    ExportSpec,
    PlanformSpec,
    SamplingSpec,
    SectionCSTSpec,
    SectionFamilySpec,
    SectionProfileRelationSpec,
    SectionedBWBModelConfig,
    SpanwiseLawSpec,
    VolumeConstraintSpec,
)
from .topology import SectionedBWBTopologySpec


@dataclass(frozen=True)
class Relation:
    target: str
    requires: tuple[str, ...]
    compute: Callable[..., float]
    description: str


@dataclass(frozen=True)
class CanonicalSectionDeclaration:
    y_0_m: float
    y_1_m: float
    y_2_m: float
    y_3_m: float
    le_x_0_m: float
    le_x_1_m: float
    le_x_2_m: float
    le_x_3_m: float
    chord_0_m: float
    chord_1_m: float
    chord_2_m: float
    chord_3_m: float


@dataclass(frozen=True)
class ExplicitSectionSpec:
    name: str
    y_m: float
    le_x_m: float
    chord_m: float


@dataclass(frozen=True)
class CaseTemplate:
    anchor_y_m: tuple[float, ...]
    anchor_le_z_m: tuple[float, ...]
    anchor_twist_deg: tuple[float, ...]
    anchor_camber_delta: tuple[float, ...]
    anchor_sections: tuple[SectionCSTSpec, ...]
    cst_n1: float
    cst_n2: float
    shared_leading_edge: bool
    profile_generation_mode: str
    body_le_fixed_points: tuple[tuple[float, float], ...]
    te_c1_y_m: float
    te_inboard_blend_y_m: float
    te_exact_segments: tuple[int, ...]
    te_inboard_blend_dx_m: float
    te_inboard_radius_factor: float
    med_3_te_helper_fraction: float
    te_outer_blend_fraction: float
    blend_fraction: float
    min_linear_core_fraction: float
    te_blend_fraction: float | None
    te_min_linear_core_fraction: float | None
    te_spline_bridge: tuple[int, int] | None
    symmetry_blend_y_m: float
    section_interpolation: str
    airfoil_distribution_mode: str
    num_airfoil_points: int
    num_base_stations: int
    section_curve_n_ctl: int
    k_span: int
    twist_interpolation: str = "pchip"
    camber_interpolation: str = "pyspline"
    vertical_interpolation: str = "pchip"


def leading_edge_sweep_from_chord_sweep(
    chord_sweep_deg: float,
    chord_fraction: float,
    dy_m: float,
    chord_in_m: float,
    chord_out_m: float,
) -> float:
    dy = float(dy_m)
    if dy <= 0.0:
        raise ValueError(f"Spanwise segment length must be positive, got {dy:.6f} m")
    tan_le = tan(radians(float(chord_sweep_deg))) - float(chord_fraction) * (
        float(chord_out_m) - float(chord_in_m)
    ) / dy
    return float(degrees(atan2(tan_le * dy, dy)))


def chord_sweep_from_leading_edge_sweep(
    le_sweep_deg: float,
    chord_fraction: float,
    dy_m: float,
    chord_in_m: float,
    chord_out_m: float,
) -> float:
    dy = float(dy_m)
    if dy <= 0.0:
        raise ValueError(f"Spanwise segment length must be positive, got {dy:.6f} m")
    tan_chord = tan(radians(float(le_sweep_deg))) + float(chord_fraction) * (
        float(chord_out_m) - float(chord_in_m)
    ) / dy
    return float(degrees(atan2(tan_chord * dy, dy)))


def solve_relations(initial_values: dict[str, float], relations: list[Relation]) -> dict[str, float]:
    values = dict(initial_values)
    pending = list(relations)
    while pending:
        progressed = False
        next_pending: list[Relation] = []
        for relation in pending:
            if relation.target in values:
                progressed = True
                continue
            if all(name in values for name in relation.requires):
                args = [values[name] for name in relation.requires]
                values[relation.target] = float(relation.compute(*args))
                progressed = True
            else:
                next_pending.append(relation)
        if not progressed:
            missing = {
                relation.target: tuple(name for name in relation.requires if name not in values)
                for relation in next_pending
            }
            raise ValueError(f"Could not resolve relations, missing dependencies: {missing}")
        pending = next_pending
    return values


def build_case_config_from_resolved(
    resolved: dict[str, float],
    template: CaseTemplate,
    anchor_y: tuple[float, ...],
    topology_name: str = "bwb_case",
) -> SectionedBWBModelConfig:
    topology = SectionedBWBTopologySpec(
        span=float(resolved["span_m"]),
        b1_span_ratio=float(resolved["b1_span_ratio"]),
        b2_span_ratio=float(resolved["b2_span_ratio"]),
        b3_span_ratio=float(resolved["b3_span_ratio"]),
        anchor_y=anchor_y,
        topology_name=topology_name,
    )
    planform = PlanformSpec(
        le_root_x=float(resolved["le_x_0_m"]),
        c1_root_chord=float(resolved["chord_0_m"]),
        c2_c1_ratio=float(resolved["c2_c1_ratio"]),
        c3_c1_ratio=float(resolved["c3_c1_ratio"]),
        c4_c1_ratio=float(resolved["c4_c1_ratio"]),
        s1_deg=float(resolved["s1_deg"]),
        s2_deg=float(resolved["s2_deg"]),
        s3_deg=float(resolved["s3_deg"]),
        med_3_te_sweep_deg=float(resolved.get("med_3_te_sweep_deg", 0.0)),
        body_le_fixed_points=template.body_le_fixed_points,
        te_exact_segments=template.te_exact_segments,
        te_c1_span_fraction=float(resolved["te_c1_span_fraction"]),
        te_inboard_blend_fraction=float(resolved["te_inboard_blend_fraction"]),
        te_inboard_blend_dx=template.te_inboard_blend_dx_m,
        te_inboard_radius_factor=template.te_inboard_radius_factor,
        med_3_te_helper_fraction=template.med_3_te_helper_fraction,
        te_outer_blend_fraction=template.te_outer_blend_fraction,
        blend_fraction=template.blend_fraction,
        min_linear_core_fraction=template.min_linear_core_fraction,
        te_blend_fraction=template.te_blend_fraction,
        te_min_linear_core_fraction=template.te_min_linear_core_fraction,
        te_spline_bridge=template.te_spline_bridge,
        symmetry_blend_y=template.symmetry_blend_y_m,
    )
    sections = SectionFamilySpec(
        cst_degree=len(template.anchor_sections[0].upper_coeffs) - 1,
        n1=template.cst_n1,
        n2=template.cst_n2,
        shared_leading_edge=template.shared_leading_edge,
        profile_generation_mode=template.profile_generation_mode,
        section_specs_override=tuple(replace(spec) for spec in template.anchor_sections),
        profile_relations=tuple(SectionProfileRelationSpec() for _ in range(len(template.anchor_sections))),
    )
    spanwise = SpanwiseLawSpec(
        twist_deg=AnchoredSpanwiseLaw(
            section_indices=tuple(range(len(anchor_y))),
            values=template.anchor_twist_deg,
            interpolation=template.twist_interpolation,
        ),
        camber_delta=AnchoredSpanwiseLaw(
            section_indices=tuple(range(len(anchor_y))),
            values=template.anchor_camber_delta,
            interpolation=template.camber_interpolation,
        ),
        vertical_offset_z=AnchoredSpanwiseLaw(
            section_indices=tuple(range(len(anchor_y))),
            values=template.anchor_le_z_m,
            interpolation=template.vertical_interpolation,
        ),
    )
    sampling = SamplingSpec(
        num_airfoil_points=template.num_airfoil_points,
        num_base_stations=template.num_base_stations,
        section_curve_n_ctl=template.section_curve_n_ctl,
        section_interpolation=template.section_interpolation,
        airfoil_distribution_mode=template.airfoil_distribution_mode,
        k_span=template.k_span,
    )
    return SectionedBWBModelConfig(
        topology=topology,
        planform=planform,
        sections=sections,
        spanwise=spanwise,
        sampling=sampling,
        export=ExportSpec(),
        volume=VolumeConstraintSpec(),
    )


def build_case_config_from_explicit_sections(
    sections_definition: tuple[ExplicitSectionSpec, ...],
    template: CaseTemplate,
    topology_name: str = "bwb_case",
    leading_edge_control_points: tuple[tuple[float, float], ...] | None = None,
    trailing_edge_control_points: tuple[tuple[float, float], ...] | None = None,
    le_exact_segments: tuple[int, ...] = (),
    te_exact_segments: tuple[int, ...] | None = None,
    le_spline_bridge: tuple[int, int] | None = None,
    te_spline_bridge: tuple[int, int] | None = None,
    le_linear_start_index: int | None = None,
) -> SectionedBWBModelConfig:
    if len(sections_definition) < 2:
        raise ValueError("sections_definition must contain at least two sections")

    section_y = tuple(float(section.y_m) for section in sections_definition)
    section_le_x = tuple(float(section.le_x_m) for section in sections_definition)
    section_chords = tuple(float(section.chord_m) for section in sections_definition)
    span = float(section_y[-1])
    anchor_y = tuple(float(value) for value in template.anchor_y_m)

    topology = SectionedBWBTopologySpec(
        span=span,
        section_y=section_y,
        anchor_y=anchor_y,
        topology_name=topology_name,
    )
    planform = PlanformSpec(
        le_root_x=section_le_x[0],
        c1_root_chord=section_chords[0],
        section_le_x=section_le_x,
        section_chords_override=section_chords,
        leading_edge_control_points=leading_edge_control_points,
        trailing_edge_control_points=trailing_edge_control_points,
        le_exact_segments=tuple(int(value) for value in le_exact_segments),
        te_exact_segments=(
            ()
            if te_exact_segments is None
            else tuple(int(value) for value in te_exact_segments)
        ),
        le_spline_bridge=le_spline_bridge,
        te_spline_bridge=te_spline_bridge,
        le_linear_start_index=le_linear_start_index,
        body_le_fixed_points=(),
        te_c1_span_fraction=0.5,
        te_inboard_blend_fraction=0.75,
        te_inboard_blend_dx=template.te_inboard_blend_dx_m,
        te_inboard_radius_factor=template.te_inboard_radius_factor,
        med_3_te_helper_fraction=template.med_3_te_helper_fraction,
        te_outer_blend_fraction=template.te_outer_blend_fraction,
        blend_fraction=template.blend_fraction,
        min_linear_core_fraction=template.min_linear_core_fraction,
        te_blend_fraction=template.te_blend_fraction,
        te_min_linear_core_fraction=template.te_min_linear_core_fraction,
        symmetry_blend_y=template.symmetry_blend_y_m,
    )
    sections = SectionFamilySpec(
        cst_degree=len(template.anchor_sections[0].upper_coeffs) - 1,
        n1=template.cst_n1,
        n2=template.cst_n2,
        shared_leading_edge=template.shared_leading_edge,
        profile_generation_mode=template.profile_generation_mode,
        section_specs_override=tuple(replace(spec) for spec in template.anchor_sections),
        profile_relations=tuple(SectionProfileRelationSpec() for _ in range(len(template.anchor_sections))),
    )
    spanwise = SpanwiseLawSpec(
        twist_deg=AnchoredSpanwiseLaw(
            section_indices=tuple(range(len(anchor_y))),
            values=template.anchor_twist_deg,
            interpolation=template.twist_interpolation,
        ),
        camber_delta=AnchoredSpanwiseLaw(
            section_indices=tuple(range(len(anchor_y))),
            values=template.anchor_camber_delta,
            interpolation=template.camber_interpolation,
        ),
        vertical_offset_z=AnchoredSpanwiseLaw(
            section_indices=tuple(range(len(anchor_y))),
            values=template.anchor_le_z_m,
            interpolation=template.vertical_interpolation,
        ),
    )
    sampling = SamplingSpec(
        num_airfoil_points=template.num_airfoil_points,
        num_base_stations=template.num_base_stations,
        section_curve_n_ctl=template.section_curve_n_ctl,
        section_interpolation=template.section_interpolation,
        airfoil_distribution_mode=template.airfoil_distribution_mode,
        k_span=template.k_span,
    )
    return SectionedBWBModelConfig(
        topology=topology,
        planform=planform,
        sections=sections,
        spanwise=spanwise,
        sampling=sampling,
        export=ExportSpec(),
        volume=VolumeConstraintSpec(),
    )
