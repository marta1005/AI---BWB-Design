from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from parametrization.CTA.reference import build_reference_design, to_cta_model_config
from parametrization.aircraft import (
    ContinuityOrder,
    IntuitiveAirfoilProfileSpec,
    InterpolationMethod,
    InterpolationSpec,
    LiftingSurfaceBuildOptions,
    ProfileCatalog,
    ScalarLawAnchor,
    ScalarLawSpec,
    WingSpec,
    WingStationSpec,
)
from parametrization.bwb.planform import build_sectioned_bwb_planform


@dataclass(frozen=True)
class CTAPlanformReference:
    span: float
    y_sections: np.ndarray
    y_samples: np.ndarray
    y_dense: np.ndarray
    eta_samples: np.ndarray
    le_x_samples: np.ndarray
    te_x_samples: np.ndarray
    chord_samples: np.ndarray
    le_x_dense: np.ndarray
    te_x_dense: np.ndarray
    section_etas: np.ndarray
    section_le_x: np.ndarray
    section_chords: np.ndarray
    section_twists_deg: np.ndarray
    continuity_order: int
    blend_fraction: float
    max_le_polyline_error_m: float
    max_te_polyline_error_m: float


def default_output_dir() -> Path:
    return SCRIPT_DIR.parent / "example_outputs" / "cta_planform_match"


def _scalar_law(
    name: str,
    etas: np.ndarray,
    values: np.ndarray,
    interpolation: InterpolationSpec,
) -> ScalarLawSpec:
    anchors = tuple(
        ScalarLawAnchor(eta=float(eta), value=float(value))
        for eta, value in zip(np.asarray(etas, dtype=float), np.asarray(values, dtype=float))
    )
    return ScalarLawSpec(name=name, anchors=anchors, interpolation=interpolation)


def build_cta_planform_reference() -> CTAPlanformReference:
    design = build_reference_design()
    config = to_cta_model_config(design)
    planform = build_sectioned_bwb_planform(config.topology, config.planform)

    span = float(config.topology.span)
    y_sections = np.asarray(config.topology.y_sections_array, dtype=float)
    te_points = np.asarray(config.planform.trailing_edge_points(config.topology), dtype=float)

    # Non-uniform spanwise sampling:
    # dense inboard and around CTA TE auxiliary points, lighter outboard.
    y_samples = np.unique(
        np.concatenate(
            [
                np.linspace(0.0, float(y_sections[1]), 101),
                np.linspace(float(y_sections[1]), float(y_sections[2]), 61),
                np.linspace(float(y_sections[2]), span, 101),
                y_sections,
                te_points[:, 1],
                np.asarray([0.0, float(config.planform.symmetry_blend_y), span], dtype=float),
            ]
        )
    )
    eta_samples = y_samples / span
    le_x_samples = np.asarray([planform.le_x(float(y)) for y in y_samples], dtype=float)
    te_x_samples = np.asarray([planform.te_x(float(y)) for y in y_samples], dtype=float)
    chord_samples = te_x_samples - le_x_samples

    dense_y = np.linspace(0.0, span, 8001)
    dense_le = np.asarray([planform.le_x(float(y)) for y in dense_y], dtype=float)
    dense_te = np.asarray([planform.te_x(float(y)) for y in dense_y], dtype=float)
    dense_le_poly = np.interp(dense_y, y_samples, le_x_samples)
    dense_te_poly = np.interp(dense_y, y_samples, te_x_samples)

    section_le_x = np.asarray([planform.le_x(float(y)) for y in y_sections], dtype=float)
    section_te_x = np.asarray([planform.te_x(float(y)) for y in y_sections], dtype=float)
    section_chords = section_te_x - section_le_x
    section_twists_deg = np.asarray(
        (
            float(design.twist_c1_deg),
            float(design.twist_c2_deg),
            float(design.twist_c3_deg),
            float(design.twist_c4_deg),
        ),
        dtype=float,
    )

    return CTAPlanformReference(
        span=span,
        y_sections=y_sections,
        y_samples=y_samples,
        y_dense=dense_y,
        eta_samples=eta_samples,
        le_x_samples=le_x_samples,
        te_x_samples=te_x_samples,
        chord_samples=chord_samples,
        le_x_dense=dense_le,
        te_x_dense=dense_te,
        section_etas=y_sections / span,
        section_le_x=section_le_x,
        section_chords=section_chords,
        section_twists_deg=section_twists_deg,
        continuity_order=int(config.planform.continuity_order),
        blend_fraction=float(config.planform.blend_fraction),
        max_le_polyline_error_m=float(np.max(np.abs(dense_le - dense_le_poly))),
        max_te_polyline_error_m=float(np.max(np.abs(dense_te - dense_te_poly))),
    )


def build_cta_planform_profiles() -> ProfileCatalog:
    return ProfileCatalog(
        (
            IntuitiveAirfoilProfileSpec(
                profile_id="cta_like_c1",
                leading_edge_radius=0.030,
                max_thickness=0.195,
                x_tmax=0.34,
                max_camber=0.012,
                x_cmax=0.43,
                trailing_edge_wedge_angle_deg=9.5,
                trailing_edge_camber_angle_deg=-0.4,
                aft_control_x=0.73,
                te_thickness=0.0018,
            ),
            IntuitiveAirfoilProfileSpec(
                profile_id="cta_like_c2",
                leading_edge_radius=0.024,
                max_thickness=0.175,
                x_tmax=0.33,
                max_camber=0.010,
                x_cmax=0.41,
                trailing_edge_wedge_angle_deg=9.0,
                trailing_edge_camber_angle_deg=-0.3,
                aft_control_x=0.72,
                te_thickness=0.0016,
            ),
            IntuitiveAirfoilProfileSpec(
                profile_id="cta_like_c3",
                leading_edge_radius=0.016,
                max_thickness=0.130,
                x_tmax=0.31,
                max_camber=0.007,
                x_cmax=0.38,
                trailing_edge_wedge_angle_deg=8.0,
                trailing_edge_camber_angle_deg=-0.2,
                aft_control_x=0.70,
                te_thickness=0.0012,
            ),
            IntuitiveAirfoilProfileSpec(
                profile_id="cta_like_c4",
                leading_edge_radius=0.008,
                max_thickness=0.090,
                x_tmax=0.28,
                max_camber=0.004,
                x_cmax=0.34,
                trailing_edge_wedge_angle_deg=7.0,
                trailing_edge_camber_angle_deg=-0.1,
                aft_control_x=0.67,
                te_thickness=0.0008,
            ),
        )
    )


def build_cta_planform_wing(
    reference: CTAPlanformReference | None = None,
) -> WingSpec:
    ref = build_cta_planform_reference() if reference is None else reference

    profile_ids = ("cta_like_c1", "cta_like_c2", "cta_like_c3", "cta_like_c4")
    station_ids = ("c1_root", "c2_break", "c3_break", "c4_tip")
    section_interp = InterpolationSpec(
        method=InterpolationMethod.SEGMENTED,
        continuity=ContinuityOrder.C2,
        blend_fraction=0.16,
    )
    twist_interp = InterpolationSpec(
        method=InterpolationMethod.SEGMENTED,
        continuity=ContinuityOrder.C2,
        blend_fraction=float(ref.blend_fraction),
    )
    planform_interp = InterpolationSpec(method=InterpolationMethod.LINEAR)

    stations = tuple(
        WingStationSpec(
            station_id=station_id,
            eta=float(eta),
            profile_id=profile_id,
            chord=float(chord),
            twist_deg=float(twist),
            x_le=float(x_le) if index > 0 else None,
        )
        for index, (station_id, eta, profile_id, chord, twist, x_le) in enumerate(
            zip(
                station_ids,
                ref.section_etas,
                profile_ids,
                ref.section_chords,
                ref.section_twists_deg,
                ref.section_le_x,
            )
        )
    )

    return WingSpec(
        wing_id="cta_planform_match",
        semispan=float(ref.span),
        stations=stations,
        root_x=float(ref.section_le_x[0]),
        root_y=0.0,
        root_z=0.0,
        section_interpolation=section_interp,
        spine_interpolation=InterpolationSpec(
            method=InterpolationMethod.SEGMENTED,
            continuity=ContinuityOrder.C2,
            blend_fraction=0.16,
        ),
        scalar_laws=(
            _scalar_law("x", ref.eta_samples, ref.le_x_samples, planform_interp),
            _scalar_law("chord", ref.eta_samples, ref.chord_samples, planform_interp),
            _scalar_law("twist_deg", ref.section_etas, ref.section_twists_deg, twist_interp),
        ),
        metadata={
            "description": "Aircraft-wing example matching the CTA reference planform",
            "reference_source": "parametrization.CTA reference planform",
            "planform_match_mode": "dense_sampled_scalar_laws",
        },
    )


def default_build_options(
    reference: CTAPlanformReference | None = None,
) -> LiftingSurfaceBuildOptions:
    ref = build_cta_planform_reference() if reference is None else reference
    return LiftingSurfaceBuildOptions(
        station_count=max(2, int(ref.eta_samples.size)),
        station_etas=tuple(float(value) for value in ref.eta_samples),
        include_anchor_sections=True,
        airfoil_sample_count=301,
        blunt_te=True,
        tip_style="rounded",
    )


def export_build_options(
    reference: CTAPlanformReference | None = None,
) -> LiftingSurfaceBuildOptions:
    ref = build_cta_planform_reference() if reference is None else reference
    critical_y = np.asarray([0.0, 2.5, 5.0, 7.0, 8.0, 12.0, float(ref.span)], dtype=float)
    export_y = np.unique(np.concatenate([ref.y_samples[::3], critical_y]))
    export_etas = export_y / float(ref.span)
    return LiftingSurfaceBuildOptions(
        station_count=max(2, int(export_etas.size)),
        station_etas=tuple(float(value) for value in export_etas),
        include_anchor_sections=True,
        airfoil_sample_count=181,
        k_span=2,
        blunt_te=True,
        tip_style="rounded",
    )
