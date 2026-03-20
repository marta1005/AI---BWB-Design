from __future__ import annotations

from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from parametrization.aircraft import (
    ContinuityOrder,
    InterpolationMethod,
    IntuitiveAirfoilProfileSpec,
    LiftingSurfaceBuildOptions,
    ProfileCatalog,
    ScalarLawAnchor,
    ScalarLawSpec,
    WingDesignerSpec,
    WingSpec,
    WingStationDesignerSpec,
    WingTransitionSpec,
)


def default_output_dir() -> Path:
    return SCRIPT_DIR.parent / "example_outputs" / "demo_main_wing"


def build_demo_profiles() -> ProfileCatalog:
    return ProfileCatalog(
        (
            IntuitiveAirfoilProfileSpec(
                profile_id="root_airfoil",
                leading_edge_radius=0.018,
                max_thickness=0.145,
                x_tmax=0.34,
                max_camber=0.024,
                x_cmax=0.44,
                trailing_edge_wedge_angle_deg=13.5,
                trailing_edge_camber_angle_deg=-0.8,
                aft_control_x=0.74,
                te_thickness=0.0018,
            ),
            IntuitiveAirfoilProfileSpec(
                profile_id="kink_airfoil",
                leading_edge_radius=0.013,
                max_thickness=0.118,
                x_tmax=0.32,
                max_camber=0.016,
                x_cmax=0.41,
                trailing_edge_wedge_angle_deg=11.0,
                trailing_edge_camber_angle_deg=-0.6,
                aft_control_x=0.72,
                te_thickness=0.0014,
            ),
            IntuitiveAirfoilProfileSpec(
                profile_id="tip_airfoil",
                leading_edge_radius=0.008,
                max_thickness=0.090,
                x_tmax=0.29,
                max_camber=0.008,
                x_cmax=0.36,
                trailing_edge_wedge_angle_deg=9.0,
                trailing_edge_camber_angle_deg=-0.3,
                aft_control_x=0.69,
                te_thickness=0.0009,
            ),
        )
    )


def build_demo_wing(
    root_x: float = 0.0,
    root_y: float = 0.0,
    root_z: float = 0.0,
    transition_continuity: ContinuityOrder = ContinuityOrder.C2,
    transition_fraction: float = 0.18,
) -> WingSpec:
    return build_demo_wing_designer(
        root_x=root_x,
        root_y=root_y,
        root_z=root_z,
        transition_continuity=transition_continuity,
        transition_fraction=transition_fraction,
    ).to_wing_spec()


def build_demo_wing_designer(
    root_x: float = 0.0,
    root_y: float = 0.0,
    root_z: float = 0.0,
    transition_continuity: ContinuityOrder = ContinuityOrder.C2,
    transition_fraction: float = 0.18,
) -> WingDesignerSpec:
    transition = WingTransitionSpec(
        method=InterpolationMethod.SEGMENTED,
        continuity=transition_continuity,
        blend_fraction=transition_fraction,
    )
    transition_interp = transition.to_interpolation_spec()
    return WingDesignerSpec(
        wing_id="demo_main_wing",
        span=29.0,
        root_chord=6.2,
        symmetric=True,
        root_le_x=float(root_x),
        root_vertical_y=float(root_y),
        root_z=float(root_z),
        spine_transition=transition,
        section_transition=transition,
        stations=(
            WingStationDesignerSpec(
                station_id="root",
                eta=0.0,
                profile_id="root_airfoil",
                chord_ratio=1.0,
                twist_deg=2.5,
            ),
            WingStationDesignerSpec(
                station_id="kink",
                eta=0.38,
                profile_id="kink_airfoil",
                chord_ratio=4.0 / 6.2,
                twist_deg=0.8,
                sweep_qc_deg=14.796803971186739,
                dihedral_deg=3.5,
            ),
            WingStationDesignerSpec(
                station_id="tip",
                eta=1.0,
                profile_id="tip_airfoil",
                chord_ratio=1.55 / 6.2,
                twist_deg=-2.4,
                sweep_qc_deg=25.92799182883875,
                dihedral_deg=5.5,
            ),
        ),
        scalar_laws=(
            ScalarLawSpec(
                name="twist_deg",
                anchors=(
                    ScalarLawAnchor(0.0, 2.5),
                    ScalarLawAnchor(0.38, 0.8),
                    ScalarLawAnchor(1.0, -2.4),
                ),
                interpolation=transition_interp,
            ),
            ScalarLawSpec(
                name="thickness_scale",
                anchors=(
                    ScalarLawAnchor(0.0, 1.00),
                    ScalarLawAnchor(0.60, 0.92),
                    ScalarLawAnchor(1.0, 0.86),
                ),
                interpolation=transition_interp,
            ),
        ),
        metadata={
            "description": "Generic demo wing for aircraft package examples",
            "transition_mode": "segmented",
            "transition_continuity": str(int(transition_continuity)),
        },
    )


def default_build_options() -> LiftingSurfaceBuildOptions:
    return LiftingSurfaceBuildOptions(
        station_count=13,
        include_anchor_sections=True,
        airfoil_sample_count=301,
        blunt_te=True,
        tip_style="rounded",
    )
