from __future__ import annotations

from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from parametrization.aircraft import (
    AircraftAssemblySpec,
    AircraftFuselageEntry,
    AircraftVerticalTailEntry,
    AircraftWingEntry,
    ContinuityOrder,
    FuselageBuildOptions,
    FuselageSectionSpec,
    FuselageSpec,
    InterpolationMethod,
    InterpolationSpec,
    IntuitiveAirfoilProfileSpec,
    LiftingSurfaceBuildOptions,
    ProfileCatalog,
    ScalarLawAnchor,
    ScalarLawSpec,
    VerticalTailSpec,
    VerticalTailStationSpec,
    WingSpec,
    WingStationSpec,
)


def default_output_dir() -> Path:
    return SCRIPT_DIR.parent / "example_outputs" / "demo_main_wing"


def default_fuselage_output_dir() -> Path:
    return SCRIPT_DIR.parent / "example_outputs" / "demo_fuselage"


def default_aircraft_output_dir() -> Path:
    return SCRIPT_DIR.parent / "example_outputs" / "demo_aircraft"


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
            IntuitiveAirfoilProfileSpec(
                profile_id="htp_root_airfoil",
                leading_edge_radius=0.010,
                max_thickness=0.105,
                x_tmax=0.31,
                max_camber=0.012,
                x_cmax=0.40,
                trailing_edge_wedge_angle_deg=10.0,
                trailing_edge_camber_angle_deg=-0.4,
                aft_control_x=0.70,
                te_thickness=0.0010,
            ),
            IntuitiveAirfoilProfileSpec(
                profile_id="htp_tip_airfoil",
                leading_edge_radius=0.006,
                max_thickness=0.080,
                x_tmax=0.29,
                max_camber=0.006,
                x_cmax=0.36,
                trailing_edge_wedge_angle_deg=8.5,
                trailing_edge_camber_angle_deg=-0.2,
                aft_control_x=0.68,
                te_thickness=0.0008,
            ),
            IntuitiveAirfoilProfileSpec(
                profile_id="vtp_root_airfoil",
                leading_edge_radius=0.011,
                max_thickness=0.112,
                x_tmax=0.31,
                max_camber=0.010,
                x_cmax=0.38,
                trailing_edge_wedge_angle_deg=10.5,
                trailing_edge_camber_angle_deg=-0.3,
                aft_control_x=0.70,
                te_thickness=0.0010,
            ),
            IntuitiveAirfoilProfileSpec(
                profile_id="vtp_tip_airfoil",
                leading_edge_radius=0.006,
                max_thickness=0.082,
                x_tmax=0.28,
                max_camber=0.004,
                x_cmax=0.34,
                trailing_edge_wedge_angle_deg=8.0,
                trailing_edge_camber_angle_deg=-0.1,
                aft_control_x=0.67,
                te_thickness=0.0008,
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
    transition_interp = InterpolationSpec(
        method=InterpolationMethod.SEGMENTED,
        continuity=transition_continuity,
        blend_fraction=transition_fraction,
    )
    return WingSpec(
        wing_id="demo_main_wing",
        semispan=14.5,
        stations=(
            WingStationSpec(
                station_id="root",
                eta=0.0,
                profile_id="root_airfoil",
                chord=6.2,
                twist_deg=2.5,
            ),
            WingStationSpec(
                station_id="kink",
                eta=0.38,
                profile_id="kink_airfoil",
                chord=4.0,
                twist_deg=0.8,
                sweep_le_deg=20.0,
                dihedral_deg=3.5,
            ),
            WingStationSpec(
                station_id="tip",
                eta=1.0,
                profile_id="tip_airfoil",
                chord=1.55,
                twist_deg=-2.4,
                sweep_le_deg=29.0,
                dihedral_deg=5.5,
            ),
        ),
        root_x=float(root_x),
        root_y=float(root_y),
        root_z=float(root_z),
        section_interpolation=transition_interp,
        spine_interpolation=transition_interp,
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


def default_tail_build_options() -> LiftingSurfaceBuildOptions:
    return LiftingSurfaceBuildOptions(
        station_count=7,
        include_anchor_sections=True,
        airfoil_sample_count=261,
        blunt_te=True,
        tip_style="rounded",
    )


def build_demo_fuselage() -> FuselageSpec:
    return FuselageSpec(
        fuselage_id="demo_fuselage",
        length=22.0,
        sections=(
            FuselageSectionSpec(
                section_id="nose_tip",
                eta=0.0,
                width=0.10,
                height=0.12,
                center_y=0.06,
                top_shape_exp=1.7,
                bottom_shape_exp=1.6,
                side_shape_exp=1.7,
            ),
            FuselageSectionSpec(
                section_id="radome_mid",
                eta=0.03,
                width=0.46,
                height=0.58,
                center_y=0.12,
                top_shape_exp=1.95,
                bottom_shape_exp=1.65,
                side_shape_exp=1.95,
            ),
            FuselageSectionSpec(
                section_id="windshield_base",
                eta=0.075,
                width=0.98,
                height=1.22,
                center_y=0.22,
                top_shape_exp=2.35,
                bottom_shape_exp=1.75,
                side_shape_exp=2.25,
            ),
            FuselageSectionSpec(
                section_id="cockpit",
                eta=0.14,
                width=1.62,
                height=1.90,
                center_y=0.25,
                top_shape_exp=2.8,
                bottom_shape_exp=1.8,
                side_shape_exp=2.6,
            ),
            FuselageSectionSpec(
                section_id="cabin_front",
                eta=0.30,
                width=2.85,
                height=3.05,
                center_y=0.08,
                top_shape_exp=2.8,
                bottom_shape_exp=1.55,
                side_shape_exp=2.5,
            ),
            FuselageSectionSpec(
                section_id="cabin_mid",
                eta=0.52,
                width=3.00,
                height=3.15,
                center_y=0.00,
                top_shape_exp=3.0,
                bottom_shape_exp=1.45,
                side_shape_exp=2.6,
            ),
            FuselageSectionSpec(
                section_id="cabin_aft",
                eta=0.74,
                width=2.55,
                height=2.70,
                center_y=-0.06,
                top_shape_exp=2.5,
                bottom_shape_exp=1.40,
                side_shape_exp=2.4,
            ),
            FuselageSectionSpec(
                section_id="tail_cone",
                eta=0.92,
                width=0.90,
                height=1.00,
                center_y=0.08,
                top_shape_exp=2.1,
                bottom_shape_exp=1.7,
                side_shape_exp=2.0,
            ),
            FuselageSectionSpec(
                section_id="tail_tip",
                eta=1.0,
                width=0.14,
                height=0.16,
                center_y=0.18,
                top_shape_exp=1.8,
                bottom_shape_exp=1.8,
                side_shape_exp=1.8,
            ),
        ),
        nose_x=0.0,
        section_interpolation=InterpolationSpec(
            method=InterpolationMethod.PYSPLINE,
            continuity=ContinuityOrder.C2,
        ),
        spine_interpolation=InterpolationSpec(
            method=InterpolationMethod.PYSPLINE,
            continuity=ContinuityOrder.C2,
        ),
        metadata={"description": "Generic demo fuselage for aircraft package examples"},
    )
def default_fuselage_build_options() -> FuselageBuildOptions:
    return FuselageBuildOptions(
        station_count=29,
        include_anchor_sections=True,
        perimeter_point_count=181,
        ku=4,
        kv=4,
        n_ctlu=19,
        n_ctlv=49,
    )


def build_demo_htp() -> WingSpec:
    return WingSpec(
        wing_id="demo_htp",
        semispan=5.6,
        stations=(
            WingStationSpec("root", 0.0, "htp_root_airfoil", chord=2.75, twist_deg=0.5),
            WingStationSpec("tip", 1.0, "htp_tip_airfoil", chord=1.10, twist_deg=-1.5, sweep_le_deg=28.0, dihedral_deg=7.0),
        ),
        root_x=17.9,
        root_y=0.92,
        root_z=0.0,
        section_interpolation=InterpolationSpec(method=InterpolationMethod.LINEAR),
        spine_interpolation=InterpolationSpec(method=InterpolationMethod.LINEAR),
        scalar_laws=(
            ScalarLawSpec(
                name="twist_deg",
                anchors=(ScalarLawAnchor(0.0, 0.5), ScalarLawAnchor(1.0, -1.5)),
                interpolation=InterpolationSpec(method=InterpolationMethod.LINEAR),
            ),
        ),
        metadata={"description": "Demo HTP"},
    )


def build_demo_vtp() -> VerticalTailSpec:
    return VerticalTailSpec(
        tail_id="demo_vtp",
        span=4.7,
        stations=(
            VerticalTailStationSpec("root", 0.0, "vtp_root_airfoil", chord=4.05),
            VerticalTailStationSpec("mid", 0.46, "vtp_root_airfoil", chord=2.55, sweep_le_deg=26.0),
            VerticalTailStationSpec("tip", 1.0, "vtp_tip_airfoil", chord=1.10, sweep_le_deg=36.0),
        ),
        root_x=17.9,
        root_y=1.20,
        root_z=0.0,
        section_interpolation=InterpolationSpec(method=InterpolationMethod.LINEAR),
        spine_interpolation=InterpolationSpec(method=InterpolationMethod.LINEAR),
        metadata={"description": "Demo VTP"},
    )


def build_demo_aircraft() -> AircraftAssemblySpec:
    return AircraftAssemblySpec(
        aircraft_id="demo_transport",
        profiles=build_demo_profiles(),
        wings=(
            AircraftWingEntry(
                wing=build_demo_wing(root_x=8.35, root_y=0.10, root_z=0.0),
                options=default_build_options(),
                mirror=True,
            ),
            AircraftWingEntry(
                wing=build_demo_htp(),
                options=default_tail_build_options(),
                mirror=True,
            ),
        ),
        vertical_tails=(
            AircraftVerticalTailEntry(
                tail=build_demo_vtp(),
                options=default_tail_build_options(),
            ),
        ),
        fuselages=(
            AircraftFuselageEntry(
                fuselage=build_demo_fuselage(),
                options=default_fuselage_build_options(),
            ),
        ),
        metadata={"description": "Demo full-aircraft assembly"},
    )
