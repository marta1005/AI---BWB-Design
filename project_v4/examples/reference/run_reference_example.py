from dataclasses import asdict
import json
from pathlib import Path
import sys
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from project_v4.builder import export_iges
from project_v4.design_variables import SectionedBWBDesignVariables


def to_jsonable(value: Any):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {key: to_jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_jsonable(item) for item in value]
    return value


def build_reference_design() -> SectionedBWBDesignVariables:
    return SectionedBWBDesignVariables.reference_seed()


def main() -> None:
    output_dir = SCRIPT_DIR.parent.parent / "example_outputs" / "reference_v4_demo"
    profiles_dir = output_dir / "profiles_reference"
    iges_path = output_dir / "reference_v4.igs"
    metadata_path = output_dir / "reference_v4.json"
    output_dir.mkdir(parents=True, exist_ok=True)

    design = build_reference_design()
    config = design.to_model_config()
    config.export.out_dir = profiles_dir
    config.export.iges_path = iges_path

    prepared = export_iges(config)

    y_sections = config.topology.y_sections_array
    le_sections = config.planform.leading_edge_x_sections(config.topology)
    te_sections = config.planform.trailing_edge_x_sections(config.topology)
    chords = te_sections - le_sections

    metadata = {
        "example_name": "reference_v4_demo",
        "design_variables": asdict(design),
        "profile_generation_mode": config.sections.profile_generation_mode,
        "topology_name": config.topology.topology_name,
        "section_names": ["root", "section_2", "section_3", "tip"],
        "section_y_m": [float(value) for value in y_sections],
        "leading_edge_x_m": [float(value) for value in le_sections],
        "trailing_edge_x_m": [float(value) for value in te_sections],
        "chord_m": [float(value) for value in chords],
        "sweeps_deg": {
            "s1": float(config.planform.s1_deg),
            "s2": float(config.planform.s2_deg),
            "s3": float(config.planform.s3_deg),
        },
        "b_segment_ratios": {
            "b1_span_ratio": float(config.topology.b1_span_ratio),
            "b2_span_ratio": float(config.topology.b2_span_ratio),
            "b3_span_ratio": float(config.topology.b3_span_ratio),
        },
        "cst_section_parameters": {
            "c1": {
                "upper_coeffs": list(design.c1_upper_cst),
                "lower_coeffs": list(design.c1_lower_cst),
                "tc_max": float(design.c1_tc_max),
                "x_tmax": float(design.c1_x_tmax),
                "te_thickness": float(design.c1_te_thickness),
            },
            "c2": {
                "upper_coeffs": list(design.c2_upper_cst),
                "lower_coeffs": list(design.c2_lower_cst),
                "tc_max": float(design.c2_tc_max),
                "x_tmax": float(design.c2_x_tmax),
                "te_thickness": float(design.c2_te_thickness),
            },
            "c3": {
                "upper_coeffs": list(design.c3_upper_cst),
                "lower_coeffs": list(design.c3_lower_cst),
                "tc_max": float(design.c3_tc_max),
                "x_tmax": float(design.c3_x_tmax),
                "te_thickness": float(design.c3_te_thickness),
            },
            "c4": {
                "upper_coeffs": list(design.c4_upper_cst),
                "lower_coeffs": list(design.c4_lower_cst),
                "tc_max": float(design.c4_tc_max),
                "x_tmax": float(design.c4_x_tmax),
                "te_thickness": float(design.c4_te_thickness),
            },
        },
        "cst_class_parameters": {
            "n1": float(design.cst_n1),
            "n2": float(design.cst_n2),
        },
        "twist_law_deg": {
            "section_indices": list(config.spanwise.twist_deg.section_indices),
            "values": list(config.spanwise.twist_deg.values),
            "interpolation": config.spanwise.twist_deg.interpolation,
        },
        "validation": {
            "min_inner_tc": prepared.validation.min_inner_tc,
            "min_inner_tc_y": prepared.validation.min_inner_tc_y,
            "min_inner_tc_xc": prepared.validation.min_inner_tc_xc,
            "num_samples": prepared.validation.num_samples,
        },
        "volume_constraint": {
            "enabled": prepared.volume.enabled,
            "satisfied": prepared.volume.satisfied,
            "enclosed_volume_m3": prepared.volume.enclosed_volume_m3,
            "required_volume_m3": prepared.volume.required_volume_m3,
            "volume_margin_m3": prepared.volume.volume_margin_m3,
            "volume_ratio": prepared.volume.volume_ratio,
            "mean_cross_section_area_m2": prepared.volume.mean_cross_section_area_m2,
            "max_cross_section_area_m2": prepared.volume.max_cross_section_area_m2,
            "max_cross_section_area_y": prepared.volume.max_cross_section_area_y,
            "span_samples": prepared.volume.span_samples,
        },
        "iges_export_success": True,
        "iges_path": str(iges_path),
        "profiles_dir": str(profiles_dir),
        "config": to_jsonable(asdict(config)),
    }
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    print(
        "Reference v4 sections: "
        f"y={', '.join(f'{value:.2f}' for value in y_sections)}"
    )
    print(
        "Reference v4 chords: "
        f"{', '.join(f'{value:.2f}' for value in chords)} m"
    )
    print(
        "Reference v4 sweeps: "
        f"S1={config.planform.s1_deg:.1f} deg, "
        f"S2={config.planform.s2_deg:.1f} deg, "
        f"S3={config.planform.s3_deg:.1f} deg"
    )
    print(
        "Reference v4 CST class: "
        f"N1={config.sections.n1:.2f}, N2={config.sections.n2:.2f} | "
        f"full 6+6 per section"
    )
    print(
        "Reference v4 profile generation mode: "
        f"{config.sections.profile_generation_mode}"
    )
    print(
        "Reference v4 x_tmax: "
        f"({design.c1_x_tmax:.3f}, {design.c2_x_tmax:.3f}, "
        f"{design.c3_x_tmax:.3f}, {design.c4_x_tmax:.3f})"
    )
    print(
        "Reference v4 nose blend: "
        f"{config.planform.symmetry_blend_y:.2f} m"
    )
    print(f"Metadata written to: {metadata_path}")


if __name__ == "__main__":
    main()
