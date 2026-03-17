import argparse
import csv
from dataclasses import asdict
import json
from pathlib import Path
import re
import sys

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from project_v4.design_space import flatten_design
from project_v4.gemseo_space import (
    available_gemseo_doe_algorithms,
    build_gemseo_design_space,
    evaluate_gemseo_sample_geometry,
    sample_gemseo_doe,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Sample the project_v4 geometric design space with GEMSEO DOE Latin Hypercube "
            "(or another available GEMSEO DOE algorithm)."
        )
    )
    parser.add_argument("--preset", default="ai_geometry_core")
    parser.add_argument("--n-samples", type=int, default=16, help="Number of DOE samples to generate.")
    parser.add_argument(
        "--algo",
        default="LHS",
        help="GEMSEO DOE algorithm name. Default: LHS. Alternatives depend on the installation, e.g. OT_LHS.",
    )
    parser.add_argument("--seed", type=int, default=7, help="Random seed passed to GEMSEO.")
    parser.add_argument(
        "--profile-generation-mode",
        default="cst_only",
        choices=("cst_only", "enforce_targets"),
        help="Section generation mode used for geometry evaluation. Default: cst_only.",
    )
    parser.add_argument(
        "--volume-required",
        type=float,
        default=None,
        help="If provided, evaluate a volume constraint with this required total enclosed volume [m^3].",
    )
    parser.add_argument(
        "--volume-span-samples",
        type=int,
        default=161,
        help="Number of spanwise samples used to integrate the enclosed volume when --volume-required is set.",
    )
    parser.add_argument(
        "--interpolation-override",
        choices=("linear", "pchip", "cubic", "pyspline"),
        default=None,
        help="Optional override for section/twist/camber spanwise interpolation during geometry evaluation.",
    )
    return parser.parse_args()


def _sanitize_name(text: str) -> str:
    return re.sub(r"[^a-zA-Z0-9]+", "_", text.strip()).strip("_").lower()


def _write_csv(path: Path, fieldnames, rows) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    args = parse_args()
    adapter = build_gemseo_design_space(args.preset)
    available_algorithms = available_gemseo_doe_algorithms()
    if args.algo not in available_algorithms:
        raise ValueError(
            "Algorithm %r is not available. Available GEMSEO DOE algorithms: %s"
            % (args.algo, ", ".join(available_algorithms))
        )

    flat_vectors = sample_gemseo_doe(
        adapter=adapter,
        n_samples=args.n_samples,
        algo_name=args.algo,
        seed=args.seed,
    )

    flat_names = adapter.flat_gemseo_variable_names()
    project_parameter_names = list(adapter.project_space.active_variables)
    output_slug = "gemseo_doe_%s_%s_n%d" % (
        _sanitize_name(args.preset),
        _sanitize_name(args.algo),
        int(args.n_samples),
    )
    output_dir = SCRIPT_DIR.parent.parent / "example_outputs" / output_slug
    output_dir.mkdir(parents=True, exist_ok=True)

    flat_rows = []
    project_rows = []
    geometry_rows = []
    grouped_rows = []
    fixed_parameters = {
        name: float(value)
        for name, value in flatten_design(adapter.project_space.reference_design).items()
        if name not in adapter.project_space.active_variables
    }

    for sample_id, flat_vector in enumerate(flat_vectors):
        sample = adapter.flat_vector_to_gemseo_sample(flat_vector)
        design = adapter.to_project_design(sample)
        project_flat = flatten_design(design)
        geometry_evaluation = evaluate_gemseo_sample_geometry(
            adapter=adapter,
            sample=sample,
            profile_generation_mode=args.profile_generation_mode,
            required_volume_m3=args.volume_required,
            volume_span_samples=args.volume_span_samples,
            interpolation_override=args.interpolation_override,
        )

        flat_row = {"sample_id": sample_id}
        flat_row.update(
            {
                variable_name: float(value)
                for variable_name, value in zip(flat_names, flat_vector)
            }
        )
        flat_rows.append(flat_row)

        project_row = {"sample_id": sample_id}
        project_row.update(
            {
                parameter_name: float(project_flat[parameter_name])
                for parameter_name in project_parameter_names
            }
        )
        project_rows.append(project_row)

        geometry_row = {"sample_id": sample_id}
        geometry_row.update(asdict(geometry_evaluation))
        geometry_rows.append(geometry_row)

        grouped_rows.append(
            {
                "sample_id": sample_id,
                "gemseo_variables": {
                    name: np_array.tolist()
                    for name, np_array in sample.items()
                },
                "project_active_parameters": {
                    parameter_name: float(project_flat[parameter_name])
                    for parameter_name in project_parameter_names
                },
                "geometry_evaluation": asdict(geometry_evaluation),
            }
        )

    flat_csv_path = output_dir / "gemseo_flat_samples.csv"
    project_csv_path = output_dir / "project_active_parameters.csv"
    geometry_csv_path = output_dir / "sample_geometry_metrics.csv"
    summary_json_path = output_dir / "sampling_summary.json"

    _write_csv(flat_csv_path, ["sample_id"] + flat_names, flat_rows)
    _write_csv(project_csv_path, ["sample_id"] + project_parameter_names, project_rows)
    _write_csv(
        geometry_csv_path,
        [
            "sample_id",
            "geometry_valid",
            "error_message",
            "profile_generation_mode",
            "volume_enabled",
            "volume_satisfied",
            "enclosed_volume_m3",
            "required_volume_m3",
            "volume_margin_m3",
            "volume_ratio",
            "mean_cross_section_area_m2",
            "max_cross_section_area_m2",
            "max_cross_section_area_y",
            "min_inner_tc",
            "min_inner_tc_y",
            "min_inner_tc_xc",
        ],
        geometry_rows,
    )

    n_geometry_valid = sum(1 for row in geometry_rows if row["geometry_valid"])
    n_volume_satisfied = sum(
        1 for row in geometry_rows if row["geometry_valid"] and row["volume_enabled"] and row["volume_satisfied"]
    )

    summary = {
        "preset": args.preset,
        "algo": args.algo,
        "seed": int(args.seed),
        "profile_generation_mode": args.profile_generation_mode,
        "interpolation_override": args.interpolation_override,
        "volume_constraint": {
            "enabled": args.volume_required is not None,
            "required_volume_m3": args.volume_required,
            "span_samples": int(args.volume_span_samples),
        },
        "n_samples_requested": int(args.n_samples),
        "n_samples_generated": len(flat_vectors),
        "n_geometry_valid": n_geometry_valid,
        "n_volume_satisfied": n_volume_satisfied,
        "gemseo_variables": [spec.name for spec in adapter.variable_specs],
        "flat_gemseo_variable_names": flat_names,
        "project_active_parameters": project_parameter_names,
        "fixed_parameters_held_constant": fixed_parameters,
        "samples": grouped_rows,
    }
    summary_json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("GEMSEO DOE algorithm: %s" % args.algo)
    print("Requested samples: %d" % int(args.n_samples))
    print("Generated samples: %d" % len(flat_vectors))
    print("Geometry-valid samples: %d" % n_geometry_valid)
    if args.volume_required is not None:
        print("Volume-satisfied samples: %d" % n_volume_satisfied)
    print("Flat GEMSEO CSV written to: %s" % flat_csv_path)
    print("Project-parameter CSV written to: %s" % project_csv_path)
    print("Geometry-metrics CSV written to: %s" % geometry_csv_path)
    print("Summary JSON written to: %s" % summary_json_path)


if __name__ == "__main__":
    main()
