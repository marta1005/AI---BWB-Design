import argparse
import csv
from dataclasses import replace
import json
import os
from pathlib import Path
import sys
import tempfile
from typing import Dict, List, Optional, Tuple

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR.parent.parent / ".mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from project_v4.builder import prepare_geometry
from project_v4.design_space import available_presets, build_design_space, flatten_design
from project_v4.design_variables import SectionedBWBDesignVariables


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Sample and visualize project_v4 design-space variants.")
    parser.add_argument("--preset", default="presentation_core", choices=available_presets())
    parser.add_argument("--count", type=int, default=4, help="Number of valid sampled variants to generate")
    parser.add_argument("--seed", type=int, default=7, help="Random seed")
    parser.add_argument(
        "--variation-scale",
        type=float,
        default=0.25,
        help="Local normalized box size around the reference design",
    )
    return parser.parse_args()


def save_figure(fig, path: Path, **kwargs) -> None:
    try:
        fig.savefig(path, **kwargs)
    except OSError as exc:
        if "Resource deadlock avoided" not in str(exc):
            raise
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir) / path.name
            fig.savefig(tmp_path, **kwargs)
            Path(path).write_bytes(tmp_path.read_bytes())


def generate_valid_variants(
    preset: str,
    count: int,
    seed: int,
    variation_scale: float,
    reference_design: SectionedBWBDesignVariables,
    root_tc_sequence: Optional[Tuple[float, ...]] = None,
) -> List[Tuple[str, SectionedBWBDesignVariables, object]]:
    design_space = build_design_space(preset, reference_design=reference_design)
    profile_generation_mode = (
        "enforce_targets" if "thickness_targets" in design_space.active_groups else "cst_only"
    )
    variants: List[Tuple[str, SectionedBWBDesignVariables, object]] = []
    accepted = 0
    batch_seed = int(seed)
    ref_tc_values = np.array(
        [
            reference_design.c1_tc_max,
            reference_design.c2_tc_max,
            reference_design.c3_tc_max,
            reference_design.c4_tc_max,
        ],
        dtype=float,
    )
    while accepted < count:
        batch = design_space.sample_designs(count=max(4, count - accepted), seed=batch_seed, variation_scale=variation_scale)
        batch_seed += 1
        for design in batch:
            candidate = design
            if root_tc_sequence is not None and accepted < len(root_tc_sequence):
                target_root_tc = float(root_tc_sequence[accepted])
                scale = target_root_tc / max(reference_design.c1_tc_max, 1e-12)
                scaled_tc = ref_tc_values * scale
                candidate = replace(
                    design,
                    c1_tc_max=float(scaled_tc[0]),
                    c2_tc_max=float(scaled_tc[1]),
                    c3_tc_max=float(scaled_tc[2]),
                    c4_tc_max=float(scaled_tc[3]),
                )
            try:
                prepared = prepare_geometry(
                    candidate.to_model_config(profile_generation_mode=profile_generation_mode)
                )
            except Exception:
                continue
            variants.append(("v%d" % (accepted + 1), candidate, prepared))
            accepted += 1
            if accepted >= count:
                break
        if batch_seed - seed > 100:
            raise RuntimeError("Could not generate enough valid design variants inside the requested local design space")
    return variants


def write_variant_table(
    csv_path: Path,
    active_variables: Tuple[str, ...],
    rows: List[Dict[str, object]],
) -> None:
    fieldnames = ["variant_id", "min_inner_tc", "min_inner_tc_y", "min_inner_tc_xc"] + list(active_variables)
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def build_variant_rows(
    active_variables: Tuple[str, ...],
    reference_design: SectionedBWBDesignVariables,
    reference_prepared,
    sampled_variants: List[Tuple[str, SectionedBWBDesignVariables, object]],
) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    variants = [("reference", reference_design, reference_prepared)] + sampled_variants
    for variant_id, design, prepared in variants:
        flat = flatten_design(design)
        row: Dict[str, object] = {
            "variant_id": variant_id,
            "min_inner_tc": prepared.validation.min_inner_tc,
            "min_inner_tc_y": prepared.validation.min_inner_tc_y,
            "min_inner_tc_xc": prepared.validation.min_inner_tc_xc,
        }
        for name in active_variables:
            row[name] = flat[name]
        rows.append(row)
    return rows


def plot_variants(
    output_path: Path,
    preset: str,
    active_groups: Tuple[str, ...],
    reference_design: SectionedBWBDesignVariables,
    reference_prepared,
    sampled_variants: List[Tuple[str, SectionedBWBDesignVariables, object]],
) -> None:
    fig = plt.figure(figsize=(14.0, 10.0), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=(2.2, 1.4))
    ax_planform = fig.add_subplot(grid[0, :])
    ax_root = fig.add_subplot(grid[1, 0])
    ax_tip = fig.add_subplot(grid[1, 1])

    colors = ["#0f172a", "#1d4ed8", "#0f766e", "#7c3aed", "#c2410c", "#be123c"]
    variants = [("reference", reference_design, reference_prepared)] + sampled_variants
    root_handles = []
    root_labels = []
    tip_handles = []
    tip_labels = []

    for idx, (variant_id, design, prepared) in enumerate(variants):
        color = colors[idx % len(colors)]
        config = design.to_model_config()
        dense_span = np.linspace(0.0, config.topology.span, 500)
        y_sections = config.topology.y_sections_array
        te_sections = config.planform.trailing_edge_x_sections(config.topology)
        le_x = np.array([prepared.planform.le_x(float(y)) for y in dense_span], dtype=float)
        te_mask = dense_span >= float(y_sections[1])
        te_span_outer = dense_span[te_mask]
        te_x_outer = np.array([prepared.planform.te_x(float(y)) for y in te_span_outer], dtype=float)
        te_exact_y = np.array([float(y_sections[0]), float(y_sections[1])], dtype=float)
        te_exact_x = np.array([float(te_sections[0]), float(te_sections[1])], dtype=float)

        label = "reference" if variant_id == "reference" else variant_id
        linewidth = 2.4 if variant_id == "reference" else 1.6
        alpha = 0.95 if variant_id == "reference" else 0.80

        ax_planform.plot(le_x, dense_span, color=color, linewidth=linewidth, alpha=alpha, label=label)
        ax_planform.plot(te_x_outer, te_span_outer, color=color, linewidth=linewidth, alpha=alpha)
        ax_planform.plot(
            te_exact_x,
            te_exact_y,
            color=color,
            linewidth=max(1.8, linewidth + 0.2),
            alpha=alpha,
            linestyle=(0, (6, 4)),
        )

        root_y = float(config.topology.y_sections_array[0])
        tip_y = float(config.topology.y_sections_array[-1])
        yu_root, yl_root, _ = prepared.section_model.coordinates_at_y(root_y)
        yu_tip, yl_tip, _ = prepared.section_model.coordinates_at_y(tip_y)
        root_metrics, _ = prepared.section_model.geometry_metrics_at_y(root_y)
        tip_metrics, _ = prepared.section_model.geometry_metrics_at_y(tip_y)
        x = prepared.section_model.x_air

        root_line = ax_root.plot(x, yu_root, color=color, linewidth=linewidth, alpha=alpha)[0]
        ax_root.plot(x, yl_root, color=color, linewidth=linewidth, alpha=alpha)
        if variant_id != "reference":
            root_handles.append(root_line)
            root_labels.append(f"t/c={root_metrics.max_tc:.2f}")

        tip_line = ax_tip.plot(x, yu_tip, color=color, linewidth=linewidth, alpha=alpha)[0]
        ax_tip.plot(x, yl_tip, color=color, linewidth=linewidth, alpha=alpha)
        if variant_id != "reference":
            tip_handles.append(tip_line)
            tip_labels.append(f"t/c={tip_metrics.max_tc:.2f}")

    ax_planform.set_title("Design-space variants: planform overlay")
    ax_planform.set_xlabel("x / chord direction")
    ax_planform.set_ylabel("semispan y")
    ax_planform.grid(True, linewidth=0.4, alpha=0.25)
    ax_planform.set_aspect("equal", adjustable="box")
    ax_planform.legend(loc="upper left", ncol=min(3, len(variants)))

    ax_root.set_title("Root section C1 (normalized)")
    ax_root.set_xlabel("x/c")
    ax_root.set_ylabel("z/c")
    ax_root.grid(True, linewidth=0.4, alpha=0.25)
    ax_root.set_aspect("equal", adjustable="box")
    ax_root.legend(root_handles, root_labels, loc="upper right", fontsize=8.0, framealpha=0.9, title="Thickness")

    ax_tip.set_title("Tip section C4 (normalized)")
    ax_tip.set_xlabel("x/c")
    ax_tip.grid(True, linewidth=0.4, alpha=0.25)
    ax_tip.set_aspect("equal", adjustable="box")
    ax_tip.legend(tip_handles, tip_labels, loc="upper right", fontsize=8.0, framealpha=0.9, title="Thickness")

    save_figure(fig, output_path, dpi=220, bbox_inches="tight")
    save_figure(fig, output_path.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    preset = args.preset
    output_dir = SCRIPT_DIR.parent.parent / "example_outputs" / ("design_space_%s" % preset)
    output_dir.mkdir(parents=True, exist_ok=True)

    reference_design = SectionedBWBDesignVariables.reference_seed()
    design_space = build_design_space(preset, reference_design=reference_design)
    profile_generation_mode = (
        "enforce_targets" if "thickness_targets" in design_space.active_groups else "cst_only"
    )
    reference_prepared = prepare_geometry(
        reference_design.to_model_config(profile_generation_mode=profile_generation_mode)
    )
    sampled_variants = generate_valid_variants(
        preset=preset,
        count=args.count,
        seed=args.seed,
        variation_scale=args.variation_scale,
        reference_design=reference_design,
        root_tc_sequence=(0.18, 0.20, 0.21, 0.22) if preset == "presentation_core" else None,
    )

    rows = build_variant_rows(
        active_variables=design_space.active_variables,
        reference_design=reference_design,
        reference_prepared=reference_prepared,
        sampled_variants=sampled_variants,
    )
    csv_path = output_dir / ("design_space_%s_variants.csv" % preset)
    json_path = output_dir / ("design_space_%s_summary.json" % preset)
    figure_path = output_dir / ("design_space_%s_overview.png" % preset)

    write_variant_table(csv_path, design_space.active_variables, rows)
    summary = {
        "preset": preset,
        "active_groups": list(design_space.active_groups),
        "active_variables": list(design_space.active_variables),
        "seed": int(args.seed),
        "variation_scale": float(args.variation_scale),
        "variant_count": 1 + len(sampled_variants),
        "variants": rows,
    }
    json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    plot_variants(
        output_path=figure_path,
        preset=preset,
        active_groups=design_space.active_groups,
        reference_design=reference_design,
        reference_prepared=reference_prepared,
        sampled_variants=sampled_variants,
    )

    print("Design-space preset: %s" % preset)
    print("Active groups: %s" % ", ".join(design_space.active_groups))
    print("Overview PNG written to: %s" % figure_path)
    print("Variant CSV written to: %s" % csv_path)
    print("Summary JSON written to: %s" % json_path)


if __name__ == "__main__":
    main()
