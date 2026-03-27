from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import sys
from typing import Dict, List

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
MPLCONFIG_DIR = REPO_ROOT / ".mplconfig"
MPLCONFIG_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG_DIR))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from parametrization.CTA.gemseo_space import build_cta_gemseo_design_space_definition
from parametrization.CTA.reference import build_reference_design, to_cta_model_config
from parametrization.bwb.planform import build_sectioned_bwb_planform


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot a random subset of CTA GEMSEO DOE samples as planform overlays."
    )
    parser.add_argument(
        "--summary-json",
        type=Path,
        default=SCRIPT_DIR.parent / "outputs" / "cta_gemseo_doe_lhs_n4" / "cta_sampling_summary.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional output PNG path. Defaults next to the summary JSON.",
    )
    parser.add_argument(
        "--count",
        type=int,
        default=9,
        help="Number of random samples to visualize.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=7,
        help="Random seed used to choose the subsamples.",
    )
    parser.add_argument("--dpi", type=int, default=220)
    return parser.parse_args()


def _dense_span(config, n_points: int = 600) -> np.ndarray:
    anchors = np.asarray(config.topology.y_sections_array, dtype=float)
    dense = np.linspace(0.0, float(config.topology.span), int(n_points))
    return np.unique(np.concatenate([dense, anchors]))


def _planform_curves(config) -> Dict[str, np.ndarray]:
    planform = build_sectioned_bwb_planform(config.topology, config.planform)
    y = _dense_span(config)
    le_x = np.array([planform.le_x(float(yy)) for yy in y], dtype=float)
    te_x = np.array([planform.te_x(float(yy)) for yy in y], dtype=float)
    return {"y": y, "le_x": le_x, "te_x": te_x}


def _sample_to_text(sample: Dict[str, object]) -> str:
    params = sample["cta_public_parameters"]
    return "\n".join(
        [
        f"Bw = {float(params['span']):.2f} m",
        f"C0 = {float(params['c1_root_chord']):.2f} m",
        f"C3 = {float(params['c2_c1_ratio']):.2f} m",
        f"C4 = {float(params['c3_c1_ratio']):.2f} m",
        f"C5 = {float(params['c4_c1_ratio']):.2f} m",
        f"B2/Bw = {float(params['b2_span_ratio']):.3f}",
        f"S1 = {float(params['s2_deg']):.1f} deg",
        f"S2 = {float(params['s3_deg']):.1f} deg",
        ]
    )


def main() -> None:
    args = parse_args()
    summary_path = args.summary_json.resolve()
    if not summary_path.exists():
        raise FileNotFoundError(f"CTA GEMSEO summary JSON not found: {summary_path}")

    output_path = (
        args.output.resolve()
        if args.output is not None
        else summary_path.with_name("cta_gemseo_random_subsamples.png")
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)

    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    samples: List[Dict[str, object]] = list(summary.get("samples", []))
    if not samples:
        raise ValueError(f"No samples found in {summary_path}")

    count = max(1, int(args.count))
    rng = np.random.default_rng(int(args.seed))
    if count >= len(samples):
        chosen_indices = np.arange(len(samples), dtype=int)
    else:
        chosen_indices = np.sort(rng.choice(len(samples), size=count, replace=False))
    selected_samples = [samples[int(idx)] for idx in chosen_indices]

    adapter = build_cta_gemseo_design_space_definition()

    reference_design = build_reference_design()
    reference_config = to_cta_model_config(reference_design)
    reference_curves = _planform_curves(reference_config)

    sample_curves: List[Dict[str, object]] = []
    x_min = float(np.min(reference_curves["le_x"]))
    x_max = float(np.max(reference_curves["te_x"]))
    y_max = float(np.max(reference_curves["y"]))

    for sample in selected_samples:
        gemseo_sample = {
            name: np.asarray(values, dtype=float)
            for name, values in sample["gemseo_variables"].items()
        }
        design = adapter.to_project_design(gemseo_sample)
        config = to_cta_model_config(design)
        curves = _planform_curves(config)
        sample_curves.append(curves)
        x_min = min(x_min, float(np.min(curves["le_x"])))
        x_max = max(x_max, float(np.max(curves["te_x"])))
        y_max = max(y_max, float(np.max(curves["y"])))

    n_samples = len(selected_samples)
    ncols = min(3, max(1, math.ceil(math.sqrt(n_samples))))
    nrows = math.ceil(n_samples / ncols)

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(6.4 * ncols, 5.2 * nrows),
        squeeze=False,
        constrained_layout=True,
    )

    n_valid = sum(1 for sample in selected_samples if sample["geometry_evaluation"].get("geometry_valid"))
    fig.suptitle(
        (
            f"CTA GEMSEO DOE Random Subsamples\n"
            f"{n_samples} shown (seed={int(args.seed)}), {n_valid} valid, {n_samples - n_valid} invalid"
        ),
        fontsize=16,
        fontweight="semibold",
    )

    for ax, sample, curves in zip(axes.flat, selected_samples, sample_curves):
        geom = sample["geometry_evaluation"]
        valid = bool(geom.get("geometry_valid"))
        outline = "#1d4ed8" if valid else "#dc2626"
        fill = "#bfdbfe" if valid else "#fecaca"
        label = "valid" if valid else "invalid"

        ax.fill_betweenx(
            reference_curves["y"],
            reference_curves["le_x"],
            reference_curves["te_x"],
            color="#e5e7eb",
            alpha=0.55,
            zorder=0,
        )
        ax.plot(reference_curves["le_x"], reference_curves["y"], color="#64748b", linestyle="--", linewidth=1.2)
        ax.plot(reference_curves["te_x"], reference_curves["y"], color="#64748b", linestyle="--", linewidth=1.2)

        ax.fill_betweenx(curves["y"], curves["le_x"], curves["te_x"], color=fill, alpha=0.45, zorder=1)
        ax.plot(curves["le_x"], curves["y"], color=outline, linewidth=2.2, zorder=2)
        ax.plot(curves["te_x"], curves["y"], color=outline, linewidth=2.2, zorder=2)
        ax.plot([curves["le_x"][0], curves["te_x"][0]], [0.0, 0.0], color=outline, linewidth=2.0, zorder=2)
        ax.plot(
            [curves["le_x"][-1], curves["te_x"][-1]],
            [curves["y"][-1], curves["y"][-1]],
            color=outline,
            linewidth=2.0,
            zorder=2,
        )

        sample_id = int(sample["sample_id"])
        ax.set_title(f"Sample {sample_id} | {label}", fontsize=12, color=outline, fontweight="semibold")
        ax.text(
            0.02,
            0.98,
            _sample_to_text(sample),
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=9,
            bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#cbd5e1", "alpha": 0.92},
        )
        ax.set_xlim(x_min - 1.5, x_max + 1.5)
        ax.set_ylim(-0.6, y_max + 1.5)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("spanwise y [m]")
        ax.grid(True, alpha=0.22)
        ax.set_aspect("equal", adjustable="box")

    for ax in axes.flat[n_samples:]:
        ax.axis("off")

    fig.savefig(output_path, dpi=int(args.dpi))
    plt.close(fig)
    print(f"CTA GEMSEO random subsamples PNG written to: {output_path}")


if __name__ == "__main__":
    main()
