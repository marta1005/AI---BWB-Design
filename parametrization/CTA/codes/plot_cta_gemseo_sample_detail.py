from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
from typing import Dict, Tuple

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
from matplotlib.lines import Line2D

from parametrization.CTA.gemseo_space import build_cta_gemseo_design_space_definition
from parametrization.CTA.reference import cta_fixed_values, to_cta_model_config
from parametrization.CTA.codes.plot_cta_reference_views import (
    annotate_public_nomenclature_scheme,
    draw_3d,
    draw_planform,
)
from parametrization.bwb.builder import prepare_geometry
from parametrization.bwb.planform import build_sectioned_bwb_planform


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate planform, profiles and 3D views for one CTA GEMSEO DOE sample."
    )
    parser.add_argument("--sample-id", type=int, required=True, help="Sample id stored in cta_sampling_summary.json.")
    parser.add_argument(
        "--summary-json",
        type=Path,
        default=SCRIPT_DIR.parent / "outputs" / "cta_gemseo_doe_lhs_n4" / "cta_sampling_summary.json",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Optional output directory. Defaults next to the summary JSON in sample_<id>/",
    )
    parser.add_argument("--dpi", type=int, default=220)
    return parser.parse_args()


def _dense_span(config, n_points: int = 700) -> np.ndarray:
    y_sections = np.asarray(config.topology.y_sections_array, dtype=float)
    dense = np.linspace(0.0, float(config.topology.span), int(n_points))
    return np.unique(np.concatenate([dense, y_sections]))


def _placeholder(ax, title: str, message: str) -> None:
    ax.axis("off")
    ax.set_title(title)
    ax.text(
        0.5,
        0.5,
        message,
        ha="center",
        va="center",
        fontsize=11,
        color="#7f1d1d",
        bbox={"boxstyle": "round,pad=0.4", "facecolor": "#fef2f2", "edgecolor": "#fecaca", "alpha": 0.96},
        transform=ax.transAxes,
    )


def _find_sample(summary: Dict[str, object], sample_id: int) -> Dict[str, object]:
    for sample in summary.get("samples", []):
        if int(sample["sample_id"]) == int(sample_id):
            return sample
    raise KeyError(f"Sample id {sample_id} not found in {summary.get('n_samples_generated', 0)} generated samples.")


def _build_sample_geometry(sample: Dict[str, object]):
    adapter = build_cta_gemseo_design_space_definition()
    gemseo_sample = {
        name: np.asarray(values, dtype=float)
        for name, values in sample["gemseo_variables"].items()
    }
    design = adapter.to_project_design(gemseo_sample)
    config = to_cta_model_config(design)
    return design, config


def _planform_curves(config) -> Tuple[object, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    planform = build_sectioned_bwb_planform(config.topology, config.planform)
    dense_span = _dense_span(config)
    leading_edge_x = np.array([planform.le_x(float(y)) for y in dense_span], dtype=float)
    trailing_edge_x = np.array([planform.te_x(float(y)) for y in dense_span], dtype=float)
    chords = trailing_edge_x - leading_edge_x
    return planform, dense_span, leading_edge_x, trailing_edge_x, chords


def _public_station_data(config, planform):
    te_points = np.asarray(config.planform.trailing_edge_points(config.topology), dtype=float)
    public_station_indices = np.asarray((0, 1, 3, 4, te_points.shape[0] - 1), dtype=int)
    c_names = ("C0", "C1", "C3", "C4", "C5")
    c_y = te_points[public_station_indices, 1].astype(float)
    c_te_x = te_points[public_station_indices, 0].astype(float)
    c_le_x = np.array([planform.le_x(float(value)) for value in c_y], dtype=float)
    c_chords = c_te_x - c_le_x
    return c_names, c_y, c_le_x, c_te_x, c_chords


def _draw_profiles_overlay(ax, prepared, section_y: np.ndarray, section_labels: Tuple[str, ...]) -> None:
    colors = ("#0f4c5c", "#0891b2", "#1d4ed8", "#7c3aed", "#c44536")
    x_air = np.asarray(prepared.section_model.x_air, dtype=float)
    z_min = 0.0
    z_max = 0.0

    for yy, color, label in zip(section_y, colors, section_labels):
        yu, yl, _ = prepared.section_model.coordinates_at_y(float(yy))
        yu = np.asarray(yu, dtype=float)
        yl = np.asarray(yl, dtype=float)
        z_min = min(z_min, float(np.min(yl)))
        z_max = max(z_max, float(np.max(yu)))

        ax.fill(
            np.concatenate([x_air, x_air[::-1]]),
            np.concatenate([yu, yl[::-1]]),
            color=color,
            alpha=0.08,
            zorder=1,
        )
        ax.plot(x_air, yu, color=color, linewidth=2.0, zorder=3)
        ax.plot(x_air, yl, color=color, linewidth=2.0, zorder=3)

    legend_handles = [
        Line2D([0], [0], color=color, lw=2.2, label=label)
        for color, label in zip(colors, section_labels)
    ]
    ax.legend(handles=legend_handles, loc="upper right", framealpha=0.95, title="Section")
    ax.text(
        0.03,
        0.04,
        "Profiles overlaid at common scale (normalized by local chord)",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9.0,
        color="#334155",
        bbox={"boxstyle": "round,pad=0.20", "facecolor": "white", "edgecolor": "#cbd5e1", "alpha": 0.94},
    )
    ax.set_title("Section profiles (overlay, x/c and z/c)")
    ax.set_xlabel("x / c [-]")
    ax.set_ylabel("z / c [-]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, linewidth=0.35, alpha=0.25)
    ax.set_xlim(-0.02, 1.02)
    z_margin = 0.08 * max(z_max - z_min, 0.12)
    ax.set_ylim(z_min - z_margin, z_max + z_margin)


def main() -> None:
    args = parse_args()
    summary_path = args.summary_json.resolve()
    if not summary_path.exists():
        raise FileNotFoundError(f"CTA GEMSEO summary JSON not found: {summary_path}")

    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    sample = _find_sample(summary, args.sample_id)
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else summary_path.parent / f"sample_{int(args.sample_id):03d}"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    design, config = _build_sample_geometry(sample)
    y_sections = np.asarray(config.topology.y_sections_array, dtype=float)
    planform, dense_span, leading_edge_x, trailing_edge_x, chord_dense = _planform_curves(config)
    fixed_values = cta_fixed_values(reference_design=design)
    c_names, c_y, c_le_x, c_te_x, c_chords = _public_station_data(config, planform)

    geometry_error = ""
    prepared = None
    try:
        prepared = prepare_geometry(config)
    except Exception as exc:
        geometry_error = str(exc)

    valid = bool(sample["geometry_evaluation"].get("geometry_valid")) and prepared is not None
    title_prefix = f"CTA sample {int(args.sample_id)}"
    status = "valid" if valid else "invalid"
    title = f"{title_prefix} | {status}"

    planform_png = output_dir / f"cta_sample_{int(args.sample_id):03d}_planform.png"
    profiles_png = output_dir / f"cta_sample_{int(args.sample_id):03d}_profiles.png"
    view3d_png = output_dir / f"cta_sample_{int(args.sample_id):03d}_3d.png"

    fig_plan, ax_plan = plt.subplots(figsize=(11.8, 8.2), constrained_layout=True)
    draw_planform(
        ax_plan,
        dense_span,
        leading_edge_x,
        trailing_edge_x,
        y_sections,
        {},
        show_cargo=False,
        half_wing=True,
        label_zones=True,
        faded_symmetry=True,
    )
    annotate_public_nomenclature_scheme(
        ax_plan,
        planform,
        y_sections,
        c_y,
        c_le_x,
        c_te_x,
        c_chords,
        config,
        fixed_values,
    )
    ax_plan.set_title(title + " | planform")
    fig_plan.savefig(planform_png, dpi=int(args.dpi), bbox_inches="tight")
    plt.close(fig_plan)

    if prepared is not None:
        fig_profiles, ax_profiles = plt.subplots(figsize=(11.4, 8.0), constrained_layout=True)
        _draw_profiles_overlay(ax_profiles, prepared, c_y, c_names)
        ax_profiles.set_title(title + " | profiles")
        fig_profiles.savefig(profiles_png, dpi=int(args.dpi), bbox_inches="tight")
        plt.close(fig_profiles)

        fig_3d = plt.figure(figsize=(12.0, 8.2), constrained_layout=True)
        ax_3d = fig_3d.add_subplot(111, projection="3d")
        draw_3d(ax_3d, prepared, dense_span, leading_edge_x, chord_dense)
        ax_3d.set_title(title + " | 3D")
        fig_3d.savefig(view3d_png, dpi=int(args.dpi), bbox_inches="tight")
        plt.close(fig_3d)
    else:
        message = "Geometry could not be prepared.\n\n" + (geometry_error or "Unknown error")
        fig_profiles, ax_profiles = plt.subplots(figsize=(10.0, 6.0), constrained_layout=True)
        _placeholder(ax_profiles, title + " | profiles", message)
        fig_profiles.savefig(profiles_png, dpi=int(args.dpi), bbox_inches="tight")
        plt.close(fig_profiles)

        fig_3d, ax_3d = plt.subplots(figsize=(10.0, 6.0), constrained_layout=True)
        _placeholder(ax_3d, title + " | 3D", message)
        fig_3d.savefig(view3d_png, dpi=int(args.dpi), bbox_inches="tight")
        plt.close(fig_3d)

    print(f"CTA GEMSEO sample planform PNG written to: {planform_png}")
    print(f"CTA GEMSEO sample profiles PNG written to: {profiles_png}")
    print(f"CTA GEMSEO sample 3D PNG written to: {view3d_png}")


if __name__ == "__main__":
    main()
