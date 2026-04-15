from __future__ import annotations

import os
from pathlib import Path
import sys
from typing import Dict, Iterable, List, Tuple

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
CTA_DIR = SCRIPT_DIR.parent.parent
REPO_ROOT = CTA_DIR.parent.parent
MPLCONFIG_DIR = REPO_ROOT / ".mplconfig"
MPLCONFIG_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG_DIR))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from parametrization.CTA.design_space import build_cta_design_space
from parametrization.CTA.codes.plotting.animate_cta_planform_parameter_sweep import _public_sample_to_cta_design
from parametrization.CTA.case import build_cta_design, to_cta_model_config
from parametrization.bwb.planform import build_sectioned_bwb_planform


OUTPUT_DIR = CTA_DIR / "outputs" / "wing"
GRID_PNG = CTA_DIR / "outputs" / "wing" / "cta_planform_parameter_variations_grid.png"


def _dense_span(config, n_points: int = 700) -> np.ndarray:
    y_sections = np.asarray(config.topology.y_sections_array, dtype=float)
    dense = np.linspace(0.0, float(config.topology.span), int(n_points))
    return np.unique(np.concatenate([dense, y_sections]))


def _planform_curves(config):
    planform = build_sectioned_bwb_planform(config.topology, config.planform)
    dense_span = _dense_span(config)
    leading_edge_x = np.array([planform.le_x(float(y)) for y in dense_span], dtype=float)
    trailing_edge_x = np.array([planform.te_x(float(y)) for y in dense_span], dtype=float)
    return dense_span, leading_edge_x, trailing_edge_x


def _base_background(ax, ref_dense_span: np.ndarray, ref_le: np.ndarray, ref_te: np.ndarray, y_sections: np.ndarray) -> None:
    zone_defs = (
        ("Center body", float(y_sections[0]), float(y_sections[1]), "#fde68a"),
        ("Transition wing", float(y_sections[1]), float(y_sections[2]), "#bfdbfe"),
        ("Outer wing", float(y_sections[2]), float(y_sections[3]), "#bbf7d0"),
    )
    for _, y0, y1, color in zone_defs:
        mask = (ref_dense_span >= y0) & (ref_dense_span <= y1)
        ax.fill_betweenx(
            ref_dense_span[mask],
            ref_le[mask],
            ref_te[mask],
            color=color,
            alpha=0.35,
            zorder=0,
        )
        ax.fill_betweenx(
            -ref_dense_span[mask],
            ref_le[mask],
            ref_te[mask],
            color=color,
            alpha=0.08,
            zorder=0,
        )

    ax.plot(ref_le, ref_dense_span, color="#475569", linewidth=1.5, linestyle="--", alpha=0.85, zorder=2)
    ax.plot(ref_te, ref_dense_span, color="#475569", linewidth=1.5, linestyle="--", alpha=0.85, zorder=2)
    ax.plot(ref_le, -ref_dense_span, color="#94a3b8", linewidth=1.0, linestyle="--", alpha=0.6, zorder=1)
    ax.plot(ref_te, -ref_dense_span, color="#94a3b8", linewidth=1.0, linestyle="--", alpha=0.6, zorder=1)

    for y_value in y_sections:
        ax.axhline(float(y_value), color="#94a3b8", linewidth=0.6, linestyle=(0, (4, 4)), alpha=0.35, zorder=0)
    ax.axhline(0.0, color="#94a3b8", linewidth=0.8, linestyle=":", alpha=0.5, zorder=0)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, linewidth=0.3, alpha=0.20)


def _style_axis(ax, xlim: Tuple[float, float], ylim: Tuple[float, float]) -> None:
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("spanwise y [m]")


def _sample_with_single_override(name: str, value: float) -> Dict[str, float]:
    space = build_cta_design_space()
    sample = dict(space.cta_flat())
    sample[name] = float(value)
    return sample


def _actual_c4_from_sample(sample: Dict[str, float]) -> float:
    return float(sample["c2_c1_ratio"]) * float(sample["c4_c3_ratio"])


def _parameter_specs(bounds: Dict[str, Tuple[float, float]], cta_view: Dict[str, float]):
    return [
        {
            "name": "c1_root_chord",
            "title": "C0 variation",
            "display": lambda v, _: f"C0 = {v:.2f} m",
            "lower": float(bounds["c1_root_chord"][0]),
            "cta": float(cta_view["c1_root_chord"]),
            "upper": float(bounds["c1_root_chord"][1]),
        },
        {
            "name": "c2_c1_ratio",
            "title": "C3 variation",
            "display": lambda v, _: f"C3 = {v:.2f} m",
            "lower": float(bounds["c2_c1_ratio"][0]),
            "cta": float(cta_view["c2_c1_ratio"]),
            "upper": float(bounds["c2_c1_ratio"][1]),
        },
        {
            "name": "c4_c3_ratio",
            "title": "C4 variation via taper ratio",
            "display": lambda v, sample: f"C4/C3 = {v:.3f}\nC4 = {_actual_c4_from_sample(sample):.2f} m",
            "lower": float(bounds["c4_c3_ratio"][0]),
            "cta": float(cta_view["c4_c3_ratio"]),
            "upper": float(bounds["c4_c3_ratio"][1]),
        },
        {
            "name": "c4_c1_ratio",
            "title": "C5 variation",
            "display": lambda v, _: f"C5 = {v:.2f} m",
            "lower": float(bounds["c4_c1_ratio"][0]),
            "cta": float(cta_view["c4_c1_ratio"]),
            "upper": float(bounds["c4_c1_ratio"][1]),
        },
        {
            "name": "s2_deg",
            "title": "S1 variation",
            "display": lambda v, _: f"S1 (50%) = {v:.2f}°",
            "lower": float(bounds["s2_deg"][0]),
            "cta": float(cta_view["s2_deg"]),
            "upper": float(bounds["s2_deg"][1]),
        },
        {
            "name": "s3_deg",
            "title": "S2 variation",
            "display": lambda v, _: f"S2 (25%) = {v:.2f}°",
            "lower": float(bounds["s3_deg"][0]),
            "cta": float(cta_view["s3_deg"]),
            "upper": float(bounds["s3_deg"][1]),
        },
        {
            "name": "b2_span_ratio",
            "title": "B2/Bw variation",
            "display": lambda v, _: f"B2/Bw = {v:.3f}",
            "lower": float(bounds["b2_span_ratio"][0]),
            "cta": float(cta_view["b2_span_ratio"]),
            "upper": float(bounds["b2_span_ratio"][1]),
        },
        {
            "name": "span",
            "title": "Bw variation",
            "display": lambda v, _: f"Bw = {v:.2f} m",
            "lower": float(bounds["span"][0]),
            "cta": float(cta_view["span"]),
            "upper": float(bounds["span"][1]),
        },
        {
            "name": "med_3_te_sweep_deg",
            "title": "med_3_TEswp",
            "display": lambda v, _: f"med_3_TEswp = {v:.2f}°",
            "lower": float(bounds["med_3_te_sweep_deg"][0]),
            "cta": float(cta_view["med_3_te_sweep_deg"]),
            "upper": float(bounds["med_3_te_sweep_deg"][1]),
        },
    ]


def _plot_variation_panel(
    ax,
    spec: Dict[str, object],
    cta_curves,
    y_sections: np.ndarray,
    xlim: Tuple[float, float],
    ylim: Tuple[float, float],
) -> None:
    ref_dense_span, ref_le, ref_te = cta_curves
    _base_background(ax, ref_dense_span, ref_le, ref_te, y_sections)

    levels = [
        ("Lower", float(spec["lower"]), "#2563eb"),
        ("CTA", float(spec["cta"]), "#0f172a"),
        ("Upper", float(spec["upper"]), "#dc2626"),
    ]

    labels = []
    for label, value, color in levels:
        sample = _sample_with_single_override(str(spec["name"]), float(value))
        design = _public_sample_to_cta_design(sample)
        config = to_cta_model_config(design)
        dense_span, le_x, te_x = _planform_curves(config)
        ax.plot(le_x, dense_span, color=color, linewidth=2.1, zorder=4)
        ax.plot(te_x, dense_span, color=color, linewidth=2.1, zorder=4)
        ax.plot(le_x, -dense_span, color=color, linewidth=1.1, alpha=0.40, zorder=3)
        ax.plot(te_x, -dense_span, color=color, linewidth=1.1, alpha=0.40, zorder=3)
        labels.append(f"{label}: {spec['display'](float(value), sample)}")

    ax.set_title(str(spec["title"]))
    ax.text(
        0.02,
        0.98,
        "\n".join(labels),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.8,
        color="#0f172a",
        bbox={"boxstyle": "round,pad=0.26", "facecolor": "white", "edgecolor": "#cbd5e1", "alpha": 0.96},
        zorder=10,
    )
    _style_axis(ax, xlim, ylim)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    GRID_PNG.parent.mkdir(parents=True, exist_ok=True)

    space = build_cta_design_space()
    cta_view = space.cta_flat()
    bounds = space.bounds
    ref_config = to_cta_model_config(build_cta_design())
    cta_curves = _planform_curves(ref_config)
    ref_dense_span, ref_le, ref_te = cta_curves
    y_sections = np.asarray(ref_config.topology.y_sections_array, dtype=float)

    x_min = float(np.min(ref_le))
    x_max = float(np.max(ref_te))
    x_margin = 0.08 * max(x_max - x_min, 1.0)
    y_max = float(np.max(ref_dense_span))
    xlim = (x_min - x_margin, x_max + x_margin)
    ylim = (-0.18 * y_max, y_max + 0.10 * y_max)

    specs = _parameter_specs(bounds, cta_view)

    fig, axes = plt.subplots(3, 3, figsize=(17.5, 15.0), constrained_layout=True)
    axes_flat = axes.ravel()
    for ax, spec in zip(axes_flat, specs):
        _plot_variation_panel(ax, spec, cta_curves, y_sections, xlim, ylim)
    for ax in axes_flat[len(specs) :]:
        ax.axis("off")
    fig.suptitle("CTA planform parameter variations", fontsize=18)
    fig.savefig(GRID_PNG, dpi=190)
    plt.close(fig)

    print(f"CTA planform variation grid written to: {GRID_PNG}")


if __name__ == "__main__":
    main()
