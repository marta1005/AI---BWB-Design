from __future__ import annotations

import os
from dataclasses import replace
from pathlib import Path
import sys
import tempfile
from typing import Dict, List

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
from PIL import Image

from parametrization.CTA.design_space import build_cta_design_space
from parametrization.CTA.codes.plotting.plot_cta_views import draw_planform
from parametrization.CTA.case import build_cta_design, to_cta_model_config
from parametrization.bwb.planform import build_sectioned_bwb_planform


OUTPUT_GIF = CTA_DIR / "outputs" / "wing" / "cta_planform_parameter_sweep.gif"


def _dense_span(config, n_points: int = 700) -> np.ndarray:
    y_sections = np.asarray(config.topology.y_sections_array, dtype=float)
    dense = np.linspace(0.0, float(config.topology.span), int(n_points))
    return np.unique(np.concatenate([dense, y_sections]))


def _planform_curves(config):
    planform = build_sectioned_bwb_planform(config.topology, config.planform)
    dense_span = _dense_span(config)
    leading_edge_x = np.array([planform.le_x(float(y)) for y in dense_span], dtype=float)
    trailing_edge_x = np.array([planform.te_x(float(y)) for y in dense_span], dtype=float)
    return planform, dense_span, leading_edge_x, trailing_edge_x


def _public_sample_to_cta_design(sample: Dict[str, float]):
    cta_design = build_cta_design()
    b1_fixed_m = float(cta_design.span * cta_design.b1_span_ratio)
    wing_span_m = float(sample["span"])
    total_span_m = float(b1_fixed_m + wing_span_m)
    b2_m = float(sample["b2_span_ratio"]) * wing_span_m
    b3_m = float(wing_span_m - b2_m)
    c0_m = float(sample["c1_root_chord"])
    c3_m = float(sample["c2_c1_ratio"])
    c4_m = float(sample["c2_c1_ratio"]) * float(sample["c4_c3_ratio"])
    c5_m = float(sample["c4_c1_ratio"])
    return replace(
        cta_design,
        span=total_span_m,
        b1_span_ratio=b1_fixed_m / total_span_m,
        b2_span_ratio=b2_m / total_span_m,
        b3_span_ratio=b3_m / total_span_m,
        c1_root_chord=c0_m,
        c2_c1_ratio=c3_m / c0_m,
        c3_c1_ratio=c4_m / c0_m,
        c4_c1_ratio=c5_m / c0_m,
        s2_deg=float(sample["s2_deg"]),
        s3_deg=float(sample["s3_deg"]),
    )


PUBLIC_PARAMETER_NAMES = (
    "span",
    "b2_span_ratio",
    "c1_root_chord",
    "c2_c1_ratio",
    "c4_c3_ratio",
    "c4_c1_ratio",
    "s2_deg",
    "s3_deg",
)


def _parameter_box_text(sample: Dict[str, float]) -> str:
    c4_m = float(sample["c2_c1_ratio"]) * float(sample["c4_c3_ratio"])
    return "\n".join(
        (
            f"Bw = {float(sample['span']):.2f} m",
            f"B2/Bw = {float(sample['b2_span_ratio']):.3f}",
            f"C0 = {float(sample['c1_root_chord']):.2f} m",
            f"C3 = {float(sample['c2_c1_ratio']):.2f} m",
            f"C4/C3 = {float(sample['c4_c3_ratio']):.3f}",
            f"C4 = {c4_m:.2f} m",
            f"C5 = {float(sample['c4_c1_ratio']):.2f} m",
            f"S1 = {float(sample['s2_deg']):.2f} deg",
            f"S2 = {float(sample['s3_deg']):.2f} deg",
        )
    )


def _latin_hypercube_samples(
    bounds: Dict[str, tuple[float, float]],
    sample_count: int,
    rng: np.random.Generator,
) -> List[Dict[str, float]]:
    if sample_count <= 0:
        return []

    samples: List[Dict[str, float]] = [
        {} for _ in range(int(sample_count))
    ]
    for name in PUBLIC_PARAMETER_NAMES:
        lower, upper = bounds[name]
        unit_values = (rng.permutation(sample_count) + rng.random(sample_count)) / float(sample_count)
        values = float(lower) + unit_values * (float(upper) - float(lower))
        for idx, value in enumerate(values):
            samples[idx][name] = float(value)
    return samples


def _normalized_distance(
    sample_a: Dict[str, float],
    sample_b: Dict[str, float],
    bounds: Dict[str, tuple[float, float]],
) -> float:
    deltas = []
    for name in PUBLIC_PARAMETER_NAMES:
        lower, upper = bounds[name]
        span = max(float(upper) - float(lower), 1.0e-12)
        deltas.append((float(sample_a[name]) - float(sample_b[name])) / span)
    return float(np.linalg.norm(np.asarray(deltas, dtype=float)))


def _order_samples_nearest_neighbor(
    samples: List[Dict[str, float]],
    start_sample: Dict[str, float],
    bounds: Dict[str, tuple[float, float]],
) -> List[Dict[str, float]]:
    remaining = [dict(sample) for sample in samples]
    ordered: List[Dict[str, float]] = []
    current = dict(start_sample)

    while remaining:
        next_index = min(
            range(len(remaining)),
            key=lambda idx: _normalized_distance(current, remaining[idx], bounds),
        )
        current = remaining.pop(next_index)
        ordered.append(current)
    return ordered


def _interpolate_sample_path(
    key_samples: List[Dict[str, float]],
    transition_steps: int,
) -> List[Dict[str, float]]:
    if not key_samples:
        return []
    if int(transition_steps) <= 1 or len(key_samples) == 1:
        return [dict(sample) for sample in key_samples]

    samples: List[Dict[str, float]] = []
    for start, end in zip(key_samples[:-1], key_samples[1:]):
        for step in range(int(transition_steps)):
            t = step / float(transition_steps)
            samples.append(
                {
                    name: float((1.0 - t) * float(start[name]) + t * float(end[name]))
                    for name in PUBLIC_PARAMETER_NAMES
                }
            )
    samples.append(dict(key_samples[-1]))
    return samples


def build_planform_sweep_samples(
    frame_count: int = 51,
    seed: int = 11,
    transition_steps: int = 1,
) -> List[Dict[str, float]]:
    space = build_cta_design_space()
    cta_view = space.cta_flat()
    bounds = space.bounds
    rng = np.random.default_rng(int(seed))
    if int(frame_count) <= 0:
        return []

    transition_steps = max(int(transition_steps), 1)
    key_frame_count = max(2, 1 + (int(frame_count) - 1) // transition_steps)
    latin_samples = _latin_hypercube_samples(bounds, key_frame_count - 1, rng)
    ordered_samples = _order_samples_nearest_neighbor(latin_samples, cta_view, bounds)
    key_samples = [dict(cta_view), *[{**cta_view, **sample} for sample in ordered_samples]]
    samples = _interpolate_sample_path(key_samples, transition_steps)
    if len(samples) < int(frame_count):
        samples.extend([dict(samples[-1])] * (int(frame_count) - len(samples)))
    return samples[: int(frame_count)]


def _render_frame(
    sample: Dict[str, float],
    cta_curves,
    frame_index: int,
    frame_count: int,
    output_path: Path,
) -> None:
    design = _public_sample_to_cta_design(sample)
    config = to_cta_model_config(design)
    _, dense_span, leading_edge_x, trailing_edge_x = _planform_curves(config)
    y_sections = np.asarray(config.topology.y_sections_array, dtype=float)

    ref_dense_span, ref_le, ref_te = cta_curves

    fig, ax = plt.subplots(figsize=(11.8, 8.2), constrained_layout=True)
    draw_planform(
        ax,
        dense_span,
        leading_edge_x,
        trailing_edge_x,
        y_sections,
        {},
        show_cargo=False,
        half_wing=True,
        label_zones=True,
        faded_symmetry=True,
        palette="dark_blue",
    )

    ax.plot(ref_le, ref_dense_span, color="#475569", linewidth=1.5, linestyle="--", alpha=0.85, zorder=6)
    ax.plot(ref_te, ref_dense_span, color="#475569", linewidth=1.5, linestyle="--", alpha=0.85, zorder=6)
    ax.text(
        0.02,
        0.98,
        _parameter_box_text(sample),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9.5,
        color="#0f172a",
        bbox={"boxstyle": "round,pad=0.34", "facecolor": "white", "edgecolor": "#cbd5e1", "alpha": 0.96},
        zorder=10,
    )
    ax.text(
        0.98,
        0.98,
        f"Frame {frame_index + 1}/{frame_count}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=10.0,
        color="#334155",
        bbox={"boxstyle": "round,pad=0.24", "facecolor": "white", "edgecolor": "#cbd5e1", "alpha": 0.96},
        zorder=10,
    )
    ax.set_title("CTA planform parameter sweep")
    ax.text(
        0.98,
        0.02,
        "Gray dashed: CTA case",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9.0,
        color="#475569",
        bbox={"boxstyle": "round,pad=0.20", "facecolor": "white", "edgecolor": "none", "alpha": 0.92},
        zorder=10,
    )
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def main() -> None:
    OUTPUT_GIF.parent.mkdir(parents=True, exist_ok=True)

    cta_design = build_cta_design()
    cta_config = to_cta_model_config(cta_design)
    _, ref_dense_span, ref_le, ref_te = _planform_curves(cta_config)
    cta_curves = (ref_dense_span, ref_le, ref_te)

    samples = build_planform_sweep_samples(frame_count=51, seed=11)
    durations = []
    images: List[Image.Image] = []

    with tempfile.TemporaryDirectory(prefix="cta_planform_gif_") as tmp_dir:
        tmp_dir_path = Path(tmp_dir)
        for idx, sample in enumerate(samples):
            frame_path = tmp_dir_path / f"frame_{idx:03d}.png"
            _render_frame(sample, cta_curves, idx, len(samples), frame_path)
            images.append(
                Image.open(frame_path).convert(
                    "P",
                    palette=Image.Palette.ADAPTIVE,
                    colors=255,
                    dither=Image.Dither.NONE,
                )
            )
            durations.append(260)
        if durations:
            durations[0] = 700
            durations[-1] = 700
        images[0].save(
            OUTPUT_GIF,
            save_all=True,
            append_images=images[1:],
            duration=durations,
            loop=0,
            disposal=2,
        )

    print(f"CTA planform GIF written to: {OUTPUT_GIF}")


if __name__ == "__main__":
    main()
