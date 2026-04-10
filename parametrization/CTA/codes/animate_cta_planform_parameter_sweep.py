from __future__ import annotations

import os
from dataclasses import replace
from pathlib import Path
import sys
import tempfile
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
from PIL import Image

from parametrization.CTA.design_space import build_cta_design_space
from parametrization.CTA.reference import build_reference_design, to_cta_model_config
from parametrization.CTA.codes.plot_cta_reference_views import draw_planform
from parametrization.bwb.planform import build_sectioned_bwb_planform


OUTPUT_GIF = SCRIPT_DIR.parent / "outputs" / "reference" / "cta_planform_parameter_sweep.gif"


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


def _public_sample_to_reference_design(sample: Dict[str, float]):
    reference_design = build_reference_design()
    b1_fixed_m = float(reference_design.span * reference_design.b1_span_ratio)
    wing_span_m = float(sample["span"])
    total_span_m = float(b1_fixed_m + wing_span_m)
    b2_m = float(sample["b2_span_ratio"]) * wing_span_m
    b3_m = float(wing_span_m - b2_m)
    c0_m = float(sample["c1_root_chord"])
    c3_m = float(sample["c2_c1_ratio"])
    c4_m = float(sample["c2_c1_ratio"]) * float(sample["c4_c3_ratio"])
    c5_m = float(sample["c4_c1_ratio"])
    return replace(
        reference_design,
        span=total_span_m,
        b1_span_ratio=b1_fixed_m / total_span_m,
        b2_span_ratio=b2_m / total_span_m,
        b3_span_ratio=b3_m / total_span_m,
        c1_root_chord=c0_m,
        c2_c1_ratio=c3_m / c0_m,
        c3_c1_ratio=c4_m / c0_m,
        c4_c1_ratio=c5_m / c0_m,
    )


PUBLIC_PARAMETER_NAMES = (
    "span",
    "b2_span_ratio",
    "c1_root_chord",
    "c2_c1_ratio",
    "c4_c3_ratio",
    "c4_c1_ratio",
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
        )
    )


def _jittered_random_series(
    lower: float,
    upper: float,
    frame_count: int,
    rng: np.random.Generator,
) -> np.ndarray:
    if frame_count <= 1:
        return np.asarray([0.5 * (float(lower) + float(upper))], dtype=float)
    edges = np.linspace(float(lower), float(upper), int(frame_count) + 1, dtype=float)
    values = np.array(
        [rng.uniform(edges[idx], edges[idx + 1]) for idx in range(int(frame_count))],
        dtype=float,
    )
    rng.shuffle(values)
    return values


def build_random_samples(
    frame_count: int = 51,
    seed: int = 11,
) -> List[Dict[str, float]]:
    space = build_cta_design_space()
    rng = np.random.default_rng(int(seed))
    bounds = space.bounds
    series = {
        name: _jittered_random_series(bounds[name][0], bounds[name][1], frame_count, rng)
        for name in PUBLIC_PARAMETER_NAMES
    }
    return [
        {name: float(series[name][frame_idx]) for name in PUBLIC_PARAMETER_NAMES}
        for frame_idx in range(int(frame_count))
    ]


def _clip_sample_to_bounds(sample: Dict[str, float], bounds: Dict[str, tuple[float, float]]) -> Dict[str, float]:
    clipped: Dict[str, float] = {}
    for name in PUBLIC_PARAMETER_NAMES:
        lower, upper = bounds[name]
        clipped[name] = float(np.clip(float(sample[name]), float(lower), float(upper)))
    return clipped


def _interpolate_samples(
    start: Dict[str, float],
    end: Dict[str, float],
    interval_frames: int,
) -> List[Dict[str, float]]:
    if interval_frames <= 0:
        return [dict(start)]
    frames: List[Dict[str, float]] = []
    for step in range(interval_frames):
        t = step / float(interval_frames)
        frames.append(
            {
                name: float((1.0 - t) * float(start[name]) + t * float(end[name]))
                for name in PUBLIC_PARAMETER_NAMES
            }
        )
    return frames


def _advance_focused_anchor_samples(
    reference: Dict[str, float],
    bounds: Dict[str, tuple[float, float]],
    rng: np.random.Generator,
) -> List[Dict[str, float]]:
    def random_anchor(local_scale: float = 0.42) -> Dict[str, float]:
        sample = dict(reference)
        for name in PUBLIC_PARAMETER_NAMES:
            lower, upper = bounds[name]
            span = float(upper) - float(lower)
            ref = float(reference[name])
            local_half = float(local_scale) * span
            sample[name] = float(
                np.clip(
                    rng.uniform(ref - local_half, ref + local_half),
                    float(lower),
                    float(upper),
                )
            )
        return sample

    forward_ref = dict(reference)
    forward_ref.update(
        {
            "span": float(reference["span"]),
            "b2_span_ratio": float(reference["b2_span_ratio"]),
            "c2_c1_ratio": float(bounds["c2_c1_ratio"][1]),
            "c4_c3_ratio": float(bounds["c4_c3_ratio"][1]),
            "c4_c1_ratio": float(bounds["c4_c1_ratio"][0]),
        }
    )

    forward_b2 = dict(reference)
    forward_b2.update(
        {
            "span": float(reference["span"]),
            "b2_span_ratio": 0.35 * float(bounds["b2_span_ratio"][0]) + 0.65 * float(bounds["b2_span_ratio"][1]),
            "c2_c1_ratio": float(bounds["c2_c1_ratio"][1]),
            "c4_c3_ratio": float(bounds["c4_c3_ratio"][1]),
            "c4_c1_ratio": float(bounds["c4_c1_ratio"][0]),
        }
    )

    forward_max = dict(reference)
    forward_max.update(
        {
            "span": float(bounds["span"][1]),
            "b2_span_ratio": float(bounds["b2_span_ratio"][1]),
            "c2_c1_ratio": float(bounds["c2_c1_ratio"][1]),
            "c4_c3_ratio": float(bounds["c4_c3_ratio"][1]),
            "c4_c1_ratio": float(bounds["c4_c1_ratio"][0]),
        }
    )

    forward_soft = dict(reference)
    forward_soft.update(
        {
            "span": float(reference["span"]),
            "b2_span_ratio": float(reference["b2_span_ratio"]),
            "c2_c1_ratio": 0.65 * float(bounds["c2_c1_ratio"][1]) + 0.35 * float(reference["c2_c1_ratio"]),
            "c4_c3_ratio": 0.75 * float(bounds["c4_c3_ratio"][1]) + 0.25 * float(reference["c4_c3_ratio"]),
            "c4_c1_ratio": float(bounds["c4_c1_ratio"][0]),
        }
    )

    aft_max = dict(reference)
    aft_max.update(
        {
            "span": float(bounds["span"][0]),
            "b2_span_ratio": float(bounds["b2_span_ratio"][0]),
            "c2_c1_ratio": float(bounds["c2_c1_ratio"][0]),
            "c4_c3_ratio": float(bounds["c4_c3_ratio"][0]),
            "c4_c1_ratio": float(bounds["c4_c1_ratio"][1]),
        }
    )

    anchors = [
        dict(reference),
        random_anchor(0.38),
        forward_ref,
        random_anchor(0.45),
        aft_max,
        random_anchor(0.40),
        forward_b2,
        random_anchor(0.46),
        forward_max,
        random_anchor(0.42),
        forward_soft,
        random_anchor(0.50),
        aft_max,
        random_anchor(0.38),
        forward_ref,
        random_anchor(0.44),
        forward_max,
        random_anchor(0.36),
        forward_b2,
        random_anchor(0.48),
        aft_max,
        random_anchor(0.40),
        forward_soft,
        random_anchor(0.44),
        forward_ref,
        dict(reference),
    ]
    return [_clip_sample_to_bounds(sample, bounds) for sample in anchors]


def build_planform_sweep_samples(frame_count: int = 51, seed: int = 11) -> List[Dict[str, float]]:
    space = build_cta_design_space()
    reference = space.reference_flat()
    bounds = space.bounds
    rng = np.random.default_rng(int(seed))

    anchors = _advance_focused_anchor_samples(reference, bounds, rng)
    interval_count = max(len(anchors) - 1, 1)
    base_interval = 2
    total_frames = interval_count * base_interval + 1
    if total_frames != int(frame_count):
        raise ValueError(
            f"Anchor layout expects {total_frames} frames with interval {base_interval}, got {frame_count}"
        )

    samples: List[Dict[str, float]] = []
    for idx in range(interval_count):
        samples.extend(_interpolate_samples(anchors[idx], anchors[idx + 1], base_interval))
    samples.append(dict(anchors[-1]))
    samples = samples[: int(frame_count)]
    if samples:
        samples[0] = dict(reference)
    return samples


def _render_frame(
    sample: Dict[str, float],
    reference_curves,
    frame_index: int,
    frame_count: int,
    output_path: Path,
) -> None:
    design = _public_sample_to_reference_design(sample)
    config = to_cta_model_config(design)
    _, dense_span, leading_edge_x, trailing_edge_x = _planform_curves(config)
    y_sections = np.asarray(config.topology.y_sections_array, dtype=float)

    ref_dense_span, ref_le, ref_te = reference_curves

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
        "Gray dashed: reference",
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

    reference_design = build_reference_design()
    reference_config = to_cta_model_config(reference_design)
    _, ref_dense_span, ref_le, ref_te = _planform_curves(reference_config)
    reference_curves = (ref_dense_span, ref_le, ref_te)

    samples = build_planform_sweep_samples(frame_count=51, seed=11)
    durations = []
    images: List[Image.Image] = []

    with tempfile.TemporaryDirectory(prefix="cta_planform_gif_") as tmp_dir:
        tmp_dir_path = Path(tmp_dir)
        for idx, sample in enumerate(samples):
            frame_path = tmp_dir_path / f"frame_{idx:03d}.png"
            _render_frame(sample, reference_curves, idx, len(samples), frame_path)
            images.append(Image.open(frame_path).convert("P", palette=Image.Palette.ADAPTIVE))
            durations.append(140)
        if durations:
            durations[0] = 420
            durations[-1] = 420
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
