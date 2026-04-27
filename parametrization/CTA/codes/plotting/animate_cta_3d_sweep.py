from __future__ import annotations

import os
from pathlib import Path
import sys
import tempfile
from typing import List, Tuple

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

from parametrization.CTA.case import to_cta_model_config
from parametrization.CTA.codes.plotting.animate_cta_planform_parameter_sweep import (
    build_planform_sweep_samples,
    _public_sample_to_cta_design,
)
from parametrization.CTA.codes.plotting.plot_cta_views import draw_3d
from parametrization.bwb.builder import prepare_geometry


OUTPUT_GIF = CTA_DIR / "outputs" / "wing" / "cta_3d_sweep.gif"
GIF_DPI = 170
FRAME_COUNT = 51
SEED = 11


def _dense_span(config, n_points: int = 900) -> np.ndarray:
    y_sections = np.asarray(config.topology.y_sections_array, dtype=float)
    anchor_y = np.asarray(config.topology.anchor_y_array, dtype=float)
    dense = np.linspace(0.0, float(config.topology.span), int(n_points))
    return np.unique(np.concatenate([dense, y_sections, anchor_y]))


def _geometry_axes_limits(payloads) -> Tuple[Tuple[float, float], Tuple[float, float], Tuple[float, float]]:
    x_min = np.inf
    x_max = -np.inf
    y_max = 0.0
    z_min = np.inf
    z_max = -np.inf

    for payload in payloads:
        config = payload["config"]
        prepared = prepare_geometry(config)
        dense_span = _dense_span(config)
        vertical_center = np.interp(dense_span, prepared.loft.span_stations, prepared.loft.vertical_y)
        twist_dense = np.array([prepared.spanwise_laws.twist_deg(float(y)) for y in dense_span], dtype=float)
        leading_edge_x = np.interp(dense_span, prepared.loft.span_stations, prepared.loft.leading_edge_x)
        chord_dense = np.interp(dense_span, prepared.loft.span_stations, prepared.loft.chord)

        x_air = prepared.section_model.x_air
        for yy, le, chord, twist, v0 in zip(dense_span, leading_edge_x, chord_dense, twist_dense, vertical_center):
            y_up, y_lo, _ = prepared.section_model.coordinates_at_y(float(yy))
            x_local = x_air * float(chord)
            z_up_local = y_up * float(chord)
            z_lo_local = y_lo * float(chord)
            tw = np.deg2rad(float(twist))
            ctw = np.cos(tw)
            stw = np.sin(tw)

            xu = float(le) + x_local * ctw - z_up_local * stw
            xl = float(le) + x_local * ctw - z_lo_local * stw
            zu = float(v0) + x_local * stw + z_up_local * ctw
            zl = float(v0) + x_local * stw + z_lo_local * ctw

            x_min = min(x_min, float(np.min(xu)), float(np.min(xl)))
            x_max = max(x_max, float(np.max(xu)), float(np.max(xl)))
            y_max = max(y_max, abs(float(yy)))
            z_min = min(z_min, float(np.min(zu)), float(np.min(zl)))
            z_max = max(z_max, float(np.max(zu)), float(np.max(zl)))

    x_pad = 0.04 * max(x_max - x_min, 1.0)
    y_pad = 0.04 * max(2.0 * y_max, 1.0)
    z_pad = 0.10 * max(z_max - z_min, 1.0)
    return (
        (x_min - x_pad, x_max + x_pad),
        (-y_max - y_pad, y_max + y_pad),
        (z_min - z_pad, z_max + z_pad),
    )


def _build_payloads() -> List[dict]:
    payloads: List[dict] = []
    for sample in build_planform_sweep_samples(frame_count=FRAME_COUNT, seed=SEED):
        design = _public_sample_to_cta_design(sample)
        config = to_cta_model_config(design, use_cta_anchor_twist=True)
        payloads.append({"sample": dict(sample), "config": config})
    return payloads


def _frame_title(frame_index: int, frame_count: int) -> str:
    return f"CTA 3D variation sweep | Frame {frame_index + 1}/{frame_count}"


def _render_frame(
    payload,
    x_limits: Tuple[float, float],
    y_limits: Tuple[float, float],
    z_limits: Tuple[float, float],
    frame_index: int,
    frame_count: int,
    output_path: Path,
) -> None:
    config = payload["config"]
    sample = payload["sample"]
    prepared = prepare_geometry(config)
    dense_span = _dense_span(config)
    leading_edge_x = np.interp(dense_span, prepared.loft.span_stations, prepared.loft.leading_edge_x)
    chord_dense = np.interp(dense_span, prepared.loft.span_stations, prepared.loft.chord)

    fig = plt.figure(figsize=(12.8, 7.4), constrained_layout=False)
    fig.patch.set_facecolor("white")
    ax = fig.add_subplot(111, projection="3d")
    ax.set_position([0.02, 0.04, 0.96, 0.90])
    draw_3d(ax, prepared, dense_span, leading_edge_x, chord_dense)
    ax.set_xlim(*x_limits)
    ax.set_ylim(*y_limits)
    ax.set_zlim(*z_limits)
    ax.set_title(_frame_title(frame_index, frame_count), pad=8.0)

    info = "\n".join(
        (
            f"Bw={float(sample['span']):.2f} m",
            f"B2/Bw={float(sample['b2_span_ratio']):.3f}",
            f"C0={float(sample['c1_root_chord']):.2f} m",
            f"C3={float(sample['c2_c1_ratio']):.2f} m",
            f"C4/C3={float(sample['c4_c3_ratio']):.3f}",
            f"C5={float(sample['c4_c1_ratio']):.2f} m",
            f"S1={float(sample['s2_deg']):.2f} deg",
            f"S2={float(sample['s3_deg']):.2f} deg",
        )
    )
    fig.text(
        0.02,
        0.96,
        info,
        ha="left",
        va="top",
        fontsize=9.2,
        color="#0f172a",
        bbox={"boxstyle": "round,pad=0.32", "facecolor": "white", "edgecolor": "#cbd5e1", "alpha": 0.95},
    )

    fig.savefig(output_path, dpi=GIF_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    OUTPUT_GIF.parent.mkdir(parents=True, exist_ok=True)

    payloads = _build_payloads()
    x_limits, y_limits, z_limits = _geometry_axes_limits(payloads)

    durations: List[int] = []
    images: List[Image.Image] = []
    with tempfile.TemporaryDirectory(prefix="cta_3d_gif_") as tmp_dir:
        tmp_dir_path = Path(tmp_dir)
        for idx, payload in enumerate(payloads):
            frame_path = tmp_dir_path / f"frame_{idx:03d}.png"
            _render_frame(
                payload,
                x_limits,
                y_limits,
                z_limits,
                idx,
                len(payloads),
                frame_path,
            )
            with Image.open(frame_path) as image:
                images.append(
                    image.convert(
                        "P",
                        palette=Image.Palette.ADAPTIVE,
                        colors=255,
                        dither=Image.Dither.NONE,
                    )
                )
            durations.append(230)

        if durations:
            durations[0] = 900
            durations[-1] = 900

        images[0].save(
            OUTPUT_GIF,
            save_all=True,
            append_images=images[1:],
            duration=durations,
            loop=0,
            disposal=2,
        )

    print(f"CTA 3D GIF written to: {OUTPUT_GIF}")


if __name__ == "__main__":
    main()
