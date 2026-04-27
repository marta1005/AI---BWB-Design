from __future__ import annotations

import os
from pathlib import Path
import sys
import tempfile
from dataclasses import replace
from typing import Dict, List, Tuple

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
from parametrization.CTA.design_space import build_cta_design_space
from parametrization.CTA.codes.plotting.animate_cta_planform_parameter_sweep import (
    build_planform_sweep_samples,
)
from parametrization.CTA.codes.plotting.animate_cta_front_thickness_sweep import (
    _apply_thickness_scaling,
    _constraint_maps,
    _latin_hypercube_samples,
    _order_samples_nearest_neighbor,
)
from parametrization.CTA.codes.plotting.plot_cta_views import draw_3d
from parametrization.bwb.builder import prepare_geometry


OUTPUT_GIF = CTA_DIR / "outputs" / "wing" / "cta_3d_planform_twist_thickness_sweep.gif"
GIF_DPI = 170
FRAME_COUNT = 51
PLANFORM_SEED = 11
AUX_RNG_SEED = 17
GIF_UPPER_COLOR = "#1e3a8a"
GIF_LOWER_COLOR = "#2563eb"
GIF_CONTOUR_COLOR = "#0f172a"
GIF_LEADING_EDGE_COLOR = "#0f172a"
GIF_TRAILING_EDGE_COLOR = "#1d4ed8"
GIF_TIP_COLOR = "#334155"
GIF_ROOT_COLOR = "#64748b"
FIG_BG = "white"
TITLE_COLOR = "#0f172a"

PLANFORM_VARIABLES: Tuple[str, ...] = (
    "span",
    "b2_span_ratio",
    "c1_root_chord",
    "c2_c1_ratio",
    "c4_c3_ratio",
    "c4_c1_ratio",
    "s2_deg",
    "s3_deg",
)

TWIST_VARIABLES: Tuple[str, ...] = (
    "twist_c1_deg",
    "twist_c3_deg",
    "twist_c4_deg",
)

THICKNESS_VARIABLES: Tuple[str, ...] = (
    "RThickness_CentreBody",
    "RThickness_MidWing",
    "RThickness_OutWing",
    "RThickness_OutWing_Tip",
)

SWEEP_VARIABLES: Tuple[str, ...] = PLANFORM_VARIABLES + TWIST_VARIABLES + THICKNESS_VARIABLES


def _dense_span(config, n_points: int = 900) -> np.ndarray:
    y_sections = np.asarray(config.topology.y_sections_array, dtype=float)
    anchor_y = np.asarray(config.topology.anchor_y_array, dtype=float)
    dense = np.linspace(0.0, float(config.topology.span), int(n_points))
    return np.unique(np.concatenate([dense, y_sections, anchor_y]))


def _reference_sample(space, constraints_by_parameter: Dict[str, Dict[str, object]]) -> Dict[str, float]:
    flat = space.cta_flat()
    sample = {name: float(flat[name]) for name in PLANFORM_VARIABLES + TWIST_VARIABLES}
    for parameter, item in constraints_by_parameter.items():
        sample[parameter] = float(item["reference"])
    return sample


def _design_from_sample(space, sample: Dict[str, float]):
    base_sample = {name: float(sample[name]) for name in PLANFORM_VARIABLES + TWIST_VARIABLES}
    design = space.to_design(base_sample)
    return _apply_thickness_scaling(design, sample)


def _build_payload(sample: Dict[str, float], space):
    design = _design_from_sample(space, sample)
    config = to_cta_model_config(design, use_cta_anchor_twist=True)
    ref_flat = space.cta_flat()
    delta_c1 = float(sample["twist_c1_deg"]) - float(ref_flat["twist_c1_deg"])
    delta_c3 = float(sample["twist_c3_deg"]) - float(ref_flat["twist_c3_deg"])
    delta_c4 = float(sample["twist_c4_deg"]) - float(ref_flat["twist_c4_deg"])
    base_twist = list(config.spanwise.twist_deg.values)
    if len(base_twist) >= 6:
        for idx in range(4):
            base_twist[idx] = float(base_twist[idx]) + delta_c1
        base_twist[4] = float(base_twist[4]) + delta_c3
        base_twist[5] = float(base_twist[5]) + delta_c4
        config.spanwise.twist_deg = replace(config.spanwise.twist_deg, values=tuple(base_twist))
    return {"sample": dict(sample), "config": config}


def _build_samples(space) -> List[Dict[str, float]]:
    constraints_by_parameter, _ = _constraint_maps()
    aux_bounds: Dict[str, Tuple[float, float]] = {
        name: tuple(float(value) for value in space.bounds[name]) for name in TWIST_VARIABLES
    }
    for parameter, item in constraints_by_parameter.items():
        aux_bounds[parameter] = (float(item["lower_bound"]), float(item["upper_bound"]))

    rng = np.random.default_rng(int(AUX_RNG_SEED))
    ref_sample = _reference_sample(space, constraints_by_parameter)
    planform_samples = build_planform_sweep_samples(frame_count=FRAME_COUNT, seed=PLANFORM_SEED)

    aux_ref = {
        name: float(ref_sample[name]) for name in TWIST_VARIABLES + THICKNESS_VARIABLES
    }
    candidate_aux = _latin_hypercube_samples(aux_bounds, FRAME_COUNT - 1, rng)
    ordered_aux = _order_samples_nearest_neighbor(
        candidate_aux,
        aux_ref,
        aux_bounds,
        TWIST_VARIABLES + THICKNESS_VARIABLES,
    )
    aux_samples = [aux_ref, *ordered_aux[: FRAME_COUNT - 1]]

    combined: List[Dict[str, float]] = []
    for planform_sample, aux_sample in zip(planform_samples, aux_samples):
        merged = dict(planform_sample)
        for name in TWIST_VARIABLES + THICKNESS_VARIABLES:
            merged[name] = float(aux_sample[name])
        combined.append(merged)
    return combined


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


def _frame_title(frame_index: int, frame_count: int) -> str:
    return f"CTA 3D | Frame {frame_index + 1}/{frame_count}"


def _frame_info(sample: Dict[str, float]) -> str:
    return "\n".join(
        (
            f"Bw={float(sample['span']):.2f} m",
            f"B2/Bw={float(sample['b2_span_ratio']):.3f}",
            f"C0={float(sample['c1_root_chord']):.2f} m",
            f"C3={float(sample['c2_c1_ratio']):.2f}",
            f"C4/C3={float(sample['c4_c3_ratio']):.3f}",
            f"C5={float(sample['c4_c1_ratio']):.2f} m",
            f"S1={float(sample['s2_deg']):.2f} deg",
            f"S2={float(sample['s3_deg']):.2f} deg",
            f"tw C0={float(sample['twist_c1_deg']):.2f} deg",
            f"tw C3={float(sample['twist_c3_deg']):.2f} deg",
            f"tw C5={float(sample['twist_c4_deg']):.2f} deg",
            f"t/c C0={float(sample['RThickness_CentreBody']):.3f}",
            f"t/c C3={float(sample['RThickness_MidWing']):.3f}",
            f"t/c C4={float(sample['RThickness_OutWing']):.3f}",
            f"t/c C5={float(sample['RThickness_OutWing_Tip']):.3f}",
        )
    )


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
    fig.patch.set_facecolor(FIG_BG)
    ax = fig.add_subplot(111, projection="3d")
    ax.set_position([0.02, 0.04, 0.96, 0.90])
    draw_3d(
        ax,
        prepared,
        dense_span,
        leading_edge_x,
        chord_dense,
        upper_color=GIF_UPPER_COLOR,
        lower_color=GIF_LOWER_COLOR,
        contour_color=GIF_CONTOUR_COLOR,
        leading_edge_color=GIF_LEADING_EDGE_COLOR,
        trailing_edge_color=GIF_TRAILING_EDGE_COLOR,
        tip_color=GIF_TIP_COLOR,
        root_color=GIF_ROOT_COLOR,
    )
    ax.set_xlim(*x_limits)
    ax.set_ylim(*y_limits)
    ax.set_zlim(*z_limits)
    ax.view_init(elev=23, azim=-138)
    try:
        ax.set_proj_type("persp", focal_length=1.12)
    except Exception:
        pass
    ax.set_title(
        _frame_title(frame_index, frame_count),
        pad=10.0,
        fontsize=18.0,
        color=TITLE_COLOR,
        fontweight="semibold",
    )
    fig.text(
        0.02,
        0.96,
        _frame_info(sample),
        ha="left",
        va="top",
        fontsize=9.2,
        color="#0f172a",
        bbox={"boxstyle": "round,pad=0.32", "facecolor": "white", "edgecolor": "#cbd5e1", "alpha": 0.95},
    )

    fig.savefig(output_path, dpi=GIF_DPI, bbox_inches="tight", facecolor=FIG_BG)
    plt.close(fig)


def main() -> None:
    OUTPUT_GIF.parent.mkdir(parents=True, exist_ok=True)
    space = build_cta_design_space()
    samples = _build_samples(space)
    payloads = [_build_payload(sample, space) for sample in samples]
    x_limits, y_limits, z_limits = _geometry_axes_limits(payloads)

    durations: List[int] = []
    images: List[Image.Image] = []
    with tempfile.TemporaryDirectory(prefix="cta_3d_twist_thickness_gif_") as tmp_dir:
        tmp_dir_path = Path(tmp_dir)
        for idx, payload in enumerate(payloads):
            frame_path = tmp_dir_path / f"frame_{idx:03d}.png"
            _render_frame(payload, x_limits, y_limits, z_limits, idx, len(payloads), frame_path)
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

    print(f"CTA 3D planform+twist+thickness GIF written to: {OUTPUT_GIF}")


if __name__ == "__main__":
    main()
