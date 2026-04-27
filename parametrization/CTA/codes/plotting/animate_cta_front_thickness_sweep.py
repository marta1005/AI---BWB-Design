from __future__ import annotations

import os
from dataclasses import replace
from pathlib import Path
import sys
import tempfile
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
from parametrization.CTA.design_space import build_cta_design_space, cta_thickness_constraints
from parametrization.CTA.codes.plotting.plot_cta_views import front_linear_envelope
from parametrization.bwb.builder import prepare_geometry
from parametrization.bwb.design_variables import SectionedBWBDesignVariables
from parametrization.shared.cst import (
    KulfanCSTAirfoil,
    cosine_spacing,
    fit_kulfan_airfoil_coefficients,
)


OUTPUT_GIF = CTA_DIR / "outputs" / "wing" / "cta_front_thickness_sweep.gif"
FRAME_COUNT = 11
RNG_SEED = 13
DENSE_SPAN_POINTS = 900
GIF_DPI = 135

REDUCED_ACTIVE_VARIABLES: Tuple[str, ...] = (
    "span",
    "c1_root_chord",
    "c2_c1_ratio",
    "c4_c3_ratio",
    "b2_span_ratio",
    "c4_c1_ratio",
    "s2_deg",
    "s3_deg",
    "med_3_te_sweep_deg",
    "twist_c1_deg",
    "twist_c3_deg",
    "twist_c4_deg",
)

ORDERING_VARIABLES: Tuple[str, ...] = REDUCED_ACTIVE_VARIABLES

THICKNESS_PREFIX_MAP: Dict[str, str] = {
    "RThickness_CentreBody": "c1",
    "RThickness_MidWing": "c2",
    "RThickness_OutWing": "c3",
    "RThickness_OutWing_Tip": "c4",
}

STATION_ORDER: Tuple[str, ...] = ("C0", "C3", "C4", "C5")
STATION_COLORS: Dict[str, str] = {
    "C0": "#b45309",
    "C3": "#2563eb",
    "C4": "#16a34a",
    "C5": "#15803d",
}

GIF_ZONE_COLORS: Tuple[str, str, str] = ("#1e3a8a", "#1d4ed8", "#2563eb")
GIF_FILL_COLOR = "#2563eb"
GIF_LINE_COLOR = "#0f172a"
GIF_REF_COLOR = "#64748b"


def _dense_span(config, n_points: int = DENSE_SPAN_POINTS) -> np.ndarray:
    y_sections = np.asarray(config.topology.y_sections_array, dtype=float)
    anchor_y = np.asarray(config.topology.anchor_y_array, dtype=float)
    dense = np.linspace(0.0, float(config.topology.span), int(n_points))
    return np.unique(np.concatenate([dense, y_sections, anchor_y]))


def _public_station_y(config) -> Tuple[Tuple[str, float], ...]:
    anchor_y = np.asarray(config.topology.anchor_y_array, dtype=float)
    return (
        ("C0", float(anchor_y[0])),
        ("C3", float(anchor_y[3])),
        ("C4", float(anchor_y[4])),
        ("C5", float(anchor_y[5])),
    )


def _actual_tc_map(prepared, config) -> Dict[str, float]:
    values: Dict[str, float] = {}
    for label, yy in _public_station_y(config):
        metrics, _ = prepared.section_model.geometry_metrics_at_y(float(yy))
        values[label] = float(metrics.max_tc)
    return values


def _constraint_maps() -> Tuple[Dict[str, Dict[str, object]], Dict[str, Dict[str, object]]]:
    by_parameter = {str(item["parameter"]): dict(item) for item in cta_thickness_constraints()}
    by_station = {str(item["station"]): dict(item) for item in cta_thickness_constraints()}
    return by_parameter, by_station


def _latin_hypercube_samples(
    bounds: Dict[str, Tuple[float, float]],
    sample_count: int,
    rng: np.random.Generator,
) -> List[Dict[str, float]]:
    if sample_count <= 0:
        return []
    samples: List[Dict[str, float]] = [{} for _ in range(int(sample_count))]
    for name, (lower, upper) in bounds.items():
        unit_values = (rng.permutation(sample_count) + rng.random(sample_count)) / float(sample_count)
        values = float(lower) + unit_values * (float(upper) - float(lower))
        for idx, value in enumerate(values):
            samples[idx][name] = float(value)
    return samples


def _normalized_distance(
    sample_a: Dict[str, float],
    sample_b: Dict[str, float],
    bounds: Dict[str, Tuple[float, float]],
    variable_names: Tuple[str, ...],
) -> float:
    deltas = []
    for name in variable_names:
        lower, upper = bounds[name]
        span = max(float(upper) - float(lower), 1.0e-12)
        deltas.append((float(sample_a[name]) - float(sample_b[name])) / span)
    return float(np.linalg.norm(np.asarray(deltas, dtype=float)))


def _order_samples_nearest_neighbor(
    samples: List[Dict[str, float]],
    start_sample: Dict[str, float],
    bounds: Dict[str, Tuple[float, float]],
    variable_names: Tuple[str, ...],
) -> List[Dict[str, float]]:
    remaining = [dict(sample) for sample in samples]
    ordered: List[Dict[str, float]] = []
    current = dict(start_sample)

    while remaining:
        next_index = min(
            range(len(remaining)),
            key=lambda idx: _normalized_distance(current, remaining[idx], bounds, variable_names),
        )
        current = remaining.pop(next_index)
        ordered.append(current)
    return ordered


def _section_tc_current(
    design: SectionedBWBDesignVariables,
    section_prefix: str,
) -> float:
    upper = np.asarray(getattr(design, f"{section_prefix}_upper_cst"), dtype=float)
    lower = np.asarray(getattr(design, f"{section_prefix}_lower_cst"), dtype=float)
    te_thickness = float(getattr(design, f"{section_prefix}_te_thickness"))
    shape = KulfanCSTAirfoil(
        degree=len(upper) - 1,
        n1=float(design.cst_n1),
        n2=float(design.cst_n2),
        x_tc_window=(0.15, 0.65),
        shared_leading_edge=False,
    )
    x_air = cosine_spacing(241)
    coeffs = np.concatenate([upper, lower], dtype=float)
    yu, yl = shape.evaluate(x_air, coeffs, te_thickness=te_thickness)
    return float(np.max(yu - yl))


def _scale_section_coefficients_to_tc(
    design: SectionedBWBDesignVariables,
    section_prefix: str,
    target_tc: float,
) -> SectionedBWBDesignVariables:
    upper = np.asarray(getattr(design, f"{section_prefix}_upper_cst"), dtype=float)
    lower = np.asarray(getattr(design, f"{section_prefix}_lower_cst"), dtype=float)
    te_thickness = float(getattr(design, f"{section_prefix}_te_thickness"))
    shape = KulfanCSTAirfoil(
        degree=len(upper) - 1,
        n1=float(design.cst_n1),
        n2=float(design.cst_n2),
        x_tc_window=(0.15, 0.65),
        shared_leading_edge=False,
    )
    x_air = cosine_spacing(241)
    coeffs = np.concatenate([upper, lower], dtype=float)
    yu, yl = shape.evaluate(x_air, coeffs, te_thickness=te_thickness)
    current_tc = float(np.max(yu - yl))
    if current_tc <= 1.0e-12:
        return design

    scale = float(target_tc) / float(current_tc)
    camber = 0.5 * (yu + yl)
    thickness = (yu - yl) * scale
    new_yu = camber + 0.5 * thickness
    new_yl = camber - 0.5 * thickness
    new_te = te_thickness * scale
    new_upper, new_lower = fit_kulfan_airfoil_coefficients(
        x_air,
        new_yu,
        new_yl,
        degree=len(upper) - 1,
        n1=float(design.cst_n1),
        n2=float(design.cst_n2),
        te_thickness=float(new_te),
        regularization=1.0e-8,
        smoothness_weight=0.0,
        shared_leading_edge=False,
    )
    return replace(
        design,
        **{
            f"{section_prefix}_upper_cst": tuple(float(value) for value in new_upper),
            f"{section_prefix}_lower_cst": tuple(float(value) for value in new_lower),
            f"{section_prefix}_tc_max": float(target_tc),
            f"{section_prefix}_te_thickness": float(new_te),
        },
    )


def _apply_thickness_scaling(
    design: SectionedBWBDesignVariables,
    sample: Dict[str, float],
) -> SectionedBWBDesignVariables:
    updated = design
    for parameter, prefix in THICKNESS_PREFIX_MAP.items():
        updated = _scale_section_coefficients_to_tc(updated, prefix, float(sample[parameter]))
    return updated


def _design_from_sample(space, sample: Dict[str, float]) -> SectionedBWBDesignVariables:
    base_sample = {name: float(sample[name]) for name in REDUCED_ACTIVE_VARIABLES}
    design = space.to_design(base_sample)
    return _apply_thickness_scaling(design, sample)


def _reference_sample(
    space,
    constraints_by_parameter: Dict[str, Dict[str, object]],
) -> Dict[str, float]:
    flat = space.cta_flat()
    sample = {name: float(flat[name]) for name in REDUCED_ACTIVE_VARIABLES}
    for parameter, item in constraints_by_parameter.items():
        sample[parameter] = float(item["reference"])
    return sample


def _build_frame_payload(
    design: SectionedBWBDesignVariables,
    sample: Dict[str, float],
) -> Dict[str, object]:
    config = to_cta_model_config(design)
    prepared = prepare_geometry(config)
    dense_span = _dense_span(config)
    anchor_y = np.asarray(config.topology.anchor_y_array, dtype=float)
    upper, lower, _, _ = front_linear_envelope(prepared, dense_span, anchor_y)
    return {
        "sample": dict(sample),
        "config": config,
        "dense_span": dense_span,
        "upper": upper,
        "lower": lower,
        "actual_tc": _actual_tc_map(prepared, config),
    }


def _build_payloads(space) -> List[Dict[str, object]]:
    constraints_by_parameter, _ = _constraint_maps()
    rng = np.random.default_rng(int(RNG_SEED))
    bounds: Dict[str, Tuple[float, float]] = {
        name: tuple(float(value) for value in space.bounds[name]) for name in REDUCED_ACTIVE_VARIABLES
    }
    for parameter, item in constraints_by_parameter.items():
        bounds[parameter] = (float(item["lower_bound"]), float(item["upper_bound"]))

    payloads: List[Dict[str, object]] = []
    ref_sample = _reference_sample(space, constraints_by_parameter)
    ref_design = _design_from_sample(space, ref_sample)
    payloads.append(_build_frame_payload(ref_design, ref_sample))

    candidate_samples = _latin_hypercube_samples(bounds, max(8 * FRAME_COUNT, 40), rng)
    if not candidate_samples:
        return payloads

    ordered_samples = _order_samples_nearest_neighbor(
        candidate_samples,
        ref_sample,
        bounds,
        ORDERING_VARIABLES,
    )
    for sample in ordered_samples:
        try:
            design = _design_from_sample(space, sample)
            payloads.append(_build_frame_payload(design, sample))
        except ValueError:
            continue
        if len(payloads) >= FRAME_COUNT:
            break
    return payloads[:FRAME_COUNT]


def _status(actual: float, lower: float, upper: float) -> str:
    if actual < lower - 1.0e-9:
        return "LOW"
    if actual > upper + 1.0e-9:
        return "HIGH"
    return "OK"


def _frame_header(frame_index: int, frame_count: int) -> str:
    return (
        f"Frame {frame_index + 1}/{frame_count}\n"
        "Latin hypercube sweep\n"
        "Profiles fixed in shape\n"
        "Thickness scales them by region"
    )


def _thickness_target_block(
    sample: Dict[str, float],
    constraints_by_station: Dict[str, Dict[str, object]],
) -> str:
    lines = ["Target thickness scaling (t/c)"]
    parameter_for_station = {
        "C0": "RThickness_CentreBody",
        "C3": "RThickness_MidWing",
        "C4": "RThickness_OutWing",
        "C5": "RThickness_OutWing_Tip",
    }
    for station in STATION_ORDER:
        parameter = parameter_for_station[station]
        lower = float(constraints_by_station[station]["lower_bound"])
        upper = float(constraints_by_station[station]["upper_bound"])
        lines.append(f"{station}: {float(sample[parameter]):.4f} [{lower:.3f}, {upper:.3f}]")
    return "\n".join(lines)


def _thickness_actual_block(
    actual_tc: Dict[str, float],
    constraints_by_station: Dict[str, Dict[str, object]],
) -> str:
    lines = ["Actual geometry thickness (t/c)"]
    for station in STATION_ORDER:
        item = constraints_by_station[station]
        actual = float(actual_tc[station])
        lines.append(
            f"{station}: {actual:.4f} | ref {float(item['reference']):.4f} | {_status(actual, float(item['lower_bound']), float(item['upper_bound']))}"
        )
    return "\n".join(lines)


def _public_value_block(sample: Dict[str, float]) -> str:
    c4_m = float(sample["c2_c1_ratio"]) * float(sample["c4_c3_ratio"])
    return "\n".join(
        (
            "Global sampled variables",
            f"Bw = {float(sample['span']):.2f} m",
            f"B2/Bw = {float(sample['b2_span_ratio']):.3f}",
            f"C0 = {float(sample['c1_root_chord']):.2f} m",
            f"C3 = {float(sample['c2_c1_ratio']):.2f} m",
            f"C4 = {c4_m:.2f} m",
            f"C5 = {float(sample['c4_c1_ratio']):.2f} m",
            f"S1 = {float(sample['s2_deg']):.2f} deg",
            f"S2 = {float(sample['s3_deg']):.2f} deg",
            f"tw_C0/C3 = {float(sample['twist_c1_deg']):.2f} deg",
            f"tw_C4 = {float(sample['twist_c3_deg']):.2f} deg",
            f"tw_C5 = {float(sample['twist_c4_deg']):.2f} deg",
        )
    )


def _station_markers(ax, dense_span, upper, lower, station_y: Tuple[Tuple[str, float], ...]) -> None:
    z_span = max(1.0e-9, float(np.max(upper) - np.min(lower)))
    for idx, (station, yy) in enumerate(station_y):
        upper_here = float(np.interp(float(yy), dense_span, upper))
        lower_here = float(np.interp(float(yy), dense_span, lower))
        color = STATION_COLORS[station]
        ax.plot([yy, yy], [lower_here, upper_here], color=color, linewidth=1.8, zorder=6)
        ax.scatter([yy, yy], [lower_here, upper_here], color=color, s=16.0, zorder=7)
        ax.text(
            float(yy) + (0.18 if idx == 0 else 0.0),
            upper_here + 0.018 * z_span,
            station,
            ha="center",
            va="bottom",
            fontsize=8.6,
            color="#0f172a",
            bbox={"boxstyle": "round,pad=0.14", "facecolor": "white", "edgecolor": "none", "alpha": 0.92},
            zorder=8,
        )


def _render_frame(
    payload: Dict[str, object],
    reference_payload: Dict[str, object],
    constraints_by_station: Dict[str, Dict[str, object]],
    z_limits: Tuple[float, float],
    frame_index: int,
    frame_count: int,
    output_path: Path,
) -> None:
    sample = dict(payload["sample"])
    config = payload["config"]
    dense_span = np.asarray(payload["dense_span"], dtype=float)
    upper = np.asarray(payload["upper"], dtype=float)
    lower = np.asarray(payload["lower"], dtype=float)
    actual_tc = dict(payload["actual_tc"])

    ref_span = np.asarray(reference_payload["dense_span"], dtype=float)
    ref_upper = np.asarray(reference_payload["upper"], dtype=float)
    ref_lower = np.asarray(reference_payload["lower"], dtype=float)

    y_sections = np.asarray(config.topology.y_sections_array, dtype=float)
    span = float(np.max(dense_span))
    b1 = float(y_sections[1])
    b2_end = float(y_sections[2])

    fig = plt.figure(figsize=(13.4, 6.1), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[2.8, 1.2])
    ax = fig.add_subplot(gs[0, 0])
    ax_info = fig.add_subplot(gs[0, 1])

    for x0, x1, color in (
        (0.0, b1, GIF_ZONE_COLORS[0]),
        (b1, b2_end, GIF_ZONE_COLORS[1]),
        (b2_end, span, GIF_ZONE_COLORS[2]),
    ):
        ax.axvspan(x0, x1, color=color, alpha=0.10, zorder=0)

    ax.fill_between(dense_span, lower, upper, color=GIF_FILL_COLOR, zorder=1, alpha=0.18)
    ax.plot(dense_span, upper, color=GIF_LINE_COLOR, linewidth=2.25, zorder=3)
    ax.plot(dense_span, lower, color=GIF_LINE_COLOR, linewidth=2.25, zorder=3)
    ax.plot(ref_span, ref_upper, color=GIF_REF_COLOR, linewidth=1.35, linestyle="--", zorder=4, alpha=0.95)
    ax.plot(ref_span, ref_lower, color=GIF_REF_COLOR, linewidth=1.35, linestyle="--", zorder=4, alpha=0.95)

    for x_value in np.asarray(config.topology.anchor_y_array, dtype=float):
        ax.axvline(float(x_value), color="#64748b", linewidth=0.85, linestyle=(0, (4, 4)), alpha=0.45, zorder=2)

    _station_markers(ax, dense_span, upper, lower, _public_station_y(config))

    ax.axhline(0.0, color="#94a3b8", linewidth=0.9, linestyle=":")
    ax.set_xlim(0.0, span + 0.8)
    ax.set_ylim(float(z_limits[0]), float(z_limits[1]))
    ax.set_xlabel("spanwise y [m]")
    ax.set_ylabel("vertical z [m]")
    ax.set_title("CTA front Latin hypercube sweep")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, linewidth=0.35, alpha=0.30)

    ax.text(
        0.98,
        0.02,
        "Gray dashed: reference geometry",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8.6,
        color="#475569",
        bbox={"boxstyle": "round,pad=0.20", "facecolor": "white", "edgecolor": "none", "alpha": 0.92},
        zorder=10,
    )

    ax_info.axis("off")
    ax_info.text(
        0.02,
        0.98,
        _frame_header(frame_index, frame_count),
        ha="left",
        va="top",
        fontsize=10.6,
        color="#0f172a",
        linespacing=1.35,
        bbox={"boxstyle": "round,pad=0.34", "facecolor": "white", "edgecolor": "#cbd5e1", "alpha": 0.96},
    )
    ax_info.text(
        0.02,
        0.67,
        _public_value_block(sample),
        ha="left",
        va="top",
        fontsize=9.4,
        color="#0f172a",
        linespacing=1.28,
        bbox={"boxstyle": "round,pad=0.32", "facecolor": "#f8fafc", "edgecolor": "#cbd5e1", "alpha": 0.97},
    )
    ax_info.text(
        0.02,
        0.34,
        _thickness_target_block(sample, constraints_by_station),
        ha="left",
        va="top",
        fontsize=9.4,
        color="#0f172a",
        linespacing=1.30,
        bbox={"boxstyle": "round,pad=0.32", "facecolor": "#eff6ff", "edgecolor": "#bfdbfe", "alpha": 0.97},
    )
    ax_info.text(
        0.02,
        0.10,
        _thickness_actual_block(actual_tc, constraints_by_station),
        ha="left",
        va="bottom",
        fontsize=9.2,
        color="#0f172a",
        linespacing=1.30,
        bbox={"boxstyle": "round,pad=0.32", "facecolor": "#ecfeff", "edgecolor": "#a5f3fc", "alpha": 0.97},
    )

    fig.savefig(output_path, dpi=GIF_DPI)
    plt.close(fig)


def main() -> None:
    OUTPUT_GIF.parent.mkdir(parents=True, exist_ok=True)

    design_space = build_cta_design_space()
    _, constraints_by_station = _constraint_maps()
    payloads = _build_payloads(design_space)
    reference_payload = payloads[0]

    all_upper = np.concatenate([np.asarray(payload["upper"], dtype=float) for payload in payloads])
    all_lower = np.concatenate([np.asarray(payload["lower"], dtype=float) for payload in payloads])
    z_span = max(1.0e-9, float(np.max(all_upper) - np.min(all_lower)))
    z_limits = (
        float(np.min(all_lower) - 0.10 * z_span),
        float(np.max(all_upper) + 0.10 * z_span),
    )

    durations: List[int] = []
    images: List[Image.Image] = []
    with tempfile.TemporaryDirectory(prefix="cta_front_thickness_gif_") as tmp_dir:
        tmp_dir_path = Path(tmp_dir)
        for idx, payload in enumerate(payloads):
            frame_path = tmp_dir_path / f"frame_{idx:03d}.png"
            _render_frame(
                payload,
                reference_payload,
                constraints_by_station,
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

    print("Reference actual geometry thickness values:")
    for station in STATION_ORDER:
        print(f"  {station}: {float(reference_payload['actual_tc'][station]):.6f}")
    print(f"CTA front thickness GIF written to: {OUTPUT_GIF}")


if __name__ == "__main__":
    main()
