from pathlib import Path
import os
import sys

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
CTA_DIR = SCRIPT_DIR.parent.parent
REPO_ROOT = CTA_DIR.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

os.environ.setdefault("MPLCONFIGDIR", str(REPO_ROOT / ".mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.signal import savgol_filter

from parametrization.CTA.case import build_cta_design, to_cta_model_config
from parametrization.CTA.internal_volume_constraints import (
    CTA_CAD_REFERENCE_FRAME,
    evaluate_cta_internal_volume_constraints,
    load_cta_internal_volume_constraint_set,
)
from parametrization.bwb.builder import prepare_geometry


OUTPUT_PNG = CTA_DIR / "outputs" / "wing" / "cta_internal_volume_constraints_views.png"

SURFACE_PALETTE = (
    "#2563eb",
    "#dc2626",
    "#059669",
    "#d97706",
    "#7c3aed",
    "#0891b2",
    "#ea580c",
    "#65a30d",
    "#db2777",
    "#4f46e5",
)


def _sorted_unique_curve(x_coords: np.ndarray, z_coords: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    order = np.argsort(np.asarray(x_coords, dtype=float))
    x_sorted = np.asarray(x_coords, dtype=float)[order]
    z_sorted = np.asarray(z_coords, dtype=float)[order]
    unique_x, inverse = np.unique(x_sorted, return_inverse=True)
    if unique_x.size == x_sorted.size:
        return x_sorted, z_sorted
    z_unique = np.zeros_like(unique_x, dtype=float)
    counts = np.zeros_like(unique_x, dtype=float)
    np.add.at(z_unique, inverse, z_sorted)
    np.add.at(counts, inverse, 1.0)
    z_unique /= np.maximum(counts, 1.0)
    return unique_x, z_unique


def _section_cad_curves(prepared, y_model_m: float):
    x_air = np.asarray(prepared.section_model.x_air, dtype=float)
    chord = float(np.interp(y_model_m, prepared.loft.span_stations, prepared.loft.chord))
    le_x = float(np.interp(y_model_m, prepared.loft.span_stations, prepared.loft.leading_edge_x))
    vertical_z = float(np.interp(y_model_m, prepared.loft.span_stations, prepared.loft.vertical_y))
    twist_deg = float(prepared.spanwise_laws.twist_deg(float(y_model_m)))

    upper, lower, _ = prepared.section_model.coordinates_at_y(float(y_model_m))
    x_local = x_air * chord
    z_upper_local = np.asarray(upper, dtype=float) * chord
    z_lower_local = np.asarray(lower, dtype=float) * chord
    twist_rad = np.deg2rad(twist_deg)
    cos_twist = float(np.cos(twist_rad))
    sin_twist = float(np.sin(twist_rad))

    x_upper = le_x + x_local * cos_twist - z_upper_local * sin_twist + CTA_CAD_REFERENCE_FRAME.offset_x_m
    x_lower = le_x + x_local * cos_twist - z_lower_local * sin_twist + CTA_CAD_REFERENCE_FRAME.offset_x_m
    z_upper = vertical_z + x_local * sin_twist + z_upper_local * cos_twist + CTA_CAD_REFERENCE_FRAME.offset_z_m
    z_lower = vertical_z + x_local * sin_twist + z_lower_local * cos_twist + CTA_CAD_REFERENCE_FRAME.offset_z_m
    y_cad = float(y_model_m + CTA_CAD_REFERENCE_FRAME.offset_y_m)

    return x_upper, z_upper, x_lower, z_lower, y_cad


def _front_envelope(prepared, dense_span: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    upper = np.empty_like(dense_span, dtype=float)
    lower = np.empty_like(dense_span, dtype=float)
    for idx, yy in enumerate(dense_span):
        _, z_upper, _, z_lower, _ = _section_cad_curves(prepared, float(yy))
        upper[idx] = float(np.max(z_upper))
        lower[idx] = float(np.min(z_lower))
    return upper, lower


def _side_envelope(prepared, dense_span: np.ndarray, x_grid: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    upper_curves = []
    lower_curves = []
    for yy in dense_span:
        x_upper, z_upper, x_lower, z_lower, _ = _section_cad_curves(prepared, float(yy))
        x_u, z_u = _sorted_unique_curve(x_upper, z_upper)
        x_l, z_l = _sorted_unique_curve(x_lower, z_lower)
        upper_curves.append(np.interp(x_grid, x_u, z_u, left=np.nan, right=np.nan))
        lower_curves.append(np.interp(x_grid, x_l, z_l, left=np.nan, right=np.nan))
    upper_stack = np.vstack(upper_curves)
    lower_stack = np.vstack(lower_curves)
    upper_env = np.nanmax(upper_stack, axis=0)
    lower_env = np.nanmin(lower_stack, axis=0)
    return _smooth_side_envelope(upper_env), _smooth_side_envelope(lower_env)


def _smooth_side_envelope(z_curve: np.ndarray) -> np.ndarray:
    z_values = np.asarray(z_curve, dtype=float).copy()
    valid = np.isfinite(z_values)
    valid_idx = np.flatnonzero(valid)
    if valid_idx.size < 9:
        return z_values

    start = int(valid_idx[0])
    stop = int(valid_idx[-1]) + 1
    segment = z_values[start:stop]
    if segment.size < 9:
        return z_values

    window = min(81, segment.size if segment.size % 2 == 1 else segment.size - 1)
    if window < 9:
        return z_values

    smooth = savgol_filter(segment, window_length=window, polyorder=3, mode="interp")
    blend = np.ones_like(segment, dtype=float)
    edge = min(20, max(3, segment.size // 10))
    if edge > 0:
        ramp = np.linspace(0.0, 1.0, edge, dtype=float)
        blend[:edge] = ramp
        blend[-edge:] = ramp[::-1]
    z_values[start:stop] = (1.0 - blend) * segment + blend * smooth
    return z_values


def _surface_color_map(constraint_set):
    return {
        surface.label: SURFACE_PALETTE[idx % len(SURFACE_PALETTE)]
        for idx, surface in enumerate(constraint_set.surfaces)
    }


def _surface_style(surface, color_map):
    color = color_map[surface.label]
    alpha = 0.34
    return color, color, alpha


def _add_surface_projection_plan(ax, surface, color_map):
    xy = np.asarray(surface.vertices_xyz_m[:, :2], dtype=float)
    face, edge, alpha = _surface_style(surface, color_map)
    ax.fill(xy[:, 0], xy[:, 1], facecolor=face, edgecolor=edge, linewidth=1.2, alpha=alpha, zorder=6)
    if np.max(np.abs(xy[:, 1])) > 1.0e-9:
        xy_m = xy.copy()
        xy_m[:, 1] *= -1.0
        ax.fill(xy_m[:, 0], xy_m[:, 1], facecolor=face, edgecolor=edge, linewidth=1.2, alpha=alpha, zorder=6)


def _add_surface_projection_front(ax, surface, color_map):
    yz = np.asarray(surface.vertices_xyz_m[:, 1:3], dtype=float)
    face, edge, alpha = _surface_style(surface, color_map)
    ax.fill(yz[:, 0], yz[:, 1], facecolor=face, edgecolor=edge, linewidth=1.2, alpha=alpha, zorder=6)
    if np.max(np.abs(yz[:, 0])) > 1.0e-9:
        yz_m = yz.copy()
        yz_m[:, 0] *= -1.0
        ax.fill(yz_m[:, 0], yz_m[:, 1], facecolor=face, edgecolor=edge, linewidth=1.2, alpha=alpha, zorder=6)


def _add_surface_projection_side(ax, surface, color_map):
    xz = np.asarray(surface.vertices_xyz_m[:, (0, 2)], dtype=float)
    face, edge, alpha = _surface_style(surface, color_map)
    ax.fill(xz[:, 0], xz[:, 1], facecolor=face, edgecolor=edge, linewidth=1.2, alpha=alpha, zorder=6)


def main() -> None:
    OUTPUT_PNG.parent.mkdir(parents=True, exist_ok=True)

    design = build_cta_design()
    config = to_cta_model_config(design, use_cta_anchor_twist=True)
    prepared = prepare_geometry(config)
    constraint_set = load_cta_internal_volume_constraint_set()
    constraint_result = evaluate_cta_internal_volume_constraints(prepared=prepared, triangle_resolution=10)
    result_by_label = {item.label: item for item in constraint_result.surface_results}
    color_map = _surface_color_map(constraint_set)

    dense_span = np.unique(
        np.concatenate(
            [
                np.linspace(0.0, config.topology.span, 700),
                prepared.loft.span_stations,
                config.topology.anchor_y_array,
            ]
        )
    )
    planform = prepared.planform
    x_le = np.array([float(planform.le_x(float(y))) for y in dense_span], dtype=float) + CTA_CAD_REFERENCE_FRAME.offset_x_m
    x_te = np.array([float(planform.te_x(float(y))) for y in dense_span], dtype=float) + CTA_CAD_REFERENCE_FRAME.offset_x_m
    y_cad = dense_span + CTA_CAD_REFERENCE_FRAME.offset_y_m

    front_upper, front_lower = _front_envelope(prepared, dense_span)
    y_front = np.concatenate([-y_cad[::-1], y_cad])
    z_front_upper = np.concatenate([front_upper[::-1], front_upper])
    z_front_lower = np.concatenate([front_lower[::-1], front_lower])

    x_side_min = min(float(np.min(x_le)), float(np.min(x_te)))
    x_side_max = max(float(np.max(x_le)), float(np.max(x_te)))
    x_grid = np.linspace(x_side_min, x_side_max, 1400, dtype=float)
    z_side_upper, z_side_lower = _side_envelope(prepared, dense_span, x_grid)

    fig, (ax_plan, ax_front, ax_side) = plt.subplots(
        1,
        3,
        figsize=(18.2, 7.2),
    )
    fig.patch.set_facecolor("white")

    wing_plan_color = "#d9e7fb"
    wing_front_color = "#d9e7fb"
    wing_side_color = "#d9e7fb"
    outline_color = "#0f4c5c"

    ax_plan.fill_betweenx(y_cad, x_le, x_te, color=wing_plan_color, alpha=0.85, zorder=1)
    ax_plan.fill_betweenx(-y_cad, x_le, x_te, color=wing_plan_color, alpha=0.85, zorder=1)
    ax_plan.plot(x_le, y_cad, color=outline_color, linewidth=1.8, zorder=3)
    ax_plan.plot(x_te, y_cad, color=outline_color, linewidth=1.8, zorder=3)
    ax_plan.plot(x_le, -y_cad, color=outline_color, linewidth=1.8, zorder=3)
    ax_plan.plot(x_te, -y_cad, color=outline_color, linewidth=1.8, zorder=3)

    ax_front.fill_between(y_front, z_front_lower, z_front_upper, color=wing_front_color, alpha=0.85, zorder=1)
    ax_front.plot(y_front, z_front_upper, color=outline_color, linewidth=1.8, zorder=3)
    ax_front.plot(y_front, z_front_lower, color=outline_color, linewidth=1.8, zorder=3)

    valid_side = np.isfinite(z_side_upper) & np.isfinite(z_side_lower)
    ax_side.fill_between(x_grid[valid_side], z_side_lower[valid_side], z_side_upper[valid_side], color=wing_side_color, alpha=0.85, zorder=1)
    ax_side.plot(x_grid[valid_side], z_side_upper[valid_side], color=outline_color, linewidth=1.8, zorder=3)
    ax_side.plot(x_grid[valid_side], z_side_lower[valid_side], color=outline_color, linewidth=1.8, zorder=3)

    for surface in constraint_set.surfaces:
        result = result_by_label[surface.label]
        _add_surface_projection_plan(ax_plan, surface, color_map)
        _add_surface_projection_front(ax_front, surface, color_map)
        _add_surface_projection_side(ax_side, surface, color_map)

    top_view_legend_lines = []
    top_view_callouts = []
    for idx, surface in enumerate(constraint_set.surfaces, start=1):
        centroid = np.mean(surface.vertices_xyz_m, axis=0)
        label = surface.sub_category.replace("_", " ")
        text_color = color_map[surface.label]
        top_view_legend_lines.append((idx, label, text_color))
        top_view_callouts.append(
            {
                "index": idx,
                "label": label,
                "x": float(centroid[0]),
                "y": float(centroid[1]),
                "color": text_color,
                "face": "#eff6ff" if surface.category == "Payload" else "#fffbeb",
            }
        )

    ax_plan.set_title("Top view (X-Y)")
    ax_front.set_title("Front view (Y-Z)")
    ax_side.set_title("Side view (X-Z)")

    ax_plan.set_xlabel("X CAD [m]")
    ax_plan.set_ylabel("Y CAD [m]")
    ax_front.set_xlabel("Y CAD [m]")
    ax_front.set_ylabel("Z CAD [m]")
    ax_side.set_xlabel("X CAD [m]")
    ax_side.set_ylabel("Z CAD [m]")

    ax_plan.set_aspect("equal", adjustable="box")
    ax_front.set_aspect("equal", adjustable="box")
    ax_side.set_aspect("equal", adjustable="box")

    x_range = max(float(np.max(x_te) - np.min(x_le)), 1.0)
    y_range = max(float(np.max(np.abs(y_cad))) * 2.0, 1.0)
    x_pad = 0.06 * x_range
    x_callout_pad = 0.08 * x_range
    y_pad = 0.08 * max(float(np.max(np.abs(y_cad))) * 2.0, 1.0)
    z_all_min = min(float(np.nanmin(z_side_lower[valid_side])), min(float(np.min(s.vertices_xyz_m[:, 2])) for s in constraint_set.surfaces))
    z_all_max = max(float(np.nanmax(z_side_upper[valid_side])), max(float(np.max(s.vertices_xyz_m[:, 2])) for s in constraint_set.surfaces))
    z_pad = 0.10 * max(z_all_max - z_all_min, 1.0)

    ax_plan.set_xlim(float(np.min(x_le)) - x_pad - x_callout_pad, float(np.max(x_te)) + x_pad + x_callout_pad)
    ax_plan.set_ylim(-float(np.max(y_cad)) - y_pad, float(np.max(y_cad)) + y_pad)
    ax_front.set_xlim(-float(np.max(y_cad)) - y_pad, float(np.max(y_cad)) + y_pad)
    ax_front.set_ylim(z_all_min - z_pad, z_all_max + z_pad)
    ax_side.set_xlim(float(np.min(x_le)) - x_pad, float(np.max(x_te)) + x_pad)
    ax_side.set_ylim(z_all_min - z_pad, z_all_max + z_pad)

    for ax in (ax_plan, ax_front, ax_side):
        ax.grid(True, linewidth=0.35, alpha=0.22)

    x_plan_min, x_plan_max = ax_plan.get_xlim()
    y_plan_min, y_plan_max = ax_plan.get_ylim()
    sorted_items = sorted(top_view_callouts, key=lambda item: item["x"])
    n_items = len(sorted_items)
    x_text_positions = np.linspace(
        x_plan_min + 0.38 * (x_plan_max - x_plan_min),
        x_plan_max - 0.06 * (x_plan_max - x_plan_min),
        n_items,
    )
    y_row_top = y_plan_max - 0.12 * y_range
    y_row_bottom = y_plan_max - 0.22 * y_range
    y_text_positions = np.array(
        [y_row_top if idx % 2 == 0 else y_row_bottom for idx in range(n_items)],
        dtype=float,
    )

    for item, x_text, y_text in zip(sorted_items, x_text_positions, y_text_positions):
        ax_plan.annotate(
            f"{item['index']}",
            xy=(item["x"], item["y"]),
            xytext=(x_text, y_text),
            textcoords="data",
            ha="center",
            va="center",
            fontsize=7.8,
            color=item["color"],
            zorder=10,
            bbox={
                "boxstyle": "circle,pad=0.24",
                "facecolor": item["face"],
                "edgecolor": item["color"],
                "linewidth": 1.1,
                "alpha": 0.98,
            },
            arrowprops={
                "arrowstyle": "-|>",
                "color": item["color"],
                "linewidth": 1.0,
                "shrinkA": 8,
                "shrinkB": 6,
                "mutation_scale": 8,
                "connectionstyle": "arc3,rad=0.0",
            },
        )
        ax_plan.plot(item["x"], item["y"], marker="o", markersize=2.8, color=item["color"], zorder=9)

    top_view_text = "Top-view labels\n" + "\n".join(
        f"{idx}. {label}" for idx, label, _ in top_view_legend_lines
    )
    ax_plan.text(
        0.02,
        0.04,
        top_view_text,
        transform=ax_plan.transAxes,
        fontsize=8.0,
        color="#0f172a",
        ha="left",
        va="bottom",
        zorder=10,
        linespacing=1.35,
        bbox={
            "boxstyle": "round,pad=0.35",
            "facecolor": "white",
            "edgecolor": "#cbd5e1",
            "alpha": 0.96,
        },
    )

    legend_items = [
        Line2D(
            [0],
            [0],
            color=color_map[surface.label],
            lw=4,
            label=surface.sub_category.replace("_", " "),
        )
        for surface in constraint_set.surfaces
    ]

    summary_text = (
        f"CTA internal volume constraints in CAD frame | "
        f"offsets = ({CTA_CAD_REFERENCE_FRAME.offset_x_m:.6f}, "
        f"{CTA_CAD_REFERENCE_FRAME.offset_y_m:.3f}, "
        f"{CTA_CAD_REFERENCE_FRAME.offset_z_m:.3f}) m\n"
        f"Enclosed volume = {constraint_result.total_enclosed_volume_m3:.2f} m³ | "
        f"minimum margin = {constraint_result.minimum_margin_m:.3f} m | "
        f"surfaces = {len(constraint_result.surface_results)}"
    )
    fig.suptitle(summary_text, y=0.97, fontsize=12.0)
    fig.legend(
        handles=legend_items,
        loc="upper center",
        ncol=5,
        frameon=True,
        bbox_to_anchor=(0.5, 0.905),
        borderpad=0.6,
        columnspacing=1.4,
        handlelength=2.6,
        facecolor="white",
        edgecolor="#cbd5e1",
        framealpha=0.96,
    )

    fig.subplots_adjust(
        left=0.05,
        right=0.99,
        bottom=0.09,
        top=0.80,
        wspace=0.18,
    )

    fig.savefig(OUTPUT_PNG, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {OUTPUT_PNG}")


if __name__ == "__main__":
    main()
