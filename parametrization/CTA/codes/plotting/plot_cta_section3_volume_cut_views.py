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

from parametrization.CTA.case import build_cta_design, to_cta_model_config
from parametrization.CTA.internal_volume_constraints import (
    CTA_CAD_REFERENCE_FRAME,
    evaluate_cta_internal_volume_constraints,
    load_cta_internal_volume_constraint_set,
)
from parametrization.bwb.builder import prepare_geometry


OUTPUT_PNG = CTA_DIR / "outputs" / "wing" / "cta_section3_volume_cut_views.png"


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
    return np.nanmax(upper_stack, axis=0), np.nanmin(lower_stack, axis=0)


def _surface_style(surface, result):
    face = "#93c5fd" if surface.category == "Payload" else "#fcd34d"
    edge = "#15803d" if result.satisfied else "#b91c1c"
    alpha = 0.50 if surface.category == "Payload" else 0.42
    return face, edge, alpha


def _surface_cut_segments(surface, y_cut_cad: float) -> list[tuple[np.ndarray, np.ndarray]]:
    vertices = np.asarray(surface.vertices_xyz_m, dtype=float)
    points = []
    tol = 1.0e-9
    n_vertices = vertices.shape[0]
    for idx in range(n_vertices):
        p0 = vertices[idx]
        p1 = vertices[(idx + 1) % n_vertices]
        y0 = float(p0[1])
        y1 = float(p1[1])
        if abs(y1 - y0) <= tol:
            if abs(y_cut_cad - y0) <= tol:
                points.append(np.asarray([p0[0], p0[2]], dtype=float))
                points.append(np.asarray([p1[0], p1[2]], dtype=float))
            continue
        if (y_cut_cad - y0) * (y_cut_cad - y1) > 0.0:
            continue
        t = (y_cut_cad - y0) / (y1 - y0)
        if -tol <= t <= 1.0 + tol:
            point = p0 + t * (p1 - p0)
            points.append(np.asarray([point[0], point[2]], dtype=float))
    if not points:
        return []
    pts = np.unique(np.round(np.vstack(points), decimals=9), axis=0)
    if pts.shape[0] < 2:
        return []
    order = np.argsort(pts[:, 0])
    pts = pts[order]
    return [(pts[0], pts[-1])]


def main() -> None:
    OUTPUT_PNG.parent.mkdir(parents=True, exist_ok=True)

    design = build_cta_design()
    config = to_cta_model_config(design, use_cta_anchor_twist=True)
    prepared = prepare_geometry(config)
    constraint_set = load_cta_internal_volume_constraint_set()
    constraint_result = evaluate_cta_internal_volume_constraints(prepared=prepared, triangle_resolution=10)
    result_by_label = {item.label: item for item in constraint_result.surface_results}

    y_cut_model = float(config.topology.anchor_y_array[3])
    x_cut_upper, z_cut_upper, x_cut_lower, z_cut_lower, y_cut_cad = _section_cad_curves(prepared, y_cut_model)
    x_cut_upper, z_cut_upper = _sorted_unique_curve(x_cut_upper, z_cut_upper)
    x_cut_lower, z_cut_lower = _sorted_unique_curve(x_cut_lower, z_cut_lower)

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

    fig, axes = plt.subplots(2, 2, figsize=(16.5, 11.2))
    ax_plan, ax_front = axes[0]
    ax_side, ax_cut = axes[1]
    fig.patch.set_facecolor("white")

    wing_color = "#d9e7fb"
    outline_color = "#0f4c5c"
    cut_color = "#7c3aed"

    ax_plan.fill_betweenx(y_cad, x_le, x_te, color=wing_color, alpha=0.85, zorder=1)
    ax_plan.fill_betweenx(-y_cad, x_le, x_te, color=wing_color, alpha=0.85, zorder=1)
    ax_plan.plot(x_le, y_cad, color=outline_color, linewidth=1.8, zorder=3)
    ax_plan.plot(x_te, y_cad, color=outline_color, linewidth=1.8, zorder=3)
    ax_plan.plot(x_le, -y_cad, color=outline_color, linewidth=1.8, zorder=3)
    ax_plan.plot(x_te, -y_cad, color=outline_color, linewidth=1.8, zorder=3)
    ax_plan.axhline(y_cut_cad, color=cut_color, linewidth=2.0, linestyle="--", zorder=7)
    ax_plan.axhline(-y_cut_cad, color=cut_color, linewidth=1.5, linestyle="--", zorder=7, alpha=0.55)

    ax_front.fill_between(y_front, z_front_lower, z_front_upper, color=wing_color, alpha=0.85, zorder=1)
    ax_front.plot(y_front, z_front_upper, color=outline_color, linewidth=1.8, zorder=3)
    ax_front.plot(y_front, z_front_lower, color=outline_color, linewidth=1.8, zorder=3)
    ax_front.axvline(y_cut_cad, color=cut_color, linewidth=2.0, linestyle="--", zorder=7)
    ax_front.axvline(-y_cut_cad, color=cut_color, linewidth=1.5, linestyle="--", zorder=7, alpha=0.55)

    valid_side = np.isfinite(z_side_upper) & np.isfinite(z_side_lower)
    ax_side.fill_between(x_grid[valid_side], z_side_lower[valid_side], z_side_upper[valid_side], color=wing_color, alpha=0.85, zorder=1)
    ax_side.plot(x_grid[valid_side], z_side_upper[valid_side], color=outline_color, linewidth=1.8, zorder=3)
    ax_side.plot(x_grid[valid_side], z_side_lower[valid_side], color=outline_color, linewidth=1.8, zorder=3)
    ax_side.plot(x_cut_upper, z_cut_upper, color=cut_color, linewidth=2.4, zorder=8)
    ax_side.plot(x_cut_lower, z_cut_lower, color=cut_color, linewidth=2.4, zorder=8)

    for surface in constraint_set.surfaces:
        result = result_by_label[surface.label]
        face, edge, alpha = _surface_style(surface, result)
        xy = np.asarray(surface.vertices_xyz_m[:, :2], dtype=float)
        yz = np.asarray(surface.vertices_xyz_m[:, 1:3], dtype=float)
        xz = np.asarray(surface.vertices_xyz_m[:, (0, 2)], dtype=float)
        ax_plan.fill(xy[:, 0], xy[:, 1], facecolor=face, edgecolor=edge, linewidth=1.0, alpha=alpha, zorder=5)
        if np.max(np.abs(xy[:, 1])) > 1.0e-9:
            ax_plan.fill(xy[:, 0], -xy[:, 1], facecolor=face, edgecolor=edge, linewidth=1.0, alpha=alpha, zorder=5)
        ax_front.fill(yz[:, 0], yz[:, 1], facecolor=face, edgecolor=edge, linewidth=1.0, alpha=alpha, zorder=5)
        if np.max(np.abs(yz[:, 0])) > 1.0e-9:
            ax_front.fill(-yz[:, 0], yz[:, 1], facecolor=face, edgecolor=edge, linewidth=1.0, alpha=alpha, zorder=5)
        ax_side.fill(xz[:, 0], xz[:, 1], facecolor=face, edgecolor=edge, linewidth=1.0, alpha=alpha, zorder=5)
        for p0, p1 in _surface_cut_segments(surface, y_cut_cad):
            ax_cut.plot([p0[0], p1[0]], [p0[1], p1[1]], color=edge, linewidth=2.2, alpha=0.9)

    ax_cut.fill_between(x_cut_upper, np.interp(x_cut_upper, x_cut_lower, z_cut_lower, left=np.nan, right=np.nan), z_cut_upper, color=wing_color, alpha=0.82, zorder=1)
    ax_cut.plot(x_cut_upper, z_cut_upper, color=outline_color, linewidth=2.2, zorder=3)
    ax_cut.plot(x_cut_lower, z_cut_lower, color=outline_color, linewidth=2.2, zorder=3)

    for ax in (ax_plan, ax_front, ax_side, ax_cut):
        ax.grid(True, linewidth=0.35, alpha=0.22)

    ax_plan.set_title(f"Top view (X-Y) | cut at y={y_cut_cad:.3f} m")
    ax_front.set_title("Front view (Y-Z)")
    ax_side.set_title("Side view (X-Z)")
    ax_cut.set_title(f"Section cut at start of section 3 | y={y_cut_cad:.3f} m")

    ax_plan.set_xlabel("X CAD [m]")
    ax_plan.set_ylabel("Y CAD [m]")
    ax_front.set_xlabel("Y CAD [m]")
    ax_front.set_ylabel("Z CAD [m]")
    ax_side.set_xlabel("X CAD [m]")
    ax_side.set_ylabel("Z CAD [m]")
    ax_cut.set_xlabel("X CAD [m]")
    ax_cut.set_ylabel("Z CAD [m]")

    ax_plan.set_aspect("equal", adjustable="box")
    ax_front.set_aspect("equal", adjustable="box")
    ax_side.set_aspect("equal", adjustable="box")
    ax_cut.set_aspect("equal", adjustable="box")

    x_range = max(float(np.max(x_te) - np.min(x_le)), 1.0)
    y_range = max(float(np.max(np.abs(y_cad))) * 2.0, 1.0)
    z_min = min(float(np.nanmin(z_side_lower[valid_side])), float(np.min(z_cut_lower)))
    z_max = max(float(np.nanmax(z_side_upper[valid_side])), float(np.max(z_cut_upper)))
    z_pad = 0.10 * max(z_max - z_min, 1.0)
    x_pad = 0.06 * x_range
    y_pad = 0.08 * y_range

    ax_plan.set_xlim(float(np.min(x_le)) - x_pad, float(np.max(x_te)) + x_pad)
    ax_plan.set_ylim(-float(np.max(y_cad)) - y_pad, float(np.max(y_cad)) + y_pad)
    ax_front.set_xlim(-float(np.max(y_cad)) - y_pad, float(np.max(y_cad)) + y_pad)
    ax_front.set_ylim(z_min - z_pad, z_max + z_pad)
    ax_side.set_xlim(float(np.min(x_le)) - x_pad, float(np.max(x_te)) + x_pad)
    ax_side.set_ylim(z_min - z_pad, z_max + z_pad)
    ax_cut.set_xlim(min(float(np.min(x_cut_lower)), float(np.min(x_cut_upper))) - 1.0, max(float(np.max(x_cut_lower)), float(np.max(x_cut_upper))) + 1.0)
    ax_cut.set_ylim(z_min - z_pad, z_max + z_pad)

    fig.suptitle(
        "CTA cut at the start of section 3\n"
        "The side view can look 'open' because it is an X-Z projection of a surface changing with Y; "
        "the cut panel shows the real section at that station.",
        fontsize=15.0,
        y=0.98,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.95))
    fig.savefig(OUTPUT_PNG, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {OUTPUT_PNG}")


if __name__ == "__main__":
    main()
