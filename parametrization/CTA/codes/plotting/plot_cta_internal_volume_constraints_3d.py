from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

SCRIPT_DIR = Path(__file__).resolve().parent
CTA_DIR = SCRIPT_DIR.parent.parent
REPO_ROOT = CTA_DIR.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from parametrization.CTA.case import build_cta_design, to_cta_model_config
from parametrization.CTA.internal_volume_constraints import (
    CTA_CAD_REFERENCE_FRAME,
    evaluate_cta_internal_volume_constraints,
    load_cta_internal_volume_constraint_set,
)
from parametrization.bwb.builder import prepare_geometry


OUTPUT_PNG = CTA_DIR / "outputs" / "wing" / "cta_internal_volume_constraints_3d.png"


def _build_wing_surface_arrays(prepared, frame, dense_span: np.ndarray) -> tuple[np.ndarray, ...]:
    x_air = np.asarray(prepared.section_model.x_air, dtype=float)
    leading_edge_x = np.interp(dense_span, prepared.loft.span_stations, prepared.loft.leading_edge_x)
    chord_dense = np.interp(dense_span, prepared.loft.span_stations, prepared.loft.chord)
    vertical_center = np.interp(dense_span, prepared.loft.span_stations, prepared.loft.vertical_y)
    twist_dense = np.array([prepared.spanwise_laws.twist_deg(float(y)) for y in dense_span], dtype=float)

    ns = dense_span.size
    nx = x_air.size
    xu = np.zeros((ns, nx), dtype=float)
    yu = np.zeros((ns, nx), dtype=float)
    zu = np.zeros((ns, nx), dtype=float)
    xl = np.zeros((ns, nx), dtype=float)
    yl = np.zeros((ns, nx), dtype=float)
    zl = np.zeros((ns, nx), dtype=float)

    for i, (yy, le, chord, twist, v0) in enumerate(
        zip(dense_span, leading_edge_x, chord_dense, twist_dense, vertical_center)
    ):
        upper, lower, _ = prepared.section_model.coordinates_at_y(float(yy))
        x_local = x_air * float(chord)
        z_up_local = np.asarray(upper, dtype=float) * float(chord)
        z_lo_local = np.asarray(lower, dtype=float) * float(chord)
        twist_rad = np.deg2rad(float(twist))
        cos_twist = float(np.cos(twist_rad))
        sin_twist = float(np.sin(twist_rad))

        x_upper_model = float(le) + x_local * cos_twist - z_up_local * sin_twist
        x_lower_model = float(le) + x_local * cos_twist - z_lo_local * sin_twist
        y_model = np.full_like(x_local, float(yy))
        z_upper_model = float(v0) + x_local * sin_twist + z_up_local * cos_twist
        z_lower_model = float(v0) + x_local * sin_twist + z_lo_local * cos_twist

        xu[i, :] = x_upper_model + float(frame.offset_x_m)
        yu[i, :] = y_model + float(frame.offset_y_m)
        zu[i, :] = z_upper_model + float(frame.offset_z_m)
        xl[i, :] = x_lower_model + float(frame.offset_x_m)
        yl[i, :] = y_model + float(frame.offset_y_m)
        zl[i, :] = z_lower_model + float(frame.offset_z_m)

    return xu, yu, zu, xl, yl, zl


def _surface_triangles(vertices_xyz_m: np.ndarray) -> list[np.ndarray]:
    vertices = np.asarray(vertices_xyz_m, dtype=float)
    return [
        np.vstack([vertices[0], vertices[idx], vertices[idx + 1]])
        for idx in range(1, vertices.shape[0] - 1)
    ]


def _mirrored_triangles(triangles: list[np.ndarray]) -> list[np.ndarray]:
    mirrored = []
    for tri in triangles:
        tri_m = np.asarray(tri, dtype=float).copy()
        tri_m[:, 1] *= -1.0
        mirrored.append(tri_m)
    return mirrored


def main() -> None:
    OUTPUT_PNG.parent.mkdir(parents=True, exist_ok=True)

    design = build_cta_design()
    config = to_cta_model_config(design, use_cta_anchor_twist=True)
    prepared = prepare_geometry(config)
    constraint_set = load_cta_internal_volume_constraint_set()
    constraint_result = evaluate_cta_internal_volume_constraints(prepared=prepared, triangle_resolution=10)
    result_by_label = {item.label: item for item in constraint_result.surface_results}

    dense_span = np.unique(
        np.concatenate(
            [
                np.linspace(0.0, config.topology.span, 320),
                prepared.loft.span_stations,
                config.topology.anchor_y_array,
            ]
        )
    )
    xu, yu, zu, xl, yl, zl = _build_wing_surface_arrays(prepared, CTA_CAD_REFERENCE_FRAME, dense_span)

    fig = plt.figure(figsize=(13.8, 8.8))
    ax = fig.add_subplot(111, projection="3d")
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    wing_kw = {
        "rstride": 1,
        "cstride": 1,
        "linewidth": 0.0,
        "antialiased": False,
        "shade": False,
        "edgecolor": "none",
        "alpha": 0.24,
    }
    ax.plot_surface(xu, yu, zu, color="#d6e6fb", **wing_kw)
    ax.plot_surface(xl, yl, zl, color="#dceaf9", **wing_kw)
    ax.plot_surface(xu, -yu, zu, color="#d6e6fb", **wing_kw)
    ax.plot_surface(xl, -yl, zl, color="#dceaf9", **wing_kw)

    category_face = {
        "Payload": "#60a5fa",
        "LandingGear": "#f59e0b",
    }

    for surface in constraint_set.surfaces:
        result = result_by_label[surface.label]
        face_color = category_face.get(surface.category, "#94a3b8")
        edge_color = "#15803d" if result.satisfied else "#b91c1c"
        triangles = _surface_triangles(surface.vertices_xyz_m)
        all_triangles = triangles + _mirrored_triangles(triangles)
        poly = Poly3DCollection(
            all_triangles,
            facecolors=face_color,
            edgecolors=edge_color,
            linewidths=1.0,
            alpha=0.62,
        )
        ax.add_collection3d(poly)

        centroid = np.mean(surface.vertices_xyz_m, axis=0)
        label = surface.sub_category.replace("_", " ")
        ax.text(
            float(centroid[0]),
            float(centroid[1]) + 0.20,
            float(centroid[2]) + 0.08,
            label,
            fontsize=8.0,
            color=edge_color,
            ha="left",
            va="bottom",
        )

    legend_items = [
        Line2D([0], [0], color="#60a5fa", lw=8, alpha=0.8, label="Payload surfaces"),
        Line2D([0], [0], color="#f59e0b", lw=8, alpha=0.8, label="Landing gear surfaces"),
        Line2D([0], [0], color="#15803d", lw=2, label="Constraint satisfied"),
        Line2D([0], [0], color="#b91c1c", lw=2, label="Constraint violated"),
    ]
    ax.legend(handles=legend_items, loc="upper left", bbox_to_anchor=(0.01, 0.98), frameon=False)

    summary_text = (
        f"CAD frame offsets: X={CTA_CAD_REFERENCE_FRAME.offset_x_m:.6f} m, "
        f"Y={CTA_CAD_REFERENCE_FRAME.offset_y_m:.3f} m, "
        f"Z={CTA_CAD_REFERENCE_FRAME.offset_z_m:.3f} m\n"
        f"Enclosed volume: {constraint_result.total_enclosed_volume_m3:.2f} m³\n"
        f"Minimum margin: {constraint_result.minimum_margin_m:.3f} m"
    )
    ax.text2D(
        0.01,
        0.02,
        summary_text,
        transform=ax.transAxes,
        fontsize=10.5,
        color="#0f172a",
        bbox={"boxstyle": "round,pad=0.40", "facecolor": "#ffffff", "edgecolor": "#cbd5e1", "alpha": 0.96},
    )

    ax.set_title("CTA geometry with internal volume constraints", pad=10.0)
    ax.set_xlabel("X CAD [m]", labelpad=8.0)
    ax.set_ylabel("Y CAD [m]", labelpad=8.0)
    ax.set_zlabel("Z CAD [m]", labelpad=8.0)
    ax.view_init(elev=21, azim=-132)

    x_min = min(float(np.min(xu)), float(np.min(xl)))
    x_max = max(float(np.max(xu)), float(np.max(xl)))
    y_max = max(float(np.max(np.abs(yu))), max(float(np.max(np.abs(s.vertices_xyz_m[:, 1]))) for s in constraint_set.surfaces))
    z_min = min(float(np.min(zl)), min(float(np.min(s.vertices_xyz_m[:, 2])) for s in constraint_set.surfaces))
    z_max = max(float(np.max(zu)), max(float(np.max(s.vertices_xyz_m[:, 2])) for s in constraint_set.surfaces))
    x_pad = 0.05 * max(x_max - x_min, 1.0)
    y_pad = 0.10 * max(2.0 * y_max, 1.0)
    z_pad = 0.10 * max(z_max - z_min, 1.0)

    ax.set_xlim(x_min - x_pad, x_max + x_pad)
    ax.set_ylim(-y_max - y_pad, y_max + y_pad)
    ax.set_zlim(z_min - z_pad, z_max + z_pad)
    try:
        ax.set_box_aspect((x_max - x_min, 2.0 * y_max, 1.25 * max(z_max - z_min, 1.0e-9)))
    except Exception:
        pass

    ax.grid(False)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        try:
            axis.pane.set_alpha(0.0)
            axis.pane.set_edgecolor((1.0, 1.0, 1.0, 0.0))
        except Exception:
            pass
    ax.tick_params(pad=2.0, labelsize=9.0)

    fig.tight_layout()
    fig.savefig(OUTPUT_PNG, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {OUTPUT_PNG}")


if __name__ == "__main__":
    main()
