from pathlib import Path
import os
import sys

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR.parent.parent / ".mplconfig"))

import matplotlib.pyplot as plt

from parametrization.CTA.reference import build_reference_design, to_cta_model_config
from parametrization.bwb.builder import prepare_geometry


def draw_3d(ax, prepared, dense_span: np.ndarray, leading_edge_x: np.ndarray, chord_dense: np.ndarray) -> None:
    x_air = prepared.section_model.x_air
    ns = dense_span.size
    nx = x_air.size

    vertical_center = np.interp(dense_span, prepared.loft.span_stations, prepared.loft.vertical_y)
    twist_dense = np.array([prepared.spanwise_laws.twist_deg(float(y)) for y in dense_span], dtype=float)

    xu = np.zeros((ns, nx), dtype=float)
    yu = np.zeros((ns, nx), dtype=float)
    zu = np.zeros((ns, nx), dtype=float)
    xl = np.zeros((ns, nx), dtype=float)
    yl = np.zeros((ns, nx), dtype=float)
    zl = np.zeros((ns, nx), dtype=float)

    for i, (yy, le, chord, twist, v0) in enumerate(
        zip(dense_span, leading_edge_x, chord_dense, twist_dense, vertical_center)
    ):
        y_up, y_lo, _ = prepared.section_model.coordinates_at_y(float(yy))
        x_local = x_air * float(chord)
        z_up_local = y_up * float(chord)
        z_lo_local = y_lo * float(chord)
        tw = np.deg2rad(float(twist))
        ctw = np.cos(tw)
        stw = np.sin(tw)

        xu[i, :] = float(le) + x_local * ctw - z_up_local * stw
        xl[i, :] = float(le) + x_local * ctw - z_lo_local * stw
        yu[i, :] = float(yy)
        yl[i, :] = float(yy)
        zu[i, :] = float(v0) + x_local * stw + z_up_local * ctw
        zl[i, :] = float(v0) + x_local * stw + z_lo_local * ctw

    stride_s = max(1, ns // 45)
    stride_x = max(1, nx // 45)
    xu_s = xu[::stride_s, ::stride_x]
    yu_s = yu[::stride_s, ::stride_x]
    zu_s = zu[::stride_s, ::stride_x]
    xl_s = xl[::stride_s, ::stride_x]
    yl_s = yl[::stride_s, ::stride_x]
    zl_s = zl[::stride_s, ::stride_x]

    ax.plot_surface(xu_s, yu_s, zu_s, rstride=1, cstride=1, linewidth=0.25, edgecolor="#3b82f6", alpha=0.82, color="#dbeafe")
    ax.plot_surface(xl_s, yl_s, zl_s, rstride=1, cstride=1, linewidth=0.25, edgecolor="#3b82f6", alpha=0.82, color="#dbeafe")
    ax.plot_surface(xu_s, -yu_s, zu_s, rstride=1, cstride=1, linewidth=0.25, edgecolor="#3b82f6", alpha=0.68, color="#e0f2fe")
    ax.plot_surface(xl_s, -yl_s, zl_s, rstride=1, cstride=1, linewidth=0.25, edgecolor="#3b82f6", alpha=0.68, color="#e0f2fe")

    ax.set_xlabel("x [m]")
    ax.set_ylabel("spanwise y [m]")
    ax.set_zlabel("vertical [m]")
    ax.set_title("CTA 3D view (interactive)")
    ax.view_init(elev=22, azim=-126)

    x_min = min(float(np.min(xu)), float(np.min(xl)))
    x_max = max(float(np.max(xu)), float(np.max(xl)))
    y_max = float(np.max(np.abs(yu)))
    z_min = min(float(np.min(zl)), float(np.min(zu)))
    z_max = max(float(np.max(zl)), float(np.max(zu)))
    ax.set_xlim(x_min - 2.0, x_max + 2.0)
    ax.set_ylim(-y_max - 1.0, y_max + 1.0)
    ax.set_zlim(z_min - 0.8, z_max + 0.8)
    try:
        ax.set_box_aspect((x_max - x_min, 2.0 * y_max, max(1e-9, z_max - z_min)))
    except Exception:
        pass


def main() -> None:
    design = build_reference_design()
    config = to_cta_model_config(design, use_reference_anchor_twist=True)
    prepared = prepare_geometry(config)
    planform = prepared.planform

    y_sections = config.topology.y_sections_array
    dense_span = np.unique(np.concatenate([np.linspace(0.0, config.topology.span, 700), y_sections]))
    leading_edge_x = np.array([planform.le_x(float(y)) for y in dense_span], dtype=float)
    trailing_edge_x = np.array([planform.te_x(float(y)) for y in dense_span], dtype=float)
    chord_dense = trailing_edge_x - leading_edge_x

    fig_3d = plt.figure(figsize=(12.0, 8.2), constrained_layout=True)
    ax_3d = fig_3d.add_subplot(111, projection="3d")
    draw_3d(ax_3d, prepared, dense_span, leading_edge_x, chord_dense)
    plt.show()


if __name__ == "__main__":
    main()
