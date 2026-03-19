from pathlib import Path
import os
import shutil
import sys
import tempfile

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR.parent.parent / ".mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Arc

from parametrization.CTA.reference import build_reference_design, to_cta_model_config
from parametrization.core.builder import prepare_geometry


def save_figure(fig, path: Path, **kwargs) -> None:
    try:
        fig.savefig(path, **kwargs)
    except OSError as exc:
        if "Resource deadlock avoided" not in str(exc):
            raise
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir) / path.name
            fig.savefig(tmp_path, **kwargs)
            shutil.copy2(tmp_path, path)


def polygon_area_xy(polygon: np.ndarray) -> float:
    pts = np.asarray(polygon, dtype=float)
    x = pts[:, 0]
    y = pts[:, 1]
    return float(0.5 * np.abs(np.dot(x[:-1], y[1:]) + x[-1] * y[0] - np.dot(y[:-1], x[1:]) - y[-1] * x[0]))


def build_cargo_and_engines(planform, y_sections: np.ndarray):
    y_cargo = float(y_sections[1])
    x_cut = float(planform.te_x(y_cargo))

    y_curve = np.linspace(0.0, y_cargo, 160)
    x_curve = np.array([float(planform.le_x(float(value))) for value in y_curve], dtype=float)

    upper = np.column_stack([x_curve, y_curve])
    top_edge = np.asarray([[x_cut, y_cargo]], dtype=float)
    rear_edge = np.asarray([[x_cut, -y_cargo]], dtype=float)
    lower_start = np.asarray([[float(planform.le_x(y_cargo)), -y_cargo]], dtype=float)
    lower = np.column_stack([x_curve[::-1], -y_curve[::-1]])
    cargo_polygon = np.vstack([upper, top_edge, rear_edge, lower_start, lower, upper[:1]])

    y_engine = np.linspace(-y_cargo, y_cargo, 800)
    x_engine_right = np.array([float(planform.te_x(abs(float(value)))) for value in y_engine], dtype=float)
    valid = x_engine_right > x_cut

    y_engine_half = np.linspace(0.0, y_cargo, 800)
    x_te_half = np.array([float(planform.te_x(float(value))) for value in y_engine_half], dtype=float)
    width_half = np.maximum(0.0, x_te_half - x_cut)
    engine_area_half = np.sum(0.5 * (width_half[1:] + width_half[:-1]) * (y_engine_half[1:] - y_engine_half[:-1]))
    engine_area = float(2.0 * engine_area_half)

    return {
        "cargo_polygon": cargo_polygon,
        "cargo_area": polygon_area_xy(cargo_polygon),
        "x_cut": x_cut,
        "y_cargo": y_cargo,
        "y_engine": y_engine,
        "x_engine_right": x_engine_right,
        "engine_valid": valid,
        "engine_area": engine_area,
    }


def front_envelope(prepared, dense_span: np.ndarray, chord_dense: np.ndarray):
    twist_dense = np.array([prepared.spanwise_laws.twist_deg(float(y)) for y in dense_span], dtype=float)
    vertical_center = np.interp(dense_span, prepared.loft.span_stations, prepared.loft.vertical_y)
    upper = np.empty_like(dense_span)
    lower = np.empty_like(dense_span)
    for idx, (yy, chord, twist_here, vertical_here) in enumerate(
        zip(dense_span, chord_dense, twist_dense, vertical_center)
    ):
        yu, yl, _ = prepared.section_model.coordinates_at_y(float(yy))
        x_local = prepared.section_model.x_air * float(chord)
        zu = yu * float(chord)
        zl = yl * float(chord)
        twist_rad = np.deg2rad(float(twist_here))
        sin_twist = np.sin(twist_rad)
        cos_twist = np.cos(twist_rad)
        upper_here = float(vertical_here) + x_local * sin_twist + zu * cos_twist
        lower_here = float(vertical_here) + x_local * sin_twist + zl * cos_twist
        upper[idx] = float(np.max(upper_here))
        lower[idx] = float(np.min(lower_here))
    return upper, lower, vertical_center


def add_span_bracket(ax, x: float, y0: float, y1: float, label: str, color: str = "#475569") -> None:
    ax.annotate(
        "",
        xy=(x, y1),
        xytext=(x, y0),
        arrowprops={"arrowstyle": "<->", "color": color, "linewidth": 1.0},
        zorder=8,
    )
    ax.text(
        x - 0.16,
        0.5 * (y0 + y1),
        label,
        rotation=90,
        ha="right",
        va="center",
        fontsize=8.0,
        color=color,
        bbox={"boxstyle": "round,pad=0.12", "facecolor": "white", "edgecolor": "none", "alpha": 0.86},
        zorder=9,
    )


def add_sweep_label(ax, x0: float, y0: float, x1: float, y1: float, text: str) -> None:
    dx = float(x1 - x0)
    dy = float(y1 - y0)
    seg_len = max(np.hypot(dx, dy), 1e-9)
    nx = -dy / seg_len
    ny = dx / seg_len
    xm = 0.5 * (x0 + x1)
    ym = 0.5 * (y0 + y1)
    ax.plot([x0, x1], [y0, y1], color="#0f4c5c", linewidth=0.8, linestyle=":", alpha=0.42, zorder=4)
    ax.text(
        xm + 1.2 * nx,
        ym + 1.2 * ny,
        text,
        fontsize=8.0,
        color="#0f4c5c",
        ha="center",
        va="center",
        bbox={"boxstyle": "round,pad=0.14", "facecolor": "white", "edgecolor": "#cbd5e1", "alpha": 0.90},
        zorder=9,
    )


def draw_planform(
    ax,
    dense_span: np.ndarray,
    leading_edge_x: np.ndarray,
    trailing_edge_x: np.ndarray,
    y_sections: np.ndarray,
    cargo_eng: dict,
    show_cargo: bool,
    half_wing: bool = False,
):
    zone_defs = (
        ("Centre body", float(y_sections[0]), float(y_sections[1]), "#fde68a"),
        ("Transition wing", float(y_sections[1]), float(y_sections[2]), "#bfdbfe"),
        ("Outer wing", float(y_sections[2]), float(y_sections[3]), "#bbf7d0"),
    )
    for zone_name, y0, y1, zone_color in zone_defs:
        mask = (dense_span >= y0) & (dense_span <= y1)
        ax.fill_betweenx(
            dense_span[mask],
            leading_edge_x[mask],
            trailing_edge_x[mask],
            color=zone_color,
            alpha=0.55,
            zorder=1,
            label=zone_name,
        )
        if not half_wing:
            ax.fill_betweenx(
                -dense_span[mask],
                leading_edge_x[mask],
                trailing_edge_x[mask],
                color=zone_color,
                alpha=0.20,
                zorder=1,
            )

    line_kwargs = {
        "solid_joinstyle": "round",
        "solid_capstyle": "round",
        "antialiased": True,
    }
    ax.plot(leading_edge_x, dense_span, color="#0f4c5c", linewidth=2.2, zorder=3, **line_kwargs)
    ax.plot(trailing_edge_x, dense_span, color="#c44536", linewidth=2.2, zorder=3, **line_kwargs)
    if not half_wing:
        ax.plot(leading_edge_x, -dense_span, color="#0f4c5c", linewidth=2.2, zorder=3, **line_kwargs)
        ax.plot(trailing_edge_x, -dense_span, color="#c44536", linewidth=2.2, zorder=3, **line_kwargs)

    # Highlight the transition region with a soft overlay so the smooth pySpline
    # evolution is visually clearer and does not look like a sharp peak.
    y1 = float(y_sections[1])
    y2 = float(y_sections[2])
    transition_mask = (dense_span >= y1) & (dense_span <= y2)
    ax.plot(
        leading_edge_x[transition_mask],
        dense_span[transition_mask],
        color="#0f4c5c",
        linewidth=5.6,
        alpha=0.16,
        zorder=2,
        **line_kwargs,
    )
    ax.plot(
        trailing_edge_x[transition_mask],
        dense_span[transition_mask],
        color="#c44536",
        linewidth=5.6,
        alpha=0.16,
        zorder=2,
        **line_kwargs,
    )
    if not half_wing:
        ax.plot(
            leading_edge_x[transition_mask],
            -dense_span[transition_mask],
            color="#0f4c5c",
            linewidth=5.6,
            alpha=0.16,
            zorder=2,
            **line_kwargs,
        )
        ax.plot(
            trailing_edge_x[transition_mask],
            -dense_span[transition_mask],
            color="#c44536",
            linewidth=5.6,
            alpha=0.16,
            zorder=2,
            **line_kwargs,
        )

    if show_cargo:
        cargo_polygon = cargo_eng["cargo_polygon"]
        ax.fill(
            cargo_polygon[:, 0],
            cargo_polygon[:, 1],
            facecolor="#bfdbfe",
            edgecolor="#1d4ed8",
            hatch="///",
            linewidth=1.2,
            alpha=0.35,
            zorder=5,
            label="Cargo",
        )
        ax.plot(cargo_polygon[:, 0], cargo_polygon[:, 1], color="#1d4ed8", linewidth=2.8, zorder=6)

        y_engine = cargo_eng["y_engine"]
        x_engine_right = cargo_eng["x_engine_right"]
        valid = cargo_eng["engine_valid"]
        x_cut = float(cargo_eng["x_cut"])
        ax.fill_betweenx(
            y_engine[valid],
            x_cut,
            x_engine_right[valid],
            color="#fde68a",
            alpha=0.62,
            hatch="\\\\",
            edgecolor="#d97706",
            linewidth=0.0,
            zorder=6,
            label="Engines (rear remainder)",
        )
        ax.plot([x_cut, x_cut], [-cargo_eng["y_cargo"], cargo_eng["y_cargo"]], color="#d97706", linewidth=1.8, zorder=7)
        ax.text(
            float(np.mean(cargo_polygon[:, 0])),
            0.0,
            f"Cargo area ≈ {cargo_eng['cargo_area']:.1f} m²",
            fontsize=8.6,
            color="#1e3a8a",
            ha="center",
            va="center",
            bbox={"boxstyle": "round,pad=0.16", "facecolor": "white", "edgecolor": "#93c5fd", "alpha": 0.95},
            zorder=8,
        )
        ax.text(
            float(np.mean(x_engine_right[valid])) if np.any(valid) else x_cut,
            0.0,
            f"Engines area ≈ {cargo_eng['engine_area']:.1f} m²",
            fontsize=8.6,
            color="#92400e",
            ha="center",
            va="center",
            bbox={"boxstyle": "round,pad=0.16", "facecolor": "white", "edgecolor": "#fcd34d", "alpha": 0.95},
            zorder=8,
        )

    x_min = float(np.min(leading_edge_x))
    x_max = float(np.max(trailing_edge_x))
    x_margin = 0.06 * (x_max - x_min)
    y_margin = 0.07 * float(y_sections[-1])
    ax.set_xlim(x_min - x_margin, x_max + x_margin)
    if half_wing:
        ax.set_ylim(0.0, float(y_sections[-1]) + y_margin)
    else:
        ax.set_ylim(-float(y_sections[-1]) - y_margin, float(y_sections[-1]) + y_margin)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("spanwise y [m]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, linewidth=0.35, alpha=0.25)


def annotate_layout_measurements(
    ax,
    planform,
    y_sections: np.ndarray,
    c_names: tuple[str, ...],
    c_y: np.ndarray,
    c_le_x: np.ndarray,
    c_te_x: np.ndarray,
    c_chords: np.ndarray,
    config,
) -> None:
    text_style = {
        "fontsize": 10.2,
        "color": "#0f172a",
        "ha": "center",
        "va": "center",
        "zorder": 10,
    }

    def outlined_text(x_pos: float, y_pos: float, label: str, rotation: float = 0.0) -> None:
        txt = ax.text(x_pos, y_pos, label, rotation=rotation, **text_style)
        txt.set_path_effects([pe.withStroke(linewidth=2.0, foreground="white", alpha=0.95)])

    # Chord measurements centered on each chord line.
    for name, yy, le_x, te_x, chord in zip(c_names, c_y, c_le_x, c_te_x, c_chords):
        y_here = float(yy)
        if y_here < 0.35:
            y_here = 0.55
        ax.annotate(
            "",
            xy=(te_x, y_here),
            xytext=(le_x, y_here),
            arrowprops={"arrowstyle": "<->", "color": "#475569", "linewidth": 1.05},
            zorder=5,
        )
        mid_x = float(0.5 * (le_x + te_x))
        outlined_text(mid_x, y_here + 0.34, f"{name}={chord:.2f} m")

    # Span segment measurements outside the wing, just beyond TE.
    span_segments = (
        ("B1", float(y_sections[0]), float(y_sections[1]), float(config.topology.b1_span_ratio)),
        ("B2", float(y_sections[1]), float(y_sections[2]), float(config.topology.b2_span_ratio)),
        ("B3", float(y_sections[2]), float(y_sections[3]), float(config.topology.b3_span_ratio)),
    )
    for idx, (name, y0, y1, ratio) in enumerate(span_segments):
        y_mid = 0.5 * (y0 + y1)
        x_te_y0 = float(planform.te_x(y0))
        x_te_y1 = float(planform.te_x(y1))
        x_pos = max(x_te_y0, x_te_y1) + 1.4 + 0.45 * idx
        ax.annotate(
            "",
            xy=(x_pos, y1),
            xytext=(x_pos, y0),
            arrowprops={"arrowstyle": "<->", "color": "#334155", "linewidth": 0.95},
            zorder=6,
        )
        outlined_text(x_pos + 0.78, y_mid, f"{name}={ratio:.3f} b/2")

    # Sweep definition at each LE segment start:
    # angle between local LE direction and the spanwise (vertical) direction.
    sweep_segments = (
        ("S", float(y_sections[0]), float(y_sections[1]), float(config.planform.s1_deg)),
        ("S1", float(y_sections[1]), float(y_sections[2]), float(config.planform.s2_deg)),
        ("S2", float(y_sections[2]), float(y_sections[3]), float(config.planform.s3_deg)),
    )
    sweep_color = "#1d4ed8"
    for idx, (label, y0, y1, angle_deg) in enumerate(sweep_segments):
        x0 = float(planform.le_x(y0))
        x1 = float(planform.le_x(y1))
        dy = float(y1 - y0)
        dx = float(x1 - x0)
        seg_len = max(np.hypot(dx, dy), 1e-9)
        # Short local LE segment to show local direction.
        local_len = 0.35 * dy
        x_local = x0 + (dx / seg_len) * local_len
        y_local = y0 + (dy / seg_len) * local_len

        # Vertical reference (spanwise direction).
        y_ref = y0 + min(0.55 * dy, 5.4)
        ax.plot(
            [x0, x0],
            [y0, y_ref],
            color=sweep_color,
            linewidth=1.8,
            linestyle=(0, (6, 5)),
            zorder=8,
        )
        ax.plot([x0, x_local], [y0, y_local], color=sweep_color, linewidth=1.9, zorder=8)

        # Arc from LE-direction angle to spanwise (90 deg).
        radius = min(0.42 * dy, 2.4)
        theta1 = max(0.0, 90.0 - angle_deg)
        theta2 = 90.0
        arc = Arc(
            (x0, y0),
            width=2.0 * radius,
            height=2.0 * radius,
            angle=0.0,
            theta1=theta1,
            theta2=theta2,
            color=sweep_color,
            linewidth=2.0,
            zorder=9,
        )
        ax.add_patch(arc)

        theta_mid = np.deg2rad(0.5 * (theta1 + theta2))
        tx = x0 + 1.18 * radius * np.cos(theta_mid)
        ty = y0 + 1.18 * radius * np.sin(theta_mid)
        outlined_text(tx + 0.15 + 0.10 * idx, ty + 0.10, f"{label}={angle_deg:.1f}°")


def draw_front(
    ax,
    dense_span: np.ndarray,
    upper: np.ndarray,
    lower: np.ndarray,
    cargo_eng: dict,
    prepared,
    section_labels: tuple[str, ...],
    section_y: np.ndarray,
    section_chords: np.ndarray,
):
    ax.fill_between(dense_span, lower, upper, color="#dce7f1", zorder=1, alpha=0.85)
    ax.fill_between(-dense_span, lower, upper, color="#dce7f1", zorder=1, alpha=0.45)
    ax.plot(dense_span, upper, color="#0f4c5c", linewidth=2.0, zorder=3)
    ax.plot(dense_span, lower, color="#0f4c5c", linewidth=2.0, zorder=3)
    ax.plot(-dense_span, upper, color="#0f4c5c", linewidth=1.6, zorder=3)
    ax.plot(-dense_span, lower, color="#0f4c5c", linewidth=1.6, zorder=3)
    ax.axvspan(-cargo_eng["y_cargo"], cargo_eng["y_cargo"], color="#bfdbfe", alpha=0.18, zorder=0)
    ax.axvline(float(cargo_eng["y_cargo"]), color="#64748b", linewidth=1.0, linestyle="--", zorder=4)
    ax.axvline(-float(cargo_eng["y_cargo"]), color="#64748b", linewidth=1.0, linestyle="--", zorder=4)
    for idx, (name, yy, chord) in enumerate(zip(section_labels, section_y, section_chords)):
        metrics, _ = prepared.section_model.geometry_metrics_at_y(float(yy))
        t_over_c = float(metrics.max_tc)
        t_abs = float(t_over_c * float(chord))
        y_here = float(yy)
        y_up = float(np.interp(y_here, dense_span, upper))
        y_lo = float(np.interp(y_here, dense_span, lower))
        ax.plot([y_here, y_here], [y_lo, y_up], color="#64748b", linewidth=0.9, linestyle="--", zorder=4)
        ax.scatter([y_here], [y_up], color="#0f4c5c", s=16, zorder=5)
        x_offset = 0.45 + 0.14 * idx
        y_offset = 0.08 + 0.10 * (idx % 2)
        ax.text(
            y_here + x_offset,
            y_up + y_offset,
            f"{name}: t/c={t_over_c:.3f}, t={t_abs:.2f} m",
            fontsize=8.0,
            color="#334155",
            ha="left",
            va="bottom",
            bbox={"boxstyle": "round,pad=0.12", "facecolor": "white", "edgecolor": "none", "alpha": 0.86},
            zorder=6,
        )
    ax.axhline(0.0, color="#94a3b8", linewidth=0.9, linestyle=":")
    ax.axvline(0.0, color="#94a3b8", linewidth=0.9, linestyle=":")
    z_margin = 0.12 * max(1e-9, float(np.max(upper) - np.min(lower)))
    span = float(np.max(dense_span))
    ax.set_xlim(-span - 1.0, span + 1.0)
    ax.set_ylim(float(np.min(lower)) - z_margin, float(np.max(upper)) + z_margin)
    ax.set_xlabel("spanwise y [m]")
    ax.set_ylabel("vertical [m]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, linewidth=0.35, alpha=0.30)


def draw_profiles(ax, prepared, y_sections: np.ndarray, chords: np.ndarray):
    colors = ("#0f4c5c", "#1d4ed8", "#7c3aed", "#c44536")
    section_labels = ("C1", "C2", "C3", "C4")
    raw_profiles = []
    thickness_ranges = []
    for yy, chord in zip(y_sections, chords):
        yu, yl, _ = prepared.section_model.coordinates_at_y(float(yy))
        x_local = prepared.section_model.x_air * float(chord)
        zu = yu * float(chord)
        zl = yl * float(chord)
        raw_profiles.append((x_local, zu, zl, float(yy), float(chord)))
        thickness_ranges.append(float(np.max(zu) - np.min(zl)))

    offsets = [0.0]
    for idx in range(1, len(raw_profiles)):
        sep = 0.60 * (thickness_ranges[idx - 1] + thickness_ranges[idx]) + 0.7
        offsets.append(offsets[-1] - sep)

    for (x_local, zu, zl, yy, chord), off, color, label in zip(raw_profiles, offsets, colors, section_labels):
        zu_s = zu + off
        zl_s = zl + off
        order = np.argsort(x_local)
        ax.fill(
            np.concatenate([x_local[order], x_local[order][::-1]]),
            np.concatenate([zu_s[order], zl_s[order][::-1]]),
            color=color,
            alpha=0.15,
            zorder=1,
        )
        ax.plot(x_local[order], zu_s[order], color=color, linewidth=2.0, zorder=2)
        ax.plot(x_local[order], zl_s[order], color=color, linewidth=2.0, zorder=2)
        ax.text(
            float(np.max(x_local)) + 0.9,
            off,
            f"{label}: y={yy:.2f} m | c={chord:.2f} m",
            fontsize=8.3,
            color="#1f2937",
            va="center",
            bbox={"boxstyle": "round,pad=0.14", "facecolor": "white", "edgecolor": "none", "alpha": 0.85},
            zorder=4,
        )
        ax.plot([0.0, float(np.max(x_local))], [off, off], color="#94a3b8", linewidth=0.8, linestyle=":", zorder=0)

    ax.text(
        0.6,
        offsets[0] + 0.55 * thickness_ranges[0],
        "Profiles at true scale (stacked only for readability, no normalization)",
        fontsize=8.6,
        color="#334155",
        va="top",
        bbox={"boxstyle": "round,pad=0.18", "facecolor": "white", "edgecolor": "#cbd5e1", "alpha": 0.92},
    )
    ax.set_title("Section profiles (true scale)")
    ax.set_xlabel("local chordwise x [m]")
    ax.set_ylabel("local vertical [m]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, linewidth=0.35, alpha=0.25)
    ax.set_yticks([])

    x_right = max(float(np.max(item[0])) for item in raw_profiles)
    x_left = min(float(np.min(item[0])) for item in raw_profiles)
    z_top = max(float(np.max(item[1] + off)) for item, off in zip(raw_profiles, offsets))
    z_bot = min(float(np.min(item[2] + off)) for item, off in zip(raw_profiles, offsets))
    ax.set_xlim(min(-0.5, x_left - 0.5), x_right + 10.0)
    ax.set_ylim(z_bot - 0.8, z_top + 0.8)


def draw_3d(ax, prepared, dense_span: np.ndarray, leading_edge_x: np.ndarray, chord_dense: np.ndarray):
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

    ax.plot_surface(xu_s, yu_s, zu_s, rstride=1, cstride=1, linewidth=0.25, edgecolor="#3b82f6", alpha=0.78, color="#dbeafe")
    ax.plot_surface(xl_s, yl_s, zl_s, rstride=1, cstride=1, linewidth=0.25, edgecolor="#3b82f6", alpha=0.78, color="#dbeafe")
    ax.plot_surface(xu_s, -yu_s, zu_s, rstride=1, cstride=1, linewidth=0.25, edgecolor="#3b82f6", alpha=0.62, color="#e0f2fe")
    ax.plot_surface(xl_s, -yl_s, zl_s, rstride=1, cstride=1, linewidth=0.25, edgecolor="#3b82f6", alpha=0.62, color="#e0f2fe")

    ax.set_xlabel("x [m]")
    ax.set_ylabel("spanwise y [m]")
    ax.set_zlabel("vertical [m]")
    ax.set_title("CTA 3D view (inclined)")
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
    config = to_cta_model_config(design)
    prepared = prepare_geometry(config)
    planform = prepared.planform

    output_dir = SCRIPT_DIR.parent / "example_outputs" / "reference"
    output_dir.mkdir(parents=True, exist_ok=True)

    layout_png = output_dir / "cta_reference_layout.png"
    layout_svg = output_dir / "cta_reference_layout.svg"
    plan_png = output_dir / "cta_reference_planform_cargo.png"
    plan_svg = output_dir / "cta_reference_planform_cargo.svg"
    front_png = output_dir / "cta_reference_front.png"
    front_svg = output_dir / "cta_reference_front.svg"
    profiles_png = output_dir / "cta_reference_profiles.png"
    profiles_svg = output_dir / "cta_reference_profiles.svg"
    view3d_png = output_dir / "cta_reference_3d.png"

    y_sections = config.topology.y_sections_array
    dense_span = np.unique(np.concatenate([np.linspace(0.0, config.topology.span, 2400), y_sections]))
    leading_edge_x = np.array([planform.le_x(float(y)) for y in dense_span], dtype=float)
    trailing_edge_x = np.array([planform.te_x(float(y)) for y in dense_span], dtype=float)
    chord_dense = trailing_edge_x - leading_edge_x
    le_sections = config.planform.leading_edge_x_sections(config.topology)
    te_sections = config.planform.trailing_edge_x_sections(config.topology)
    chords = te_sections - le_sections
    te_points = config.planform.trailing_edge_points(config.topology)
    c_names = ("C0", "C1", "C2", "C3", "C4", "C5")
    c_y = te_points[:, 1].astype(float)
    c_te_x = te_points[:, 0].astype(float)
    c_le_x = np.array([planform.le_x(float(value)) for value in c_y], dtype=float)
    c_chords = c_te_x - c_le_x
    cargo_eng = build_cargo_and_engines(planform, y_sections)
    upper, lower, _ = front_envelope(prepared, dense_span, chord_dense)

    # Layout figure (top view only, half-wing)
    fig_layout, ax_l_plan = plt.subplots(figsize=(12.6, 8.0), constrained_layout=True)
    draw_planform(
        ax_l_plan,
        dense_span,
        leading_edge_x,
        trailing_edge_x,
        y_sections,
        cargo_eng,
        show_cargo=False,
        half_wing=True,
    )
    annotate_layout_measurements(
        ax_l_plan,
        planform,
        y_sections,
        c_names,
        c_y,
        c_le_x,
        c_te_x,
        c_chords,
        config,
    )
    ax_l_plan.set_title("CTA reference planform (top view, half-wing)")
    ax_l_plan.legend(loc="upper left", ncol=1, fontsize=10)

    save_figure(fig_layout, layout_png, dpi=230, bbox_inches="tight")
    save_figure(fig_layout, layout_svg, bbox_inches="tight")
    plt.close(fig_layout)

    # Planform with cargo/engines
    fig_plan, ax_plan = plt.subplots(figsize=(12.0, 7.0), constrained_layout=True)
    draw_planform(ax_plan, dense_span, leading_edge_x, trailing_edge_x, y_sections, cargo_eng, show_cargo=True)
    ax_plan.set_title("CTA planform (true scale) with cargo and engines")
    ax_plan.legend(loc="upper left", ncol=2, fontsize=8)
    save_figure(fig_plan, plan_png, dpi=230, bbox_inches="tight")
    save_figure(fig_plan, plan_svg, bbox_inches="tight")
    plt.close(fig_plan)

    # Front only
    fig_front, ax_front = plt.subplots(figsize=(12.0, 5.2), constrained_layout=True)
    draw_front(
        ax_front,
        dense_span,
        upper,
        lower,
        cargo_eng,
        prepared,
        c_names,
        c_y,
        c_chords,
    )
    ax_front.set_title("CTA front view (true scale)")
    save_figure(fig_front, front_png, dpi=230, bbox_inches="tight")
    save_figure(fig_front, front_svg, bbox_inches="tight")
    plt.close(fig_front)

    # Profiles only
    fig_profiles, ax_profiles = plt.subplots(figsize=(11.6, 8.1), constrained_layout=True)
    draw_profiles(ax_profiles, prepared, y_sections, chords)
    save_figure(fig_profiles, profiles_png, dpi=230, bbox_inches="tight")
    save_figure(fig_profiles, profiles_svg, bbox_inches="tight")
    plt.close(fig_profiles)

    # 3D inclined view
    fig_3d = plt.figure(figsize=(12.0, 8.2), constrained_layout=True)
    ax_3d = fig_3d.add_subplot(111, projection="3d")
    draw_3d(ax_3d, prepared, dense_span, leading_edge_x, chord_dense)
    save_figure(fig_3d, view3d_png, dpi=240, bbox_inches="tight")
    plt.close(fig_3d)

    print(f"CTA layout PNG written to: {layout_png}")
    print(f"CTA layout SVG written to: {layout_svg}")
    print(f"CTA planform PNG written to: {plan_png}")
    print(f"CTA planform SVG written to: {plan_svg}")
    print(f"CTA front PNG written to: {front_png}")
    print(f"CTA front SVG written to: {front_svg}")
    print(f"CTA profiles PNG written to: {profiles_png}")
    print(f"CTA profiles SVG written to: {profiles_svg}")
    print(f"CTA 3D PNG written to: {view3d_png}")
    print(f"Cargo area: {cargo_eng['cargo_area']:.3f} m^2")
    print(f"Engines area (rear remainder): {cargo_eng['engine_area']:.3f} m^2")
    print(f"Section interpolation: {prepared.section_model.interpolation_name}")
    print(f"Twist interpolation: {config.spanwise.twist_deg.interpolation}")
    print(f"Camber interpolation: {config.spanwise.camber_delta.interpolation}")


if __name__ == "__main__":
    main()
