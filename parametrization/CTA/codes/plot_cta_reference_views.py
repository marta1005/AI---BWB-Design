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

from parametrization.CTA.reference import (
    build_reference_design,
    cta_fixed_values,
    to_cta_model_config,
)
from parametrization.CTA.design_space import build_cta_design_space
from parametrization.bwb.builder import prepare_geometry


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


def glider_reference_z_limits(padding: float = 0.12):
    geo_path = SCRIPT_DIR.parent / "airfoils" / "bwb_glider.geo"
    if not geo_path.exists():
        return None
    z_values = []
    for line in geo_path.read_text().splitlines():
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            z_values.append(float(parts[2]))
        except ValueError:
            continue
    if not z_values:
        return None
    z_array = np.asarray(z_values, dtype=float)
    z_pad = float(padding) * max(float(np.ptp(z_array)), 1.0)
    return float(np.min(z_array) - z_pad), float(np.max(z_array) + z_pad)


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


def front_anchor_sections(prepared, anchor_y: np.ndarray):
    twist_anchor = np.array([prepared.spanwise_laws.twist_deg(float(y)) for y in anchor_y], dtype=float)
    vertical_center = np.interp(anchor_y, prepared.loft.span_stations, prepared.loft.vertical_y)
    chords = np.interp(anchor_y, prepared.loft.span_stations, prepared.loft.chord)
    upper = np.empty_like(anchor_y)
    lower = np.empty_like(anchor_y)
    for idx, (yy, chord, twist_here, vertical_here) in enumerate(
        zip(anchor_y, chords, twist_anchor, vertical_center)
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


def front_linear_envelope(prepared, dense_y: np.ndarray, anchor_y: np.ndarray):
    upper_anchor, lower_anchor, _ = front_anchor_sections(prepared, anchor_y)
    upper_dense = blend_linear_sections(dense_y, anchor_y, upper_anchor)
    lower_dense = blend_linear_sections(dense_y, anchor_y, lower_anchor)
    return upper_dense, lower_dense, upper_anchor, lower_anchor


def front_edge_traces(prepared, spanwise_y: np.ndarray):
    vertical_center = np.interp(spanwise_y, prepared.loft.span_stations, prepared.loft.vertical_y)
    chord_dense = np.interp(spanwise_y, prepared.loft.span_stations, prepared.loft.chord)
    twist_dense = np.array([prepared.spanwise_laws.twist_deg(float(y)) for y in spanwise_y], dtype=float)

    z_le = np.empty_like(spanwise_y, dtype=float)
    z_te_upper = np.empty_like(spanwise_y, dtype=float)
    z_te_lower = np.empty_like(spanwise_y, dtype=float)
    z_te_mid = np.empty_like(spanwise_y, dtype=float)

    for idx, (yy, chord, twist_here, vertical_here) in enumerate(
        zip(spanwise_y, chord_dense, twist_dense, vertical_center)
    ):
        yu, yl, _ = prepared.section_model.coordinates_at_y(float(yy))
        x_local = prepared.section_model.x_air * float(chord)
        zu = yu * float(chord)
        zl = yl * float(chord)
        twist_rad = np.deg2rad(float(twist_here))
        sin_twist = np.sin(twist_rad)
        cos_twist = np.cos(twist_rad)

        z_le[idx] = float(vertical_here + x_local[0] * sin_twist + zu[0] * cos_twist)
        z_te_upper[idx] = float(vertical_here + x_local[-1] * sin_twist + zu[-1] * cos_twist)
        z_te_lower[idx] = float(vertical_here + x_local[-1] * sin_twist + zl[-1] * cos_twist)
        z_te_mid[idx] = 0.5 * (z_te_upper[idx] + z_te_lower[idx])

    return z_le, z_te_upper, z_te_lower, z_te_mid


def front_linear_edge_traces(prepared, dense_y: np.ndarray, anchor_y: np.ndarray):
    z_le_anchor, z_te_upper_anchor, z_te_lower_anchor, z_te_mid_anchor = front_edge_traces(prepared, anchor_y)
    z_le = blend_linear_sections(dense_y, anchor_y, z_le_anchor)
    z_te_upper = blend_linear_sections(dense_y, anchor_y, z_te_upper_anchor)
    z_te_lower = blend_linear_sections(dense_y, anchor_y, z_te_lower_anchor)
    z_te_mid = blend_linear_sections(dense_y, anchor_y, z_te_mid_anchor)
    return z_le, z_te_upper, z_te_lower, z_te_mid


def blend_linear_sections(
    dense_y: np.ndarray,
    anchor_y: np.ndarray,
    anchor_values: np.ndarray,
    blend_fraction: float = 0.07,
    root_blend_fraction: float = 0.45,
) -> np.ndarray:
    dense_y = np.asarray(dense_y, dtype=float)
    anchor_y = np.asarray(anchor_y, dtype=float)
    anchor_values = np.asarray(anchor_values, dtype=float)
    blended = np.interp(dense_y, anchor_y, anchor_values)

    if anchor_y.size < 3:
        return blended

    for idx in range(1, anchor_y.size - 1):
        y_prev = float(anchor_y[idx - 1])
        y_here = float(anchor_y[idx])
        y_next = float(anchor_y[idx + 1])
        v_prev = float(anchor_values[idx - 1])
        v_here = float(anchor_values[idx])
        v_next = float(anchor_values[idx + 1])

        left_span = max(y_here - y_prev, 1e-12)
        right_span = max(y_next - y_here, 1e-12)
        half_width = blend_fraction * min(left_span, right_span)
        if half_width <= 0.0:
            continue

        left_slope = (v_here - v_prev) / left_span
        right_slope = (v_next - v_here) / right_span
        y_left = y_here - half_width
        y_right = y_here + half_width
        v_left = v_here - left_slope * half_width
        v_right = v_here + right_slope * half_width

        blend_mask = (dense_y >= y_left) & (dense_y <= y_right)
        if np.any(blend_mask):
            local_y = dense_y[blend_mask]
            left_line = v_here + left_slope * (local_y - y_here)
            right_line = v_here + right_slope * (local_y - y_here)
            t = (local_y - y_left) / max(y_right - y_left, 1e-12)
            # Quintic smoothstep: value, slope, and curvature match the linear
            # segments at both ends, giving a short local C2 blend.
            w = t * t * t * (10.0 + t * (-15.0 + 6.0 * t))
            blended[blend_mask] = (1.0 - w) * left_line + w * right_line

    first_span = float(anchor_y[1] - anchor_y[0])
    root_blend_y = max(0.0, root_blend_fraction * first_span)
    if root_blend_y > 1e-12:
        root_value = float(anchor_values[0])
        root_slope = float((anchor_values[1] - anchor_values[0]) / max(first_span, 1e-12))
        join_value = root_value + root_slope * root_blend_y
        root_mask = dense_y <= root_blend_y
        if np.any(root_mask):
            y_local = dense_y[root_mask]
            t = np.clip(y_local / root_blend_y, 0.0, 1.0)
            # Quintic Hermite segment:
            # f(0)=root_value, f'(0)=0, f''(0)=0
            # f(root_blend_y)=join_value, f'(root_blend_y)=root_slope, f''(root_blend_y)=0
            delta = join_value - root_value
            h00 = 1.0 - 10.0 * t**3 + 15.0 * t**4 - 6.0 * t**5
            h10 = t - 6.0 * t**3 + 8.0 * t**4 - 3.0 * t**5
            h20 = 0.5 * t**2 - 1.5 * t**3 + 1.5 * t**4 - 0.5 * t**5
            h01 = 10.0 * t**3 - 15.0 * t**4 + 6.0 * t**5
            h11 = -4.0 * t**3 + 7.0 * t**4 - 3.0 * t**5
            h21 = 0.5 * t**3 - t**4 + 0.5 * t**5
            blended[root_mask] = (
                h00 * root_value
                + h10 * root_blend_y * 0.0
                + h20 * (root_blend_y**2) * 0.0
                + h01 * join_value
                + h11 * root_blend_y * root_slope
                + h21 * (root_blend_y**2) * 0.0
            )

    return blended


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


def add_horizontal_bracket(ax, x0: float, x1: float, y: float, label: str, color: str = "#475569") -> None:
    ax.annotate(
        "",
        xy=(x1, y),
        xytext=(x0, y),
        arrowprops={"arrowstyle": "<->", "color": color, "linewidth": 1.0},
        zorder=8,
    )
    ax.text(
        0.5 * (x0 + x1),
        y + 0.12,
        label,
        ha="center",
        va="bottom",
        fontsize=8.4,
        color=color,
        bbox={"boxstyle": "round,pad=0.12", "facecolor": "white", "edgecolor": "none", "alpha": 0.88},
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
    label_zones: bool = False,
    faded_symmetry: bool = False,
):
    zone_defs = (
        ("Center body", float(y_sections[0]), float(y_sections[1]), "#fde68a"),
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
        if not half_wing or faded_symmetry:
            ax.fill_betweenx(
                -dense_span[mask],
                leading_edge_x[mask],
                trailing_edge_x[mask],
                color=zone_color,
                alpha=0.12 if faded_symmetry else 0.20,
                zorder=1,
            )
        if label_zones:
            y_mid = 0.5 * (y0 + y1)
            x_mid = 0.5 * (
                float(np.interp(y_mid, dense_span, leading_edge_x))
                + float(np.interp(y_mid, dense_span, trailing_edge_x))
            )
            ax.text(
                x_mid,
                y_mid,
                zone_name,
                fontsize=11.4,
                color="#0f172a",
                ha="center",
                va="center",
                bbox={"boxstyle": "round,pad=0.22", "facecolor": "white", "edgecolor": "#cbd5e1", "alpha": 0.92},
                zorder=7,
            )

    line_kwargs = {
        "solid_joinstyle": "round",
        "solid_capstyle": "round",
        "antialiased": True,
    }
    ax.plot(leading_edge_x, dense_span, color="#0f4c5c", linewidth=2.2, zorder=3, **line_kwargs)
    ax.plot(trailing_edge_x, dense_span, color="#c44536", linewidth=2.2, zorder=3, **line_kwargs)
    if not half_wing or faded_symmetry:
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
    if not half_wing or faded_symmetry:
        ax.plot(
            leading_edge_x[transition_mask],
            -dense_span[transition_mask],
            color="#0f4c5c",
            linewidth=5.6,
            alpha=0.09 if faded_symmetry else 0.16,
            zorder=2,
            **line_kwargs,
        )
        ax.plot(
            trailing_edge_x[transition_mask],
            -dense_span[transition_mask],
            color="#c44536",
            linewidth=5.6,
            alpha=0.09 if faded_symmetry else 0.16,
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
        if faded_symmetry:
            ax.set_ylim(-2.0, float(y_sections[-1]) + y_margin)
        else:
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
    public_reference: dict,
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

    taper_ratio = float(public_reference["c4_c3_ratio"])

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
        if name == "C4":
            outlined_text(mid_x, y_here + 0.46, f"C4/C3={taper_ratio:.3f}\nC4={chord:.2f} m")
        else:
            outlined_text(mid_x, y_here + 0.34, f"{name}={chord:.2f} m")

    # Span segment measurements outside the wing, just beyond TE.
    wing_span = float(y_sections[3] - y_sections[1])
    b2_length = float(y_sections[2] - y_sections[1])
    b3_length = float(y_sections[3] - y_sections[2])
    span_segments = (
        ("B1", float(y_sections[0]), float(y_sections[1]), f"{float(y_sections[1] - y_sections[0]):.2f} m"),
        ("B2", float(y_sections[1]), float(y_sections[2]), f"{100.0 * b2_length / max(wing_span, 1e-12):.1f}% wing span"),
        ("B3", float(y_sections[2]), float(y_sections[3]), f"{100.0 * b3_length / max(wing_span, 1e-12):.1f}% wing span"),
    )
    for idx, (name, y0, y1, label_text) in enumerate(span_segments):
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
        outlined_text(x_pos + 0.92, y_mid, f"{name}={label_text}")

    # Sweep definition at each LE segment start:
    # angle between local LE direction and the spanwise (vertical) direction.
    le_helper_points = np.asarray(config.planform.body_le_fixed_points, dtype=float)
    sweep_segments = [
        (
            "S",
            float(y_sections[0]),
            float(le_helper_points[0, 1]),
            float(planform.le_x(float(y_sections[0]))),
            float(planform.le_x(float(le_helper_points[0, 1]))),
            66.87,
        ),
        (
            "S1 (50%)",
            float(c_y[2]),
            float(c_y[3]),
            float(c_le_x[2] + 0.50 * c_chords[2]),
            float(c_le_x[3] + 0.50 * c_chords[3]),
            float(public_reference["s2_deg"]),
        ),
        (
            "S2 (25%)",
            float(c_y[3]),
            float(c_y[4]),
            float(c_le_x[3] + 0.25 * c_chords[3]),
            float(c_le_x[4] + 0.25 * c_chords[4]),
            float(public_reference["s3_deg"]),
        ),
    ]
    sweep_color = "#1d4ed8"
    for idx, (label, y0, y1, x0, x1, angle_deg) in enumerate(sweep_segments):
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
        txt = ax.text(
            tx + 0.15 + 0.10 * idx,
            ty + 0.10,
            f"{label}={angle_deg:.1f}°",
            fontsize=10.9,
            color=sweep_color,
            ha="center",
            va="center",
            zorder=10,
        )
        txt.set_path_effects([pe.withStroke(linewidth=2.0, foreground="white", alpha=0.95)])

def annotate_public_nomenclature_scheme(
    ax,
    planform,
    y_sections: np.ndarray,
    c_y: np.ndarray,
    c_le_x: np.ndarray,
    c_te_x: np.ndarray,
    c_chords: np.ndarray,
    config,
    fixed_values: dict,
    public_reference: dict,
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

    def add_section_guides(x_end: float) -> None:
        guide_color = "#64748b"
        helper_points = np.asarray(config.planform.body_le_fixed_points, dtype=float)
        y_values = [float(y_sections[0])]
        y_values.extend(float(value) for value in helper_points[:, 1].tolist())
        y_values.extend(float(value) for value in y_sections[1:].tolist())
        unique_y_values = []
        for value in y_values:
            if not unique_y_values or abs(value - unique_y_values[-1]) > 1e-9:
                unique_y_values.append(value)

        for idx, y_value in enumerate(unique_y_values):
            ax.plot(
                [0.0, x_end],
                [y_value, y_value],
                color=guide_color,
                linewidth=0.95,
                linestyle=(0, (5, 4)),
                alpha=0.72,
                zorder=4,
            )
            x_label = 0.85
            y_offset = 0.38 if idx == 0 else 0.22
            outlined_text(
                x_label,
                y_value + y_offset,
                f"{idx}",
            )

    # Mark the public profile/chord stations shown in the CTA reference
    # sketch: C0, C1, C3, C4 and C5. Each station is centered above its own
    # LE->TE chord arrow to avoid piling labels up on the trailing edge.
    public_chord_labels = {
        0: ("Body chord", float(fixed_values["c0_body_chord_m"]), 1.18),
        3: ("Transition taper ratio", float(public_reference["c4_c3_ratio"]), 1.18),
        4: ("Wing tip", float(fixed_values["c5_wing_tip_m"]), 1.05),
    }
    for idx, name in enumerate(("C0", "C1", "C3", "C4", "C5")):
        yy = float(c_y[idx])
        le_x = float(c_le_x[idx])
        te_x = float(c_te_x[idx])
        y_here = 0.55 if yy < 0.35 else yy
        ax.annotate(
            "",
            xy=(te_x, y_here),
            xytext=(le_x, y_here),
            arrowprops={"arrowstyle": "<->", "color": "#475569", "linewidth": 1.05},
            zorder=5,
        )
        mid_x = float(0.5 * (le_x + te_x))
        label_y = y_here - (0.72 if idx == 0 else 0.58)
        outlined_text(mid_x, label_y, name)
        if idx in public_chord_labels:
            label, label_value, y_offset = public_chord_labels[idx]
            if idx == 3:
                outlined_text(
                    mid_x,
                    y_here + y_offset,
                    f"{label}\nC4/C3 = {label_value:.3f}\nC4 = {float(fixed_values['c4_wing_chord_m']):.3f} m",
                )
            else:
                outlined_text(mid_x, y_here + y_offset, f"{label}\n{label_value:.3f} m")

    wing_span = float(y_sections[3] - y_sections[1])
    b2_length = float(y_sections[2] - y_sections[1])

    span_annotations = (
        ("B1", float(y_sections[0]), float(y_sections[1]), f"{float(y_sections[1] - y_sections[0]):.2f} m"),
        ("B2", float(y_sections[1]), float(y_sections[2]), f"{100.0 * b2_length / max(wing_span, 1e-12):.1f}% wing span"),
        ("Wing span", float(y_sections[1]), float(y_sections[3]), f"{wing_span:.2f} m"),
    )
    span_arrow_positions: list[float] = []
    for idx, (name, y0, y1, label_text) in enumerate(span_annotations):
        y_mid = 0.5 * (y0 + y1)
        x_te_y0 = float(planform.te_x(y0))
        x_te_y1 = float(planform.te_x(y1))
        x_pos = max(x_te_y0, x_te_y1) + 1.5 + 0.62 * idx
        span_arrow_positions.append(x_pos)
        ax.annotate(
            "",
            xy=(x_pos, y1),
            xytext=(x_pos, y0),
            arrowprops={"arrowstyle": "<->", "color": "#334155", "linewidth": 0.95},
            zorder=6,
        )
        x_text = x_pos + 1.15 if name == "B2" else x_pos + 1.00
        outlined_text(x_text, y_mid, f"{name}={label_text}")

    x_guide_end = max(span_arrow_positions) if span_arrow_positions else float(np.max(c_te_x)) + 2.0
    add_section_guides(x_guide_end)

    sweep_segments = (
        (
            "S1 (50%)",
            float(c_y[2]),
            float(c_y[3]),
            float(c_le_x[2] + 0.50 * c_chords[2]),
            float(c_le_x[3] + 0.50 * c_chords[3]),
            float(public_reference["s2_deg"]),
        ),
        (
            "S2 (25%)",
            float(c_y[3]),
            float(c_y[4]),
            float(c_le_x[3] + 0.25 * c_chords[3]),
            float(c_le_x[4] + 0.25 * c_chords[4]),
            float(public_reference["s3_deg"]),
        ),
    )
    sweep_color = "#1d4ed8"
    for idx, (label, y0, y1, x0, x1, angle_deg) in enumerate(sweep_segments):
        dy = float(y1 - y0)
        dx = float(x1 - x0)
        seg_len = max(np.hypot(dx, dy), 1e-9)
        local_len = 0.35 * dy
        x_local = x0 + (dx / seg_len) * local_len
        y_local = y0 + (dy / seg_len) * local_len
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


def annotate_helper_geometry_scheme(
    ax,
    planform,
    y_sections: np.ndarray,
    le_points: np.ndarray,
    te_points: np.ndarray,
    c_y: np.ndarray,
    c_le_x: np.ndarray,
    c_te_x: np.ndarray,
    c_chords: np.ndarray,
    config,
    fixed_values: dict,
    public_reference: dict,
) -> None:
    annotate_public_nomenclature_scheme(
        ax,
        planform,
        y_sections,
        c_y,
        c_le_x,
        c_te_x,
        c_chords,
        config,
        fixed_values,
        public_reference,
    )

    helper_color = "#0891b2"
    helper_edge = "#0f172a"

    def add_local_sweep_marker(
        y0: float,
        label: str,
        radius: float,
        text_dx: float = 0.12,
        angle_override: float | None = None,
    ) -> None:
        x0 = float(planform.le_x(y0))
        if angle_override is None:
            dy_probe = 0.85
            y1 = min(float(y_sections[-1]), y0 + dy_probe)
            x1 = float(planform.le_x(y1))
            dy = max(y1 - y0, 1e-9)
            dx = x1 - x0
            angle_deg = float(np.degrees(np.arctan2(dx, dy)))
        else:
            angle_deg = float(angle_override)
            dy = 1.0
            dx = float(np.tan(np.deg2rad(angle_deg))) * dy

        ax.plot(
            [x0, x0],
            [y0, y0 + radius + 0.35],
            color="#1d4ed8",
            linewidth=1.7,
            linestyle=(0, (6, 5)),
            zorder=8,
        )
        seg_len = max(np.hypot(dx, dy), 1e-9)
        x_local = x0 + (dx / seg_len) * (radius + 0.12)
        y_local = y0 + (dy / seg_len) * (radius + 0.12)
        ax.plot([x0, x_local], [y0, y_local], color="#1d4ed8", linewidth=1.9, zorder=8)

        theta1 = max(0.0, 90.0 - angle_deg)
        theta2 = 90.0
        helper_arc = Arc(
            (x0, y0),
            width=2.0 * radius,
            height=2.0 * radius,
            angle=0.0,
            theta1=theta1,
            theta2=theta2,
            color="#1d4ed8",
            linewidth=2.0,
            zorder=9,
        )
        ax.add_patch(helper_arc)
        theta_mid = np.deg2rad(0.5 * (theta1 + theta2))
        tx = x0 + 1.12 * radius * np.cos(theta_mid)
        ty = y0 + 1.12 * radius * np.sin(theta_mid)
        ax.text(
            tx + text_dx,
            ty + 0.05,
            f"{label}={angle_deg:.2f}°",
            fontsize=10.0,
            color="#1d4ed8",
            ha="left",
            va="bottom",
            zorder=10,
        )

    le_helper_1_x = float(le_points[1, 0])
    le_helper_1_y = float(le_points[1, 1])
    le_helper_2_x = float(le_points[2, 0])
    le_helper_2_y = float(le_points[2, 1])

    for idx, (hx, hy) in enumerate(
        (
            (le_helper_1_x, le_helper_1_y),
            (le_helper_2_x, le_helper_2_y),
        )
    ):
        ax.scatter(
            [hx],
            [hy],
            s=56,
            marker="s",
            color=helper_color,
            edgecolors=helper_edge,
            linewidths=0.9,
            zorder=11,
        )
        ax.plot(
            [hx, hx],
            [0.0, hy],
            color=helper_color,
            linewidth=1.0,
            linestyle=(0, (5, 4)),
            zorder=9,
        )
        ax.plot(
            [0.0, hx],
            [hy, hy],
            color=helper_color,
            linewidth=1.0,
            linestyle=(0, (5, 4)),
            zorder=9,
        )

    c1_x = float(c_te_x[1])
    c1_y = float(c_y[1])
    add_local_sweep_marker(le_helper_1_y, "S", radius=1.18, text_dx=0.10, angle_override=66.87)
    add_local_sweep_marker(le_helper_2_y, "C1 sweep", radius=1.05, text_dx=0.14, angle_override=64.85)

    # Public C1 TE station used for the straight C0->C1 segment.
    ax.plot(
        [0.0, c1_x],
        [c1_y, c1_y],
        color="#94a3b8",
        linewidth=0.9,
        linestyle=(0, (2, 3)),
        zorder=4,
    )

    # Hidden TE helper that shapes the spline from C1 to C3.
    te_helper_x = float(te_points[2, 0])
    te_helper_y = float(te_points[2, 1])
    c3_x = float(c_te_x[2])
    c3_y = float(c_y[2])
    ax.scatter(
        [te_helper_x],
        [te_helper_y],
        s=56,
        marker="s",
        color=helper_color,
        edgecolors=helper_edge,
        linewidths=0.9,
        zorder=11,
    )
    ax.plot(
        [te_helper_x, te_helper_x],
        [te_helper_y, c3_y],
        color=helper_color,
        linewidth=1.0,
        linestyle=(0, (5, 4)),
        zorder=9,
    )
    ax.plot(
        [c3_x, te_helper_x],
        [c3_y + 0.02, c3_y + 0.02],
        color=helper_color,
        linewidth=1.0,
        linestyle=(0, (5, 4)),
        zorder=9,
    )


def draw_front(
    ax,
    dense_span: np.ndarray,
    upper: np.ndarray,
    lower: np.ndarray,
    anchor_y: np.ndarray,
    cargo_eng: dict,
    y_sections: np.ndarray,
    z_limits=None,
):
    del cargo_eng
    span = float(np.max(dense_span))
    b1 = float(y_sections[1])
    b2_end = float(y_sections[2])
    wing_span = float(y_sections[3] - y_sections[1])
    b2 = float(y_sections[2] - y_sections[1])

    zone_spans = (
        (-b1, b1, "#fde68a"),
        (-b2_end, -b1, "#bfdbfe"),
        (b1, b2_end, "#bfdbfe"),
        (-span, -b2_end, "#bbf7d0"),
        (b2_end, span, "#bbf7d0"),
    )
    for x0, x1, color in zone_spans:
        ax.axvspan(x0, x1, color=color, alpha=0.12, zorder=0)

    ax.fill_between(dense_span, lower, upper, color="#dce7f1", zorder=1, alpha=0.85)
    ax.fill_between(-dense_span, lower, upper, color="#dce7f1", zorder=1, alpha=0.45)
    ax.plot(dense_span, upper, color="#0f4c5c", linewidth=2.2, zorder=3)
    ax.plot(dense_span, lower, color="#0f4c5c", linewidth=2.2, zorder=3)
    ax.plot(-dense_span, upper, color="#0f4c5c", linewidth=1.7, zorder=3)
    ax.plot(-dense_span, lower, color="#0f4c5c", linewidth=1.7, zorder=3)

    for x_value in np.unique(np.concatenate((-anchor_y[1:], anchor_y[1:]))):
        ax.axvline(float(x_value), color="#64748b", linewidth=0.95, linestyle=(0, (4, 4)), zorder=4)

    ax.axhline(0.0, color="#94a3b8", linewidth=0.9, linestyle=":")
    ax.axvline(0.0, color="#94a3b8", linewidth=0.9, linestyle=":")
    z_span = max(1e-9, float(np.max(upper) - np.min(lower)))
    z_trace_max = float(np.max(upper))
    z_trace_min = float(np.min(lower))
    z_margin_top = 0.10 * z_span
    z_margin_bottom = 0.20 * z_span
    y_bracket_0 = z_trace_min - z_margin_bottom
    y_bracket_1 = y_bracket_0 + 0.55
    y_bracket_2 = y_bracket_0 + 1.10

    add_horizontal_bracket(ax, 0.0, b1, y_bracket_0, f"B1={b1:.3f} m")
    add_horizontal_bracket(ax, b1, b2_end, y_bracket_1, f"B2={b2:.3f} m")
    add_horizontal_bracket(ax, b1, span, y_bracket_2, f"Wing span={wing_span:.3f} m")

    label_y_base = z_trace_max + 0.025 * z_span
    for idx, x_value in enumerate(anchor_y):
        ax.text(
            float(x_value),
            label_y_base,
            f"{idx}",
            fontsize=11.0,
            color="#334155",
            ha="center",
            va="bottom",
            bbox={"boxstyle": "round,pad=0.14", "facecolor": "white", "edgecolor": "none", "alpha": 0.94},
            zorder=6,
        )

    ax.set_xlim(-span - 1.0, span + 1.0)
    if z_limits is None:
        ax.set_ylim(y_bracket_0 - 0.10, z_trace_max + z_margin_top)
    else:
        ax.set_ylim(float(z_limits[0]), float(z_limits[1]))
    ax.set_xlabel("spanwise y [m]")
    ax.set_ylabel("vertical [m]")
    ax.grid(True, linewidth=0.35, alpha=0.30)


def draw_front_half_le_te(
    ax,
    dense_span: np.ndarray,
    upper: np.ndarray,
    lower: np.ndarray,
    z_le: np.ndarray,
    z_te_upper: np.ndarray,
    z_te_lower: np.ndarray,
    z_te_mid: np.ndarray,
    anchor_y: np.ndarray,
    y_sections: np.ndarray,
):
    span = float(np.max(dense_span))
    b1 = float(y_sections[1])
    b2_end = float(y_sections[2])
    wing_span = float(y_sections[3] - y_sections[1])
    b2 = float(y_sections[2] - y_sections[1])

    zone_spans = (
        (0.0, b1, "#fde68a"),
        (b1, b2_end, "#bfdbfe"),
        (b2_end, span, "#bbf7d0"),
    )
    for x0, x1, color in zone_spans:
        ax.axvspan(x0, x1, color=color, alpha=0.14, zorder=0)

    ax.fill_between(dense_span, lower, upper, color="#dce7f1", zorder=1, alpha=0.88)
    ax.plot(dense_span, upper, color="#0f4c5c", linewidth=2.25, zorder=3)
    ax.plot(dense_span, lower, color="#0f4c5c", linewidth=2.25, zorder=3)

    ax.plot(dense_span, z_le, color="#7c3aed", linewidth=2.0, linestyle=":", zorder=4, label="LE z(y)")
    ax.plot(
        dense_span,
        z_te_upper,
        color="#ea580c",
        linewidth=1.8,
        linestyle="--",
        zorder=4,
        label="TE upper z(y)",
    )
    ax.plot(
        dense_span,
        z_te_lower,
        color="#ea580c",
        linewidth=1.8,
        linestyle="--",
        zorder=4,
        label="TE lower z(y)",
    )
    ax.plot(
        dense_span,
        z_te_mid,
        color="#475569",
        linewidth=1.5,
        linestyle="-.",
        zorder=4,
        label="TE mid z(y)",
    )

    for idx, x_value in enumerate(anchor_y):
        ax.axvline(float(x_value), color="#64748b", linewidth=0.95, linestyle=(0, (4, 4)), zorder=2, alpha=0.8)
        ax.text(
            float(x_value),
            float(np.max(upper)) + 0.025 * max(1e-9, float(np.max(upper) - np.min(lower))),
            f"{float(x_value):.3f} m",
            fontsize=9.0,
            color="#334155",
            ha="center",
            va="bottom",
            bbox={"boxstyle": "round,pad=0.10", "facecolor": "white", "edgecolor": "none", "alpha": 0.92},
            zorder=6,
        )

    ax.axhline(0.0, color="#94a3b8", linewidth=0.9, linestyle=":")
    z_span = max(1e-9, float(np.max(upper) - np.min(lower)))
    z_trace_max = float(np.max(np.concatenate([upper, z_le, z_te_upper, z_te_mid])))
    z_trace_min = float(np.min(np.concatenate([lower, z_te_lower, z_te_mid])))
    y_bracket_0 = z_trace_min - 0.22 * z_span
    y_bracket_1 = y_bracket_0
    y_bracket_2 = y_bracket_0 + 0.68
    add_horizontal_bracket(ax, 0.0, b1, y_bracket_0, f"B1={b1:.3f} m")
    add_horizontal_bracket(ax, b1, b2_end, y_bracket_1, f"B2={b2:.3f} m")
    add_horizontal_bracket(ax, b1, span, y_bracket_2, f"Wing span={wing_span:.3f} m")

    ax.set_xlim(0.0, span + 1.0)
    ax.set_ylim(y_bracket_0 - 0.08, z_trace_max + 0.08 * z_span)
    ax.set_xlabel("spanwise y [m]")
    ax.set_ylabel("vertical z [m]")
    ax.grid(True, linewidth=0.35, alpha=0.30)
    ax.legend(loc="upper right", ncol=2, fontsize=8, framealpha=0.94)


def draw_profiles(
    ax,
    prepared,
    y_sections: np.ndarray,
    chords: np.ndarray,
    section_labels: tuple[str, ...],
):
    colors = (
        "#0f4c5c",
        "#2563eb",
        "#7c3aed",
        "#c026d3",
        "#ea580c",
        "#c2410c",
    )
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

    color_cycle = [colors[idx % len(colors)] for idx in range(len(raw_profiles))]
    for (x_local, zu, zl, yy, chord), off, color, label in zip(raw_profiles, offsets, color_cycle, section_labels):
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


def draw_definition_panel(ax, config, c_names: tuple[str, ...], c_chords: np.ndarray, fixed_values: dict) -> None:
    ax.axis("off")

    fixed_lines = [
        "Fixed CTA geometry",
        f"S = {float(fixed_values.get('s_deg', config.planform.s1_deg)):.1f} deg",
        f"B1 = {float(config.topology.y_sections_array[1]):.2f} m",
        f"Body chord = {float(fixed_values['c0_body_chord_m']):.3f} m",
        f"Wing chord = {float(fixed_values['c4_wing_chord_m']):.3f} m",
        f"Wing tip = {float(fixed_values['c5_wing_tip_m']):.3f} m",
    ]

    active_lines = [
        "Active AI design variables",
        "span",
        "B2/(B2+B3)",
        "S1 = s2_deg",
        "S2 = s3_deg",
        "twist at C0, C3, C4, C5",
        "CST upper/lower coefficients for C0, C3, C4, C5",
    ]

    def write_block(x0: float, y0: float, title_lines: list[str], face: str, edge: str) -> float:
        text = "\n".join(title_lines)
        ax.text(
            x0,
            y0,
            text,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=10.0,
            color="#0f172a",
            linespacing=1.45,
            bbox={"boxstyle": "round,pad=0.42", "facecolor": face, "edgecolor": edge, "alpha": 0.96},
        )
        return y0

    write_block(0.02, 0.98, fixed_lines, "#fff7ed", "#fdba74")
    write_block(0.02, 0.52, active_lines, "#eff6ff", "#93c5fd")


def main() -> None:
    design = build_reference_design()
    fixed_values = cta_fixed_values(reference_design=design)
    public_reference = build_cta_design_space(reference_design=design).reference_flat()
    config = to_cta_model_config(design, use_reference_anchor_twist=True)
    prepared = prepare_geometry(config)
    planform = prepared.planform

    output_dir = SCRIPT_DIR.parent / "outputs" / "reference"
    output_dir.mkdir(parents=True, exist_ok=True)

    layout_png = output_dir / "cta_reference_layout.png"
    plan_png = output_dir / "cta_reference_planform_cargo.png"
    front_png = output_dir / "cta_reference_front.png"
    front_half_png = output_dir / "cta_reference_front_half_le_te.png"
    profiles_png = output_dir / "cta_reference_profiles.png"
    view3d_png = output_dir / "cta_reference_3d.png"
    scheme_png = output_dir / "cta_reference_definition_scheme.png"
    helper_scheme_png = output_dir / "cta_reference_helper_scheme.png"

    for stale_svg in (
        output_dir / "cta_reference_layout.svg",
        output_dir / "cta_reference_planform_cargo.svg",
        output_dir / "cta_reference_front.svg",
        output_dir / "cta_reference_profiles.svg",
    ):
        if stale_svg.exists():
            stale_svg.unlink()

    y_sections = config.topology.y_sections_array
    profile_anchor_y = config.topology.anchor_y_array
    dense_span = np.unique(
        np.concatenate([np.linspace(0.0, config.topology.span, 2400), y_sections, profile_anchor_y])
    )
    leading_edge_x = np.array([planform.le_x(float(y)) for y in dense_span], dtype=float)
    trailing_edge_x = np.array([planform.te_x(float(y)) for y in dense_span], dtype=float)
    chord_dense = trailing_edge_x - leading_edge_x
    le_sections = config.planform.leading_edge_x_sections(config.topology)
    te_sections = config.planform.trailing_edge_x_sections(config.topology)
    chords = te_sections - le_sections
    te_points = config.planform.trailing_edge_points(config.topology)
    public_station_indices = np.asarray((0, 1, 3, 4, te_points.shape[0] - 1), dtype=int)
    c_names = ("C0", "C1", "C3", "C4", "C5")
    c_y = te_points[public_station_indices, 1].astype(float)
    c_te_x = te_points[public_station_indices, 0].astype(float)
    c_le_x = np.array([planform.le_x(float(value)) for value in c_y], dtype=float)
    c_chords = c_te_x - c_le_x
    profile_labels = (
        "Section 0 / C0",
        "Section 1",
        "Section 2 / C1",
        "Section 3 / C3",
        "Section 4 / C4",
        "Section 5 / C5",
    )
    profile_chords = np.array(
        [float(planform.te_x(float(yy)) - planform.le_x(float(yy))) for yy in profile_anchor_y],
        dtype=float,
    )
    cargo_eng = build_cargo_and_engines(planform, y_sections)
    front_anchor_y = profile_anchor_y
    upper, lower, _, _ = front_linear_envelope(prepared, dense_span, front_anchor_y)
    z_le, z_te_upper, z_te_lower, z_te_mid = front_linear_edge_traces(prepared, dense_span, front_anchor_y)
    glider_z_limits = glider_reference_z_limits()

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
        public_reference,
    )
    ax_l_plan.set_title("CTA reference planform (top view, half-wing)")
    ax_l_plan.legend(loc="upper left", ncol=1, fontsize=10)

    save_figure(fig_layout, layout_png, dpi=230, bbox_inches="tight")
    plt.close(fig_layout)

    # Planform with cargo/engines
    fig_plan, ax_plan = plt.subplots(figsize=(12.0, 7.0), constrained_layout=True)
    draw_planform(ax_plan, dense_span, leading_edge_x, trailing_edge_x, y_sections, cargo_eng, show_cargo=True)
    ax_plan.set_title("CTA planform (true scale) with cargo and engines")
    ax_plan.legend(loc="upper left", ncol=2, fontsize=8)
    save_figure(fig_plan, plan_png, dpi=230, bbox_inches="tight")
    plt.close(fig_plan)

    # Front only
    fig_front, ax_front = plt.subplots(figsize=(14.0, 5.4), constrained_layout=True)
    draw_front(
        ax_front,
        dense_span,
        upper,
        lower,
        front_anchor_y,
        cargo_eng,
        y_sections,
        z_limits=glider_z_limits,
    )
    ax_front.set_title("CTA front view with reference span measures")
    save_figure(fig_front, front_png, dpi=220)
    plt.close(fig_front)

    fig_front_half, ax_front_half = plt.subplots(figsize=(14.0, 5.4), constrained_layout=True)
    draw_front_half_le_te(
        ax_front_half,
        dense_span,
        upper,
        lower,
        z_le,
        z_te_upper,
        z_te_lower,
        z_te_mid,
        front_anchor_y,
        y_sections,
    )
    ax_front_half.set_title("CTA half-wing front view with LE/TE z evolution")
    save_figure(fig_front_half, front_half_png, dpi=220)
    plt.close(fig_front_half)

    # Profiles only
    fig_profiles, ax_profiles = plt.subplots(figsize=(11.6, 8.1), constrained_layout=True)
    draw_profiles(ax_profiles, prepared, profile_anchor_y, profile_chords, profile_labels)
    save_figure(fig_profiles, profiles_png, dpi=230, bbox_inches="tight")
    plt.close(fig_profiles)

    # Final definition scheme: annotated planform only.
    fig_scheme, ax_scheme_plan = plt.subplots(figsize=(12.8, 8.6), constrained_layout=True)
    draw_planform(
        ax_scheme_plan,
        dense_span,
        leading_edge_x,
        trailing_edge_x,
        y_sections,
        cargo_eng,
        show_cargo=False,
        half_wing=True,
        label_zones=True,
        faded_symmetry=True,
    )
    annotate_public_nomenclature_scheme(
        ax_scheme_plan,
        planform,
        y_sections,
        c_y,
        c_le_x,
        c_te_x,
        c_chords,
        config,
        fixed_values,
        public_reference,
    )
    ax_scheme_plan.set_title("CTA final geometric definition scheme")
    save_figure(fig_scheme, scheme_png, dpi=230, bbox_inches="tight")
    plt.close(fig_scheme)

    fig_helper, ax_helper = plt.subplots(figsize=(13.1, 9.0), constrained_layout=True)
    draw_planform(
        ax_helper,
        dense_span,
        leading_edge_x,
        trailing_edge_x,
        y_sections,
        cargo_eng,
        show_cargo=False,
        half_wing=True,
        label_zones=True,
        faded_symmetry=True,
    )
    annotate_helper_geometry_scheme(
        ax_helper,
        planform,
        y_sections,
        config.planform.leading_edge_points(config.topology),
        te_points,
        c_y,
        c_le_x,
        c_te_x,
        c_chords,
        config,
        fixed_values,
        public_reference,
    )
    ax_helper.set_title("CTA helper geometry definition scheme")
    save_figure(fig_helper, helper_scheme_png, dpi=230, bbox_inches="tight")
    plt.close(fig_helper)

    # 3D inclined view
    fig_3d = plt.figure(figsize=(12.0, 8.2), constrained_layout=True)
    ax_3d = fig_3d.add_subplot(111, projection="3d")
    draw_3d(ax_3d, prepared, dense_span, leading_edge_x, chord_dense)
    save_figure(fig_3d, view3d_png, dpi=240, bbox_inches="tight")
    plt.close(fig_3d)

    print(f"CTA layout PNG written to: {layout_png}")
    print(f"CTA planform PNG written to: {plan_png}")
    print(f"CTA front PNG written to: {front_png}")
    print(f"CTA half-wing front PNG written to: {front_half_png}")
    print(f"CTA profiles PNG written to: {profiles_png}")
    print(f"CTA final definition scheme PNG written to: {scheme_png}")
    print(f"CTA helper geometry scheme PNG written to: {helper_scheme_png}")
    print(f"CTA 3D PNG written to: {view3d_png}")
    print(f"Cargo area: {cargo_eng['cargo_area']:.3f} m^2")
    print(f"Engines area (rear remainder): {cargo_eng['engine_area']:.3f} m^2")
    print(f"Section interpolation: {prepared.section_model.interpolation_name}")
    print(f"Twist interpolation: {config.spanwise.twist_deg.interpolation}")
    print(f"Camber interpolation: {config.spanwise.camber_delta.interpolation}")


if __name__ == "__main__":
    main()
