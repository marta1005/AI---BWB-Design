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

from project_v4.builder import prepare_geometry
from project_v4.examples.reference.run_reference_example import build_reference_design


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


def section_max_tc(section_model, y: float) -> float:
    yu, yl, _ = section_model.coordinates_at_y(float(y))
    return float(np.max(yu - yl))


def section_profile_geometry(section_model, y: float, chord: float):
    yu, yl, params = section_model.coordinates_at_y(float(y))
    x = section_model.x_air * chord
    zu = yu * chord
    zl = yl * chord
    return x, zu, zl, params


def rotate_profile_with_twist(x: np.ndarray, zu: np.ndarray, zl: np.ndarray, twist_deg: float):
    twist_rad = np.deg2rad(float(twist_deg))
    cos_twist = np.cos(twist_rad)
    sin_twist = np.sin(twist_rad)
    xu = x * cos_twist - zu * sin_twist
    zu_rot = x * sin_twist + zu * cos_twist
    xl = x * cos_twist - zl * sin_twist
    zl_rot = x * sin_twist + zl * cos_twist
    return xu, zu_rot, xl, zl_rot


def front_projected_profile(
    section_model,
    y: float,
    chord: float,
    twist_deg: float,
    vertical_offset: float,
):
    x, zu, zl, _ = section_profile_geometry(section_model, y, chord)
    twist_rad = np.deg2rad(float(twist_deg))
    sin_twist = np.sin(twist_rad)
    cos_twist = np.cos(twist_rad)
    upper = vertical_offset + x * sin_twist + zu * cos_twist
    lower = vertical_offset + x * sin_twist + zl * cos_twist
    te_center = vertical_offset + chord * sin_twist
    return upper, lower, te_center


def add_span_bracket(ax, x: float, y0: float, y1: float, label: str, color: str = "#475569") -> None:
    ax.annotate(
        "",
        xy=(x, y1),
        xytext=(x, y0),
        arrowprops={"arrowstyle": "<->", "color": color, "linewidth": 1.1},
        zorder=5,
    )
    ax.text(
        x - 0.18,
        0.5 * (y0 + y1),
        label,
        rotation=90,
        ha="right",
        va="center",
        fontsize=8.2,
        color=color,
        bbox={"boxstyle": "round,pad=0.14", "facecolor": "white", "edgecolor": "none", "alpha": 0.82},
        zorder=6,
    )


def add_sweep_label(ax, x0: float, y0: float, x1: float, y1: float, text: str) -> None:
    dx = float(x1 - x0)
    dy = float(y1 - y0)
    seg_len = max(np.hypot(dx, dy), 1e-9)
    nx = -dy / seg_len
    ny = dx / seg_len
    xm = 0.5 * (x0 + x1)
    ym = 0.5 * (y0 + y1)
    ax.plot([x0, x1], [y0, y1], color="#0f4c5c", linewidth=0.8, linestyle=":", alpha=0.45, zorder=2)
    ax.text(
        xm + 1.4 * nx,
        ym + 1.4 * ny,
        text,
        fontsize=8.2,
        color="#0f4c5c",
        ha="center",
        va="center",
        bbox={"boxstyle": "round,pad=0.16", "facecolor": "white", "edgecolor": "#cbd5e1", "alpha": 0.9},
        zorder=6,
    )


def main() -> None:
    design = build_reference_design()
    config = design.to_model_config()

    output_dir = SCRIPT_DIR.parent.parent / "example_outputs" / "reference_v4_demo"
    output_dir.mkdir(parents=True, exist_ok=True)
    png_path = output_dir / "reference_v4_layout.png"
    svg_path = output_dir / "reference_v4_layout.svg"
    profiles_png_path = output_dir / "reference_v4_profiles.png"
    profiles_svg_path = output_dir / "reference_v4_profiles.svg"

    prepared = prepare_geometry(config)
    planform = prepared.planform
    dense_span = np.unique(
        np.concatenate(
            [
                np.linspace(0.0, config.topology.span, 600),
                config.topology.y_sections_array,
            ]
        )
    )
    leading_edge_x = np.array([planform.le_x(float(y)) for y in dense_span], dtype=float)
    trailing_edge_x = np.array([planform.te_x(float(y)) for y in dense_span], dtype=float)
    twist_dense = np.array([prepared.spanwise_laws.twist_deg(float(y)) for y in dense_span], dtype=float)
    vertical_center = np.interp(dense_span, prepared.loft.span_stations, prepared.loft.vertical_y)

    fig, (ax, ax_front) = plt.subplots(
        2,
        1,
        figsize=(12.2, 12.8),
        constrained_layout=True,
        gridspec_kw={"height_ratios": (3.2, 2.0)},
    )
    y_sections = config.topology.y_sections_array
    le_sections = config.planform.leading_edge_x_sections(config.topology)
    te_sections = config.planform.trailing_edge_x_sections(config.topology)
    te_points = config.planform.trailing_edge_points(config.topology)
    chords = te_sections - le_sections

    section_labels = ("C1", "C2", "C3", "C4")
    zone_defs = (
        ("Centre Body", float(y_sections[0]), float(y_sections[1]), "#fde68a"),
        ("Transition Wing", float(y_sections[1]), float(y_sections[2]), "#bfdbfe"),
        ("Outer Wing", float(y_sections[2]), float(y_sections[3]), "#bbf7d0"),
    )
    for zone_name, y0, y1, zone_color in zone_defs:
        mask = (dense_span >= y0) & (dense_span <= y1)
        ax.fill_betweenx(
            dense_span[mask],
            leading_edge_x[mask],
            trailing_edge_x[mask],
            color=zone_color,
            alpha=0.68,
            zorder=1,
            label=zone_name,
        )
        ax.fill_betweenx(
            -dense_span[mask],
            leading_edge_x[mask],
            trailing_edge_x[mask],
            color=zone_color,
            alpha=0.25,
            zorder=1,
        )
    ax.plot(leading_edge_x, dense_span, color="#0f4c5c", linewidth=2.2, label="LE", zorder=3)
    ax.plot(trailing_edge_x, dense_span, color="#c44536", linewidth=2.2, label="TE", zorder=3)
    ax.plot(leading_edge_x, -dense_span, color="#0f4c5c", linewidth=2.2, zorder=3)
    ax.plot(trailing_edge_x, -dense_span, color="#c44536", linewidth=2.2, zorder=3)

    c_names = ("C0", "C1", "C2", "C3", "C4", "C5")
    c_y = te_points[:, 1].astype(float)
    c_te_x = te_points[:, 0].astype(float)
    c_le_x = np.array([planform.le_x(float(value)) for value in c_y], dtype=float)
    c_chords = c_te_x - c_le_x
    te_exact_segments = (
        (
            np.array([float(te_points[0, 0]), float(te_points[1, 0])], dtype=float),
            np.array([float(te_points[0, 1]), float(te_points[1, 1])], dtype=float),
            "TE exact straight segment C0->C1",
            "TE C0->C1\nexact straight",
            (1.2, 0.9),
        ),
        (
            np.array([float(te_points[3, 0]), float(te_points[4, 0])], dtype=float),
            np.array([float(te_points[3, 1]), float(te_points[4, 1])], dtype=float),
            "TE exact straight segment C3->C4",
            "TE C3->C4\nexact straight",
            (1.2, 0.7),
        ),
    )

    for name, y, le_x, te_x, chord in zip(c_names, c_y, c_le_x, c_te_x, c_chords):
        ax.plot([le_x, te_x], [y, y], color="#64748b", linewidth=1.0, alpha=0.9, zorder=2)
        ax.plot([le_x, te_x], [-y, -y], color="#64748b", linewidth=1.0, alpha=0.22, zorder=2)
        ax.scatter([le_x, te_x], [y, y], color=["#0f4c5c", "#c44536"], s=16, zorder=4)
        ax.text(
            te_x + 0.35,
            y + 0.18,
            f"{name}: y={y:.2f} m, c={chord:.2f} m",
            fontsize=8.0,
            color="#334155",
            ha="left",
            va="center",
            bbox={"boxstyle": "round,pad=0.14", "facecolor": "white", "edgecolor": "none", "alpha": 0.80},
            zorder=6,
        )

    for exact_x, exact_y, legend_label, text_label, text_offset in te_exact_segments:
        ax.plot(
            exact_x,
            exact_y,
            color="#111827",
            linewidth=3.0,
            linestyle=(0, (6, 4)),
            label=legend_label,
            zorder=5,
        )
        ax.plot(
            exact_x,
            -exact_y,
            color="#111827",
            linewidth=3.0,
            linestyle=(0, (6, 4)),
            zorder=5,
        )
        te_mid_x = float(np.mean(exact_x))
        te_mid_y = float(np.mean(exact_y))
        ax.text(
            te_mid_x + float(text_offset[0]),
            te_mid_y + float(text_offset[1]),
            text_label,
            fontsize=8.2,
            color="#111827",
            ha="left",
            va="bottom",
            bbox={"boxstyle": "round,pad=0.16", "facecolor": "white", "edgecolor": "#d1d5db", "alpha": 0.92},
            zorder=6,
        )

    add_sweep_label(
        ax,
        le_sections[0],
        y_sections[0],
        le_sections[1],
        y_sections[1],
        f"S1 fixed = {config.planform.s1_deg:.1f} deg",
    )
    add_sweep_label(ax, le_sections[1], y_sections[1], le_sections[2], y_sections[2], f"S2 = {config.planform.s2_deg:.1f} deg")
    add_sweep_label(ax, le_sections[2], y_sections[2], le_sections[3], y_sections[3], f"S3 = {config.planform.s3_deg:.1f} deg")

    bracket_x = float(np.min(leading_edge_x) - 3.8)
    add_span_bracket(
        ax,
        bracket_x,
        float(y_sections[0]),
        float(y_sections[1]),
        f"B1 = {config.topology.b1_span_ratio:.2f} b/2\n({float(y_sections[1] - y_sections[0]):.2f} m)",
    )
    add_span_bracket(
        ax,
        bracket_x,
        float(y_sections[1]),
        float(y_sections[2]),
        f"B2 = {config.topology.b2_span_ratio:.2f} b/2\n({float(y_sections[2] - y_sections[1]):.2f} m)",
    )
    add_span_bracket(
        ax,
        bracket_x,
        float(y_sections[2]),
        float(y_sections[3]),
        f"B3 = {config.topology.b3_span_ratio:.2f} b/2\n({float(y_sections[3] - y_sections[2]):.2f} m)",
    )

    ax.axhline(0.0, color="#9fb3c8", linewidth=0.9, linestyle="--", zorder=0)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x / chord direction")
    ax.set_ylabel("spanwise y")
    ax.set_title("CTA reference planform (platform parameters only)")
    ax.legend(loc="upper left")
    ax.grid(True, linewidth=0.4, alpha=0.25)

    span_margin = 0.06 * config.topology.span
    chord_min = float(np.min(leading_edge_x))
    chord_max = float(np.max(trailing_edge_x))
    chord_margin = 0.06 * (chord_max - chord_min)
    ax.set_xlim(bracket_x - 1.2, float(np.max(c_te_x) + 9.0))
    ax.set_ylim(-(config.topology.span + span_margin), config.topology.span + span_margin)
    upper = np.empty_like(dense_span)
    lower = np.empty_like(dense_span)
    te_center = np.empty_like(dense_span)
    for idx, (yy, chord, twist_here, vertical_here) in enumerate(
        zip(dense_span, trailing_edge_x - leading_edge_x, twist_dense, vertical_center)
    ):
        upper_here, lower_here, te_center_here = front_projected_profile(
            prepared.section_model,
            float(yy),
            float(chord),
            float(twist_here),
            float(vertical_here),
        )
        upper[idx] = float(np.max(upper_here))
        lower[idx] = float(np.min(lower_here))
        te_center[idx] = float(te_center_here)

    ax_front.fill_between(dense_span, lower, upper, color="#dce7f1", zorder=1)
    ax_front.plot(dense_span, upper, color="#0f4c5c", linewidth=2.0, zorder=3)
    ax_front.plot(dense_span, lower, color="#0f4c5c", linewidth=2.0, zorder=3)
    ax_front.plot(dense_span, vertical_center, color="#475569", linewidth=1.1, linestyle="--", label="LE height", zorder=4)
    ax_front.plot(dense_span, te_center, color="#c44536", linewidth=1.4, label="TE height with twist", zorder=4)

    c_upper = np.interp(c_y, dense_span, upper)
    c_lower = np.interp(c_y, dense_span, lower)
    c_te = np.interp(c_y, dense_span, te_center)
    for name, y, y_top, y_bot, y_te, y_center_here in zip(
        c_names,
        c_y,
        c_upper,
        c_lower,
        c_te,
        np.interp(c_y, dense_span, vertical_center),
    ):
        ax_front.plot([float(y), float(y)], [float(y_bot), float(y_top)], color="#64748b", linewidth=0.85, linestyle="--", zorder=2)
        ax_front.scatter([float(y), float(y)], [float(y_center_here), float(y_te)], color=["#475569", "#c44536"], s=14, zorder=5)
        ax_front.text(
            float(y),
            float(y_top + 0.20),
            name,
            fontsize=9.0,
            color="#334155",
            ha="center",
            va="bottom",
            bbox={"boxstyle": "round,pad=0.12", "facecolor": "white", "edgecolor": "none", "alpha": 0.80},
            zorder=6,
        )

    ax_front.axhline(0.0, color="#94a3b8", linewidth=0.9, linestyle=":")
    ax_front.axvline(0.0, color="#94a3b8", linewidth=0.9, linestyle=":")
    ax_front.set_title("Front view (CTA station nomenclature)", pad=12)
    ax_front.set_xlabel("semispan y")
    ax_front.set_ylabel("vertical")
    ax_front.grid(True, linewidth=0.4, alpha=0.3)
    ax_front.set_xlim(-0.8, config.topology.span + 1.5)
    max_upper = float(np.max(upper))
    min_lower = float(np.min(lower))
    pad = 0.10 * max(1e-9, (max_upper - min_lower))
    ax_front.set_ylim(min_lower - 1.1 * pad, max_upper + 1.1 * pad)
    ax_front.set_aspect("equal", adjustable="box")
    ax_front.legend(loc="upper right")

    save_figure(fig, png_path, dpi=220, bbox_inches="tight")
    save_figure(fig, svg_path, bbox_inches="tight")
    plt.close(fig)

    fig_profiles, ax_profiles = plt.subplots(figsize=(11.5, 8.0), constrained_layout=True)
    profile_colors = ("#0f4c5c", "#1d4ed8", "#7c3aed", "#c44536")
    max_profile_thickness = []
    raw_profiles = []
    profile_params = []
    profile_metrics = []
    for y, chord, twist_deg in zip(y_sections, chords, config.spanwise.twist_deg.values):
        x_profile, zu_profile, zl_profile, params = section_profile_geometry(
            prepared.section_model, float(y), float(chord)
        )
        metrics, _ = prepared.section_model.geometry_metrics_at_y(float(y))
        xu, zu_rot, xl, zl_rot = rotate_profile_with_twist(x_profile, zu_profile, zl_profile, float(twist_deg))
        raw_profiles.append((xu, zu_rot, xl, zl_rot))
        profile_params.append(params)
        profile_metrics.append(metrics)
        max_profile_thickness.append(float(max(np.max(zu_rot), np.max(zl_rot)) - min(np.min(zu_rot), np.min(zl_rot))))

    profile_offsets = [0.0]
    for idx in range(1, len(raw_profiles)):
        clearance = 0.65 * (max_profile_thickness[idx - 1] + max_profile_thickness[idx]) + 1.0
        profile_offsets.append(profile_offsets[-1] - clearance)

    label_x = float(np.max(chords) + 3.0)
    for idx, ((label, y, chord), (xu, zu_profile, xl, zl_profile), params, metrics, offset, color) in enumerate(
        zip(zip(section_labels, y_sections, chords), raw_profiles, profile_params, profile_metrics, profile_offsets, profile_colors)
    ):
        zu_shift = zu_profile + offset
        zl_shift = zl_profile + offset
        camber_line = 0.5 * (zu_shift + zl_shift)
        x_mid = 0.5 * (xu + xl)
        upper_order = np.argsort(xu)
        lower_order = np.argsort(xl)
        camber_order = np.argsort(x_mid)
        ax_profiles.fill(
            np.concatenate([xu[upper_order], xl[lower_order][::-1]]),
            np.concatenate([zu_shift[upper_order], zl_shift[lower_order][::-1]]),
            color=color,
            alpha=0.12,
            zorder=1,
        )
        ax_profiles.plot(xu[upper_order], zu_shift[upper_order], color=color, linewidth=2.0, zorder=3)
        ax_profiles.plot(xl[lower_order], zl_shift[lower_order], color=color, linewidth=2.0, zorder=3)
        ax_profiles.plot(
            x_mid[camber_order],
            camber_line[camber_order],
            color=color,
            linewidth=1.0,
            linestyle="--",
            alpha=0.9,
            zorder=4,
        )
        x_ref_max = float(max(np.max(xu), np.max(xl)))
        ax_profiles.plot([0.0, x_ref_max], [offset, offset], color="#94a3b8", linewidth=0.8, linestyle=":", zorder=0)
        ax_profiles.plot([x_ref_max + 0.2, label_x - 0.2], [offset, offset], color="#94a3b8", linewidth=0.8, linestyle=":", zorder=0)
        max_camber = float(np.max(np.abs(camber_line - offset)))
        ax_profiles.text(
            label_x,
            offset,
            f"{label}: c={float(chord):.2f} m | y={float(y):.2f} m | twist={float(config.spanwise.twist_deg.values[idx]):.1f} deg\n"
            f"tc_measured={metrics.max_tc:.3f} | x_tmax_measured={metrics.x_tmax:.3f} | te={metrics.te_thickness:.4f} | max camber={max_camber / float(chord):.4f}",
            fontsize=8.3,
            color="#243b53",
            va="center",
            bbox={"boxstyle": "round,pad=0.16", "facecolor": "white", "edgecolor": "none", "alpha": 0.82},
            zorder=5,
        )

    ax_profiles.text(
        0.6,
        profile_offsets[0] + 0.45 * max_profile_thickness[0],
        "Profiles at true scale\nLE aligned at x = 0 | twist applied | dashed line = camber",
        fontsize=8.6,
        color="#263238",
        va="top",
        bbox={"boxstyle": "round,pad=0.20", "facecolor": "white", "edgecolor": "#cfd8dc", "alpha": 0.90},
    )
    ax_profiles.set_title("Section profiles (true scale)")
    ax_profiles.set_xlabel("x / chord direction")
    ax_profiles.set_ylabel("profile vertical")
    ax_profiles.grid(True, linewidth=0.4, alpha=0.25)
    ax_profiles.set_aspect("equal", adjustable="box")
    z_upper_max = max(float(max(np.max(zu), np.max(zl)) + offset) for (_, zu, _, zl), offset in zip(raw_profiles, profile_offsets))
    z_lower_min = min(float(min(np.min(zu), np.min(zl)) + offset) for (_, zu, _, zl), offset in zip(raw_profiles, profile_offsets))
    x_right_max = max(float(max(np.max(xu), np.max(xl))) for (xu, _, xl, _) in raw_profiles)
    x_left_min = min(float(min(np.min(xu), np.min(xl))) for (xu, _, xl, _) in raw_profiles)
    ax_profiles.set_xlim(min(-0.8, x_left_min - 0.8), max(label_x + 10.0, x_right_max + 10.0))
    ax_profiles.set_ylim(z_lower_min - 1.0, z_upper_max + 1.0)
    ax_profiles.set_yticks([])

    save_figure(fig_profiles, profiles_png_path, dpi=220, bbox_inches="tight")
    save_figure(fig_profiles, profiles_svg_path, bbox_inches="tight")
    plt.close(fig_profiles)

    print(f"Combined layout PNG written to: {png_path}")
    print(f"Combined layout SVG written to: {svg_path}")
    print(f"Profiles PNG written to: {profiles_png_path}")
    print(f"Profiles SVG written to: {profiles_svg_path}")


if __name__ == "__main__":
    main()
