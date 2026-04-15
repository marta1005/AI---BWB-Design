from __future__ import annotations

from collections import OrderedDict
import os
from pathlib import Path
import sys
from typing import Dict

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
from scipy.interpolate import PchipInterpolator

from parametrization.shared.cst import cosine_spacing
from parametrization.shared.airfoil_fit import fit_airfoil_section_cst
from parametrization.CTA.codes.airfoils.fit_bwb_glider_cst import FIT_OPTIONS


def _load_sections(path: Path) -> "OrderedDict[float, np.ndarray]":
    sections: "OrderedDict[float, list[list[float]]]" = OrderedDict()
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        parts = raw_line.split()
        if len(parts) != 3:
            continue
        x_val, y_val, z_val = (float(parts[0]), float(parts[1]), float(parts[2]))
        sections.setdefault(y_val, []).append([x_val, z_val])
    return OrderedDict((yy, np.asarray(points, dtype=float)) for yy, points in sections.items())


def _section_metrics(section: np.ndarray) -> Dict[str, float]:
    x_values = section[:, 0]
    z_values = section[:, 1]
    x_le = float(np.min(x_values))
    x_te = float(np.max(x_values))
    chord = x_te - x_le
    thickness = float(np.max(z_values) - np.min(z_values))
    return {
        "x_le": x_le,
        "x_te": x_te,
        "chord": chord,
        "thickness": thickness,
    }


def _normalized_section(section: np.ndarray) -> np.ndarray:
    metrics = _section_metrics(section)
    chord = max(metrics["chord"], 1e-12)
    normalized = np.empty_like(section, dtype=float)
    normalized[:, 0] = (section[:, 0] - metrics["x_le"]) / chord
    normalized[:, 1] = section[:, 1] / chord
    return normalized


def _normalized_local_surface_points(
    section: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    le_idx = int(np.argmin(section[:, 0]))
    upper_raw = section[: le_idx + 1]
    lower_raw = section[le_idx:]

    te_upper = np.asarray(section[0], dtype=float)
    te_lower = np.asarray(section[-1], dtype=float)
    te_mid = 0.5 * (te_upper + te_lower)
    le_point = np.asarray(section[le_idx], dtype=float)
    chord_vector = te_mid - le_point
    chord = float(np.linalg.norm(chord_vector))
    if chord <= 0.0:
        raise ValueError(f"Invalid chord length {chord}")

    ex = chord_vector / chord
    ez = np.array([-ex[1], ex[0]], dtype=float)

    # Keep the original section points in the local chord frame instead of
    # interpolating in x/c. Reinterpolation was visually flattening the leading
    # edge, especially in the grid/overlay inspection plots.
    upper_local = np.column_stack(((upper_raw - le_point) @ ex, (upper_raw - le_point) @ ez))[::-1]
    lower_local = np.column_stack(((lower_raw - le_point) @ ex, (lower_raw - le_point) @ ez))

    upper_x = upper_local[:, 0] / chord
    upper_z = upper_local[:, 1] / chord
    lower_x = lower_local[:, 0] / chord
    lower_z = lower_local[:, 1] / chord
    return upper_x, upper_z, lower_x, lower_z


def _section_local_surface_points(
    section: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    upper_raw, lower_raw = _split_raw_section_surfaces(section)

    te_upper = np.asarray(section[0], dtype=float)
    te_lower = np.asarray(section[-1], dtype=float)
    te_mid = 0.5 * (te_upper + te_lower)
    le_point = np.asarray(section[int(np.argmin(section[:, 0]))], dtype=float)
    chord_vector = te_mid - le_point
    chord = float(np.linalg.norm(chord_vector))
    if chord <= 0.0:
        raise ValueError(f"Invalid chord length {chord}")

    ex = chord_vector / chord
    ez = np.array([-ex[1], ex[0]], dtype=float)
    upper_local = np.column_stack(((upper_raw - le_point) @ ex, (upper_raw - le_point) @ ez))
    lower_local = np.column_stack(((lower_raw - le_point) @ ex, (lower_raw - le_point) @ ez))
    return upper_local, lower_local


def _section_local_closed_curve(section: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    te_upper = np.asarray(section[0], dtype=float)
    te_lower = np.asarray(section[-1], dtype=float)
    te_mid = 0.5 * (te_upper + te_lower)
    le_point = np.asarray(section[int(np.argmin(section[:, 0]))], dtype=float)
    chord_vector = te_mid - le_point
    chord = float(np.linalg.norm(chord_vector))
    if chord <= 0.0:
        raise ValueError(f"Invalid chord length {chord}")

    ex = chord_vector / chord
    ez = np.array([-ex[1], ex[0]], dtype=float)
    section_local = np.column_stack(((section - le_point) @ ex, (section - le_point) @ ez))

    # Follow the original contour order TE_upper -> LE -> TE_lower and close the
    # profile only at the trailing edge, using the raw points exactly as they
    # appear in the source file.
    closed_local = np.vstack([section_local, section_local[:1]])
    return closed_local[:, 0], closed_local[:, 1]


def _split_raw_section_surfaces(section: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    le_idx = int(np.argmin(section[:, 0]))
    upper_raw = section[: le_idx + 1][::-1]
    lower_raw = section[le_idx:]
    return upper_raw, lower_raw


def _smooth_parametric_surface(
    surface_points: np.ndarray,
    sample_count: int = 301,
) -> tuple[np.ndarray, np.ndarray]:
    points = np.asarray(surface_points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError("surface_points must have shape (n, 2)")
    if len(points) < 2:
        return points[:, 0], points[:, 1]

    deltas = np.diff(points, axis=0)
    ds = np.hypot(deltas[:, 0], deltas[:, 1])
    arc = np.concatenate(([0.0], np.cumsum(ds)))
    unique_arc, unique_idx = np.unique(arc, return_index=True)
    unique_points = points[unique_idx]

    if unique_arc.size < 2:
        return unique_points[:, 0], unique_points[:, 1]

    s_dense = np.linspace(float(unique_arc[0]), float(unique_arc[-1]), sample_count)
    x_interp = PchipInterpolator(unique_arc, unique_points[:, 0])
    z_interp = PchipInterpolator(unique_arc, unique_points[:, 1])
    return np.asarray(x_interp(s_dense), dtype=float), np.asarray(z_interp(s_dense), dtype=float)


def _smooth_raw_surface_global_curves(
    section: np.ndarray,
    sample_count: int = 301,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    upper_raw, lower_raw = _split_raw_section_surfaces(section)
    upper_x, upper_z = _smooth_parametric_surface(upper_raw, sample_count=sample_count)
    lower_x, lower_z = _smooth_parametric_surface(lower_raw, sample_count=sample_count)
    return upper_x, upper_z, lower_x, lower_z


def _fit_surface_global_curves(
    fit,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    theta = np.deg2rad(float(fit.twist_deg))
    ex = np.array([np.cos(theta), np.sin(theta)], dtype=float)
    ez = np.array([-np.sin(theta), np.cos(theta)], dtype=float)
    le_point = np.array([float(fit.le_x), float(fit.le_z)], dtype=float)
    x_fit = np.asarray(fit.x_fit, dtype=float)
    upper_local = (
        np.outer(x_fit * float(fit.upper_scale), ex)
        + np.outer(np.asarray(fit.fit_y_upper, dtype=float) * float(fit.upper_scale), ez)
    )
    lower_local = (
        np.outer(x_fit * float(fit.lower_scale), ex)
        + np.outer(np.asarray(fit.fit_y_lower, dtype=float) * float(fit.lower_scale), ez)
    )
    upper_global = le_point + upper_local
    lower_global = le_point + lower_local
    return (
        upper_global[:, 0],
        upper_global[:, 1],
        lower_global[:, 0],
        lower_global[:, 1],
    )


def _centered_z_limits(z_values: np.ndarray, padding: float = 0.18) -> tuple[float, float]:
    z_min = float(np.min(z_values))
    z_max = float(np.max(z_values))
    z_center = 0.5 * (z_min + z_max)
    half_span = max(0.5 * (z_max - z_min), 1e-6)
    half_span *= 1.0 + float(padding)
    return z_center - half_span, z_center + half_span


def _section_global_surface_curves(
    section: np.ndarray,
    sample_count: int = 181,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    le_idx = int(np.argmin(section[:, 0]))
    upper_raw = section[: le_idx + 1]
    lower_raw = section[le_idx:]

    te_upper = np.asarray(section[0], dtype=float)
    te_lower = np.asarray(section[-1], dtype=float)
    te_mid = 0.5 * (te_upper + te_lower)
    le_point = np.asarray(section[le_idx], dtype=float)
    chord_vector = te_mid - le_point
    chord = float(np.linalg.norm(chord_vector))
    if chord <= 0.0:
        raise ValueError(f"Invalid chord length {chord}")

    ex = chord_vector / chord
    ez = np.array([-ex[1], ex[0]], dtype=float)

    upper_local = np.column_stack(((upper_raw - le_point) @ ex, (upper_raw - le_point) @ ez))[::-1]
    lower_local = np.column_stack(((lower_raw - le_point) @ ex, (lower_raw - le_point) @ ez))

    upper_x = upper_local[:, 0] / chord
    upper_z = upper_local[:, 1] / chord
    lower_x = lower_local[:, 0] / chord
    lower_z = lower_local[:, 1] / chord

    x_fit = cosine_spacing(sample_count)
    z_upper_local = np.interp(x_fit, upper_x, upper_z)
    z_lower_local = np.interp(x_fit, lower_x, lower_z)

    upper_global = le_point + chord * (
        np.outer(x_fit, ex) + np.outer(z_upper_local, ez)
    )
    lower_global = le_point + chord * (
        np.outer(x_fit, ex) + np.outer(z_lower_local, ez)
    )
    return x_fit, upper_global[:, 1], lower_global[:, 1]


def _save_grid_plot(sections: "OrderedDict[float, np.ndarray]", output_path: Path) -> None:
    n_sections = len(sections)
    ncols = 2
    nrows = int(np.ceil(n_sections / ncols))
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12.0, 3.6 * nrows), squeeze=False, constrained_layout=True)
    colors = plt.cm.viridis(np.linspace(0.12, 0.90, n_sections))

    for ax, (yy, section), color in zip(axes.flat, sections.items(), colors):
        metrics = _section_metrics(section)
        contour_x, contour_z = _section_local_closed_curve(section)

        ax.fill(contour_x, contour_z, color=color, alpha=0.10)
        ax.plot(
            contour_x,
            contour_z,
            color="#2563eb",
            linewidth=1.2,
            alpha=0.90,
            solid_joinstyle="round",
            solid_capstyle="round",
            marker="o",
            markersize=2.8,
            markerfacecolor="#2563eb",
            markeredgewidth=0.0,
        )
        ax.set_title(
            f"y = {yy:.3f} m | c = {metrics['chord']:.3f} m"
        )
        x_all = contour_x
        z_all = contour_z
        x_pad = 0.05 * max(float(np.ptp(x_all)), 1.0)
        z_pad = 0.08 * max(float(np.ptp(z_all)), 1.0)
        ax.set_xlim(float(np.min(x_all)) - x_pad, float(np.max(x_all)) + x_pad)
        ax.set_ylim(float(np.min(z_all)) - z_pad, float(np.max(z_all)) + z_pad)
        ax.set_aspect("equal", adjustable="box")
        ax.grid(True, linewidth=0.35, alpha=0.25)
        ax.set_xlabel("local x [m]")
        ax.set_ylabel("local z [m]")

    for ax in axes.flat[n_sections:]:
        ax.axis("off")

    fig.suptitle(
        "BWB glider airfoils from bwb_glider.geo (raw file points, local section frame)",
        fontsize=16,
        fontweight="semibold",
    )
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _save_zy_evolution_plot(sections: "OrderedDict[float, np.ndarray]", output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(14.0, 5.4), constrained_layout=True)
    section_y = np.asarray(list(sections.keys()), dtype=float)

    x_fit_list = []
    upper_list = []
    lower_list = []
    z_max = []
    z_min = []
    for section in sections.values():
        x_fit, z_upper, z_lower = _section_global_surface_curves(section)
        x_fit_list.append(x_fit)
        upper_list.append(z_upper)
        lower_list.append(z_lower)
        z_max.append(float(np.max(section[:, 1])))
        z_min.append(float(np.min(section[:, 1])))

    x_fit_ref = x_fit_list[0]
    upper_array = np.asarray(upper_list, dtype=float)
    lower_array = np.asarray(lower_list, dtype=float)

    ax.fill_between(
        section_y,
        np.asarray(z_min, dtype=float),
        np.asarray(z_max, dtype=float),
        color="#93c5fd",
        alpha=0.18,
        label="Section thickness envelope",
    )

    ax.plot(section_y, np.asarray(z_max, dtype=float), color="#1d4ed8", linewidth=2.2, alpha=0.9)
    ax.plot(section_y, np.asarray(z_min, dtype=float), color="#1d4ed8", linewidth=2.2, alpha=0.9)

    slice_indices = np.unique(
        np.round(np.linspace(0, len(x_fit_ref) - 1, 9)).astype(int)
    )
    slice_colors = plt.cm.viridis(np.linspace(0.10, 0.90, len(slice_indices)))
    for color, idx in zip(slice_colors, slice_indices):
        x_over_c = float(x_fit_ref[idx])
        label = f"x/c = {x_over_c:.2f}"
        ax.plot(
            section_y,
            upper_array[:, idx],
            color=color,
            linewidth=1.7,
            alpha=0.95,
            label=label,
        )
        ax.plot(
            section_y,
            lower_array[:, idx],
            color=color,
            linewidth=1.7,
            alpha=0.95,
        )

    for yy, section in sections.items():
        ax.axvline(yy, color="#64748b", linestyle="--", linewidth=0.9, alpha=0.35)
        ax.plot(
            np.full(section.shape[0], yy),
            section[:, 1],
            color="#0f172a",
            linewidth=0.55,
            alpha=0.08,
        )
        ax.text(
            yy,
            max(section[:, 1]) + 0.12,
            f"{yy:.3f} m",
            ha="center",
            va="bottom",
            fontsize=9,
            color="#334155",
            rotation=0,
        )

    ax.set_title("BWB glider section evolution in the z-y plane (real z)")
    ax.set_xlabel("spanwise y [m]")
    ax.set_ylabel("vertical z [m]")
    ax.grid(True, linewidth=0.35, alpha=0.25)
    ax.legend(loc="upper right", ncol=2, framealpha=0.95, title="Chordwise slice")
    z_all = np.concatenate([np.asarray(z_min, dtype=float), np.asarray(z_max, dtype=float)])
    z_pad = 0.12 * max(float(np.ptp(z_all)), 1.0)
    ax.set_xlim(section_y.min() - 0.8, section_y.max() + 0.8)
    ax.set_ylim(float(z_all.min()) - z_pad, float(z_all.max()) + z_pad)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    geo_path = CTA_DIR / "airfoils" / "bwb_glider.geo"
    if not geo_path.exists():
        raise FileNotFoundError(f"Could not find {geo_path}")

    sections = _load_sections(geo_path)
    if not sections:
        raise ValueError(f"No airfoil sections found in {geo_path}")

    airfoil_dir = geo_path.parent
    grid_path = airfoil_dir / "bwb_glider_sections_grid.png"
    zy_path = airfoil_dir / "bwb_glider_sections_zy.png"

    _save_grid_plot(sections, grid_path)
    _save_zy_evolution_plot(sections, zy_path)

    print(f"BWB glider section grid written to: {grid_path}")
    print(f"BWB glider z-y evolution written to: {zy_path}")


if __name__ == "__main__":
    main()
