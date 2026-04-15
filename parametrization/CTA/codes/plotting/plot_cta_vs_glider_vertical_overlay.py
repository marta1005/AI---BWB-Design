from __future__ import annotations

import os
from pathlib import Path
import sys

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

from parametrization.CTA.case import build_cta_design, to_cta_model_config
from parametrization.shared.airfoil_fit import load_xyz_sections_by_span
from parametrization.CTA.codes.plotting.plot_cta_views import (
    front_linear_edge_traces,
    front_linear_envelope,
)
from parametrization.bwb.builder import prepare_geometry


def _glider_envelope(sections: dict[float, np.ndarray]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    y_sections = np.asarray(sorted(sections.keys()), dtype=float)
    upper = np.asarray([float(np.max(sections[yy][:, 1])) for yy in y_sections], dtype=float)
    lower = np.asarray([float(np.min(sections[yy][:, 1])) for yy in y_sections], dtype=float)
    return y_sections, upper, lower


def _cta_envelope():
    design = build_cta_design()
    config = to_cta_model_config(design, use_cta_anchor_twist=True)
    prepared = prepare_geometry(config)
    anchor_y = np.asarray(config.topology.anchor_y_array, dtype=float)
    dense_span = np.unique(
        np.concatenate([np.linspace(0.0, config.topology.span, 2400), config.topology.y_sections_array, anchor_y])
    )
    upper, lower, _, _ = front_linear_envelope(prepared, dense_span, anchor_y)
    z_le, z_te_upper, z_te_lower, z_te_mid = front_linear_edge_traces(prepared, dense_span, anchor_y)
    return dense_span, upper, lower, z_le, z_te_upper, z_te_lower, z_te_mid, anchor_y


def _anchor_mismatch_summary(
    sections: dict[float, np.ndarray],
    anchor_y: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    design = build_cta_design()
    config = to_cta_model_config(design, use_cta_anchor_twist=True)
    prepared = prepare_geometry(config)

    d_upper = []
    d_lower = []
    for yy in anchor_y:
        raw_key = min(sections.keys(), key=lambda value: abs(value - float(yy)))
        raw = sections[raw_key]
        raw_upper = float(np.max(raw[:, 1]))
        raw_lower = float(np.min(raw[:, 1]))

        yu, yl, _ = prepared.section_model.coordinates_at_y(float(yy))
        chord = float(prepared.planform.te_x(float(yy)) - prepared.planform.le_x(float(yy)))
        vertical = float(np.interp(float(yy), prepared.loft.span_stations, prepared.loft.vertical_y))
        twist = float(prepared.spanwise_laws.twist_deg(float(yy)))
        x_local = prepared.section_model.x_air * chord
        twist_rad = np.deg2rad(twist)
        z_up = vertical + x_local * np.sin(twist_rad) + yu * chord * np.cos(twist_rad)
        z_lo = vertical + x_local * np.sin(twist_rad) + yl * chord * np.cos(twist_rad)
        d_upper.append(float(np.max(z_up)) - raw_upper)
        d_lower.append(float(np.min(z_lo)) - raw_lower)

    return np.asarray(d_upper, dtype=float), np.asarray(d_lower, dtype=float)


def main() -> None:
    geo_path = CTA_DIR / "airfoils" / "bwb_glider.geo"
    sections = load_xyz_sections_by_span(geo_path)
    glider_y, glider_upper, glider_lower = _glider_envelope(sections)
    (
        cta_y,
        cta_upper,
        cta_lower,
        cta_z_le,
        cta_z_te_upper,
        cta_z_te_lower,
        cta_z_te_mid,
        anchor_y,
    ) = _cta_envelope()
    d_upper, d_lower = _anchor_mismatch_summary(sections, anchor_y)

    fig, ax = plt.subplots(figsize=(14.0, 5.4), constrained_layout=True)

    ax.fill_between(cta_y, cta_lower, cta_upper, color="#fecaca", alpha=0.18, zorder=1)
    ax.plot(cta_y, cta_upper, color="#dc2626", linewidth=2.2, zorder=3, label="CTA upper envelope")
    ax.plot(cta_y, cta_lower, color="#dc2626", linewidth=2.2, zorder=3, label="CTA lower envelope")

    ax.plot(glider_y, glider_upper, color="#2563eb", linewidth=2.0, marker="o", markersize=4.0, zorder=4, label="Glider upper envelope")
    ax.plot(glider_y, glider_lower, color="#2563eb", linewidth=2.0, linestyle="--", marker="o", markersize=4.0, zorder=4, label="Glider lower envelope")

    ax.plot(cta_y, cta_z_le, color="#7c3aed", linewidth=1.8, linestyle=":", alpha=0.8, zorder=2, label="CTA LE z(y)")
    ax.plot(cta_y, cta_z_te_upper, color="#ea580c", linewidth=1.5, linestyle="--", alpha=0.75, zorder=2, label="CTA TE upper z(y)")
    ax.plot(cta_y, cta_z_te_lower, color="#ea580c", linewidth=1.5, linestyle="--", alpha=0.75, zorder=2, label="CTA TE lower z(y)")
    ax.plot(cta_y, cta_z_te_mid, color="#475569", linewidth=1.4, linestyle="-.", alpha=0.75, zorder=2, label="CTA TE mid z(y)")

    for yy, du, dl in zip(anchor_y, d_upper, d_lower):
        ax.axvline(float(yy), color="#64748b", linewidth=0.8, linestyle=(0, (4, 4)), alpha=0.45, zorder=0)
        ax.text(
            float(yy),
            float(np.max(np.concatenate([cta_upper, glider_upper]))) + 0.08,
            f"{float(yy):.3f} m",
            fontsize=8.3,
            color="#334155",
            ha="center",
            va="bottom",
            bbox={"boxstyle": "round,pad=0.10", "facecolor": "white", "edgecolor": "none", "alpha": 0.92},
            zorder=6,
        )
        ax.plot([yy], [np.interp(yy, cta_y, cta_upper)], marker="x", color="#dc2626", markersize=6.0, zorder=5)
        ax.plot([yy], [np.interp(yy, cta_y, cta_lower)], marker="x", color="#dc2626", markersize=6.0, zorder=5)

    summary = (
        f"Anchor mismatch CTA - glider\n"
        f"Upper max: {float(np.max(np.abs(d_upper))):.3f} m\n"
        f"Lower max: {float(np.max(np.abs(d_lower))):.3f} m"
    )
    ax.text(
        0.015,
        0.98,
        summary,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9.0,
        color="#1f2937",
        bbox={"boxstyle": "round,pad=0.22", "facecolor": "white", "edgecolor": "#cbd5e1", "alpha": 0.95},
    )

    ax.set_title("CTA vertical overlay vs bwb_glider raw sections")
    ax.set_xlabel("spanwise y [m]")
    ax.set_ylabel("vertical z [m]")
    ax.grid(True, linewidth=0.35, alpha=0.30)
    ax.set_xlim(0.0, max(float(np.max(cta_y)), float(np.max(glider_y))) + 0.8)
    z_all = np.concatenate([cta_upper, cta_lower, glider_upper, glider_lower])
    z_pad = 0.10 * max(float(np.ptp(z_all)), 1.0)
    ax.set_ylim(float(np.min(z_all)) - z_pad, float(np.max(z_all)) + z_pad)
    ax.legend(loc="lower right", ncol=2, fontsize=8, framealpha=0.94)

    output_path = CTA_DIR / "outputs" / "wing" / "cta_vs_bwb_glider_vertical_overlay.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)

    print(f"CTA vs glider vertical overlay PNG written to: {output_path}")


if __name__ == "__main__":
    main()
