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
from parametrization.shared.cst import cosine_spacing
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
    root_dense = (
        float(config.planform.symmetry_blend_y) * cosine_spacing(240)
        if float(config.planform.symmetry_blend_y) > 1.0e-12
        else np.array([], dtype=float)
    )
    dense_span = np.unique(
        np.concatenate(
            [
                np.linspace(0.0, config.topology.span, 2400),
                root_dense,
                config.topology.y_sections_array,
                anchor_y,
            ]
        )
    )
    upper, lower, _, _ = front_linear_envelope(prepared, dense_span, anchor_y)
    z_le, z_te_upper, z_te_lower, z_te_mid = front_linear_edge_traces(prepared, dense_span, anchor_y)
    return dense_span, upper, lower, z_le, z_te_upper, z_te_lower, z_te_mid, anchor_y


def _save_clean_compare_figure(
    output_path: Path,
    cta_y: np.ndarray,
    cta_upper: np.ndarray,
    cta_lower: np.ndarray,
    glider_y: np.ndarray,
    glider_upper: np.ndarray,
    glider_lower: np.ndarray,
    anchor_y: np.ndarray,
) -> None:
    fig, ax = plt.subplots(figsize=(14.0, 5.6), constrained_layout=True)

    ax.fill_between(cta_y, cta_lower, cta_upper, color="#dce7f1", zorder=1, alpha=0.92)
    ax.fill_between(-cta_y, cta_lower, cta_upper, color="#dce7f1", zorder=1, alpha=0.92)
    ax.plot(cta_y, cta_upper, color="#0f4c5c", linewidth=2.3, zorder=3)
    ax.plot(cta_y, cta_lower, color="#0f4c5c", linewidth=2.3, zorder=3)
    ax.plot(-cta_y, cta_upper, color="#0f4c5c", linewidth=2.3, zorder=3)
    ax.plot(-cta_y, cta_lower, color="#0f4c5c", linewidth=2.3, zorder=3, label="CTA envelope")

    glider_style = {
        "color": "#c2410c",
        "linewidth": 2.0,
        "linestyle": "--",
        "zorder": 4,
        "alpha": 0.95,
    }
    ax.plot(glider_y, glider_upper, **glider_style)
    ax.plot(glider_y, glider_lower, **glider_style)
    ax.plot(-glider_y, glider_upper, **glider_style)
    ax.plot(-glider_y, glider_lower, label="bwb_glider envelope", **glider_style)

    for x_value in np.unique(np.concatenate((-anchor_y[1:], anchor_y[1:]))):
        ax.axvline(float(x_value), color="#64748b", linewidth=0.8, linestyle=(0, (4, 4)), alpha=0.30, zorder=2)

    ax.axhline(0.0, color="#94a3b8", linewidth=0.9, linestyle=":")
    ax.axvline(0.0, color="#94a3b8", linewidth=0.9, linestyle=":")
    ax.set_title("CTA envelope vs bwb_glider envelope")
    ax.set_xlabel("spanwise y [m]")
    ax.set_ylabel("vertical z [m]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, linewidth=0.35, alpha=0.30)

    span = max(float(np.max(cta_y)), float(np.max(glider_y)))
    z_all = np.concatenate([cta_upper, cta_lower, glider_upper, glider_lower])
    z_span = max(1e-9, float(np.ptp(z_all)))
    ax.set_xlim(-span - 1.0, span + 1.0)
    ax.set_ylim(float(np.min(z_all)) - 0.08 * z_span, float(np.max(z_all)) + 0.08 * z_span)
    ax.legend(loc="lower right", fontsize=8.6, framealpha=0.94)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


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

    for yy in anchor_y:
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

    clean_compare_path = CTA_DIR / "outputs" / "wing" / "cta_envelope_vs_bwb_glider.png"
    _save_clean_compare_figure(
        clean_compare_path,
        cta_y,
        cta_upper,
        cta_lower,
        glider_y,
        glider_upper,
        glider_lower,
        anchor_y,
    )

    print(f"CTA vs glider vertical overlay PNG written to: {output_path}")
    print(f"CTA envelope vs glider PNG written to: {clean_compare_path}")


if __name__ == "__main__":
    main()
