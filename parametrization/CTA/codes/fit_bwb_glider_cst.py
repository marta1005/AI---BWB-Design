from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import sys
from typing import List

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
CTA_DIR = SCRIPT_DIR.parent
REPO_ROOT = CTA_DIR.parent.parent
MPLCONFIG_DIR = REPO_ROOT / ".mplconfig"
MPLCONFIG_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG_DIR))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from parametrization.shared.airfoil_fit import (
    CSTAirfoilFitOptions,
    fit_airfoil_section_cst,
    load_xyz_sections_by_span,
)


DEGREE = 5
N1 = 0.5
N2 = 1.0
SHARED_LEADING_EDGE = False
SMOOTHNESS_WEIGHT = 1.0e-2
# Match the CTA order-5 parametrization, but avoid forcing the first
# few tenths of chord where some lower surfaces have a very sharp local cusp.
# This follows the spirit of cst-modeling's normalize_foil/fit_curve flow:
# normalize in local chord coordinates, keep TE thickness separately,
# and fit the CST shape on the smooth part of the profile.
FIT_XMIN = 5.0e-3
FIT_XMAX = 0.99

FIT_OPTIONS = CSTAirfoilFitOptions(
    degree=DEGREE,
    n1=N1,
    n2=N2,
    shared_leading_edge=SHARED_LEADING_EDGE,
    smoothness_weight=SMOOTHNESS_WEIGHT,
    fit_xmin=FIT_XMIN,
    fit_xmax=FIT_XMAX,
)


def _write_csv(rows: List[dict[str, object]], path: Path) -> None:
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _plot_fit_grid(results: List[dict[str, object]], output_path: Path) -> None:
    n_sections = len(results)
    ncols = 2
    nrows = int(np.ceil(n_sections / ncols))
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(12.0, 3.8 * nrows),
        squeeze=False,
        constrained_layout=True,
    )
    fit_color = "#ef4444"
    original_color = "#2563eb"

    for ax, result in zip(axes.flat, results):
        x_fit = np.asarray(result["x_fit"], dtype=float)
        yu = np.asarray(result["y_upper"], dtype=float)
        yl = np.asarray(result["y_lower"], dtype=float)
        fit_yu = np.asarray(result["fit_y_upper"], dtype=float)
        fit_yl = np.asarray(result["fit_y_lower"], dtype=float)
        yy = float(result["y_section"])

        ax.plot(
            x_fit,
            yu,
            color=original_color,
            linewidth=2.2,
            alpha=0.55,
            label="Original upper",
        )
        ax.plot(
            x_fit,
            yl,
            color=original_color,
            linewidth=2.2,
            linestyle="--",
            alpha=0.55,
            label="Original lower",
        )
        ax.plot(
            x_fit,
            fit_yu,
            color=fit_color,
            linewidth=2.2,
            linestyle="--",
            label="CST fit upper",
        )
        ax.plot(
            x_fit,
            fit_yl,
            color=fit_color,
            linewidth=2.2,
            linestyle="--",
            label="CST fit lower",
        )
        ax.set_title(
            f"y = {yy:.3f} m | c = {float(result['chord']):.3f} m | twist = {float(result['twist_deg']):.3f} deg\n"
            f"max err = {float(result['max_abs_error']):.5f}"
        )
        ax.set_xlim(-0.02, 1.02)
        ax.set_aspect("equal", adjustable="box")
        ax.grid(True, linewidth=0.35, alpha=0.25)
        ax.set_xlabel("x / c [-]")
        ax.set_ylabel("z / c [-]")

    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, framealpha=0.95)
    fig.suptitle(
        "BWB glider CST fits by section "
        f"(degree {DEGREE}, smoothness={SMOOTHNESS_WEIGHT:g}, fit x/c >= {FIT_XMIN:g})",
        fontsize=16,
        fontweight="semibold",
    )
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    geo_path = CTA_DIR / "airfoils" / "bwb_glider.geo"
    if not geo_path.exists():
        raise FileNotFoundError(f"Could not find {geo_path}")

    sections = load_xyz_sections_by_span(geo_path)
    results: List[dict[str, object]] = []
    csv_rows: List[dict[str, object]] = []
    json_rows: List[dict[str, object]] = []

    for y_section in sorted(sections):
        fit = fit_airfoil_section_cst(sections[y_section], FIT_OPTIONS)
        fit_result = {
            "x_fit": fit.x_fit,
            "y_upper": fit.y_upper,
            "y_lower": fit.y_lower,
            "fit_y_upper": fit.fit_y_upper,
            "fit_y_lower": fit.fit_y_lower,
            "upper_cst": fit.upper_cst,
            "lower_cst": fit.lower_cst,
            "chord": fit.chord,
            "twist_deg": fit.twist_deg,
            "te_thickness": fit.te_thickness,
            "te_mid_x": fit.te_mid_x,
            "te_mid_z": fit.te_mid_z,
            "le_x": fit.le_x,
            "le_z": fit.le_z,
            "fit_xmin": fit.fit_xmin,
            "fit_xmax": fit.fit_xmax,
            "rmse_upper": fit.rmse_upper,
            "rmse_lower": fit.rmse_lower,
            "max_abs_error": fit.max_abs_error,
            "y_section": float(y_section),
        }
        results.append(fit_result)

        csv_row: dict[str, object] = {
            "y_section": float(y_section),
            "chord": float(fit.chord),
            "twist_deg": float(fit.twist_deg),
            "le_x": float(fit.le_x),
            "le_z": float(fit.le_z),
            "te_mid_x": float(fit.te_mid_x),
            "te_mid_z": float(fit.te_mid_z),
            "te_thickness": float(fit.te_thickness),
            "fit_xmin": float(fit.fit_xmin),
            "fit_xmax": float(fit.fit_xmax),
            "rmse_upper": float(fit.rmse_upper),
            "rmse_lower": float(fit.rmse_lower),
            "max_abs_error": float(fit.max_abs_error),
        }
        for idx, value in enumerate(np.asarray(fit.upper_cst, dtype=float)):
            csv_row[f"upper_cst_{idx}"] = float(value)
        for idx, value in enumerate(np.asarray(fit.lower_cst, dtype=float)):
            csv_row[f"lower_cst_{idx}"] = float(value)
        csv_rows.append(csv_row)

        json_rows.append(
            {
                "y_section": float(y_section),
                "chord": float(fit.chord),
                "twist_deg": float(fit.twist_deg),
                "te_thickness": float(fit.te_thickness),
                "shared_leading_edge": SHARED_LEADING_EDGE,
                "degree": DEGREE,
                "n1": N1,
                "n2": N2,
                "smoothness_weight": SMOOTHNESS_WEIGHT,
                "fit_xmin": float(fit.fit_xmin),
                "fit_xmax": float(fit.fit_xmax),
                "upper_cst": np.asarray(fit.upper_cst, dtype=float).tolist(),
                "lower_cst": np.asarray(fit.lower_cst, dtype=float).tolist(),
                "rmse_upper": float(fit.rmse_upper),
                "rmse_lower": float(fit.rmse_lower),
                "max_abs_error": float(fit.max_abs_error),
            }
        )

    airfoil_dir = geo_path.parent
    csv_path = airfoil_dir / "bwb_glider_cst_coefficients.csv"
    json_path = airfoil_dir / "bwb_glider_cst_coefficients.json"
    plot_path = airfoil_dir / "bwb_glider_cst_fit_grid.png"

    _write_csv(csv_rows, csv_path)
    json_path.write_text(json.dumps(json_rows, indent=2), encoding="utf-8")
    _plot_fit_grid(results, plot_path)

    print(f"BWB glider CST coefficients CSV written to: {csv_path}")
    print(f"BWB glider CST coefficients JSON written to: {json_path}")
    print(f"BWB glider CST fit grid written to: {plot_path}")


if __name__ == "__main__":
    main()
