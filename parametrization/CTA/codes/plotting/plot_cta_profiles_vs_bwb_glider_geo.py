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
from parametrization.bwb.builder import prepare_geometry
from parametrization.shared.airfoil_fit import load_xyz_sections_by_span, normalize_airfoil_section


OUTPUT_PNG = CTA_DIR / "outputs" / "wing" / "cta_profiles_vs_bwb_glider_geo.png"
GEO_PATH = CTA_DIR / "airfoils" / "bwb_glider.geo"


def _cta_profile_at_y(prepared, y_station_m: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x_air = np.asarray(prepared.section_model.x_air, dtype=float)
    upper, lower, _ = prepared.section_model.coordinates_at_y(float(y_station_m))
    return x_air, np.asarray(upper, dtype=float), np.asarray(lower, dtype=float)


def _cta_chord_twist(prepared, y_station_m: float) -> tuple[float, float]:
    chord = float(np.interp(y_station_m, prepared.loft.span_stations, prepared.loft.chord))
    twist = float(prepared.spanwise_laws.twist_deg(float(y_station_m)))
    return chord, twist


def _interp(values_x: np.ndarray, values_y: np.ndarray, x_query: np.ndarray) -> np.ndarray:
    return np.interp(x_query, np.asarray(values_x, dtype=float), np.asarray(values_y, dtype=float))


def main() -> None:
    OUTPUT_PNG.parent.mkdir(parents=True, exist_ok=True)

    design = build_cta_design()
    config = to_cta_model_config(design, use_cta_anchor_twist=True)
    prepared = prepare_geometry(config)
    glider_sections = load_xyz_sections_by_span(GEO_PATH)
    y_stations = sorted(glider_sections.keys())

    n_sections = len(y_stations)
    ncols = 2
    nrows = int(np.ceil(n_sections / ncols))
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(12.8, 3.8 * nrows),
        squeeze=False,
        constrained_layout=True,
    )

    cta_color = "#0f4c5c"
    glider_color = "#c2410c"

    for ax, y_station in zip(axes.flat, y_stations):
        glider_norm = normalize_airfoil_section(glider_sections[y_station], sample_count=301)
        x_geo = np.asarray(glider_norm.x_fit, dtype=float)
        yu_geo = np.asarray(glider_norm.y_upper, dtype=float)
        yl_geo = np.asarray(glider_norm.y_lower, dtype=float)

        x_cta, yu_cta, yl_cta = _cta_profile_at_y(prepared, float(y_station))
        chord_cta, twist_cta = _cta_chord_twist(prepared, float(y_station))

        yu_cta_on_geo = _interp(x_cta, yu_cta, x_geo)
        yl_cta_on_geo = _interp(x_cta, yl_cta, x_geo)
        rmse_upper = float(np.sqrt(np.mean((yu_cta_on_geo - yu_geo) ** 2)))
        rmse_lower = float(np.sqrt(np.mean((yl_cta_on_geo - yl_geo) ** 2)))

        ax.fill_between(x_geo, yl_geo, yu_geo, color=glider_color, alpha=0.10, zorder=1)
        ax.plot(x_geo, yu_geo, color=glider_color, linewidth=2.0, linestyle="--", zorder=3, label="bwb_glider.geo")
        ax.plot(x_geo, yl_geo, color=glider_color, linewidth=2.0, linestyle="--", zorder=3)

        ax.fill_between(x_cta, yl_cta, yu_cta, color=cta_color, alpha=0.10, zorder=1)
        ax.plot(x_cta, yu_cta, color=cta_color, linewidth=2.0, zorder=4, label="CTA current")
        ax.plot(x_cta, yl_cta, color=cta_color, linewidth=2.0, zorder=4)

        ax.axhline(0.0, color="#94a3b8", linewidth=0.8, linestyle=":")
        ax.axvline(0.0, color="#94a3b8", linewidth=0.8, linestyle=":")
        ax.set_xlim(-0.02, 1.02)
        ax.set_aspect("equal", adjustable="box")
        ax.grid(True, linewidth=0.35, alpha=0.22)
        ax.set_xlabel("x/c [-]")
        ax.set_ylabel("z/c [-]")
        ax.set_title(
            (
                f"y = {float(y_station):.3f} m | "
                f"c_geo={float(glider_norm.chord):.3f} m, c_cta={chord_cta:.3f} m\n"
                f"twist_geo={float(glider_norm.twist_deg):.3f} deg, twist_cta={twist_cta:.3f} deg | "
                f"RMSE u/l = {rmse_upper:.4f}/{rmse_lower:.4f}"
            ),
            fontsize=10.0,
        )

    for ax in axes.flat[n_sections:]:
        ax.axis("off")

    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.995),
        ncol=2,
        framealpha=0.95,
    )
    fig.suptitle("CTA vs bwb_glider.geo profiles by defined station", fontsize=15, y=1.02)
    fig.savefig(OUTPUT_PNG, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {OUTPUT_PNG}")


if __name__ == "__main__":
    main()
