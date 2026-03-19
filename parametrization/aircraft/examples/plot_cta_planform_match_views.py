from __future__ import annotations

from pathlib import Path
import json
import os
import shutil
import sys
import tempfile

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR.parent / ".mplconfig"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from cta_planform_match import (
    build_cta_planform_profiles,
    build_cta_planform_reference,
    build_cta_planform_wing,
    default_build_options,
    default_output_dir,
)

from parametrization.aircraft import prepare_lifting_surface
from parametrization.aircraft.plotting import (
    build_lifting_surface_mesh,
    plot_front_view,
    plot_side_view,
    save_lifting_surface_3d,
)


def _save_figure(fig, path: Path, **kwargs) -> None:
    try:
        fig.savefig(path, **kwargs)
    except OSError as exc:
        if "Resource deadlock avoided" not in str(exc):
            raise
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir) / path.name
            fig.savefig(tmp_path, **kwargs)
            shutil.copy2(tmp_path, path)


def _nearest_station(prepared, eta: float):
    return min(prepared.stations, key=lambda station: abs(float(station.eta) - float(eta)))


def _plot_selected_profiles(ax, prepared, reference) -> None:
    labels = ("C1", "C2", "C3", "C4")
    selected = [_nearest_station(prepared, eta) for eta in reference.section_etas]
    spans = [float(np.max(station.yu) - np.min(station.yl)) * float(station.chord) for station in selected]
    base_gap = max(max(spans, default=0.0), 1.0) * 0.60
    offsets = []
    current = 0.0
    for span in spans:
        offsets.append(current)
        current += span + base_gap

    colors = ("#0f4c5c", "#0f766e", "#2563eb", "#7c3aed")
    for label, station, offset, color in zip(labels, selected, offsets, colors):
        x_vals = float(station.chord) * np.asarray(station.x_air, dtype=float)
        yu = float(station.chord) * np.asarray(station.yu, dtype=float)
        yl = float(station.chord) * np.asarray(station.yl, dtype=float)
        ax.plot(x_vals, yu + offset, color=color, linewidth=1.7)
        ax.plot(x_vals, yl + offset, color=color, linewidth=1.7)
        ax.fill_between(x_vals, yl + offset, yu + offset, color=color, alpha=0.12)
        ax.text(
            float(x_vals[-1]) + 0.35,
            offset,
            f"{label} | eta={station.eta:.3f} | c={station.chord:.2f} m",
            fontsize=8.5,
            va="center",
            color="#334155",
        )

    ax.set_title("Selected section profiles")
    ax.set_xlabel("local chordwise x [m]")
    ax.set_ylabel("stacked local thickness y [m]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, linewidth=0.35, alpha=0.25)
    ax.set_yticks([])


def _planform_summary(prepared, reference) -> dict[str, float | int]:
    model_le = np.interp(reference.y_dense, prepared.z, prepared.x)
    model_te = np.interp(reference.y_dense, prepared.z, prepared.x + prepared.chord)
    le_error = model_le - reference.le_x_dense
    te_error = model_te - reference.te_x_dense
    chord_error = (model_te - model_le) - (reference.te_x_dense - reference.le_x_dense)
    return {
        "span_m": float(reference.span),
        "sample_station_count": int(len(prepared.stations)),
        "max_le_error_m": float(np.max(np.abs(le_error))),
        "max_te_error_m": float(np.max(np.abs(te_error))),
        "max_chord_error_m": float(np.max(np.abs(chord_error))),
        "reference_polyline_le_error_m": float(reference.max_le_polyline_error_m),
        "reference_polyline_te_error_m": float(reference.max_te_polyline_error_m),
    }


def _create_overview_figure(prepared, reference) -> plt.Figure:
    mesh = build_lifting_surface_mesh(prepared)
    fig, axes = plt.subplots(2, 2, figsize=(16.0, 11.0), constrained_layout=True)

    ax_top = axes[0, 0]
    ax_top.plot(reference.le_x_dense, reference.y_dense, color="#0f172a", linewidth=2.0, linestyle="--", label="CTA LE")
    ax_top.plot(reference.te_x_dense, reference.y_dense, color="#475569", linewidth=2.0, linestyle="--", label="CTA TE")
    ax_top.plot(prepared.x, prepared.z, color="#0f4c5c", linewidth=1.9, label="aircraft LE")
    ax_top.plot(prepared.x + prepared.chord, prepared.z, color="#c44536", linewidth=1.9, label="aircraft TE")
    ax_top.scatter(reference.section_le_x, reference.y_sections, color="#1d4ed8", s=26, zorder=6, label="CTA sections")
    ax_top.set_title("Top view | CTA planform vs aircraft wing")
    ax_top.set_xlabel("x [m]")
    ax_top.set_ylabel("z spanwise [m]")
    ax_top.grid(True, linewidth=0.35, alpha=0.25)
    ax_top.set_aspect("equal", adjustable="box")
    ax_top.legend(loc="best", fontsize=8.5)

    plot_side_view(axes[0, 1], mesh, title="Side view")
    plot_front_view(axes[1, 0], mesh, title="Front view")
    _plot_selected_profiles(axes[1, 1], prepared, reference)
    return fig


def _create_comparison_figure(prepared, reference) -> tuple[plt.Figure, dict[str, float | int]]:
    summary = _planform_summary(prepared, reference)
    model_le = np.interp(reference.y_dense, prepared.z, prepared.x)
    model_te = np.interp(reference.y_dense, prepared.z, prepared.x + prepared.chord)
    le_error = model_le - reference.le_x_dense
    te_error = model_te - reference.te_x_dense

    fig, axes = plt.subplots(2, 1, figsize=(14.0, 10.0), constrained_layout=True)
    ax0, ax1 = axes

    ax0.plot(reference.le_x_dense, reference.y_dense, color="#0f172a", linewidth=2.1, linestyle="--", label="CTA LE")
    ax0.plot(reference.te_x_dense, reference.y_dense, color="#475569", linewidth=2.1, linestyle="--", label="CTA TE")
    ax0.plot(model_le, reference.y_dense, color="#0f4c5c", linewidth=1.8, label="aircraft LE")
    ax0.plot(model_te, reference.y_dense, color="#c44536", linewidth=1.8, label="aircraft TE")
    ax0.set_title("CTA reference planform match")
    ax0.set_xlabel("x [m]")
    ax0.set_ylabel("z spanwise [m]")
    ax0.grid(True, linewidth=0.35, alpha=0.25)
    ax0.set_aspect("equal", adjustable="box")
    ax0.legend(loc="best", fontsize=8.8)

    ax1.plot(reference.y_dense, le_error, color="#0f4c5c", linewidth=1.6, label=f"LE error | max={summary['max_le_error_m']:.4f} m")
    ax1.plot(reference.y_dense, te_error, color="#c44536", linewidth=1.6, label=f"TE error | max={summary['max_te_error_m']:.4f} m")
    ax1.axhline(0.0, color="#0f172a", linewidth=1.0, alpha=0.7)
    ax1.set_title("Planform deviation relative to CTA reference")
    ax1.set_xlabel("z spanwise [m]")
    ax1.set_ylabel("x error [m]")
    ax1.grid(True, linewidth=0.35, alpha=0.25)
    ax1.legend(loc="best", fontsize=8.8)

    return fig, summary


def main() -> None:
    reference = build_cta_planform_reference()
    profiles = build_cta_planform_profiles()
    wing = build_cta_planform_wing(reference=reference)
    component = wing.to_component_spec()
    prepared = prepare_lifting_surface(component, profiles, options=default_build_options(reference=reference))
    output_dir = default_output_dir()
    output_dir.mkdir(parents=True, exist_ok=True)

    overview_fig = _create_overview_figure(prepared, reference)
    overview_png = output_dir / "cta_planform_match_views.png"
    overview_svg = output_dir / "cta_planform_match_views.svg"
    _save_figure(overview_fig, overview_png, dpi=220, bbox_inches="tight")
    _save_figure(overview_fig, overview_svg, bbox_inches="tight")
    plt.close(overview_fig)

    comparison_fig, summary = _create_comparison_figure(prepared, reference)
    comparison_png = output_dir / "cta_planform_match_comparison.png"
    comparison_svg = output_dir / "cta_planform_match_comparison.svg"
    _save_figure(comparison_fig, comparison_png, dpi=220, bbox_inches="tight")
    _save_figure(comparison_fig, comparison_svg, bbox_inches="tight")
    plt.close(comparison_fig)

    view3d_png = save_lifting_surface_3d(
        prepared,
        output_dir,
        stem="cta_planform_match",
        title="CTA planform match | 3D view",
    )

    summary_path = output_dir / "cta_planform_match_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(f"Overview PNG written to: {overview_png}")
    print(f"Overview SVG written to: {overview_svg}")
    print(f"Comparison PNG written to: {comparison_png}")
    print(f"Comparison SVG written to: {comparison_svg}")
    print(f"3D PNG written to: {view3d_png}")
    print(f"Summary JSON written to: {summary_path}")


if __name__ == "__main__":
    main()
