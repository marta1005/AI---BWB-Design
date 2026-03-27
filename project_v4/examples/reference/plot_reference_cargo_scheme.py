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


def build_cargo_polygon(planform, y_sections: np.ndarray):
    # Cargo occupies the full centre-body contour up to the start of
    # transition-wing trailing edge: x = TE(y=B1).
    y_cargo = float(y_sections[1])
    x_cut = float(planform.te_x(y_cargo))

    y_curve = np.linspace(0.0, y_cargo, 120)
    x_curve = np.array([float(planform.le_x(float(value))) for value in y_curve], dtype=float)

    upper = np.column_stack([x_curve, y_curve])
    top_edge = np.asarray([[x_cut, y_cargo]], dtype=float)
    rear_edge = np.asarray([[x_cut, -y_cargo]], dtype=float)
    lower_start = np.asarray([[float(planform.le_x(y_cargo)), -y_cargo]], dtype=float)
    lower = np.column_stack([x_curve[::-1], -y_curve[::-1]])
    cargo_polygon = np.vstack([upper, top_edge, rear_edge, lower_start, lower, upper[:1]])

    return cargo_polygon, x_cut, y_cargo


def polygon_area_xy(polygon: np.ndarray) -> float:
    pts = np.asarray(polygon, dtype=float)
    if pts.shape[0] < 3:
        return 0.0
    x = pts[:, 0]
    y = pts[:, 1]
    return float(0.5 * np.abs(np.dot(x[:-1], y[1:]) + x[-1] * y[0] - np.dot(y[:-1], x[1:]) - y[-1] * x[0]))


def main() -> None:
    design = build_reference_design()
    config = design.to_model_config()
    prepared = prepare_geometry(config)
    planform = prepared.planform

    output_dir = SCRIPT_DIR.parent.parent / "example_outputs" / "reference_v4_demo"
    output_dir.mkdir(parents=True, exist_ok=True)
    png_path = output_dir / "reference_v4_cargo_scheme.png"
    svg_path = output_dir / "reference_v4_cargo_scheme.svg"

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

    y_sections = config.topology.y_sections_array
    cargo, x_cut, y_cargo = build_cargo_polygon(planform, y_sections)
    cargo_half_width = float(np.max(np.abs(cargo[:, 1])))
    cargo_length = float(np.max(cargo[:, 0]) - np.min(cargo[:, 0]))
    cargo_area = polygon_area_xy(cargo)

    y_engine_half = np.linspace(0.0, y_cargo, 600)
    x_te_half = np.array([float(planform.te_x(float(value))) for value in y_engine_half], dtype=float)
    engine_width_half = np.maximum(0.0, x_te_half - x_cut)
    engine_area_half = np.sum(
        0.5 * (engine_width_half[1:] + engine_width_half[:-1]) * (y_engine_half[1:] - y_engine_half[:-1])
    )
    engine_area = float(2.0 * engine_area_half)

    fig, ax = plt.subplots(figsize=(12.0, 7.0), constrained_layout=True)

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
            alpha=0.55,
            zorder=1,
            label=zone_name,
        )
        ax.fill_betweenx(
            -dense_span[mask],
            leading_edge_x[mask],
            trailing_edge_x[mask],
            color=zone_color,
            alpha=0.22,
            zorder=1,
        )

    ax.plot(leading_edge_x, dense_span, color="#0f4c5c", linewidth=2.1, label="LE", zorder=3)
    ax.plot(trailing_edge_x, dense_span, color="#c44536", linewidth=2.1, label="TE", zorder=3)
    ax.plot(leading_edge_x, -dense_span, color="#0f4c5c", linewidth=2.1, zorder=3)
    ax.plot(trailing_edge_x, -dense_span, color="#c44536", linewidth=2.1, zorder=3)

    ax.fill(
        cargo[:, 0],
        cargo[:, 1],
        facecolor="#bfdbfe",
        edgecolor="#2563eb",
        hatch="///",
        linewidth=1.4,
        alpha=0.35,
        zorder=4,
        label="Cargo area",
    )
    ax.plot(cargo[:, 0], cargo[:, 1], color="#2563eb", linewidth=3.2, zorder=5)

    y_engine = np.linspace(-y_cargo, y_cargo, 600)
    x_engine_right = np.array([float(planform.te_x(abs(float(value)))) for value in y_engine], dtype=float)
    valid = x_engine_right > x_cut
    ax.fill_betweenx(
        y_engine[valid],
        x_cut,
        x_engine_right[valid],
        color="#fde68a",
        alpha=0.55,
        zorder=5,
        hatch="\\\\",
        edgecolor="#d97706",
        linewidth=0.0,
        label="Engines area (rear remainder)",
    )
    ax.plot([x_cut, x_cut], [-y_cargo, y_cargo], color="#d97706", linewidth=2.0, zorder=6)
    ax.text(
        float(np.mean(cargo[:, 0])) + 0.2,
        0.35,
        f"Cargo length ≈ {cargo_length:.2f} m\n"
        f"Cargo width ≈ {2.0 * cargo_half_width:.2f} m\n"
        f"Cargo area ≈ {cargo_area:.2f} m²",
        fontsize=9.2,
        color="#1e3a8a",
        ha="center",
        va="center",
        bbox={"boxstyle": "round,pad=0.18", "facecolor": "white", "edgecolor": "#93c5fd", "alpha": 0.93},
        zorder=7,
    )
    x_engine_label = float(np.mean(x_engine_right[valid])) if np.any(valid) else float(x_cut + 0.2)
    ax.text(
        x_engine_label,
        0.0,
        f"Engines area ≈ {engine_area:.2f} m²",
        fontsize=8.9,
        color="#92400e",
        ha="center",
        va="center",
        bbox={"boxstyle": "round,pad=0.16", "facecolor": "white", "edgecolor": "#fcd34d", "alpha": 0.94},
        zorder=7,
    )

    ax.axhline(0.0, color="#94a3b8", linewidth=0.9, linestyle="--", zorder=0)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x / chord direction")
    ax.set_ylabel("spanwise y")
    ax.set_title("CTA reference planform with cargo envelope")
    ax.grid(True, linewidth=0.4, alpha=0.25)
    ax.legend(loc="upper left", ncol=2)

    span_margin = 0.08 * config.topology.span
    x_min = float(np.min(leading_edge_x))
    x_max = float(np.max(trailing_edge_x))
    x_margin = 0.08 * (x_max - x_min)
    ax.set_xlim(x_min - x_margin, x_max + x_margin)
    ax.set_ylim(-(config.topology.span + span_margin), config.topology.span + span_margin)

    save_figure(fig, png_path, dpi=240, bbox_inches="tight")
    save_figure(fig, svg_path, bbox_inches="tight")
    plt.close(fig)

    print(f"Cargo scheme PNG written to: {png_path}")
    print(f"Cargo scheme SVG written to: {svg_path}")
    print(f"Cargo area (planform): {cargo_area:.3f} m^2")
    print(f"Engines area (rear remainder): {engine_area:.3f} m^2")


if __name__ == "__main__":
    main()
