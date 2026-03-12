from pathlib import Path
import os
import sys
import tempfile

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR.parent.parent / ".mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from project_v4.airfoil_io import write_airfoil_dat
from project_v4.builder import prepare_geometry
from project_v4.dependency_setup import ensure_local_dependency_paths
from project_v4.examples.reference.run_reference_example import build_reference_design


def te_gap_from_contour(x: np.ndarray, y: np.ndarray) -> float:
    return float(abs(y[0] - y[-1]))


def main() -> None:
    ensure_local_dependency_paths()
    from pygeo.geo_utils.file_io import readAirfoilFile

    output_dir = SCRIPT_DIR.parent.parent / "example_outputs" / "reference_v4_demo"
    output_dir.mkdir(parents=True, exist_ok=True)
    png_path = output_dir / "reference_v4_blunt_te_comparison.png"
    svg_path = output_dir / "reference_v4_blunt_te_comparison.svg"

    design = build_reference_design()
    config = design.to_model_config()
    prepared = prepare_geometry(config)

    section_y = float(config.topology.y_sections_array[0])
    yu, yl, params = prepared.section_model.coordinates_at_y(section_y)
    x = prepared.section_model.x_air

    with tempfile.TemporaryDirectory() as tmp_dir:
        dat_path = Path(tmp_dir) / "root_section.dat"
        write_airfoil_dat(str(dat_path), x, yu, yl, name="ROOT_SECTION")

        x_sharp, y_sharp = readAirfoilFile(str(dat_path), bluntTe=False)
        x_blunt, y_blunt = readAirfoilFile(
            str(dat_path),
            bluntTe=True,
            bluntThickness=float(params.te_thickness),
        )

    xu_in = x[::-1]
    yu_in = yu[::-1]
    xl_in = x[1:]
    yl_in = yl[1:]
    x_in = np.concatenate([xu_in, xl_in])
    y_in = np.concatenate([yu_in, yl_in])

    fig, (ax_full, ax_zoom) = plt.subplots(1, 2, figsize=(13.0, 5.6), constrained_layout=True)
    curves = (
        ("Input profile", x_in, y_in, "#0f4c5c", 2.2),
        ("pyGeo with bluntTe=False", x_sharp, y_sharp, "#c44536", 1.8),
        ("pyGeo with bluntTe=True", x_blunt, y_blunt, "#1d4ed8", 1.8),
    )

    for label, xx, yy, color, width in curves:
        ax_full.plot(xx, yy, color=color, linewidth=width, label=label)
        ax_zoom.plot(xx, yy, color=color, linewidth=width, label=label)

    ax_full.set_title("Full section contour")
    ax_full.set_xlabel("x/c")
    ax_full.set_ylabel("z/c")
    ax_full.grid(True, linewidth=0.4, alpha=0.3)
    ax_full.set_aspect("equal", adjustable="box")
    ax_full.legend(loc="upper right")

    ax_zoom.set_title("Trailing-edge zoom")
    ax_zoom.set_xlabel("x/c")
    ax_zoom.grid(True, linewidth=0.4, alpha=0.3)
    ax_zoom.set_aspect("equal", adjustable="box")
    ax_zoom.set_xlim(0.88, 1.005)
    y_zoom = np.concatenate([y_in[x_in >= 0.88], y_sharp[x_sharp >= 0.88], y_blunt[x_blunt >= 0.88]])
    y_pad = 0.15 * max(1e-6, float(np.max(y_zoom) - np.min(y_zoom)))
    ax_zoom.set_ylim(float(np.min(y_zoom) - y_pad), float(np.max(y_zoom) + y_pad))

    text = (
        f"Section: C1 / root at y={section_y:.2f} m\n"
        f"Target TE thickness = {params.te_thickness:.4f} c\n"
        f"Input contour TE gap = {te_gap_from_contour(x_in, y_in):.4f} c\n"
        f"pyGeo sharp mode TE gap = {te_gap_from_contour(x_sharp, y_sharp):.4f} c\n"
        f"pyGeo blunt mode TE gap = {te_gap_from_contour(x_blunt, y_blunt):.4f} c"
    )
    ax_zoom.text(
        0.885,
        float(np.max(y_zoom)),
        text,
        ha="left",
        va="top",
        fontsize=8.5,
        color="#243b53",
        bbox={"boxstyle": "round,pad=0.20", "facecolor": "white", "edgecolor": "#cbd5e1", "alpha": 0.92},
    )

    fig.suptitle("pyGeo trailing-edge handling: sharp vs blunt", fontsize=14)
    fig.savefig(png_path, dpi=220, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)

    print(f"Blunt-TE comparison PNG written to: {png_path}")
    print(f"Blunt-TE comparison SVG written to: {svg_path}")


if __name__ == "__main__":
    main()
