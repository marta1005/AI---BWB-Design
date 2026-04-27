from __future__ import annotations

from collections import Counter, defaultdict
import json
import os
from pathlib import Path
import sys
from typing import Iterable

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
CTA_DIR = SCRIPT_DIR.parent.parent
REPO_ROOT = CTA_DIR.parent.parent
os.environ.setdefault("MPLCONFIGDIR", str(REPO_ROOT / ".mplconfig"))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from parametrization.shared.airfoil_fit import (
    CSTAirfoilFitOptions,
    fit_airfoil_section_cst,
    normalize_airfoil_section,
)
from parametrization.shared.airfoil_io import write_airfoil_dat


STL_PATH = Path("/Users/martaarnabatmartin/Downloads/BWB_LH_V0_GLIDER_CORRECTED2_ALLCATPART.stl")
SECTION_Y_MM = 7000.0
FIT_OPTIONS = CSTAirfoilFitOptions(
    degree=5,
    n1=0.35,
    n2=1.1,
    smoothness_weight=1.0e-2,
    fit_xmin=0.0,
    fit_xmax=0.99,
)


def _load_ascii_stl_triangles(path: Path) -> np.ndarray:
    triangles: list[list[tuple[float, float, float]]] = []
    current: list[tuple[float, float, float]] = []
    with path.open("r", encoding="utf-8", errors="ignore") as stream:
        for raw_line in stream:
            line = raw_line.strip()
            if not line.startswith("vertex"):
                continue
            _, x_val, y_val, z_val = line.split()
            current.append((float(x_val), float(y_val), float(z_val)))
            if len(current) == 3:
                triangles.append(current)
                current = []
    if not triangles:
        raise ValueError(f"No triangles found in ASCII STL: {path}")
    return np.asarray(triangles, dtype=float)


def _boundary_edges_on_plane(
    triangles: np.ndarray,
    plane_y_mm: float,
    tol_mm: float = 1.0e-6,
) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    mask = np.all(np.abs(triangles[:, :, 1] - float(plane_y_mm)) <= float(tol_mm), axis=1)
    on_plane = np.asarray(triangles[mask], dtype=float)
    if on_plane.size == 0:
        raise ValueError(f"No coplanar triangles found at y={plane_y_mm:.6f} mm")

    edge_counter: Counter[tuple[tuple[float, float], tuple[float, float]]] = Counter()
    for tri in on_plane:
        pts = [(round(float(x_val), 6), round(float(z_val), 6)) for x_val, _, z_val in tri]
        for point_a, point_b in ((pts[0], pts[1]), (pts[1], pts[2]), (pts[2], pts[0])):
            edge = tuple(sorted((point_a, point_b)))
            edge_counter[edge] += 1

    boundary = [edge for edge, count in edge_counter.items() if count == 1]
    if not boundary:
        raise ValueError(f"No boundary edges extracted from y={plane_y_mm:.6f} mm plane")
    return boundary


def _ordered_loop_from_boundary(
    boundary_edges: Iterable[tuple[tuple[float, float], tuple[float, float]]],
) -> np.ndarray:
    adjacency: defaultdict[tuple[float, float], list[tuple[float, float]]] = defaultdict(list)
    for point_a, point_b in boundary_edges:
        adjacency[point_a].append(point_b)
        adjacency[point_b].append(point_a)

    invalid = {point: neighbors for point, neighbors in adjacency.items() if len(neighbors) != 2}
    if invalid:
        raise ValueError(f"Boundary is not a single closed loop: {invalid}")

    start = max(adjacency, key=lambda point: (point[0], point[1]))
    loop = [start]
    previous = None
    current = start
    while True:
        neighbors = adjacency[current]
        next_point = neighbors[0] if neighbors[0] != previous else neighbors[1]
        if next_point == start:
            break
        loop.append(next_point)
        previous, current = current, next_point

    points = np.asarray(loop, dtype=float)
    leading_edge_idx = int(np.argmin(points[:, 0]))
    if np.mean(points[: leading_edge_idx + 1, 1]) < np.mean(points[leading_edge_idx:, 1]):
        points = points[::-1]

    max_x = float(np.max(points[:, 0]))
    candidate_idx = np.where(np.isclose(points[:, 0], max_x, atol=5.0))[0]
    start_idx = int(candidate_idx[np.argmax(points[candidate_idx, 1])])
    points = np.vstack([points[start_idx:], points[:start_idx], points[start_idx : start_idx + 1]])

    leading_edge_idx = int(np.argmin(points[:, 0]))
    if np.mean(points[: leading_edge_idx + 1, 1]) < np.mean(points[leading_edge_idx:, 1]):
        points = points[::-1]
        points = points[:-1]
        max_x = float(np.max(points[:, 0]))
        candidate_idx = np.where(np.isclose(points[:, 0], max_x, atol=5.0))[0]
        start_idx = int(candidate_idx[np.argmax(points[candidate_idx, 1])])
        points = np.vstack([points[start_idx:], points[:start_idx], points[start_idx : start_idx + 1]])

    return points


def _section_anchor_payload(fit_result) -> dict[str, object]:
    thickness = np.asarray(fit_result.fit_y_upper - fit_result.fit_y_lower, dtype=float)
    thickness_idx = int(np.argmax(thickness))
    return {
        "anchor_y_m": float(SECTION_Y_MM / 1000.0),
        "anchor_le_z_m": float(fit_result.le_z / 1000.0),
        "anchor_twist_deg": float(fit_result.twist_deg),
        "upper_coeffs": np.asarray(fit_result.upper_cst, dtype=float).tolist(),
        # SectionCSTSpec stores the lower surface coefficients as positive
        # magnitudes; KulfanCSTAirfoil applies the minus sign internally.
        "lower_coeffs": np.asarray(-fit_result.lower_cst, dtype=float).tolist(),
        "tc_max": float(thickness[thickness_idx]),
        "x_tmax": float(fit_result.x_fit[thickness_idx]),
        "te_thickness": float(fit_result.te_thickness),
        "chord_m": float(fit_result.chord / 1000.0),
        "le_x_m": float(fit_result.le_x / 1000.0),
        "max_abs_error": float(fit_result.max_abs_error),
    }


def main() -> None:
    if not STL_PATH.exists():
        raise FileNotFoundError(f"Could not find STL at {STL_PATH}")

    triangles = _load_ascii_stl_triangles(STL_PATH)
    boundary = _boundary_edges_on_plane(triangles, SECTION_Y_MM)
    section_points = _ordered_loop_from_boundary(boundary)

    normalized = normalize_airfoil_section(section_points)
    fit_result = fit_airfoil_section_cst(section_points, FIT_OPTIONS)
    anchor_payload = _section_anchor_payload(fit_result)

    airfoil_dir = CTA_DIR / "airfoils"
    output_dir = CTA_DIR / "outputs" / "wing"
    airfoil_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    raw_dat = airfoil_dir / "stl_section_y7_raw_local.dat"
    fit_dat = airfoil_dir / "stl_section_y7_cst_fit.dat"
    summary_json = airfoil_dir / "stl_section_y7_anchor_summary.json"
    plot_png = output_dir / "stl_section_y7_fit_check.png"

    write_airfoil_dat(str(raw_dat), normalized.x_fit, normalized.y_upper, normalized.y_lower, name="STL_Y7_RAW")
    write_airfoil_dat(
        str(fit_dat),
        fit_result.x_fit,
        fit_result.fit_y_upper,
        fit_result.fit_y_lower,
        name="STL_Y7_CST_FIT",
    )
    summary_json.write_text(json.dumps(anchor_payload, indent=2), encoding="utf-8")

    fig, ax = plt.subplots(figsize=(12.0, 4.2), constrained_layout=True)
    ax.plot(section_points[:, 0], section_points[:, 1], color="#2563eb", linewidth=2.0, label="STL section raw")

    x_local = fit_result.x_fit * fit_result.chord
    z_upper = fit_result.fit_y_upper * fit_result.chord
    z_lower = fit_result.fit_y_lower * fit_result.chord
    x_closed = np.concatenate([x_local[::-1], x_local[1:]])
    z_closed = np.concatenate([z_upper[::-1], z_lower[1:]])

    twist_rad = np.deg2rad(float(fit_result.twist_deg))
    ex = np.array([np.cos(twist_rad), np.sin(twist_rad)], dtype=float)
    ez = np.array([-np.sin(twist_rad), np.cos(twist_rad)], dtype=float)
    fit_global = np.column_stack(
        [
            float(fit_result.le_x) + x_closed * ex[0] + z_closed * ez[0],
            float(fit_result.le_z) + x_closed * ex[1] + z_closed * ez[1],
        ]
    )
    ax.plot(fit_global[:, 0], fit_global[:, 1], color="#c9706b", linewidth=2.0, label="CST fit")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.25)
    ax.legend()
    ax.set_title(
        "STL y=7000 mm section | "
        f"chord={fit_result.chord / 1000.0:.3f} m | "
        f"twist={fit_result.twist_deg:.3f} deg | "
        f"LE_z={fit_result.le_z / 1000.0:.3f} m"
    )
    fig.savefig(plot_png, dpi=220, bbox_inches="tight")
    plt.close(fig)

    print(f"STL raw local DAT written to: {raw_dat}")
    print(f"STL CST-fit DAT written to: {fit_dat}")
    print(f"STL anchor summary written to: {summary_json}")
    print(f"STL fit check plot written to: {plot_png}")


if __name__ == "__main__":
    main()
