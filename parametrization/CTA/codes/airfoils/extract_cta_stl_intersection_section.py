from __future__ import annotations

import argparse
from collections import defaultdict
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

from parametrization.shared.airfoil_fit import CSTAirfoilFitOptions, fit_airfoil_section_cst, normalize_airfoil_section
from parametrization.shared.airfoil_io import write_airfoil_dat


DEFAULT_STL_PATH = Path("/Users/martaarnabatmartin/Downloads/CTA.stl")
DEFAULT_SECTION_Y_M = 4.0
FIT_OPTIONS = CSTAirfoilFitOptions(
    degree=5,
    n1=0.35,
    n2=1.1,
    smoothness_weight=1.0e-3,
    fit_xmin=0.0,
    fit_xmax=0.95,
)


def _safe_y_tag(section_y_m: float) -> str:
    text = f"{section_y_m:.3f}".rstrip("0").rstrip(".")
    return text.replace("-", "m").replace(".", "p")


def _load_ascii_stl_triangles(path: Path) -> np.ndarray:
    triangles: list[list[tuple[float, float, float]]] = []
    current: list[tuple[float, float, float]] = []
    with path.open("r", encoding="utf-8", errors="ignore") as stream:
        for raw_line in stream:
            line = raw_line.strip().split()
            if len(line) != 4 or line[0] != "vertex":
                continue
            current.append((float(line[1]) / 1000.0, float(line[2]) / 1000.0, float(line[3]) / 1000.0))
            if len(current) == 3:
                triangles.append(current)
                current = []
    if not triangles:
        raise ValueError(f"No triangles found in ASCII STL: {path}")
    return np.asarray(triangles, dtype=float)


def _deduplicate_points(points: Iterable[np.ndarray], tol: float = 1.0e-9) -> list[np.ndarray]:
    unique: list[np.ndarray] = []
    for point in points:
        arr = np.asarray(point, dtype=float)
        if any(np.linalg.norm(arr - existing) <= tol for existing in unique):
            continue
        unique.append(arr)
    return unique


def _triangle_plane_segment(triangle: np.ndarray, section_y_m: float, tol: float = 1.0e-9) -> tuple[np.ndarray, np.ndarray] | None:
    distances = triangle[:, 1] - float(section_y_m)
    candidates: list[np.ndarray] = []
    for i0, i1 in ((0, 1), (1, 2), (2, 0)):
        p0 = triangle[i0]
        p1 = triangle[i1]
        d0 = distances[i0]
        d1 = distances[i1]

        if abs(d0) <= tol and abs(d1) <= tol:
            candidates.extend([p0[[0, 2]], p1[[0, 2]]])
            continue
        if abs(d0) <= tol:
            candidates.append(p0[[0, 2]])
            continue
        if abs(d1) <= tol:
            candidates.append(p1[[0, 2]])
            continue
        if d0 * d1 < 0.0:
            t = -d0 / (d1 - d0)
            point = p0 + t * (p1 - p0)
            candidates.append(point[[0, 2]])

    unique = _deduplicate_points(candidates, tol=1.0e-8)
    if len(unique) < 2:
        return None
    if len(unique) > 2:
        max_pair = None
        max_distance = -1.0
        for idx in range(len(unique)):
            for jdx in range(idx + 1, len(unique)):
                dist = float(np.linalg.norm(unique[idx] - unique[jdx]))
                if dist > max_distance:
                    max_distance = dist
                    max_pair = (unique[idx], unique[jdx])
        if max_pair is None:
            return None
        return max_pair
    return unique[0], unique[1]


def _slice_stl_at_y(triangles: np.ndarray, section_y_m: float) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    segments: list[tuple[tuple[float, float], tuple[float, float]]] = []
    for triangle in triangles:
        segment = _triangle_plane_segment(triangle, section_y_m)
        if segment is None:
            continue
        point_a, point_b = segment
        if np.linalg.norm(point_a - point_b) <= 1.0e-9:
            continue
        a = (round(float(point_a[0]), 8), round(float(point_a[1]), 8))
        b = (round(float(point_b[0]), 8), round(float(point_b[1]), 8))
        if a == b:
            continue
        segments.append((a, b))
    if not segments:
        raise ValueError(f"No intersection segments found at y={section_y_m:.6f} m")
    return segments


def _ordered_loop_from_segments(
    segments: Iterable[tuple[tuple[float, float], tuple[float, float]]],
) -> np.ndarray:
    adjacency: defaultdict[tuple[float, float], list[tuple[float, float]]] = defaultdict(list)
    for point_a, point_b in segments:
        adjacency[point_a].append(point_b)
        adjacency[point_b].append(point_a)

    visited_edges: set[tuple[tuple[float, float], tuple[float, float]]] = set()
    loops: list[np.ndarray] = []
    for start in adjacency:
        for neighbor in adjacency[start]:
            edge = tuple(sorted((start, neighbor)))
            if edge in visited_edges:
                continue
            loop = [start]
            previous = None
            current = start
            while True:
                neighbors = adjacency[current]
                next_point = None
                for candidate in neighbors:
                    candidate_edge = tuple(sorted((current, candidate)))
                    if candidate_edge in visited_edges:
                        continue
                    if candidate == previous and len(neighbors) > 1:
                        continue
                    next_point = candidate
                    break
                if next_point is None:
                    break
                visited_edges.add(tuple(sorted((current, next_point))))
                if next_point == start:
                    break
                loop.append(next_point)
                previous, current = current, next_point
            if len(loop) >= 3:
                loops.append(np.asarray(loop, dtype=float))

    if not loops:
        raise ValueError("Could not reconstruct a closed loop from the STL slice")

    def loop_score(points: np.ndarray) -> float:
        closed = np.vstack([points, points[0]])
        seg = np.diff(closed, axis=0)
        return float(np.sum(np.linalg.norm(seg, axis=1)))

    points = max(loops, key=loop_score)
    leading_edge_idx = int(np.argmin(points[:, 0]))
    if np.mean(points[: leading_edge_idx + 1, 1]) < np.mean(points[leading_edge_idx:, 1]):
        points = points[::-1]

    max_x = float(np.max(points[:, 0]))
    candidate_idx = np.where(np.isclose(points[:, 0], max_x, atol=5.0e-5))[0]
    start_idx = int(candidate_idx[np.argmax(points[candidate_idx, 1])])
    points = np.vstack([points[start_idx:], points[:start_idx], points[start_idx : start_idx + 1]])

    leading_edge_idx = int(np.argmin(points[:, 0]))
    if np.mean(points[: leading_edge_idx + 1, 1]) < np.mean(points[leading_edge_idx:, 1]):
        points = points[::-1]
        points = points[:-1]
        max_x = float(np.max(points[:, 0]))
        candidate_idx = np.where(np.isclose(points[:, 0], max_x, atol=5.0e-5))[0]
        start_idx = int(candidate_idx[np.argmax(points[candidate_idx, 1])])
        points = np.vstack([points[start_idx:], points[:start_idx], points[start_idx : start_idx + 1]])
    return points


def _section_anchor_payload(section_y_m: float, fit_result) -> dict[str, object]:
    thickness = np.asarray(fit_result.fit_y_upper - fit_result.fit_y_lower, dtype=float)
    thickness_idx = int(np.argmax(thickness))
    return {
        "anchor_y_m": float(section_y_m),
        "anchor_le_z_m": float(fit_result.le_z),
        "anchor_twist_deg": float(fit_result.twist_deg),
        "upper_coeffs": np.asarray(fit_result.upper_cst, dtype=float).tolist(),
        "lower_coeffs": np.asarray(-fit_result.lower_cst, dtype=float).tolist(),
        "tc_max": float(thickness[thickness_idx]),
        "x_tmax": float(fit_result.x_fit[thickness_idx]),
        "te_thickness": float(fit_result.te_thickness),
        "chord_m": float(fit_result.chord),
        "le_x_m": float(fit_result.le_x),
        "max_abs_error": float(fit_result.max_abs_error),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract a CST-fit airfoil section from a CTA STL by plane slicing.")
    parser.add_argument("--stl", type=Path, default=DEFAULT_STL_PATH, help="Path to the ASCII STL file.")
    parser.add_argument("--section-y", type=float, default=DEFAULT_SECTION_Y_M, help="Spanwise section y in meters.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    stl_path = Path(args.stl)
    section_y_m = float(args.section_y)
    if not stl_path.exists():
        raise FileNotFoundError(f"Could not find STL at {stl_path}")

    triangles = _load_ascii_stl_triangles(stl_path)
    segments = _slice_stl_at_y(triangles, section_y_m)
    section_points = _ordered_loop_from_segments(segments)

    normalized = normalize_airfoil_section(section_points)
    fit_result = fit_airfoil_section_cst(section_points, FIT_OPTIONS)
    anchor_payload = _section_anchor_payload(section_y_m, fit_result)

    tag = _safe_y_tag(section_y_m)
    airfoil_dir = CTA_DIR / "airfoils"
    output_dir = CTA_DIR / "outputs" / "wing"
    airfoil_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    raw_dat = airfoil_dir / f"stl_section_y{tag}_raw_local.dat"
    fit_dat = airfoil_dir / f"stl_section_y{tag}_cst_fit.dat"
    summary_json = airfoil_dir / f"stl_section_y{tag}_anchor_summary.json"
    plot_png = output_dir / f"stl_section_y{tag}_fit_check.png"

    write_airfoil_dat(str(raw_dat), normalized.x_fit, normalized.y_upper, normalized.y_lower, name=f"STL_Y{tag}_RAW")
    write_airfoil_dat(
        str(fit_dat),
        fit_result.x_fit,
        fit_result.fit_y_upper,
        fit_result.fit_y_lower,
        name=f"STL_Y{tag}_CST_FIT",
    )
    summary_json.write_text(json.dumps(anchor_payload, indent=2), encoding="utf-8")

    fig, ax = plt.subplots(figsize=(12.0, 4.2), constrained_layout=True)
    ax.plot(section_points[:, 0], section_points[:, 1], color="#2563eb", linewidth=2.0, label="STL slice raw")

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
        f"STL y={section_y_m:.3f} m section | "
        f"chord={fit_result.chord:.3f} m | "
        f"twist={fit_result.twist_deg:.3f} deg | "
        f"LE_z={fit_result.le_z:.3f} m"
    )
    fig.savefig(plot_png, dpi=220, bbox_inches="tight")
    plt.close(fig)

    print(f"STL raw local DAT written to: {raw_dat}")
    print(f"STL CST-fit DAT written to: {fit_dat}")
    print(f"STL anchor summary written to: {summary_json}")
    print(f"STL fit check plot written to: {plot_png}")


if __name__ == "__main__":
    main()
