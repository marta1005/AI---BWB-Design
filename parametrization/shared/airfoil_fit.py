from __future__ import annotations

from dataclasses import dataclass
from math import atan2, degrees
from pathlib import Path
from typing import Dict

import numpy as np

from .cst import KulfanCSTAirfoil, cosine_spacing, fit_kulfan_airfoil_coefficients


@dataclass(frozen=True)
class CSTAirfoilFitOptions:
    degree: int = 5
    n1: float = 0.5
    n2: float = 1.0
    shared_leading_edge: bool = False
    smoothness_weight: float = 0.0
    fit_xmin: float = 1.0e-4
    fit_xmax: float = 0.99
    sample_count: int = 301


@dataclass(frozen=True)
class NormalizedAirfoilSection:
    x_fit: np.ndarray
    y_upper: np.ndarray
    y_lower: np.ndarray
    chord: float
    twist_deg: float
    te_thickness: float
    te_mid_x: float
    te_mid_z: float
    le_x: float
    le_z: float


@dataclass(frozen=True)
class AirfoilSectionFitResult:
    x_fit: np.ndarray
    y_upper: np.ndarray
    y_lower: np.ndarray
    fit_y_upper: np.ndarray
    fit_y_lower: np.ndarray
    upper_cst: np.ndarray
    lower_cst: np.ndarray
    chord: float
    twist_deg: float
    te_thickness: float
    te_mid_x: float
    te_mid_z: float
    le_x: float
    le_z: float
    fit_xmin: float
    fit_xmax: float
    rmse_upper: float
    rmse_lower: float
    max_abs_error: float


def load_xyz_sections_by_span(path: str | Path) -> Dict[float, np.ndarray]:
    sections: Dict[float, list[list[float]]] = {}
    for raw_line in Path(path).read_text(encoding="utf-8").splitlines():
        parts = raw_line.split()
        if len(parts) != 3:
            continue
        x_val, y_val, z_val = (float(parts[0]), float(parts[1]), float(parts[2]))
        sections.setdefault(y_val, []).append([x_val, z_val])
    return {yy: np.asarray(points, dtype=float) for yy, points in sections.items()}


def split_airfoil_upper_lower(section: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    le_idx = int(np.argmin(section[:, 0]))
    upper = section[: le_idx + 1]
    lower = section[le_idx:]
    return upper, lower


def normalize_airfoil_section(section: np.ndarray, sample_count: int = 301) -> NormalizedAirfoilSection:
    upper_raw, lower_raw = split_airfoil_upper_lower(section)

    te_upper = np.asarray(section[0], dtype=float)
    te_lower = np.asarray(section[-1], dtype=float)
    te_mid = 0.5 * (te_upper + te_lower)
    le_point = np.asarray(section[int(np.argmin(section[:, 0]))], dtype=float)
    chord_vector = te_mid - le_point
    chord = float(np.linalg.norm(chord_vector))
    if chord <= 0.0:
        raise ValueError(f"invalid chord length {chord}")
    twist_deg = float(degrees(atan2(float(chord_vector[1]), float(chord_vector[0]))))

    ex = chord_vector / chord
    ez = np.array([-ex[1], ex[0]], dtype=float)
    upper_local = np.column_stack(((upper_raw - le_point) @ ex, (upper_raw - le_point) @ ez))
    lower_local = np.column_stack(((lower_raw - le_point) @ ex, (lower_raw - le_point) @ ez))

    # Upper arrives TE->LE in the source order, so reverse it to interpolate LE->TE.
    upper_local = upper_local[::-1]

    upper_x = upper_local[:, 0] / chord
    upper_z = upper_local[:, 1] / chord
    lower_x = lower_local[:, 0] / chord
    lower_z = lower_local[:, 1] / chord

    x_fit = cosine_spacing(sample_count)
    y_upper = np.interp(x_fit, upper_x, upper_z)
    y_lower = np.interp(x_fit, lower_x, lower_z)
    te_thickness = float(y_upper[-1] - y_lower[-1])

    return NormalizedAirfoilSection(
        x_fit=x_fit,
        y_upper=y_upper,
        y_lower=y_lower,
        chord=chord,
        twist_deg=twist_deg,
        te_thickness=te_thickness,
        te_mid_x=float(te_mid[0]),
        te_mid_z=float(te_mid[1]),
        le_x=float(le_point[0]),
        le_z=float(le_point[1]),
    )


def fit_airfoil_section_cst(
    section: np.ndarray,
    options: CSTAirfoilFitOptions | None = None,
) -> AirfoilSectionFitResult:
    fit_options = options or CSTAirfoilFitOptions()
    normalized = normalize_airfoil_section(section, sample_count=fit_options.sample_count)

    x_fit = np.asarray(normalized.x_fit, dtype=float)
    y_upper = np.asarray(normalized.y_upper, dtype=float)
    y_lower = np.asarray(normalized.y_lower, dtype=float)
    te_thickness = float(normalized.te_thickness)

    fit_mask = (x_fit >= float(fit_options.fit_xmin)) & (x_fit < float(fit_options.fit_xmax))
    if np.count_nonzero(fit_mask) < fit_options.degree + 1:
        raise ValueError(
            f"not enough fit samples for degree {fit_options.degree}: "
            f"{np.count_nonzero(fit_mask)} available"
        )

    upper_cst, lower_cst = fit_kulfan_airfoil_coefficients(
        x=x_fit[fit_mask],
        y_upper=y_upper[fit_mask],
        y_lower=y_lower[fit_mask],
        degree=fit_options.degree,
        n1=fit_options.n1,
        n2=fit_options.n2,
        te_thickness=te_thickness,
        smoothness_weight=fit_options.smoothness_weight,
        shared_leading_edge=fit_options.shared_leading_edge,
    )

    airfoil = KulfanCSTAirfoil(
        degree=fit_options.degree,
        n1=fit_options.n1,
        n2=fit_options.n2,
        shared_leading_edge=fit_options.shared_leading_edge,
    )
    coeffs = np.concatenate([upper_cst, lower_cst], dtype=float)
    fit_y_upper, fit_y_lower = airfoil.evaluate(x_fit, coeffs, te_thickness=te_thickness)

    error_upper = fit_y_upper - y_upper
    error_lower = fit_y_lower - y_lower

    return AirfoilSectionFitResult(
        x_fit=x_fit,
        y_upper=y_upper,
        y_lower=y_lower,
        fit_y_upper=fit_y_upper,
        fit_y_lower=fit_y_lower,
        upper_cst=upper_cst,
        lower_cst=lower_cst,
        chord=float(normalized.chord),
        twist_deg=float(normalized.twist_deg),
        te_thickness=te_thickness,
        te_mid_x=float(normalized.te_mid_x),
        te_mid_z=float(normalized.te_mid_z),
        le_x=float(normalized.le_x),
        le_z=float(normalized.le_z),
        fit_xmin=float(fit_options.fit_xmin),
        fit_xmax=float(fit_options.fit_xmax),
        rmse_upper=float(np.sqrt(np.mean(error_upper**2))),
        rmse_lower=float(np.sqrt(np.mean(error_lower**2))),
        max_abs_error=float(max(np.max(np.abs(error_upper)), np.max(np.abs(error_lower)))),
    )
