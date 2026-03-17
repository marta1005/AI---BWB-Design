from dataclasses import dataclass, field, fields
from functools import lru_cache
from math import comb
from typing import Dict, Optional, Tuple

import numpy as np

from .specs import (
    AnchoredSpanwiseLaw,
    PlanformSpec,
    SamplingSpec,
    SectionCSTSpec,
    SectionFamilySpec,
    SectionProfileRelationSpec,
    SectionedBWBModelConfig,
    SpanwiseLawSpec,
)
from .topology import SectionedBWBTopologySpec


def _full_surface_seed(
    degree: int,
    leading: float,
    peak: float,
    width: float,
    center: float,
    tail: float,
) -> Tuple[float, ...]:
    modes = np.arange(degree + 1, dtype=float)
    gaussian = np.exp(-0.5 * ((modes - center) / max(width, 1e-6)) ** 2)
    gaussian /= max(float(np.max(gaussian)), 1.0)
    taper = np.exp(-0.55 * modes)
    coeffs = leading * taper + peak * gaussian
    coeffs[-1] = tail
    coeffs[0] = max(coeffs[0], leading)
    return tuple(float(value) for value in coeffs)


def _flatten_tuple_names(prefix: str, values: Tuple[float, ...]) -> Tuple[str, ...]:
    return tuple(f"{prefix}_{idx}" for idx in range(len(values)))


def _cst_bound_window(reference: float, coeff_idx: int) -> Tuple[float, float]:
    if coeff_idx <= 1:
        delta = max(0.012, 0.18 * float(reference))
    elif coeff_idx <= 3:
        delta = max(0.010, 0.20 * float(reference))
    elif coeff_idx == 4:
        delta = max(0.008, 0.22 * float(reference))
    else:
        delta = max(0.004, 0.25 * float(reference))
    lower = max(0.0, float(reference) - delta)
    upper = float(reference) + delta
    return lower, upper


# The paper does not provide direct Bernstein-coefficient bounds.
# Instead, it fits fifth-order CST coefficients to a family of NACA 4-digit
# airfoils with N1=0.5 and N2=1.0. We reproduce that logic here and derive
# coefficient-wise bounds from a section-specific NACA family centered on the
# intended airfoil thickness level of each station.
_SECTION_PAPER_TC_RANGES: Dict[str, Tuple[float, float]] = {
    "c1": (0.14, 0.19),
    "c2": (0.11, 0.16),
    "c3": (0.07, 0.10),
    "c4": (0.05, 0.08),
}


def _naca4_surfaces(x: np.ndarray, m: float, p: float, t: float) -> Tuple[np.ndarray, np.ndarray]:
    x = np.asarray(x, dtype=float)
    yt = 5.0 * float(t) * (
        0.2969 * np.sqrt(np.clip(x, 0.0, 1.0))
        - 0.1260 * x
        - 0.3516 * x**2
        + 0.2843 * x**3
        - 0.1015 * x**4
    )

    if m <= 0.0 or p <= 0.0:
        yc = np.zeros_like(x)
        dyc_dx = np.zeros_like(x)
    else:
        yc = np.where(
            x < p,
            m / (p**2) * (2.0 * p * x - x**2),
            m / ((1.0 - p) ** 2) * ((1.0 - 2.0 * p) + 2.0 * p * x - x**2),
        )
        dyc_dx = np.where(
            x < p,
            2.0 * m / (p**2) * (p - x),
            2.0 * m / ((1.0 - p) ** 2) * (p - x),
        )

    del dyc_dx
    return yc + yt, yc - yt


def _cosine_spacing(n: int) -> np.ndarray:
    beta = np.linspace(0.0, np.pi, int(n))
    return 0.5 * (1.0 - np.cos(beta))


def _bernstein_matrix(x: np.ndarray, degree: int) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    basis = np.zeros((x.size, degree + 1), dtype=float)
    for idx in range(degree + 1):
        basis[:, idx] = float(comb(degree, idx)) * (x**idx) * ((1.0 - x) ** (degree - idx))
    return basis


def _fit_cst_surface_coeffs(
    x: np.ndarray,
    y_upper: np.ndarray,
    y_lower: np.ndarray,
    degree: int,
    n1: float,
    n2: float,
) -> Tuple[np.ndarray, np.ndarray]:
    x = np.asarray(x, dtype=float)
    y_upper = np.asarray(y_upper, dtype=float)
    y_lower = np.asarray(y_lower, dtype=float)

    interior = (x > 1.0e-4) & (x < 0.99)
    x_fit = x[interior]
    class_fun = (x_fit**float(n1)) * ((1.0 - x_fit) ** float(n2))
    basis = _bernstein_matrix(x_fit, degree)

    upper_rhs = y_upper[interior] / class_fun
    lower_rhs = -y_lower[interior] / class_fun

    regularization = 1.0e-6 * np.eye(degree + 1)
    upper_coeffs = np.linalg.solve(basis.T @ basis + regularization, basis.T @ upper_rhs)
    lower_coeffs = np.linalg.solve(basis.T @ basis + regularization, basis.T @ lower_rhs)
    return upper_coeffs.astype(float), lower_coeffs.astype(float)


@lru_cache(maxsize=None)
def _paper_section_cst_bounds(
    section_name: str,
    degree: int = 5,
    n1: float = 0.5,
    n2: float = 1.0,
) -> Dict[str, Tuple[float, ...]]:
    if section_name not in _SECTION_PAPER_TC_RANGES:
        raise KeyError(f"Unsupported section name for paper CST bounds: {section_name}")

    # Paper-compatible family:
    # - NACA 4-digit
    # - max camber 0% to 6%
    # - camber position 0% to 50% chord
    # - fifth-order CST with N1=0.5, N2=1.0
    m_values = (0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06)
    p_values = (0.10, 0.20, 0.30, 0.40, 0.50)
    t_low, t_high = _SECTION_PAPER_TC_RANGES[section_name]
    t_values = tuple(float(value) for value in np.linspace(t_low, t_high, 6))

    x = _cosine_spacing(241)
    upper_coeff_samples = []
    lower_coeff_samples = []

    for t_value in t_values:
        for m_value in m_values:
            if m_value == 0.0:
                yu, yl = _naca4_surfaces(x, m=0.0, p=0.0, t=t_value)
                upper, lower = _fit_cst_surface_coeffs(x, yu, yl, degree, n1, n2)
                upper_coeff_samples.append(upper)
                lower_coeff_samples.append(lower)
                continue

            for p_value in p_values:
                yu, yl = _naca4_surfaces(x, m=m_value, p=p_value, t=t_value)
                upper, lower = _fit_cst_surface_coeffs(x, yu, yl, degree, n1, n2)
                upper_coeff_samples.append(upper)
                lower_coeff_samples.append(lower)

    upper_arr = np.asarray(upper_coeff_samples, dtype=float)
    lower_arr = np.asarray(lower_coeff_samples, dtype=float)

    upper_q10 = np.quantile(upper_arr, 0.10, axis=0)
    upper_q90 = np.quantile(upper_arr, 0.90, axis=0)
    lower_q10 = np.quantile(lower_arr, 0.10, axis=0)
    lower_q90 = np.quantile(lower_arr, 0.90, axis=0)

    def _halfwidth(low: np.ndarray, high: np.ndarray) -> Tuple[float, ...]:
        span = np.maximum(high - low, 0.0)
        floor = np.array((0.010, 0.010, 0.008, 0.008, 0.006, 0.004), dtype=float)
        width = np.maximum(floor, 0.25 * span)
        return tuple(float(value) for value in width)

    return {
        "upper_halfwidth": _halfwidth(upper_q10, upper_q90),
        "lower_halfwidth": _halfwidth(lower_q10, lower_q90),
    }


@dataclass
class SectionedBWBDesignVariables:
    span: float = 39.5
    b1_span_ratio: float = 0.18
    b2_span_ratio: float = 0.12
    b3_span_ratio: float = 0.70

    le_root_x: float = 0.0
    c1_root_chord: float = 40.0
    c2_c1_ratio: float = 0.70
    c3_c1_ratio: float = 0.23
    c4_c1_ratio: float = 0.08
    s1_deg: float = 58.0
    s2_deg: float = 50.0
    s3_deg: float = 32.0
    nose_blend_y: float = 2.50

    cst_n1: float = 0.50
    cst_n2: float = 1.00

    dihedral_deg: float = 0.0
    twist_c1_deg: float = 0.0
    twist_c2_deg: float = 0.0
    twist_c3_deg: float = -1.0
    twist_c4_deg: float = -3.0
    camber_c1: float = 0.0
    camber_c2: float = 0.0
    camber_c3: float = 0.0
    camber_c4: float = 0.0

    c1_tc_max: float = 0.22
    c2_tc_max: float = 0.18
    c3_tc_max: float = 0.11
    c4_tc_max: float = 0.08

    c1_x_tmax: float = 0.33
    c2_x_tmax: float = 0.31
    c3_x_tmax: float = 0.28
    c4_x_tmax: float = 0.24

    c1_te_thickness: float = 0.0020
    c2_te_thickness: float = 0.0020
    c3_te_thickness: float = 0.0015
    c4_te_thickness: float = 0.0010

    c1_upper_cst: Tuple[float, ...] = field(
        default_factory=lambda: _full_surface_seed(5, 0.24, 0.26, 1.30, 1.7, 0.05)
    )
    c1_lower_cst: Tuple[float, ...] = field(
        default_factory=lambda: _full_surface_seed(5, 0.13, 0.06, 1.45, 2.7, 0.010)
    )
    c2_upper_cst: Tuple[float, ...] = field(
        default_factory=lambda: _full_surface_seed(5, 0.20, 0.21, 1.25, 1.8, 0.04)
    )
    c2_lower_cst: Tuple[float, ...] = field(
        default_factory=lambda: _full_surface_seed(5, 0.11, 0.05, 1.40, 2.7, 0.008)
    )
    c3_upper_cst: Tuple[float, ...] = field(
        default_factory=lambda: _full_surface_seed(5, 0.14, 0.13, 1.20, 1.9, 0.025)
    )
    c3_lower_cst: Tuple[float, ...] = field(
        default_factory=lambda: _full_surface_seed(5, 0.075, 0.030, 1.30, 2.8, 0.005)
    )
    c4_upper_cst: Tuple[float, ...] = field(
        default_factory=lambda: _full_surface_seed(5, 0.11, 0.10, 1.15, 2.0, 0.018)
    )
    c4_lower_cst: Tuple[float, ...] = field(
        default_factory=lambda: _full_surface_seed(5, 0.055, 0.018, 1.25, 3.0, 0.003)
    )

    @classmethod
    def reference_seed(cls) -> "SectionedBWBDesignVariables":
        return cls()

    @classmethod
    def _tuple_field_map(cls) -> Dict[str, int]:
        default = cls.reference_seed()
        mapping: Dict[str, int] = {}
        for item in fields(default):
            value = getattr(default, item.name)
            if isinstance(value, tuple):
                mapping[item.name] = len(value)
        return mapping

    @classmethod
    def variable_names(cls) -> Tuple[str, ...]:
        default = cls.reference_seed()
        names = []
        for item in fields(default):
            value = getattr(default, item.name)
            if isinstance(value, tuple):
                names.extend(_flatten_tuple_names(item.name, value))
            else:
                names.append(item.name)
        return tuple(names)

    @classmethod
    def default_bounds(cls) -> Dict[str, Tuple[float, float]]:
        bounds = {
            "span": (35.0, 45.0),
            "b1_span_ratio": (0.16, 0.24),
            "b2_span_ratio": (0.09, 0.16),
            "b3_span_ratio": (0.60, 0.75),
            "le_root_x": (-5.0, 8.0),
            "c1_root_chord": (36.0, 44.0),
            "c2_c1_ratio": (0.55, 0.85),
            "c3_c1_ratio": (0.18, 0.28),
            "c4_c1_ratio": (0.06, 0.09),
            "s1_deg": (40.0, 60.0),
            "s2_deg": (40.0, 60.0),
            "s3_deg": (24.0, 40.0),
            "nose_blend_y": (1.5, 4.0),
            "cst_n1": (0.45, 0.55),
            "cst_n2": (0.95, 1.10),
            "dihedral_deg": (0.0, 12.0),
            "twist_c1_deg": (-1.0, 1.0),
            "twist_c2_deg": (-1.5, 1.0),
            "twist_c3_deg": (-4.0, 0.0),
            "twist_c4_deg": (-7.0, -0.5),
            "camber_c1": (-0.03, 0.03),
            "camber_c2": (-0.03, 0.03),
            "camber_c3": (-0.03, 0.03),
            "camber_c4": (-0.03, 0.03),
            "c1_tc_max": (0.14, 0.28),
            "c2_tc_max": (0.12, 0.24),
            "c3_tc_max": (0.08, 0.18),
            "c4_tc_max": (0.05, 0.12),
            "c1_x_tmax": (0.22, 0.48),
            "c2_x_tmax": (0.20, 0.45),
            "c3_x_tmax": (0.18, 0.42),
            "c4_x_tmax": (0.14, 0.35),
            "c1_te_thickness": (0.0, 0.010),
            "c2_te_thickness": (0.0, 0.010),
            "c3_te_thickness": (0.0, 0.008),
            "c4_te_thickness": (0.0, 0.006),
        }
        default = cls.reference_seed()
        paper_cst_bounds = {
            section_name: _paper_section_cst_bounds(section_name, degree=5, n1=default.cst_n1, n2=default.cst_n2)
            for section_name in ("c1", "c2", "c3", "c4")
        }
        for tuple_name, tuple_size in cls._tuple_field_map().items():
            tuple_values = getattr(default, tuple_name)
            for idx in range(tuple_size):
                key = f"{tuple_name}_{idx}"
                reference = float(tuple_values[idx])
                tuple_lower = None
                tuple_upper = None
                if tuple_name.endswith("_upper_cst"):
                    section_name = tuple_name.split("_")[0]
                    halfwidth = paper_cst_bounds[section_name]["upper_halfwidth"][idx]
                    tuple_lower = max(0.0, reference - halfwidth)
                    tuple_upper = reference + halfwidth
                elif tuple_name.endswith("_lower_cst"):
                    section_name = tuple_name.split("_")[0]
                    halfwidth = paper_cst_bounds[section_name]["lower_halfwidth"][idx]
                    tuple_lower = max(0.0, reference - halfwidth)
                    tuple_upper = reference + halfwidth

                if tuple_lower is None or tuple_upper is None:
                    bounds[key] = _cst_bound_window(reference, idx)
                else:
                    bounds[key] = (float(tuple_lower), float(tuple_upper))
        return bounds

    def as_vector(self) -> np.ndarray:
        values = []
        for item in fields(self):
            value = getattr(self, item.name)
            if isinstance(value, tuple):
                values.extend(float(entry) for entry in value)
            else:
                values.append(float(value))
        return np.asarray(values, dtype=float)

    @classmethod
    def from_vector(cls, vector: np.ndarray) -> "SectionedBWBDesignVariables":
        values = np.asarray(vector, dtype=float).ravel()
        template = cls.reference_seed()
        kwargs = {}
        offset = 0
        for item in fields(template):
            current = getattr(template, item.name)
            if isinstance(current, tuple):
                size = len(current)
                kwargs[item.name] = tuple(float(value) for value in values[offset : offset + size])
                offset += size
            else:
                kwargs[item.name] = float(values[offset])
                offset += 1
        if offset != values.size:
            raise ValueError(f"Expected {offset} design variables, received {values.size}")
        return cls(**kwargs)

    def as_normalized_vector(self) -> np.ndarray:
        bounds = self.default_bounds()
        lower = np.array([bounds[name][0] for name in self.variable_names()], dtype=float)
        upper = np.array([bounds[name][1] for name in self.variable_names()], dtype=float)
        return (self.as_vector() - lower) / (upper - lower)

    @classmethod
    def from_normalized_vector(cls, vector: np.ndarray) -> "SectionedBWBDesignVariables":
        values = np.asarray(vector, dtype=float).ravel()
        names = cls.variable_names()
        if values.size != len(names):
            raise ValueError(f"Expected {len(names)} normalized design variables, received {values.size}")
        bounds = cls.default_bounds()
        lower = np.array([bounds[name][0] for name in names], dtype=float)
        upper = np.array([bounds[name][1] for name in names], dtype=float)
        return cls.from_vector(lower + values * (upper - lower))

    def validate(self) -> None:
        if self.span <= 0.0:
            raise ValueError(f"span must be positive, got {self.span:.6f}")

        b_ratios = (self.b1_span_ratio, self.b2_span_ratio, self.b3_span_ratio)
        if not all(value > 0.0 for value in b_ratios):
            raise ValueError(f"B ratios must be positive, got {b_ratios}")
        if abs(sum(b_ratios) - 1.0) > 1e-9:
            raise ValueError(
                "b1_span_ratio + b2_span_ratio + b3_span_ratio must equal 1.0 "
                f"for semispan-based segment ratios, got {sum(b_ratios):.12f}"
            )

        chord_ratios = (self.c2_c1_ratio, self.c3_c1_ratio, self.c4_c1_ratio)
        if self.c1_root_chord <= 0.0 or not all(value > 0.0 for value in chord_ratios):
            raise ValueError(
                "root chord and chord ratios must be positive, "
                f"got c1={self.c1_root_chord}, ratios={chord_ratios}"
            )

        for label, value in (("s1_deg", self.s1_deg), ("s2_deg", self.s2_deg), ("s3_deg", self.s3_deg)):
            if not (0.0 < value < 85.0):
                raise ValueError(f"{label} must lie in (0, 85), got {value}")

        if self.nose_blend_y <= 0.0:
            raise ValueError(f"nose_blend_y must be positive, got {self.nose_blend_y}")
        if self.cst_n1 <= 0.0 or self.cst_n2 <= 0.0:
            raise ValueError(f"cst_n1 and cst_n2 must be positive, got {(self.cst_n1, self.cst_n2)}")

        for label, tc_value in (
            ("c1_tc_max", self.c1_tc_max),
            ("c2_tc_max", self.c2_tc_max),
            ("c3_tc_max", self.c3_tc_max),
            ("c4_tc_max", self.c4_tc_max),
        ):
            if tc_value <= 0.0:
                raise ValueError(f"{label} must be positive, got {tc_value}")

        for label, x_tmax_value in (
            ("c1_x_tmax", self.c1_x_tmax),
            ("c2_x_tmax", self.c2_x_tmax),
            ("c3_x_tmax", self.c3_x_tmax),
            ("c4_x_tmax", self.c4_x_tmax),
        ):
            if not (0.15 <= x_tmax_value <= 0.65):
                raise ValueError(f"{label} must lie in [0.15, 0.65], got {x_tmax_value}")

        for label, te_value in (
            ("c1_te_thickness", self.c1_te_thickness),
            ("c2_te_thickness", self.c2_te_thickness),
            ("c3_te_thickness", self.c3_te_thickness),
            ("c4_te_thickness", self.c4_te_thickness),
        ):
            if te_value < 0.0:
                raise ValueError(f"{label} must be non-negative, got {te_value}")

        for label, coeffs in (
            ("c1_upper_cst", self.c1_upper_cst),
            ("c1_lower_cst", self.c1_lower_cst),
            ("c2_upper_cst", self.c2_upper_cst),
            ("c2_lower_cst", self.c2_lower_cst),
            ("c3_upper_cst", self.c3_upper_cst),
            ("c3_lower_cst", self.c3_lower_cst),
            ("c4_upper_cst", self.c4_upper_cst),
            ("c4_lower_cst", self.c4_lower_cst),
        ):
            if len(coeffs) != 6:
                raise ValueError(f"{label} must have 6 CST coefficients, got {len(coeffs)}")
            if any(value < 0.0 for value in coeffs):
                raise ValueError(f"{label} must be non-negative")

    def to_model_config(
        self,
        profile_generation_mode: str = "cst_only",
        profile_relations: Optional[Tuple[SectionProfileRelationSpec, ...]] = None,
    ) -> SectionedBWBModelConfig:
        self.validate()
        topology = SectionedBWBTopologySpec(
            span=self.span,
            b1_span_ratio=self.b1_span_ratio,
            b2_span_ratio=self.b2_span_ratio,
            b3_span_ratio=self.b3_span_ratio,
        )
        topology.anchor_y = topology.y_sections

        return SectionedBWBModelConfig(
            topology=topology,
            planform=PlanformSpec(
                le_root_x=self.le_root_x,
                c1_root_chord=self.c1_root_chord,
                c2_c1_ratio=self.c2_c1_ratio,
                c3_c1_ratio=self.c3_c1_ratio,
                c4_c1_ratio=self.c4_c1_ratio,
                s1_deg=self.s1_deg,
                s2_deg=self.s2_deg,
                s3_deg=self.s3_deg,
                continuity_order=2,
                blend_fraction=0.10,
                min_linear_core_fraction=0.75,
                symmetry_blend_y=self.nose_blend_y,
            ),
            sections=SectionFamilySpec(
                cst_degree=5,
                n1=self.cst_n1,
                n2=self.cst_n2,
                profile_generation_mode=profile_generation_mode,
                c1_spec=SectionCSTSpec(
                    upper_coeffs=self.c1_upper_cst,
                    lower_coeffs=self.c1_lower_cst,
                    tc_max=self.c1_tc_max,
                    x_tmax=self.c1_x_tmax,
                    te_thickness=self.c1_te_thickness,
                ),
                c2_spec=SectionCSTSpec(
                    upper_coeffs=self.c2_upper_cst,
                    lower_coeffs=self.c2_lower_cst,
                    tc_max=self.c2_tc_max,
                    x_tmax=self.c2_x_tmax,
                    te_thickness=self.c2_te_thickness,
                ),
                c3_spec=SectionCSTSpec(
                    upper_coeffs=self.c3_upper_cst,
                    lower_coeffs=self.c3_lower_cst,
                    tc_max=self.c3_tc_max,
                    x_tmax=self.c3_x_tmax,
                    te_thickness=self.c3_te_thickness,
                ),
                c4_spec=SectionCSTSpec(
                    upper_coeffs=self.c4_upper_cst,
                    lower_coeffs=self.c4_lower_cst,
                    tc_max=self.c4_tc_max,
                    x_tmax=self.c4_x_tmax,
                    te_thickness=self.c4_te_thickness,
                ),
                profile_relations=profile_relations
                if profile_relations is not None
                else (
                    SectionProfileRelationSpec(),
                    SectionProfileRelationSpec(),
                    SectionProfileRelationSpec(),
                    SectionProfileRelationSpec(),
                ),
            ),
            spanwise=SpanwiseLawSpec(
                dihedral_deg=self.dihedral_deg,
                twist_deg=AnchoredSpanwiseLaw(
                    section_indices=(0, 1, 2, 3),
                    values=(self.twist_c1_deg, self.twist_c2_deg, self.twist_c3_deg, self.twist_c4_deg),
                    interpolation="pyspline",
                ),
                camber_delta=AnchoredSpanwiseLaw(
                    section_indices=(0, 1, 2, 3),
                    values=(self.camber_c1, self.camber_c2, self.camber_c3, self.camber_c4),
                    interpolation="pyspline",
                ),
            ),
            sampling=SamplingSpec(
                section_interpolation="pyspline",
            ),
        )
