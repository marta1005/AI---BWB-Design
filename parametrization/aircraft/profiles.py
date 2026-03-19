from __future__ import annotations

from dataclasses import dataclass
from math import comb, gamma, isclose, radians, sqrt, tan
from typing import Dict, Tuple

import numpy as np

from parametrization.core.cst_sharedle import KulfanCSTAirfoil, cosine_spacing


@dataclass(frozen=True)
class ProfileSample:
    x: np.ndarray
    yu: np.ndarray
    yl: np.ndarray


@dataclass(frozen=True)
class ICSTConstraint:
    psi: float
    value: float
    order: int = 0
    label: str = ""

    @classmethod
    def point(cls, psi: float, value: float, label: str = "") -> "ICSTConstraint":
        return cls(psi=float(psi), value=float(value), order=0, label=label)

    @classmethod
    def slope(cls, psi: float, value: float, label: str = "") -> "ICSTConstraint":
        return cls(psi=float(psi), value=float(value), order=1, label=label)

    @classmethod
    def derivative(
        cls,
        psi: float,
        value: float,
        order: int,
        label: str = "",
    ) -> "ICSTConstraint":
        return cls(psi=float(psi), value=float(value), order=int(order), label=label)

    def validate(self) -> None:
        if not (0.0 <= self.psi <= 1.0):
            raise ValueError(f"constraint psi must lie in [0, 1], got {self.psi}")
        if self.order < 0:
            raise ValueError(f"constraint order must be >= 0, got {self.order}")
        if self.psi <= 0.0 and self.order > 0:
            raise ValueError(
                "generic iCST derivative constraints at the leading edge are not supported; "
                "use leading-edge radius instead"
            )
        if self.psi >= 1.0 and self.order > 1:
            raise ValueError(
                "generic iCST derivative constraints at the trailing edge are only supported "
                "up to first derivative"
            )


def _falling_factorial(exponent: float, order: int) -> float:
    if order == 0:
        return 1.0
    if order < 0:
        raise ValueError(f"order must be >= 0, got {order}")

    rounded = round(exponent)
    if isclose(exponent, rounded, abs_tol=1e-12) and order > int(rounded):
        return 0.0
    return float(gamma(exponent + 1.0) / gamma(exponent - order + 1.0))


def _class_shape_bernstein_derivative(
    psi: float,
    degree: int,
    coeff_idx: int,
    n1: float,
    n2: float,
    order: int,
) -> float:
    psi_value = float(psi)
    if psi_value < 0.0 or psi_value > 1.0:
        raise ValueError(f"psi must lie in [0, 1], got {psi}")

    if psi_value == 1.0:
        if order == 0:
            return 0.0
        if order == 1:
            return -1.0 if coeff_idx == degree else 0.0
        raise ValueError("derivatives of order > 1 at the trailing edge are not supported")

    bernstein_prefactor = float(comb(degree, coeff_idx))
    exponent_psi = n1 + coeff_idx
    exponent_one_minus = n2 + degree - coeff_idx

    total = 0.0
    for j in range(order + 1):
        left = _falling_factorial(exponent_psi, j)
        right = _falling_factorial(exponent_one_minus, order - j)
        term = (
            comb(order, j)
            * left
            * ((-1.0) ** (order - j))
            * right
            * (psi_value ** (exponent_psi - j))
            * ((1.0 - psi_value) ** (exponent_one_minus - (order - j)))
        )
        total += float(term)
    return float(bernstein_prefactor * total)


def _te_term_derivative(te_y: float, psi: float, order: int) -> float:
    if order == 0:
        return float(te_y) * float(psi)
    if order == 1:
        return float(te_y)
    return 0.0


def _build_icst_matrix(
    constraints: Tuple[ICSTConstraint, ...],
    degree: int,
    n1: float,
    n2: float,
    te_y: float,
) -> tuple[np.ndarray, np.ndarray]:
    matrix = np.zeros((len(constraints), degree + 1), dtype=float)
    rhs = np.zeros(len(constraints), dtype=float)

    for row_idx, constraint in enumerate(constraints):
        constraint.validate()
        for coeff_idx in range(degree + 1):
            matrix[row_idx, coeff_idx] = _class_shape_bernstein_derivative(
                psi=constraint.psi,
                degree=degree,
                coeff_idx=coeff_idx,
                n1=n1,
                n2=n2,
                order=constraint.order,
            )
        rhs[row_idx] = float(constraint.value) - _te_term_derivative(te_y, constraint.psi, constraint.order)
    return matrix, rhs


def _solve_icst_surface_coeffs(
    degree: int,
    leading_edge_radius: float,
    constraints: Tuple[ICSTConstraint, ...],
    te_y: float,
    sign: float,
    n1: float,
    n2: float,
) -> np.ndarray:
    matrix_constraints, rhs_constraints = _build_icst_matrix(
        constraints=constraints,
        degree=degree,
        n1=n1,
        n2=n2,
        te_y=te_y,
    )

    if matrix_constraints.shape[0] != degree:
        raise ValueError(
            f"iCST surface of degree {degree} requires exactly {degree} user constraints plus "
            f"the leading-edge-radius condition, got {matrix_constraints.shape[0]} constraints"
        )

    matrix = np.zeros((degree + 1, degree + 1), dtype=float)
    rhs = np.zeros(degree + 1, dtype=float)
    matrix[0, 0] = 1.0
    rhs[0] = float(sign) * sqrt(2.0 * float(leading_edge_radius))
    matrix[1:, :] = matrix_constraints
    rhs[1:] = rhs_constraints
    return np.linalg.solve(matrix, rhs).astype(float)


def _cubic_hermite_value_slope(
    x0: float,
    y0: float,
    dy0: float,
    x1: float,
    y1: float,
    dy1: float,
    x: float,
) -> tuple[float, float]:
    x0_value = float(x0)
    x1_value = float(x1)
    x_value = float(x)
    if not x0_value < x1_value:
        raise ValueError(f"expected x0 < x1, got x0={x0_value}, x1={x1_value}")
    if x_value < x0_value - 1e-12 or x_value > x1_value + 1e-12:
        raise ValueError(f"x={x_value} must lie in [{x0_value}, {x1_value}]")

    h = x1_value - x0_value
    t = (x_value - x0_value) / h

    h00 = 2.0 * t**3 - 3.0 * t**2 + 1.0
    h10 = t**3 - 2.0 * t**2 + t
    h01 = -2.0 * t**3 + 3.0 * t**2
    h11 = t**3 - t**2

    y_value = (
        h00 * float(y0)
        + h10 * h * float(dy0)
        + h01 * float(y1)
        + h11 * h * float(dy1)
    )

    dh00 = 6.0 * t**2 - 6.0 * t
    dh10 = 3.0 * t**2 - 4.0 * t + 1.0
    dh01 = -6.0 * t**2 + 6.0 * t
    dh11 = 3.0 * t**2 - 2.0 * t

    dy_dt = (
        dh00 * float(y0)
        + dh10 * h * float(dy0)
        + dh01 * float(y1)
        + dh11 * h * float(dy1)
    )
    return float(y_value), float(dy_dt / h)


def _canonical_camber_value_slope(
    max_camber: float,
    x_cmax: float,
    trailing_edge_camber_slope: float,
    psi: float,
) -> tuple[float, float]:
    x_value = float(psi)
    if abs(float(max_camber)) <= 1e-15:
        return _cubic_hermite_value_slope(
            x0=0.0,
            y0=0.0,
            dy0=0.0,
            x1=1.0,
            y1=0.0,
            dy1=float(trailing_edge_camber_slope),
            x=x_value,
        )

    crest_x = float(x_cmax)
    crest_y = float(max_camber)
    leading_slope = 2.0 * crest_y / crest_x

    if x_value <= crest_x:
        return _cubic_hermite_value_slope(
            x0=0.0,
            y0=0.0,
            dy0=leading_slope,
            x1=crest_x,
            y1=crest_y,
            dy1=0.0,
            x=x_value,
        )

    return _cubic_hermite_value_slope(
        x0=crest_x,
        y0=crest_y,
        dy0=0.0,
        x1=1.0,
        y1=0.0,
        dy1=float(trailing_edge_camber_slope),
        x=x_value,
    )


def _surface_leading_edge_radius_from_airfoil_radius(leading_edge_radius: float) -> float:
    return 0.25 * float(leading_edge_radius)


@dataclass(frozen=True)
class CSTAirfoilProfileSpec:
    profile_id: str
    degree: int
    upper_coeffs: Tuple[float, ...]
    lower_coeffs: Tuple[float, ...]
    n1: float = 0.5
    n2: float = 1.0
    te_thickness: float = 0.0
    x_tc_window: Tuple[float, float] = (0.15, 0.65)

    def validate(self) -> None:
        if not self.profile_id:
            raise ValueError("profile_id must be non-empty")
        if self.degree < 1:
            raise ValueError(f"degree must be >= 1, got {self.degree}")
        expected = self.degree + 1
        if len(self.upper_coeffs) != expected or len(self.lower_coeffs) != expected:
            raise ValueError(
                f"profile {self.profile_id!r} expects {expected} CST coefficients per side, "
                f"got upper={len(self.upper_coeffs)}, lower={len(self.lower_coeffs)}"
            )
        if self.n1 <= 0.0 or self.n2 <= 0.0:
            raise ValueError(f"profile {self.profile_id!r} requires positive n1 and n2")
        if self.te_thickness < 0.0:
            raise ValueError(f"profile {self.profile_id!r} te_thickness must be non-negative")
        x0, x1 = self.x_tc_window
        if not (0.0 <= x0 < x1 <= 1.0):
            raise ValueError(
                f"profile {self.profile_id!r} x_tc_window must satisfy 0 <= x0 < x1 <= 1, "
                f"got {self.x_tc_window}"
            )

    def to_cst_profile(self) -> "CSTAirfoilProfileSpec":
        self.validate()
        return self

    def evaluate(
        self,
        sample_count: int = 241,
        tc_target: float | None = None,
        x_tmax: float | None = None,
    ) -> ProfileSample:
        self.validate()
        x = cosine_spacing(int(sample_count))
        airfoil = KulfanCSTAirfoil(
            degree=self.degree,
            n1=self.n1,
            n2=self.n2,
            x_tc_window=self.x_tc_window,
        )
        coeffs = np.concatenate(
            [
                np.asarray(self.upper_coeffs, dtype=float),
                np.asarray(self.lower_coeffs, dtype=float),
            ]
        )
        kwargs = {"te_thickness": float(self.te_thickness)}
        if tc_target is not None:
            kwargs["tc_target"] = float(tc_target)
        if x_tmax is not None:
            kwargs["x_tmax"] = float(x_tmax)
        yu, yl = airfoil.evaluate(x, coeffs, **kwargs)
        return ProfileSample(x=x, yu=yu, yl=yl)


@dataclass(frozen=True)
class ICSTAirfoilProfileSpec:
    profile_id: str
    degree: int
    upper_leading_edge_radius: float
    lower_leading_edge_radius: float
    upper_constraints: Tuple[ICSTConstraint, ...]
    lower_constraints: Tuple[ICSTConstraint, ...]
    n1: float = 0.5
    n2: float = 1.0
    te_thickness: float = 0.0
    x_tc_window: Tuple[float, float] = (0.15, 0.65)

    def validate(self) -> None:
        if not self.profile_id:
            raise ValueError("profile_id must be non-empty")
        if self.degree < 1:
            raise ValueError(f"degree must be >= 1, got {self.degree}")
        if not isclose(self.n1, 0.5, abs_tol=1e-12) or not isclose(self.n2, 1.0, abs_tol=1e-12):
            raise ValueError(
                "the current airfoil iCST solver supports the classical airfoil class exponents "
                "n1=0.5 and n2=1.0 only"
            )
        if self.upper_leading_edge_radius <= 0.0 or self.lower_leading_edge_radius <= 0.0:
            raise ValueError("leading-edge radii must be positive")
        if self.te_thickness < 0.0:
            raise ValueError(f"profile {self.profile_id!r} te_thickness must be non-negative")
        x0, x1 = self.x_tc_window
        if not (0.0 <= x0 < x1 <= 1.0):
            raise ValueError(
                f"profile {self.profile_id!r} x_tc_window must satisfy 0 <= x0 < x1 <= 1, "
                f"got {self.x_tc_window}"
            )
        if len(self.upper_constraints) != self.degree or len(self.lower_constraints) != self.degree:
            raise ValueError(
                f"iCST profile {self.profile_id!r} expects degree={self.degree} user constraints per side "
                "plus the leading-edge-radius condition, "
                f"got upper={len(self.upper_constraints)}, lower={len(self.lower_constraints)}"
            )
        for constraint in self.upper_constraints + self.lower_constraints:
            constraint.validate()

    def to_cst_profile(self) -> CSTAirfoilProfileSpec:
        self.validate()
        upper_signed = _solve_icst_surface_coeffs(
            degree=self.degree,
            leading_edge_radius=self.upper_leading_edge_radius,
            constraints=self.upper_constraints,
            te_y=0.5 * float(self.te_thickness),
            sign=+1.0,
            n1=self.n1,
            n2=self.n2,
        )
        lower_signed = _solve_icst_surface_coeffs(
            degree=self.degree,
            leading_edge_radius=self.lower_leading_edge_radius,
            constraints=self.lower_constraints,
            te_y=-0.5 * float(self.te_thickness),
            sign=-1.0,
            n1=self.n1,
            n2=self.n2,
        )
        return CSTAirfoilProfileSpec(
            profile_id=self.profile_id,
            degree=self.degree,
            upper_coeffs=tuple(float(value) for value in upper_signed),
            lower_coeffs=tuple(float(-value) for value in lower_signed),
            n1=self.n1,
            n2=self.n2,
            te_thickness=self.te_thickness,
            x_tc_window=self.x_tc_window,
        )

    def evaluate(
        self,
        sample_count: int = 241,
        tc_target: float | None = None,
        x_tmax: float | None = None,
    ) -> ProfileSample:
        return self.to_cst_profile().evaluate(
            sample_count=sample_count,
            tc_target=tc_target,
            x_tmax=x_tmax,
        )


@dataclass(frozen=True)
class IntuitiveAirfoilProfileSpec:
    profile_id: str
    leading_edge_radius: float
    max_thickness: float
    x_tmax: float
    max_camber: float = 0.0
    x_cmax: float = 0.4
    trailing_edge_wedge_angle_deg: float = 14.0
    trailing_edge_camber_angle_deg: float = 0.0
    aft_control_x: float = 0.70
    degree: int = 5
    n1: float = 0.5
    n2: float = 1.0
    te_thickness: float = 0.0
    x_tc_window: Tuple[float, float] = (0.15, 0.65)

    def validate(self) -> None:
        if not self.profile_id:
            raise ValueError("profile_id must be non-empty")
        if self.degree != 5:
            raise ValueError(
                "IntuitiveAirfoilProfileSpec currently supports degree=5 only, "
                f"got {self.degree}"
            )
        if not isclose(self.n1, 0.5, abs_tol=1e-12) or not isclose(self.n2, 1.0, abs_tol=1e-12):
            raise ValueError(
                "IntuitiveAirfoilProfileSpec currently supports the classical airfoil class exponents "
                "n1=0.5 and n2=1.0 only"
            )
        if self.leading_edge_radius <= 0.0:
            raise ValueError(f"leading_edge_radius must be positive, got {self.leading_edge_radius}")
        if self.max_thickness <= 0.0:
            raise ValueError(f"max_thickness must be positive, got {self.max_thickness}")
        if self.te_thickness < 0.0:
            raise ValueError(f"te_thickness must be non-negative, got {self.te_thickness}")
        if self.max_thickness <= self.te_thickness:
            raise ValueError(
                f"max_thickness must exceed te_thickness, got max_thickness={self.max_thickness} "
                f"and te_thickness={self.te_thickness}"
            )
        if not (0.0 < self.x_tmax < 1.0):
            raise ValueError(f"x_tmax must lie in (0, 1), got {self.x_tmax}")
        if not (0.0 < self.x_cmax < 1.0):
            raise ValueError(f"x_cmax must lie in (0, 1), got {self.x_cmax}")
        if not (max(self.x_tmax, self.x_cmax) < self.aft_control_x < 1.0):
            raise ValueError(
                "aft_control_x must lie strictly aft of both x_tmax and x_cmax and before the trailing edge, "
                f"got aft_control_x={self.aft_control_x}, x_tmax={self.x_tmax}, x_cmax={self.x_cmax}"
            )
        if self.trailing_edge_wedge_angle_deg < 0.0:
            raise ValueError(
                f"trailing_edge_wedge_angle_deg must be non-negative, got {self.trailing_edge_wedge_angle_deg}"
            )
        x0, x1 = self.x_tc_window
        if not (0.0 <= x0 < x1 <= 1.0):
            raise ValueError(
                f"profile {self.profile_id!r} x_tc_window must satisfy 0 <= x0 < x1 <= 1, "
                f"got {self.x_tc_window}"
            )

    def surface_te_angles_deg(self) -> tuple[float, float]:
        self.validate()
        upper = float(self.trailing_edge_camber_angle_deg) - 0.5 * float(self.trailing_edge_wedge_angle_deg)
        lower = float(self.trailing_edge_camber_angle_deg) + 0.5 * float(self.trailing_edge_wedge_angle_deg)
        return upper, lower

    def surface_te_slopes(self) -> tuple[float, float]:
        upper_angle_deg, lower_angle_deg = self.surface_te_angles_deg()
        return float(tan(radians(upper_angle_deg))), float(tan(radians(lower_angle_deg)))

    def to_icst_profile(self) -> ICSTAirfoilProfileSpec:
        self.validate()
        upper_te_slope, lower_te_slope = self.surface_te_slopes()
        camber_te_slope = 0.5 * (upper_te_slope + lower_te_slope)
        half_thickness_te_slope = 0.5 * (upper_te_slope - lower_te_slope)

        camber_t, camber_slope_t = _canonical_camber_value_slope(
            max_camber=self.max_camber,
            x_cmax=self.x_cmax,
            trailing_edge_camber_slope=camber_te_slope,
            psi=self.x_tmax,
        )
        camber_aft, camber_slope_aft = _canonical_camber_value_slope(
            max_camber=self.max_camber,
            x_cmax=self.x_cmax,
            trailing_edge_camber_slope=camber_te_slope,
            psi=self.aft_control_x,
        )

        half_thickness_t = 0.5 * float(self.max_thickness)
        half_thickness_aft, half_thickness_slope_aft = _cubic_hermite_value_slope(
            x0=float(self.x_tmax),
            y0=half_thickness_t,
            dy0=0.0,
            x1=1.0,
            y1=0.5 * float(self.te_thickness),
            dy1=half_thickness_te_slope,
            x=float(self.aft_control_x),
        )

        upper_constraints = (
            ICSTConstraint.point(float(self.x_tmax), camber_t + half_thickness_t, "upper_crest"),
            ICSTConstraint.slope(float(self.x_tmax), camber_slope_t, "upper_crest_slope"),
            ICSTConstraint.point(float(self.aft_control_x), camber_aft + half_thickness_aft, "upper_aft_point"),
            ICSTConstraint.slope(
                float(self.aft_control_x),
                camber_slope_aft + half_thickness_slope_aft,
                "upper_aft_slope",
            ),
            ICSTConstraint.slope(1.0, upper_te_slope, "upper_te_angle"),
        )
        lower_constraints = (
            ICSTConstraint.point(float(self.x_tmax), camber_t - half_thickness_t, "lower_crest"),
            ICSTConstraint.slope(float(self.x_tmax), camber_slope_t, "lower_crest_slope"),
            ICSTConstraint.point(float(self.aft_control_x), camber_aft - half_thickness_aft, "lower_aft_point"),
            ICSTConstraint.slope(
                float(self.aft_control_x),
                camber_slope_aft - half_thickness_slope_aft,
                "lower_aft_slope",
            ),
            ICSTConstraint.slope(1.0, lower_te_slope, "lower_te_angle"),
        )

        surface_radius = _surface_leading_edge_radius_from_airfoil_radius(self.leading_edge_radius)
        return ICSTAirfoilProfileSpec(
            profile_id=self.profile_id,
            degree=self.degree,
            upper_leading_edge_radius=surface_radius,
            lower_leading_edge_radius=surface_radius,
            upper_constraints=upper_constraints,
            lower_constraints=lower_constraints,
            n1=self.n1,
            n2=self.n2,
            te_thickness=self.te_thickness,
            x_tc_window=self.x_tc_window,
        )

    def to_cst_profile(self) -> CSTAirfoilProfileSpec:
        return self.to_icst_profile().to_cst_profile()

    def evaluate(
        self,
        sample_count: int = 241,
        tc_target: float | None = None,
        x_tmax: float | None = None,
    ) -> ProfileSample:
        return self.to_cst_profile().evaluate(
            sample_count=sample_count,
            tc_target=tc_target,
            x_tmax=x_tmax,
        )


AirfoilProfileSpec = CSTAirfoilProfileSpec | ICSTAirfoilProfileSpec | IntuitiveAirfoilProfileSpec


@dataclass(frozen=True)
class ProfileCatalog:
    profiles: Tuple[AirfoilProfileSpec, ...]

    def validate(self) -> None:
        if not self.profiles:
            raise ValueError("profile catalog must contain at least one profile")
        seen = set()
        for profile in self.profiles:
            profile.validate()
            if profile.profile_id in seen:
                raise ValueError(f"duplicate profile_id detected: {profile.profile_id!r}")
            seen.add(profile.profile_id)

    def as_dict(self) -> Dict[str, AirfoilProfileSpec]:
        self.validate()
        return {profile.profile_id: profile for profile in self.profiles}

    def get(self, profile_id: str) -> AirfoilProfileSpec:
        catalog = self.as_dict()
        if profile_id not in catalog:
            raise KeyError(f"unknown profile_id: {profile_id!r}")
        return catalog[profile_id]

    def resolve(self, profile_id: str) -> CSTAirfoilProfileSpec:
        return self.get(profile_id).to_cst_profile()
