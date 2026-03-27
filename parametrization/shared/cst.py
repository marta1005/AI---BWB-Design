from math import comb
from typing import Optional, Tuple

import numpy as np


def cosine_spacing(n: int) -> np.ndarray:
    beta = np.linspace(0.0, np.pi, n)
    return 0.5 * (1.0 - np.cos(beta))


def bernstein_matrix(x: np.ndarray, degree: int) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    basis = np.zeros((x.size, degree + 1), dtype=float)
    for idx in range(degree + 1):
        basis[:, idx] = comb(degree, idx) * (x**idx) * ((1.0 - x) ** (degree - idx))
    return basis


def fit_kulfan_airfoil_coefficients(
    x: np.ndarray,
    y_upper: np.ndarray,
    y_lower: np.ndarray,
    degree: int,
    n1: float = 0.5,
    n2: float = 1.0,
    te_thickness: float = 0.0,
    regularization: float = 1.0e-6,
    smoothness_weight: float = 0.0,
    shared_leading_edge: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(x, dtype=float)
    y_upper = np.asarray(y_upper, dtype=float)
    y_lower = np.asarray(y_lower, dtype=float)
    if x.ndim != 1 or y_upper.ndim != 1 or y_lower.ndim != 1:
        raise ValueError("x, y_upper and y_lower must be 1D arrays")
    if not (x.size == y_upper.size == y_lower.size):
        raise ValueError("x, y_upper and y_lower must have the same length")

    interior = (x > 1.0e-4) & (x < 0.99)
    if np.count_nonzero(interior) < degree + 1:
        raise ValueError(
            f"not enough interior samples to fit degree {degree} CST coefficients; "
            f"need at least {degree + 1}, got {np.count_nonzero(interior)}"
        )

    x_fit = x[interior]
    class_fun = (x_fit**float(n1)) * ((1.0 - x_fit) ** float(n2))
    basis = bernstein_matrix(x_fit, degree)

    # Kulfan evaluation uses a symmetric trailing-edge thickness term:
    # yu = yu_zero + 0.5 * x * te_thickness
    # yl = yl_zero - 0.5 * x * te_thickness
    upper_rhs = (y_upper[interior] - 0.5 * x_fit * float(te_thickness)) / class_fun
    lower_rhs = -(y_lower[interior] + 0.5 * x_fit * float(te_thickness)) / class_fun

    eye = float(regularization) * np.eye(degree + 1)
    normal_matrix = basis.T @ basis + eye
    smoothness = float(smoothness_weight)
    if smoothness > 0.0 and degree >= 2:
        d2 = np.zeros((degree - 1, degree + 1), dtype=float)
        for row in range(degree - 1):
            d2[row, row] = 1.0
            d2[row, row + 1] = -2.0
            d2[row, row + 2] = 1.0
        normal_matrix = normal_matrix + smoothness * (d2.T @ d2)
    upper = np.linalg.solve(normal_matrix, basis.T @ upper_rhs)
    lower = np.linalg.solve(normal_matrix, basis.T @ lower_rhs)

    if shared_leading_edge and upper.size > 0:
        a0_shared = 0.5 * (float(upper[0]) + float(lower[0]))
        upper[0] = a0_shared
        lower[0] = a0_shared

    return upper.astype(float), lower.astype(float)


class CST:
    def __init__(self, coeffs, n1: float = 0.5, n2: float = 1.0, delta_te: float = 0.0):
        self.coeffs = np.asarray(coeffs, dtype=float)
        self.degree = len(self.coeffs) - 1
        self.n1 = float(n1)
        self.n2 = float(n2)
        self.delta_te = float(delta_te)

    def class_function(self, x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        return (x**self.n1) * ((1.0 - x) ** self.n2)

    def evaluate(self, x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        shape = bernstein_matrix(x, self.degree) @ self.coeffs
        return self.class_function(x) * shape + x * self.delta_te


class KulfanCSTAirfoil:
    def __init__(
        self,
        degree: int,
        n1: float = 0.5,
        n2: float = 1.0,
        x_tc_window: Tuple[float, float] = (0.15, 0.65),
        shared_leading_edge: bool = False,
    ):
        self.degree = int(degree)
        self.n1 = float(n1)
        self.n2 = float(n2)
        self.x_tc_window = (float(x_tc_window[0]), float(x_tc_window[1]))
        self.shared_leading_edge = bool(shared_leading_edge)

    @property
    def ncoeff(self) -> int:
        return self.degree + 1

    @property
    def total_coeff_count(self) -> int:
        return 2 * self.ncoeff

    def unpack(self, coeffs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        values = np.asarray(coeffs, dtype=float).ravel()
        expected = self.total_coeff_count
        if values.size != expected:
            raise ValueError(f"CST coeff vector must have size {expected}, received {values.size}")
        split = self.ncoeff
        upper = values[:split].copy()
        lower = values[split:].copy()
        if self.shared_leading_edge and upper.size > 0:
            shared_a0 = 0.5 * (float(upper[0]) + float(lower[0]))
            upper[0] = shared_a0
            lower[0] = shared_a0
        return upper, lower

    def _evaluate_zero_te(
        self,
        x: np.ndarray,
        coeffs: np.ndarray,
        gamma: float = 1.0,
    ) -> tuple[np.ndarray, np.ndarray]:
        upper, lower = self.unpack(coeffs)
        x_eval = np.clip(np.asarray(x, dtype=float), 0.0, 1.0) ** float(gamma)
        yu = CST(upper, n1=self.n1, n2=self.n2, delta_te=0.0).evaluate(x_eval)
        yl = -CST(lower, n1=self.n1, n2=self.n2, delta_te=0.0).evaluate(x_eval)
        return yu, yl

    def _tc_mask(self, x: np.ndarray) -> np.ndarray:
        return (x >= self.x_tc_window[0]) & (x <= self.x_tc_window[1])

    def _thickness_scale(
        self,
        x: np.ndarray,
        thickness_shape: np.ndarray,
        te_thickness: float,
        tc_target: Optional[float],
    ) -> float:
        if tc_target is None:
            return 1.0
        if tc_target <= 0.0:
            raise ValueError(f"tc_target must be positive, got {tc_target}")

        x = np.asarray(x, dtype=float)
        thickness_shape = np.asarray(thickness_shape, dtype=float)
        mask = self._tc_mask(x)
        if not np.any(mask):
            raise ValueError(f"x_tc_window={self.x_tc_window} leaves no valid samples")

        def max_tc_for(scale: float) -> float:
            thickness = scale * thickness_shape + x * float(te_thickness)
            return float(np.max(thickness[mask]))

        current_tc = max_tc_for(1.0)
        if current_tc <= 0.0:
            raise ValueError(f"invalid raw CST thickness {current_tc} for tc scaling")
        if abs(current_tc - tc_target) <= 1e-10:
            return 1.0

        low = 0.0
        high = max(2.0, 2.0 * tc_target / current_tc)
        while max_tc_for(high) < tc_target:
            high *= 2.0
            if high > 1e6:
                raise ValueError("could not bracket tc scaling factor")

        for _ in range(60):
            mid = 0.5 * (low + high)
            if max_tc_for(mid) < tc_target:
                low = mid
            else:
                high = mid
        return 0.5 * (low + high)

    def _peak_x_for_gamma(
        self,
        coeffs: np.ndarray,
        te_thickness: float,
        tc_target: Optional[float],
        gamma: float,
        sample_count: int = 401,
    ) -> float:
        x_dense = cosine_spacing(sample_count)
        yu_zero, yl_zero = self._evaluate_zero_te(x_dense, coeffs, gamma=gamma)
        thickness_shape = yu_zero - yl_zero
        scale = self._thickness_scale(x_dense, thickness_shape, te_thickness, tc_target)
        thickness = scale * thickness_shape + x_dense * float(te_thickness)
        mask = self._tc_mask(x_dense)
        x_window = x_dense[mask]
        thickness_window = thickness[mask]
        peak_idx = int(np.argmax(thickness_window))
        return float(x_window[peak_idx])

    def _solve_gamma_for_xtmax(
        self,
        coeffs: np.ndarray,
        te_thickness: float,
        tc_target: Optional[float],
        x_tmax: Optional[float],
    ) -> float:
        if x_tmax is None:
            return 1.0

        lo_target = max(self.x_tc_window[0], 1e-4)
        hi_target = min(self.x_tc_window[1], 1.0 - 1e-4)
        target = float(np.clip(x_tmax, lo_target, hi_target))
        base_peak = self._peak_x_for_gamma(coeffs, te_thickness, tc_target, gamma=1.0)
        if abs(base_peak - target) <= 5e-4:
            return 1.0

        # Monotone chordwise warp x_eval = x**gamma.
        # gamma > 1 moves the peak aft, gamma < 1 moves it forward.
        if target > base_peak:
            low = 1.0
            high = min(8.0, max(1.05, np.log(max(base_peak, 1e-6)) / np.log(max(target, 1e-6))))
            if high <= low:
                high = 1.25
            f_low = base_peak - target
            f_high = self._peak_x_for_gamma(coeffs, te_thickness, tc_target, gamma=high) - target
            while f_high < 0.0 and high < 8.0:
                high *= 1.5
                f_high = self._peak_x_for_gamma(coeffs, te_thickness, tc_target, gamma=high) - target
        else:
            high = 1.0
            low = max(0.12, np.log(max(base_peak, 1e-6)) / np.log(max(target, 1e-6)))
            if low >= high:
                low = 0.8
            f_high = base_peak - target
            f_low = self._peak_x_for_gamma(coeffs, te_thickness, tc_target, gamma=low) - target
            while f_low > 0.0 and low > 0.12:
                low = max(0.12, low / 1.5)
                f_low = self._peak_x_for_gamma(coeffs, te_thickness, tc_target, gamma=low) - target

        if f_low == 0.0:
            return low
        if f_high == 0.0:
            return high
        if f_low * f_high > 0.0:
            candidates = np.geomspace(0.12, 8.0, 31)
            peaks = np.array(
                [self._peak_x_for_gamma(coeffs, te_thickness, tc_target, gamma=value) for value in candidates],
                dtype=float,
            )
            best_idx = int(np.argmin(np.abs(peaks - target)))
            return float(candidates[best_idx])

        for _ in range(24):
            mid = 0.5 * (low + high)
            f_mid = self._peak_x_for_gamma(coeffs, te_thickness, tc_target, gamma=mid) - target
            if abs(f_mid) <= 5e-4:
                return mid
            if f_low * f_mid <= 0.0:
                high = mid
                f_high = f_mid
            else:
                low = mid
                f_low = f_mid
        return 0.5 * (low + high)

    def evaluate(
        self,
        x: np.ndarray,
        coeffs: np.ndarray,
        te_thickness: float = 0.0,
        tc_target: Optional[float] = None,
        x_tmax: Optional[float] = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        x = np.clip(np.asarray(x, dtype=float), 0.0, 1.0)
        gamma = self._solve_gamma_for_xtmax(coeffs, te_thickness, tc_target, x_tmax)
        yu_zero, yl_zero = self._evaluate_zero_te(x, coeffs, gamma=gamma)
        camber = 0.5 * (yu_zero + yl_zero)
        thickness_shape = yu_zero - yl_zero
        scale = self._thickness_scale(x, thickness_shape, te_thickness, tc_target)
        thickness = scale * thickness_shape + x * float(te_thickness)
        yu = camber + 0.5 * thickness
        yl = camber - 0.5 * thickness
        return yu, yl
