import numpy as np
from cst.math_utils import binomial_coefficient


def bernstein(n: int, i: int, psi: np.ndarray) -> np.ndarray:
    """
    Bernstein polynomial B_i^n(psi).
    """
    return binomial_coefficient(n, i) * psi**i * (1 - psi)**(n - i)


def bernstein_derivative(n: int, i: int, psi: np.ndarray) -> np.ndarray:
    """
    Derivative of Bernstein polynomial using recursive relation.
    d/dψ B_i^n = n (B_{i-1}^{n-1} - B_i^{n-1})
    """
    if n == 0:
        return np.zeros_like(psi)

    if i == 0:
        return -n * bernstein(n - 1, 0, psi)

    if i == n:
        return n * bernstein(n - 1, n - 1, psi)

    return n * (
        bernstein(n - 1, i - 1, psi)
        - bernstein(n - 1, i, psi)
    )
