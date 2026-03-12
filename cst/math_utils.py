import numpy as np


def binomial_coefficient(n: int, k: int) -> float:
    """
    Compute binomial coefficient C(n, k).
    """
    from math import comb
    return comb(n, k)


def fixed_point_iteration(func, x0, tol=1e-8, max_iter=100):
    """
    Solve x = func(x) using fixed-point iteration.
    """
    x = x0
    for _ in range(max_iter):
        x_new = func(x)
        if abs(x_new - x) < tol:
            return x_new
        x = x_new
    raise RuntimeError("Fixed point iteration did not converge.")
