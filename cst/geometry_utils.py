import numpy as np


def spar_height(airfoil, psi):
    """
    Compute dimensional spar height at psi.
    """
    z_u = airfoil.upper_surface(psi)
    z_l = airfoil.lower_surface(psi)
    return airfoil.chord * abs(z_u - z_l)


def tangent_vector(airfoil, psi, surface="upper"):
    dz_dpsi = airfoil.derivative(psi, surface)
    norm = np.sqrt(1 + dz_dpsi**2)
    return np.array([1 / norm, dz_dpsi / norm])
