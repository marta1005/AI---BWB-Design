import numpy as np
from cst.bernstein import bernstein, bernstein_derivative


class CSTAirfoil:
    """
    Class Shape Transformation airfoil representation.
    """

    def __init__(
        self,
        n: int,
        upper_coeffs: np.ndarray,
        lower_coeffs: np.ndarray,
        delta_z_te: float = 0.0,
        N1: float = 0.5,
        N2: float = 1.0,
        chord: float = 1.0,
    ):
        self.n = n
        self.upper_coeffs = upper_coeffs
        self.lower_coeffs = lower_coeffs
        self.delta_z_te = delta_z_te
        self.N1 = N1
        self.N2 = N2
        self.chord = chord

    # --------------------------------------------------
    # Class function
    # --------------------------------------------------
    def class_function(self, psi):
        return psi**self.N1 * (1 - psi)**self.N2

    def class_derivative(self, psi):
        return (
            self.N1 * psi**(self.N1 - 1) * (1 - psi)**self.N2
            - self.N2 * psi**self.N1 * (1 - psi)**(self.N2 - 1)
        )

    # --------------------------------------------------
    # Shape function
    # --------------------------------------------------
    def shape_function(self, psi, coeffs):
        S = np.zeros_like(psi)
        for i in range(self.n + 1):
            S += coeffs[i] * bernstein(self.n, i, psi)
        return S

    def shape_derivative(self, psi, coeffs):
        S = np.zeros_like(psi)
        for i in range(self.n + 1):
            S += coeffs[i] * bernstein_derivative(self.n, i, psi)
        return S

    # --------------------------------------------------
    # Surfaces
    # --------------------------------------------------
    def upper_surface(self, psi):
        C = self.class_function(psi)
        S = self.shape_function(psi, self.upper_coeffs)
        return C * S + psi * self.delta_z_te / 2

    def lower_surface(self, psi):
        C = self.class_function(psi)
        S = self.shape_function(psi, self.lower_coeffs)
        return -C * S - psi * self.delta_z_te / 2

    # --------------------------------------------------
    # Derivatives
    # --------------------------------------------------
    def derivative(self, psi, surface="upper"):
        if surface == "upper":
            coeffs = self.upper_coeffs
            sign = 1.0
        else:
            coeffs = self.lower_coeffs
            sign = -1.0

        C = self.class_function(psi)
        C_prime = self.class_derivative(psi)
        S = self.shape_function(psi, coeffs)
        S_prime = self.shape_derivative(psi, coeffs)

        return sign * (C_prime * S + C * S_prime) + self.delta_z_te / 2

    # --------------------------------------------------
    # Arc length
    # --------------------------------------------------
    def arc_length(self, surface="upper", n_points=500):
        eps = 1e-8
        psi = np.linspace(eps, 1 - eps, n_points)
        dz_dpsi = self.derivative(psi, surface=surface)
        integrand = np.sqrt(1 + dz_dpsi**2)
        return self.chord * np.trapezoid(integrand, psi)
    
    def derivative_dimensional(self, psi, surface="upper"):
        return self.derivative(psi, surface)
    
    # --------------------------------------------------
    # Radius
    # --------------------------------------------------    
    def set_thickness_le_radius(self, R_le):
        # current sum of A0
        A0_u = self.upper_coeffs[0]
        A0_l = self.lower_coeffs[0]

        current_sum = A0_u + A0_l

        target_sum = np.sqrt(2.0 * R_le)

        scale = target_sum / current_sum

        self.upper_coeffs[0] *= scale
        self.lower_coeffs[0] *= scale

    def get_thickness_le_radius(self):
        A0_u = self.upper_coeffs[0]
        A0_l = self.lower_coeffs[0]

        return 0.5 * (A0_u + A0_l)**2



