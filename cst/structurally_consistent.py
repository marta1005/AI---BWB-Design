import numpy as np
from cst.bernstein import bernstein


class StructurallyConsistentMorpher:
    """
    Structurally consistent CST morphing
    Lower surface = passive
    Upper surface = active
    """

    def __init__(self, parent_airfoil, spar_locations):
        self.parent = parent_airfoil
        self.spar_locations = np.array(spar_locations)

    # --------------------------------------------------
    # Parent spar heights (dimensional)
    # --------------------------------------------------
    def compute_parent_spar_heights(self):

        heights = []

        for psi in self.spar_locations:
            z_u = self.parent.upper_surface(psi)
            z_l = self.parent.lower_surface(psi)
            h = self.parent.chord * (z_u - z_l)
            heights.append(h)

        return np.array(heights)

    # --------------------------------------------------
    # Compute chord from passive surface length conservation
    # --------------------------------------------------
    def compute_child_chord(self, child_airfoil, tol=1e-10, max_iter=50):

        L_parent = self.parent.arc_length(surface="lower")

        c = self.parent.chord

        for _ in range(max_iter):

            child_airfoil.chord = c
            L_child = child_airfoil.arc_length(surface="lower")

            f = c * L_child - self.parent.chord * L_parent

            eps = 1e-6
            child_airfoil.chord = c + eps
            L_eps = child_airfoil.arc_length(surface="lower")
            f_eps = (c + eps) * L_eps - self.parent.chord * L_parent

            df = (f_eps - f) / eps

            c_new = c - f / df

            if abs(c_new - c) < tol:
                return c_new

            c = c_new

        raise RuntimeError("Chord solver did not converge.")

    # --------------------------------------------------
    # Forward morph
    # --------------------------------------------------
    def forward_morph(self, child_airfoil):

        n = child_airfoil.n

        # Step 1: update chord
        child_airfoil.chord = self.compute_child_chord(child_airfoil)

        c_p = self.parent.chord
        c_c = child_airfoil.chord

        # Step 2: update TE thickness
        child_airfoil.delta_z_te = (
            self.parent.delta_z_te * c_p / c_c
        )

        # Step 3: compute spar heights (parent)
        h_parent = self.compute_parent_spar_heights()

        # Scale spar heights to child
        h_child = (c_p / c_c) * h_parent

        # -------------------------------------------
        # FOR TESTING deformation:
        # Uncomment this line if you want deformation
        h_child *= 1.5
        # -------------------------------------------

        # Step 4: build linear system
        A0 = child_airfoil.upper_coeffs[0]

        F = []
        f_vec = []

        for j, psi in enumerate(self.spar_locations):

            # New psi in child due to chord scaling
            psi_c = psi * (c_p / c_c)

            C = child_airfoil.class_function(psi_c)

            z_lower = child_airfoil.lower_surface(psi_c)

            rhs = (
                z_lower
                + h_child[j] / c_c
                - psi_c * child_airfoil.delta_z_te / 2
            ) / C - A0 * (1 - psi_c) ** n

            row = []
            for i in range(1, n + 1):
                row.append(bernstein(n, i, psi_c))

            F.append(row)
            f_vec.append(rhs)

        F = np.array(F)
        f_vec = np.array(f_vec)

        A_active = np.linalg.lstsq(F, f_vec, rcond=None)[0]

        new_coeffs = np.zeros(n + 1)
        new_coeffs[0] = A0
        new_coeffs[1:] = A_active

        child_airfoil.upper_coeffs = new_coeffs

        return child_airfoil

    # --------------------------------------------------
    # Backward morph
    # --------------------------------------------------
    def backward_morph(self, child_airfoil):

        n = child_airfoil.n

        # Clone initial guess
        parent_guess = child_airfoil

        c_c = child_airfoil.chord

        # Step 1: passive length conservation
        L_child = child_airfoil.arc_length(surface="lower")

        c_p = c_c  # initial guess

        for _ in range(50):

            parent_guess.chord = c_p
            L_parent = parent_guess.arc_length(surface="lower")

            f = c_p * L_parent - c_c * L_child

            eps = 1e-6
            parent_guess.chord = c_p + eps
            L_eps = parent_guess.arc_length(surface="lower")
            f_eps = (c_p + eps) * L_eps - c_c * L_child

            df = (f_eps - f) / eps

            c_new = c_p - f / df

            if abs(c_new - c_p) < 1e-10:
                break

            c_p = c_new

        parent_guess.chord = c_p

        # Step 2: update TE thickness
        parent_guess.delta_z_te = (
            child_airfoil.delta_z_te * c_c / c_p
        )

        # Step 3: spar heights from child
        h_child = []

        for psi in self.spar_locations:
            z_u = child_airfoil.upper_surface(psi)
            z_l = child_airfoil.lower_surface(psi)
            h_child.append(c_c * (z_u - z_l))

        h_child = np.array(h_child)

        # Scale heights back to parent
        h_parent = (c_c / c_p) * h_child

        # Step 4: solve linear system for parent
        A0 = parent_guess.upper_coeffs[0]

        F = []
        f_vec = []

        for j, psi in enumerate(self.spar_locations):

            psi_p = psi * (c_c / c_p)

            C = parent_guess.class_function(psi_p)

            z_lower = parent_guess.lower_surface(psi_p)

            rhs = (
                z_lower
                + h_parent[j] / c_p
                - psi_p * parent_guess.delta_z_te / 2
            ) / C - A0 * (1 - psi_p) ** n

            row = []
            for i in range(1, n + 1):
                row.append(bernstein(n, i, psi_p))

            F.append(row)
            f_vec.append(rhs)

        F = np.array(F)
        f_vec = np.array(f_vec)

        lambda_reg = 1e-6
        A_active = np.linalg.solve(
            F.T @ F + lambda_reg * np.eye(F.shape[1]),
            F.T @ f_vec
        )


        new_coeffs = np.zeros(n + 1)
        new_coeffs[0] = A0
        new_coeffs[1:] = A_active

        parent_guess.upper_coeffs = new_coeffs

        return parent_guess
