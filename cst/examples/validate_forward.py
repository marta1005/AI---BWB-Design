import numpy as np
import matplotlib.pyplot as plt

from cst.cst_airfoil import CSTAirfoil
from cst.structurally_consistent import StructurallyConsistentMorpher


# --------------------------------------------------
# Utility: approximate NACA0008 with CST
# --------------------------------------------------

def naca_0008_thickness(x):
    t = 0.08
    return (
        5 * t * (
            0.2969 * np.sqrt(x)
            - 0.1260 * x
            - 0.3516 * x**2
            + 0.2843 * x**3
            - 0.1015 * x**4
        )
    )


def approximate_naca0008_with_cst(n=8):

    psi = np.linspace(1e-6, 1-1e-6, 300)
    thickness = naca_0008_thickness(psi)

    from cst.bernstein import bernstein

    B = np.zeros((len(psi), n + 1))
    for i in range(n + 1):
        B[:, i] = bernstein(n, i, psi)

    C = psi**0.5 * (1 - psi)

    A = np.linalg.lstsq(C[:, None] * B, thickness, rcond=None)[0]

    return A


# --------------------------------------------------
# 1️ Parent airfoil
# --------------------------------------------------

n = 8

upper_coeffs = approximate_naca0008_with_cst(n)
lower_coeffs = upper_coeffs.copy()

parent = CSTAirfoil(
    n=n,
    upper_coeffs=upper_coeffs.copy(),
    lower_coeffs=lower_coeffs.copy(),
    chord=1.0,
)

# --------------------------------------------------
# 2️ Spars
# --------------------------------------------------

spar_locations = np.linspace(0.1, 0.75, n)

morpher = StructurallyConsistentMorpher(parent, spar_locations)

# --------------------------------------------------
# 3️ Create Child (initially identical)
# --------------------------------------------------

child = CSTAirfoil(
    n=n,
    upper_coeffs=upper_coeffs.copy(),
    lower_coeffs=lower_coeffs.copy(),
    chord=1.0,
)

# --------------------------------------------------
# 4️ FORWARD
# --------------------------------------------------

child_morphed = morpher.forward_morph(child)

# --------------------------------------------------
# 5️ BACKWARD
# --------------------------------------------------

parent_recovered = morpher.backward_morph(child_morphed)

# --------------------------------------------------
# 6️ Plot
# --------------------------------------------------

psi = np.linspace(1e-6, 1-1e-6, 400)

plt.figure(figsize=(10,5))

# Original parent
plt.plot(psi, parent.upper_surface(psi), label="Parent upper")
plt.plot(psi, parent.lower_surface(psi), label="Parent lower")

# Child
plt.plot(psi, child_morphed.upper_surface(psi), "--", label="Child upper")
plt.plot(psi, child_morphed.lower_surface(psi), "--", label="Child lower")

# Recovered parent
plt.plot(psi, parent_recovered.upper_surface(psi), ":", label="Recovered upper")
plt.plot(psi, parent_recovered.lower_surface(psi), ":", label="Recovered lower")

# Spars (vertical lines)
for psi_s in spar_locations:

    # Parent
    z_l = parent.lower_surface(psi_s)
    z_u = parent.upper_surface(psi_s)
    plt.plot([psi_s, psi_s], [z_l, z_u], color="black", linewidth=2)

    # posición child (puede cambiar si cambia chord)
    psi_child = psi_s * (parent.chord / child_morphed.chord)

    z_l_c = child_morphed.lower_surface(psi_child)
    z_u_c = child_morphed.upper_surface(psi_child)

    plt.plot(
        [psi_child, psi_child],
        [z_l_c, z_u_c],
        color="red",
        linewidth=2,
        linestyle="--"
    )

plt.axis("equal")
plt.legend()
plt.title("Forward + Backward CST Morph Validation")
plt.show()
