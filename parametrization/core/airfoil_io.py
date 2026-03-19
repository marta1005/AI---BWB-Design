import numpy as np


def write_airfoil_dat(
    path: str,
    x: np.ndarray,
    yu: np.ndarray,
    yl: np.ndarray,
    name: str = "CST_AIRFOIL",
) -> None:
    x = np.asarray(x, float)
    yu = np.asarray(yu, float)
    yl = np.asarray(yl, float)

    xu = x[::-1]
    zu = yu[::-1]
    xl = x[1:]
    zl = yl[1:]

    pts = np.vstack([np.column_stack([xu, zu]), np.column_stack([xl, zl])])

    with open(path, "w", encoding="utf-8") as stream:
        stream.write(f"{name}\n")
        for xx, zz in pts:
            stream.write(f"{xx:.10f} {zz:.10f}\n")
