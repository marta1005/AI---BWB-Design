from dataclasses import dataclass

import numpy as np
from scipy.interpolate import CubicSpline

from .dependency_setup import load_pyspline_curve
from .specs import AnchoredSpanwiseLaw, SectionedBWBModelConfig, SpanwiseLawSpec
from .topology import SectionedBWBTopologySpec

try:
    from scipy.interpolate import PchipInterpolator

    HAS_PCHIP = True
except Exception:
    HAS_PCHIP = False


@dataclass
class ResolvedSpanwiseLaws:
    twist_deg: callable
    camber_delta: callable


def build_scalar_interpolant(
    y_sections: np.ndarray,
    values: np.ndarray,
    interpolation: str,
):
    y_sections = np.asarray(y_sections, dtype=float)
    values = np.asarray(values, dtype=float)

    if interpolation == "linear" or values.size == 2:
        return lambda yy: float(np.interp(float(yy), y_sections, values))

    if interpolation == "pchip":
        if not HAS_PCHIP:
            raise ValueError("interpolation='pchip' requested, but PchipInterpolator is unavailable")
        interpolant = PchipInterpolator(y_sections, values)
        return lambda yy: float(interpolant(float(yy)))

    if interpolation == "pyspline":
        Curve = load_pyspline_curve(rebuild_if_needed=True)
        order = min(4, int(values.size))
        curve = Curve(x=values, s=y_sections, k=order)
        return lambda yy: float(np.asarray(curve(float(yy))).reshape(-1)[0])

    interpolant = CubicSpline(y_sections, values, bc_type="natural")
    return lambda yy: float(interpolant(float(yy)))


def build_anchored_interpolant(
    topology: SectionedBWBTopologySpec,
    law: AnchoredSpanwiseLaw,
):
    y_sections = topology.y_sections_array[np.asarray(law.section_indices, dtype=int)]
    values = np.asarray(law.values, dtype=float)
    return build_scalar_interpolant(y_sections, values, law.interpolation)


def resolve_spanwise_laws(config: SectionedBWBModelConfig) -> ResolvedSpanwiseLaws:
    return ResolvedSpanwiseLaws(
        twist_deg=build_anchored_interpolant(config.topology, config.spanwise.twist_deg),
        camber_delta=build_anchored_interpolant(config.topology, config.spanwise.camber_delta),
    )
def vertical_offsets(
    topology: SectionedBWBTopologySpec,
    spanwise: SpanwiseLawSpec,
    span_stations: np.ndarray,
) -> np.ndarray:
    span_z = np.asarray(span_stations, dtype=float)
    return np.tan(np.deg2rad(spanwise.dihedral_deg)) * span_z
