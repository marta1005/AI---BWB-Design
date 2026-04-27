from dataclasses import dataclass

import numpy as np
from scipy.interpolate import CubicSpline

from parametrization.shared.dependency_setup import load_pyspline_curve
from .specs import AnchoredSpanwiseLaw, SectionedBWBModelConfig, SpanwiseLawSpec
from .topology import SectionedBWBTopologySpec

try:
    from scipy.interpolate import PchipInterpolator

    HAS_PCHIP = True
except Exception:
    HAS_PCHIP = False


INTERPOLATION_DISPLAY_NAMES = {
    "linear": "Linear",
    "blended_linear": "Blended linear",
    "pchip": "PCHIP",
    "cubic": "CubicSpline(natural, C2)",
    "pyspline": "pySpline Curve(k=4)",
}


def interpolation_display_name(
    interpolation: str,
    interpolation_factory=None,
) -> str:
    if interpolation_factory is not None:
        return str(interpolation)
    return INTERPOLATION_DISPLAY_NAMES.get(str(interpolation), str(interpolation))


def quintic_c2_transition(
    y: float,
    y0: float,
    y1: float,
    x0: float,
    x1: float,
    dx0: float,
    dx1: float,
    ddx0: float = 0.0,
    ddx1: float = 0.0,
) -> float:
    length = max(y1 - y0, 1e-12)
    t = np.clip((y - y0) / length, 0.0, 1.0)
    rhs = np.array(
        [x0, dx0 * length, ddx0 * length * length, x1, dx1 * length, ddx1 * length * length],
        dtype=float,
    )
    matrix = np.array(
        [
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 2.0, 0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
            [0.0, 0.0, 2.0, 6.0, 12.0, 20.0],
        ],
        dtype=float,
    )
    coeffs = np.linalg.solve(matrix, rhs)
    powers = np.array([1.0, t, t * t, t**3, t**4, t**5], dtype=float)
    return float(coeffs @ powers)


def with_root_blend(
    interpolant,
    root_value: float,
    root_blend_y: float,
):
    root_blend_y = float(max(root_blend_y, 0.0))
    if root_blend_y <= 1.0e-12:
        return interpolant

    x1 = float(interpolant(root_blend_y))
    eps = max(root_blend_y * 1.0e-4, 1.0e-6)
    y_left = max(root_blend_y - eps, 0.0)
    y_right = root_blend_y + eps
    if y_right <= y_left:
        dx1 = 0.0
    else:
        dx1 = (float(interpolant(y_right)) - float(interpolant(y_left))) / (y_right - y_left)

    def wrapped(yy: float) -> float:
        yy = float(yy)
        if yy <= 0.0:
            return float(root_value)
        if yy >= root_blend_y:
            return float(interpolant(yy))
        return quintic_c2_transition(
            y=yy,
            y0=0.0,
            y1=root_blend_y,
            x0=float(root_value),
            x1=x1,
            dx0=0.0,
            dx1=float(dx1),
        )

    return wrapped


@dataclass
class ResolvedSpanwiseLaws:
    twist_deg: callable
    camber_delta: callable


def build_scalar_interpolant(
    y_sections: np.ndarray,
    values: np.ndarray,
    interpolation: str,
    linear_start_index: int | None = None,
    interpolation_factory=None,
):
    y_sections = np.asarray(y_sections, dtype=float)
    values = np.asarray(values, dtype=float)

    if linear_start_index is not None and values.size > 2:
        linear_start_index = int(linear_start_index)
        if not (1 <= linear_start_index < values.size - 1):
            raise ValueError(
                "linear_start_index must lie inside [1, point_count - 2], "
                f"got {linear_start_index} for {values.size} points"
            )
        prefix_interpolant = build_scalar_interpolant(
            y_sections[: linear_start_index + 1],
            values[: linear_start_index + 1],
            interpolation,
            linear_start_index=None,
            interpolation_factory=interpolation_factory,
        )
        y_linear = y_sections[linear_start_index:]
        values_linear = values[linear_start_index:]

        def hybrid_interpolant(yy: float) -> float:
            yy = float(yy)
            if yy <= float(y_sections[linear_start_index]):
                return float(prefix_interpolant(yy))
            return float(np.interp(yy, y_linear, values_linear))

        return hybrid_interpolant

    if interpolation_factory is not None:
        interpolant = interpolation_factory(np.asarray(y_sections, dtype=float), np.asarray(values, dtype=float))
        if not callable(interpolant):
            raise ValueError("interpolation_factory must return a callable interpolant")
        return lambda yy: float(interpolant(float(yy)))

    if interpolation == "linear" or values.size == 2:
        return lambda yy: float(np.interp(float(yy), y_sections, values))

    if interpolation == "blended_linear":
        y_sections = np.asarray(y_sections, dtype=float)
        values = np.asarray(values, dtype=float)
        base_linear = lambda yy: float(np.interp(float(yy), y_sections, values))

        def blended_interpolant(yy: float) -> float:
            yy = float(yy)
            blended = base_linear(yy)
            if y_sections.size < 3:
                return blended

            for idx in range(1, y_sections.size - 1):
                y_prev = float(y_sections[idx - 1])
                y_here = float(y_sections[idx])
                y_next = float(y_sections[idx + 1])
                v_prev = float(values[idx - 1])
                v_here = float(values[idx])
                v_next = float(values[idx + 1])

                left_span = max(y_here - y_prev, 1.0e-12)
                right_span = max(y_next - y_here, 1.0e-12)
                half_width = 0.07 * min(left_span, right_span)
                if half_width <= 1.0e-12:
                    continue
                y_left = y_here - half_width
                y_right = y_here + half_width
                if not (y_left <= yy <= y_right):
                    continue

                left_slope = (v_here - v_prev) / left_span
                right_slope = (v_next - v_here) / right_span
                left_line = v_here + left_slope * (yy - y_here)
                right_line = v_here + right_slope * (yy - y_here)
                t = (yy - y_left) / max(y_right - y_left, 1.0e-12)
                w = t * t * t * (10.0 + t * (-15.0 + 6.0 * t))
                return float((1.0 - w) * left_line + w * right_line)

            first_span = float(y_sections[1] - y_sections[0])
            root_blend_y = max(0.0, 0.45 * first_span)
            if yy <= root_blend_y and root_blend_y > 1.0e-12:
                root_value = float(values[0])
                root_slope = float((values[1] - values[0]) / max(first_span, 1.0e-12))
                join_value = root_value + root_slope * root_blend_y
                t = np.clip(yy / root_blend_y, 0.0, 1.0)
                h00 = 1.0 - 10.0 * t**3 + 15.0 * t**4 - 6.0 * t**5
                h01 = 10.0 * t**3 - 15.0 * t**4 + 6.0 * t**5
                h11 = -4.0 * t**3 + 7.0 * t**4 - 3.0 * t**5
                return float(
                    h00 * root_value
                    + h01 * join_value
                    + h11 * root_blend_y * root_slope
                )

            return blended

        return blended_interpolant

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
    anchor_y = topology.anchor_y_array[np.asarray(law.section_indices, dtype=int)]
    values = np.asarray(law.values, dtype=float)
    interpolant = build_scalar_interpolant(
        anchor_y,
        values,
        law.interpolation,
        linear_start_index=law.linear_start_index,
        interpolation_factory=law.interpolation_factory,
    )
    if law.root_blend_y is not None:
        interpolant = with_root_blend(
            interpolant,
            root_value=float(values[0]),
            root_blend_y=float(law.root_blend_y),
        )
    return interpolant


def resolve_spanwise_laws(config: SectionedBWBModelConfig) -> ResolvedSpanwiseLaws:
    return ResolvedSpanwiseLaws(
        twist_deg=(
            config.spanwise.twist_callable
            if config.spanwise.twist_callable is not None
            else build_anchored_interpolant(config.topology, config.spanwise.twist_deg)
        ),
        camber_delta=build_anchored_interpolant(config.topology, config.spanwise.camber_delta),
    )
def vertical_offsets(
    topology: SectionedBWBTopologySpec,
    spanwise: SpanwiseLawSpec,
    span_stations: np.ndarray,
) -> np.ndarray:
    span_z = np.asarray(span_stations, dtype=float)
    if spanwise.vertical_offset_callable is not None:
        return np.asarray([spanwise.vertical_offset_callable(float(yy)) for yy in span_z], dtype=float)
    if spanwise.vertical_offset_z is not None:
        interpolant = build_anchored_interpolant(topology, spanwise.vertical_offset_z)
        return np.asarray([interpolant(float(yy)) for yy in span_z], dtype=float)
    return float(spanwise.vertical_offset_m) + np.tan(np.deg2rad(spanwise.dihedral_deg)) * span_z
