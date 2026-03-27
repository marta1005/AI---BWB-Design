from dataclasses import dataclass
from typing import Callable, List

import numpy as np

from parametrization.shared.cst import KulfanCSTAirfoil, cosine_spacing
from .specs import SectionedBWBModelConfig
from .spanwise_laws import (
    ResolvedSpanwiseLaws,
    build_scalar_interpolant,
)


@dataclass
class SectionShapeParameters:
    coeffs: np.ndarray
    tc_max: float
    x_tmax: float
    te_thickness: float


@dataclass
class SectionGeometryMetrics:
    max_tc: float
    x_tmax: float
    te_thickness: float
    min_inner_tc: float
    min_inner_tc_xc: float
    max_camber: float


@dataclass
class SectionModel:
    shape: KulfanCSTAirfoil
    x_air: np.ndarray
    y_sections: np.ndarray
    coeff_sections: np.ndarray
    tc_sections: np.ndarray
    x_tmax_sections: np.ndarray
    te_sections: np.ndarray
    mask_tc_window: np.ndarray
    mask_valid_window: np.ndarray
    coeff_interpolants: List[Callable[[float], float]]
    tc_interpolant: Callable[[float], float]
    x_tmax_interpolant: Callable[[float], float]
    camber_interpolant: Callable[[float], float]
    te_interpolant: Callable[[float], float]
    camber_mode_center: float
    camber_mode_width: float
    profile_generation_mode: str
    interpolation_name: str

    def params_at_y(self, y: float) -> SectionShapeParameters:
        coeffs = np.array([interp(y) for interp in self.coeff_interpolants], dtype=float)
        coeffs = apply_camber_delta(
            coeffs,
            camber_delta=float(self.camber_interpolant(y)),
            degree=self.shape.degree,
            center=self.camber_mode_center,
            width=self.camber_mode_width,
        )
        return SectionShapeParameters(
            coeffs=coeffs,
            tc_max=float(self.tc_interpolant(y)),
            x_tmax=float(self.x_tmax_interpolant(y)),
            te_thickness=float(self.te_interpolant(y)),
        )

    def coordinates_at_y(self, y: float) -> tuple[np.ndarray, np.ndarray, SectionShapeParameters]:
        params = self.params_at_y(float(y))
        evaluate_kwargs = {
            "te_thickness": params.te_thickness,
        }
        if self.profile_generation_mode == "enforce_targets":
            evaluate_kwargs["tc_target"] = params.tc_max
            evaluate_kwargs["x_tmax"] = params.x_tmax
        yu, yl = self.shape.evaluate(self.x_air, params.coeffs, **evaluate_kwargs)
        return yu, yl, params

    def geometry_metrics_at_y(
        self,
        y: float,
    ) -> tuple[SectionGeometryMetrics, SectionShapeParameters]:
        yu, yl, params = self.coordinates_at_y(float(y))
        thickness = yu - yl
        camber = 0.5 * (yu + yl)

        tc_window_x = self.x_air[self.mask_tc_window]
        tc_window_thickness = thickness[self.mask_tc_window]
        tc_idx = int(np.argmax(tc_window_thickness))

        valid_window_x = self.x_air[self.mask_valid_window]
        valid_window_thickness = thickness[self.mask_valid_window]
        valid_idx = int(np.argmin(valid_window_thickness))

        metrics = SectionGeometryMetrics(
            max_tc=float(tc_window_thickness[tc_idx]),
            x_tmax=float(tc_window_x[tc_idx]),
            te_thickness=float(thickness[-1]),
            min_inner_tc=float(valid_window_thickness[valid_idx]),
            min_inner_tc_xc=float(valid_window_x[valid_idx]),
            max_camber=float(np.max(np.abs(camber))),
        )
        return metrics, params


def camber_mode_vector(
    degree: int,
    center: float,
    width: float,
) -> np.ndarray:
    modes = np.arange(degree + 1, dtype=float)
    weights = np.exp(-0.5 * ((modes - center) / max(width, 1e-6)) ** 2)
    weights[0] = 0.0
    peak = max(float(weights.max()), 1e-12)
    return weights / peak


def apply_camber_delta(
    coeffs: np.ndarray,
    camber_delta: float,
    degree: int,
    center: float,
    width: float,
) -> np.ndarray:
    values = np.array(coeffs, dtype=float)
    if abs(camber_delta) < 1e-15:
        return values

    ncoeff = degree + 1
    weights = camber_mode_vector(degree, center=center, width=width)
    upper_slice = slice(0, ncoeff)
    lower_slice = slice(ncoeff, 2 * ncoeff)
    values[upper_slice] += camber_delta * weights
    values[lower_slice] -= camber_delta * weights
    return values


def build_positive_interpolant(
    y_sections: np.ndarray,
    values: np.ndarray,
    interpolation: str,
    floor: float = 1e-8,
):
    values = np.maximum(np.asarray(values, dtype=float), floor)
    log_interpolant = build_scalar_interpolant(y_sections, np.log(values), interpolation)
    return lambda yy: float(np.exp(log_interpolant(float(yy))))


def build_section_model(
    config: SectionedBWBModelConfig,
    laws: ResolvedSpanwiseLaws,
) -> SectionModel:
    topology = config.topology
    sections = config.sections
    sampling = config.sampling

    shape = KulfanCSTAirfoil(
        degree=sections.cst_degree,
        n1=sections.n1,
        n2=sections.n2,
        x_tc_window=sections.x_tc_window,
        shared_leading_edge=sections.shared_leading_edge,
    )
    x_air = cosine_spacing(sampling.num_airfoil_points)
    mask_tc_window = (x_air >= sections.x_tc_window[0]) & (x_air <= sections.x_tc_window[1])
    mask_valid_window = (x_air >= sections.x_valid_window[0]) & (x_air <= sections.x_valid_window[1])

    y_sections = topology.anchor_y_array
    section_specs = sections.section_specs
    coeff_sections = np.asarray(
        [np.concatenate([spec.upper_coeffs, spec.lower_coeffs]) for spec in section_specs],
        dtype=float,
    )
    tc_sections = np.asarray([spec.tc_max for spec in section_specs], dtype=float)
    x_tmax_sections = np.asarray([spec.x_tmax for spec in section_specs], dtype=float)
    te_sections = np.asarray([spec.te_thickness for spec in section_specs], dtype=float)

    scalar_target_interpolation = sampling.section_interpolation
    if sections.profile_generation_mode == "enforce_targets" and sampling.section_interpolation == "pyspline":
        # Shape-preserving interpolation avoids pySpline overshoots in scalar
        # targets such as x_tmax, which can otherwise move the thickness peak
        # unrealistically close to the leading edge between control sections.
        scalar_target_interpolation = "pchip"

    # CST coefficients can legitimately be signed for realistic cambered
    # airfoils, so they must use a signed interpolant instead of the
    # log-space positive interpolator used for thickness-like scalars.
    coeff_interpolants = [
        build_scalar_interpolant(y_sections, coeff_sections[:, idx], sampling.section_interpolation)
        for idx in range(sections.total_coeff_count)
    ]
    tc_interpolant = build_positive_interpolant(y_sections, tc_sections, scalar_target_interpolation)
    x_tmax_interpolant = build_scalar_interpolant(y_sections, x_tmax_sections, scalar_target_interpolation)
    te_interpolant = build_positive_interpolant(y_sections, te_sections, scalar_target_interpolation)

    interpolation_name_map = {
        "linear": "Linear",
        "pchip": "PCHIP",
        "cubic": "CubicSpline(natural, C2)",
        "pyspline": "pySpline Curve(k=4)",
    }
    interpolation_name = interpolation_name_map[sampling.section_interpolation]

    return SectionModel(
        shape=shape,
        x_air=x_air,
        y_sections=y_sections,
        coeff_sections=coeff_sections,
        tc_sections=tc_sections,
        x_tmax_sections=x_tmax_sections,
        te_sections=te_sections,
        mask_tc_window=mask_tc_window,
        mask_valid_window=mask_valid_window,
        coeff_interpolants=coeff_interpolants,
        tc_interpolant=tc_interpolant,
        x_tmax_interpolant=x_tmax_interpolant,
        camber_interpolant=laws.camber_delta,
        te_interpolant=te_interpolant,
        camber_mode_center=sections.camber_mode_center,
        camber_mode_width=sections.camber_mode_width,
        profile_generation_mode=sections.profile_generation_mode,
        interpolation_name=interpolation_name,
    )
