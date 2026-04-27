from dataclasses import dataclass
from typing import Callable, List, Optional

import numpy as np
from scipy.interpolate import PchipInterpolator

from parametrization.shared.cst import KulfanCSTAirfoil, cosine_spacing
from .specs import SectionedBWBModelConfig
from .spanwise_laws import (
    ResolvedSpanwiseLaws,
    build_scalar_interpolant,
    interpolation_display_name,
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
    camber_curve_interpolants: Optional[List[Callable[[float], float]]]
    thickness_curve_interpolants: Optional[List[Callable[[float], float]]]
    camber_mode_center: float
    camber_mode_width: float
    profile_generation_mode: str
    interpolation_name: str
    shape_hold_ranges: tuple[tuple[float, float, int], ...]
    shape_hold_tc_ramp_ranges: tuple[tuple[float, float, int, float], ...]
    exact_profile_overrides: tuple[
        tuple[float, Optional[float], bool, str, np.ndarray, np.ndarray, SectionShapeParameters],
        ...
    ]
    profile_post_transform: Optional[Callable[..., tuple[np.ndarray, np.ndarray]]] = None

    def _exact_profile_override(
        self,
        y: float,
    ) -> Optional[tuple[float, Optional[float], bool, str, np.ndarray, np.ndarray, SectionShapeParameters]]:
        for override_y, override_y_end, blend_to_base, blend_curve, upper, lower, params in self.exact_profile_overrides:
            if override_y_end is None:
                matches = abs(float(y) - override_y) <= 1.0e-12
            else:
                matches = (float(y) >= override_y - 1.0e-12) and (float(y) < float(override_y_end) - 1.0e-12)
            if matches:
                return override_y, override_y_end, blend_to_base, blend_curve, upper.copy(), lower.copy(), params
        return None

    def _base_params_at_y(self, y: float) -> SectionShapeParameters:
        source_index = self._shape_hold_source_index(float(y))
        if source_index is not None:
            coeffs = np.array(self.coeff_sections[source_index], dtype=float)
            coeffs = apply_camber_delta(
                coeffs,
                camber_delta=float(self.camber_interpolant(y)),
                degree=self.shape.degree,
                center=self.camber_mode_center,
                width=self.camber_mode_width,
            )
            return SectionShapeParameters(
                coeffs=coeffs,
                tc_max=float(self.tc_sections[source_index]),
                x_tmax=float(self.x_tmax_sections[source_index]),
                te_thickness=float(self.te_sections[source_index]),
            )
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

    def _base_coordinates_at_y(self, y: float) -> tuple[np.ndarray, np.ndarray, SectionShapeParameters]:
        source_index = self._shape_hold_source_index(y)
        params = self._base_params_at_y(y)
        if self.profile_generation_mode == "camber_thickness_interp" and source_index is None:
            if self.camber_curve_interpolants is None or self.thickness_curve_interpolants is None:
                raise RuntimeError(
                    "camber_thickness_interp mode requires pointwise camber/thickness interpolants"
                )
            camber = np.array([interp(y) for interp in self.camber_curve_interpolants], dtype=float)
            thickness = np.array([interp(y) for interp in self.thickness_curve_interpolants], dtype=float)
            thickness = thickness + self.x_air * (params.te_thickness - float(thickness[-1]))
            thickness = np.maximum(thickness, 0.0)
            yu = camber + 0.5 * thickness
            yl = camber - 0.5 * thickness
            return yu, yl, params
        evaluate_kwargs = {
            "te_thickness": params.te_thickness,
        }
        if self.profile_generation_mode == "enforce_targets":
            evaluate_kwargs["tc_target"] = params.tc_max
            evaluate_kwargs["x_tmax"] = params.x_tmax
        yu, yl = self.shape.evaluate(self.x_air, params.coeffs, **evaluate_kwargs)
        tc_ramp = self._shape_hold_tc_ramp(y)
        if source_index is not None and tc_ramp is not None and source_index == tc_ramp[2]:
            y_start, y_end, _, target_tc_end = tc_ramp
            t = np.clip((float(y) - float(y_start)) / max(float(y_end) - float(y_start), 1.0e-12), 0.0, 1.0)
            start_tc = float(self.tc_sections[source_index])
            target_tc = (1.0 - t) * start_tc + t * float(target_tc_end)
            camber = 0.5 * (yu + yl)
            thickness = yu - yl
            current_tc = float(np.max(thickness))
            if current_tc > 1.0e-12:
                thickness = thickness * (target_tc / current_tc)
                yu = camber + 0.5 * thickness
                yl = camber - 0.5 * thickness
                params = SectionShapeParameters(
                    coeffs=params.coeffs,
                    tc_max=float(np.max(thickness)),
                    x_tmax=float(self.x_air[int(np.argmax(thickness))]),
                    te_thickness=float(thickness[-1]),
                )
        return yu, yl, params

    def _shape_hold_source_index(self, y: float) -> Optional[int]:
        for y_start, y_end, source_index in self.shape_hold_ranges:
            if y >= y_start - 1.0e-12 and y < y_end - 1.0e-12:
                return source_index
        return None

    def _shape_hold_tc_ramp(self, y: float) -> Optional[tuple[float, float, int, float]]:
        for y_start, y_end, source_index, target_tc in self.shape_hold_tc_ramp_ranges:
            if y >= y_start - 1.0e-12 and y < y_end - 1.0e-12:
                return y_start, y_end, source_index, target_tc
        return None

    def params_at_y(self, y: float) -> SectionShapeParameters:
        exact_override = self._exact_profile_override(float(y))
        if exact_override is not None:
            override_y, override_y_end, blend_to_base, _, upper, lower, params = exact_override
            if not blend_to_base or override_y_end is None:
                return params
            base_params = self._base_params_at_y(float(y))
            thickness = upper - lower
            return SectionShapeParameters(
                coeffs=base_params.coeffs,
                tc_max=float(np.max(thickness)),
                x_tmax=float(self.x_air[int(np.argmax(thickness))]),
                te_thickness=float(thickness[-1]),
            )
        return self._base_params_at_y(float(y))

    def coordinates_at_y(self, y: float) -> tuple[np.ndarray, np.ndarray, SectionShapeParameters]:
        y = float(y)
        exact_override = self._exact_profile_override(y)
        if exact_override is not None:
            override_y, override_y_end, blend_to_base, blend_curve, upper, lower, params = exact_override
            if not blend_to_base or override_y_end is None:
                if self.profile_post_transform is None:
                    return upper, lower, params
                out_upper, out_lower = self.profile_post_transform(y, self.x_air, upper.copy(), lower.copy(), params)
                thickness = out_upper - out_lower
                return out_upper, out_lower, SectionShapeParameters(
                    coeffs=params.coeffs,
                    tc_max=float(np.max(thickness)),
                    x_tmax=float(self.x_air[int(np.argmax(thickness))]),
                    te_thickness=float(thickness[-1]),
                )
            base_upper, base_lower, base_params = self._base_coordinates_at_y(y)
            t = np.clip((y - float(override_y)) / max(float(override_y_end) - float(override_y), 1e-12), 0.0, 1.0)
            smooth_t = t * t * (3.0 - 2.0 * t)
            if blend_curve == "ellipse":
                w = 1.0 - np.sqrt(max(1.0 - smooth_t * smooth_t, 0.0))
            else:
                w = smooth_t
            blend_upper = (1.0 - w) * upper + w * base_upper
            blend_lower = (1.0 - w) * lower + w * base_lower
            thickness = blend_upper - blend_lower
            params_out = SectionShapeParameters(
                coeffs=base_params.coeffs,
                tc_max=float(np.max(thickness)),
                x_tmax=float(self.x_air[int(np.argmax(thickness))]),
                te_thickness=float(thickness[-1]),
            )
            if self.profile_post_transform is None:
                return blend_upper, blend_lower, params_out
            out_upper, out_lower = self.profile_post_transform(
                y,
                self.x_air,
                blend_upper.copy(),
                blend_lower.copy(),
                params_out,
            )
            thickness = out_upper - out_lower
            return out_upper, out_lower, SectionShapeParameters(
                coeffs=params_out.coeffs,
                tc_max=float(np.max(thickness)),
                x_tmax=float(self.x_air[int(np.argmax(thickness))]),
                te_thickness=float(thickness[-1]),
            )
        upper, lower, params = self._base_coordinates_at_y(y)
        if self.profile_post_transform is None:
            return upper, lower, params
        out_upper, out_lower = self.profile_post_transform(y, self.x_air, upper.copy(), lower.copy(), params)
        thickness = out_upper - out_lower
        return out_upper, out_lower, SectionShapeParameters(
            coeffs=params.coeffs,
            tc_max=float(np.max(thickness)),
            x_tmax=float(self.x_air[int(np.argmax(thickness))]),
            te_thickness=float(thickness[-1]),
        )

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
    linear_start_index: int | None = None,
    interpolation_factory=None,
):
    values = np.maximum(np.asarray(values, dtype=float), floor)
    log_interpolant = build_scalar_interpolant(
        y_sections,
        np.log(values),
        interpolation,
        linear_start_index=linear_start_index,
        interpolation_factory=interpolation_factory,
    )
    return lambda yy: float(np.exp(log_interpolant(float(yy))))


def build_anchor_profile_arrays(
    shape: KulfanCSTAirfoil,
    x_air: np.ndarray,
    y_sections: np.ndarray,
    coeff_sections: np.ndarray,
    tc_sections: np.ndarray,
    x_tmax_sections: np.ndarray,
    te_sections: np.ndarray,
    camber_interpolant: Callable[[float], float],
    degree: int,
    center: float,
    width: float,
) -> tuple[np.ndarray, np.ndarray]:
    camber_sections = np.empty((y_sections.size, x_air.size), dtype=float)
    thickness_sections = np.empty((y_sections.size, x_air.size), dtype=float)

    for idx, yy in enumerate(y_sections):
        coeffs = apply_camber_delta(
            coeff_sections[idx],
            camber_delta=float(camber_interpolant(float(yy))),
            degree=degree,
            center=center,
            width=width,
        )
        yu, yl = shape.evaluate(
            x_air,
            coeffs,
            te_thickness=float(te_sections[idx]),
            tc_target=float(tc_sections[idx]),
            x_tmax=float(x_tmax_sections[idx]),
        )
        camber_sections[idx] = 0.5 * (yu + yl)
        thickness_sections[idx] = yu - yl

    return camber_sections, thickness_sections


def _interp_override_curve(x_src: np.ndarray, y_src: np.ndarray, x_dst: np.ndarray) -> np.ndarray:
    if x_src.shape == x_dst.shape and np.allclose(x_src, x_dst, rtol=0.0, atol=1.0e-12):
        return np.asarray(y_src, dtype=float)
    interpolant = PchipInterpolator(np.asarray(x_src, dtype=float), np.asarray(y_src, dtype=float))
    return np.asarray(interpolant(np.asarray(x_dst, dtype=float)), dtype=float)


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
        build_scalar_interpolant(
            y_sections,
            coeff_sections[:, idx],
            sampling.section_interpolation,
            linear_start_index=sampling.section_linear_start_index,
            interpolation_factory=sampling.section_interpolation_factory,
        )
        for idx in range(sections.total_coeff_count)
    ]
    tc_interpolant = build_positive_interpolant(
        y_sections,
        tc_sections,
        scalar_target_interpolation,
        linear_start_index=sampling.section_linear_start_index,
        interpolation_factory=sampling.section_interpolation_factory,
    )
    x_tmax_interpolant = build_scalar_interpolant(
        y_sections,
        x_tmax_sections,
        scalar_target_interpolation,
        linear_start_index=sampling.section_linear_start_index,
        interpolation_factory=sampling.section_interpolation_factory,
    )
    te_interpolant = build_positive_interpolant(
        y_sections,
        te_sections,
        scalar_target_interpolation,
        linear_start_index=sampling.section_linear_start_index,
        interpolation_factory=sampling.section_interpolation_factory,
    )
    camber_curve_interpolants = None
    thickness_curve_interpolants = None

    if sections.profile_generation_mode == "camber_thickness_interp":
        camber_sections, thickness_sections = build_anchor_profile_arrays(
            shape=shape,
            x_air=x_air,
            y_sections=y_sections,
            coeff_sections=coeff_sections,
            tc_sections=tc_sections,
            x_tmax_sections=x_tmax_sections,
            te_sections=te_sections,
            camber_interpolant=laws.camber_delta,
            degree=sections.cst_degree,
            center=sections.camber_mode_center,
            width=sections.camber_mode_width,
        )
        camber_curve_interpolants = [
            build_scalar_interpolant(
                y_sections,
                camber_sections[:, idx],
                sampling.section_interpolation,
                linear_start_index=sampling.section_linear_start_index,
                interpolation_factory=sampling.section_interpolation_factory,
            )
            for idx in range(x_air.size)
        ]
        thickness_curve_interpolants = [
            build_scalar_interpolant(
                y_sections,
                thickness_sections[:, idx],
                sampling.section_interpolation,
                linear_start_index=sampling.section_linear_start_index,
                interpolation_factory=sampling.section_interpolation_factory,
            )
            for idx in range(x_air.size)
        ]

    interpolation_name = interpolation_display_name(
        sampling.section_interpolation,
        interpolation_factory=sampling.section_interpolation_factory,
    )
    if sampling.section_linear_start_index is not None:
        interpolation_name = (
            f"{interpolation_name} -> Linear from section {int(sampling.section_linear_start_index)}"
        )
    shape_hold_ranges = tuple(
        (
            float(y_sections[int(start_idx)]),
            float(y_sections[int(start_idx) + 1]),
            int(start_idx),
        )
        for start_idx in sections.shape_hold_segments
    )
    shape_hold_tc_ramp_ranges = tuple(
        (
            float(y_sections[int(start_idx)]),
            float(y_sections[int(start_idx) + 1]),
            int(start_idx),
            float(tc_sections[int(start_idx) + 1]),
        )
        for start_idx in sections.shape_hold_tc_ramp_segments
    )
    exact_profile_overrides = tuple(
            (
                float(y_sections[int(override.section_index)]),
                (
                    float(y_sections[int(override.y_end_section_index)])
                    if override.y_end_section_index is not None
                    else (None if override.y_end is None else float(override.y_end))
                ),
                bool(override.blend_to_base),
                str(override.blend_curve),
                _interp_override_curve(np.asarray(override.x, dtype=float), np.asarray(override.upper, dtype=float), x_air),
                _interp_override_curve(np.asarray(override.x, dtype=float), np.asarray(override.lower, dtype=float), x_air),
                SectionShapeParameters(
                coeffs=np.zeros(sections.total_coeff_count, dtype=float),
                tc_max=float(override.tc_max),
                x_tmax=float(override.x_tmax),
                te_thickness=float(override.te_thickness),
            ),
        )
        for override in sections.exact_profile_overrides
    )

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
        camber_curve_interpolants=camber_curve_interpolants,
        thickness_curve_interpolants=thickness_curve_interpolants,
        camber_mode_center=sections.camber_mode_center,
        camber_mode_width=sections.camber_mode_width,
        profile_generation_mode=sections.profile_generation_mode,
        interpolation_name=interpolation_name,
        shape_hold_ranges=shape_hold_ranges,
        shape_hold_tc_ramp_ranges=shape_hold_tc_ramp_ranges,
        exact_profile_overrides=exact_profile_overrides,
        profile_post_transform=sections.profile_post_transform,
    )
