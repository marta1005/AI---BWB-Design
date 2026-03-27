from dataclasses import dataclass
import re
from typing import Dict, List, Optional, Tuple

import numpy as np

from .design_variables import SectionedBWBDesignVariables


def parameter_groups() -> Dict[str, Tuple[str, ...]]:
    return {
        "spans": (
            "span",
            "b1_span_ratio",
            "b2_span_ratio",
            "b3_span_ratio",
        ),
        "nose_blend": (
            "nose_blend_y",
        ),
        "chords": (
            "c1_root_chord",
            "c2_c1_ratio",
            "c3_c1_ratio",
            "c4_c1_ratio",
        ),
        "sweeps": (
            "s2_deg",
            "s3_deg",
        ),
        "sweeps_full": (
            "s1_deg",
            "s2_deg",
            "s3_deg",
        ),
        "topology": (
            "span",
            "b1_span_ratio",
            "b2_span_ratio",
            "b3_span_ratio",
        ),
        "planform": (
            "le_root_x",
            "c1_root_chord",
            "c2_c1_ratio",
            "c3_c1_ratio",
            "c4_c1_ratio",
            "s1_deg",
            "s2_deg",
            "s3_deg",
            "nose_blend_y",
        ),
        "class_function": (
            "cst_n1",
            "cst_n2",
        ),
        "attitude": (
            "dihedral_deg",
        ),
        "twist": (
            "twist_c1_deg",
            "twist_c2_deg",
            "twist_c3_deg",
            "twist_c4_deg",
        ),
        "camber_mode": (
            "camber_c1",
            "camber_c2",
            "camber_c3",
            "camber_c4",
        ),
        "thickness_targets": (
            "c1_tc_max",
            "c2_tc_max",
            "c3_tc_max",
            "c4_tc_max",
            "c1_x_tmax",
            "c2_x_tmax",
            "c3_x_tmax",
            "c4_x_tmax",
            "c1_te_thickness",
            "c2_te_thickness",
            "c3_te_thickness",
            "c4_te_thickness",
        ),
        "cst_c1": tuple(
            ["c1_upper_cst_%d" % idx for idx in range(6)] + ["c1_lower_cst_%d" % idx for idx in range(6)]
        ),
        "cst_c2": tuple(
            ["c2_upper_cst_%d" % idx for idx in range(6)] + ["c2_lower_cst_%d" % idx for idx in range(6)]
        ),
        "cst_c3": tuple(
            ["c3_upper_cst_%d" % idx for idx in range(6)] + ["c3_lower_cst_%d" % idx for idx in range(6)]
        ),
        "cst_c4": tuple(
            ["c4_upper_cst_%d" % idx for idx in range(6)] + ["c4_lower_cst_%d" % idx for idx in range(6)]
        ),
    }


def preset_groups() -> Dict[str, Tuple[str, ...]]:
    return {
        "presentation_core": ("planform", "twist", "thickness_targets"),
        "planform_only": ("topology", "planform"),
        "aero_sections": ("twist", "camber_mode", "thickness_targets"),
        "ai_geometry_core": (
            "spans",
            "chords",
            "sweeps",
            "nose_blend",
            "twist",
            "cst_c1",
            "cst_c2",
            "cst_c3",
            "cst_c4",
        ),
        "gemseo_geometry_core": (
            "spans",
            "chords",
            "sweeps",
            "twist",
            "cst_c1",
            "cst_c2",
            "cst_c3",
            "cst_c4",
        ),
        "section_shapes": ("thickness_targets", "cst_c1", "cst_c2", "cst_c3", "cst_c4"),
        "full_geometry": (
            "topology",
            "planform",
            "class_function",
            "attitude",
            "twist",
            "camber_mode",
            "thickness_targets",
            "cst_c1",
            "cst_c2",
            "cst_c3",
            "cst_c4",
        ),
    }


def flatten_design(design: SectionedBWBDesignVariables) -> Dict[str, float]:
    values = design.as_vector()
    return dict(zip(SectionedBWBDesignVariables.variable_names(), values))


def _scalar_parameter_catalog() -> Dict[str, Dict[str, str]]:
    return {
        "span": {
            "display_name": "Semi-span",
            "symbol": "b/2",
            "units": "m",
            "normalization": "absolute",
            "description": "Half-span of the BWB planform.",
        },
        "b1_span_ratio": {
            "display_name": "Span segment B1",
            "symbol": "B1/(b/2)",
            "units": "-",
            "normalization": "fraction of semi-span",
            "description": "Inboard span partition ratio.",
        },
        "b2_span_ratio": {
            "display_name": "Span segment B2",
            "symbol": "B2/(b/2)",
            "units": "-",
            "normalization": "fraction of semi-span",
            "description": "Transition span partition ratio.",
        },
        "b3_span_ratio": {
            "display_name": "Span segment B3",
            "symbol": "B3/(b/2)",
            "units": "-",
            "normalization": "fraction of semi-span",
            "description": "Outboard span partition ratio.",
        },
        "le_root_x": {
            "display_name": "Root LE x-position",
            "symbol": "x_LE,root",
            "units": "m",
            "normalization": "absolute",
            "description": "Leading-edge x-location at the root section.",
        },
        "c1_root_chord": {
            "display_name": "Root chord",
            "symbol": "C1",
            "units": "m",
            "normalization": "absolute",
            "description": "Chord at the root control section.",
        },
        "c2_c1_ratio": {
            "display_name": "Section C2 chord ratio",
            "symbol": "C2/C1",
            "units": "-",
            "normalization": "ratio to C1",
            "description": "Chord ratio at section C2 relative to the root chord.",
        },
        "c3_c1_ratio": {
            "display_name": "Section C3 chord ratio",
            "symbol": "C3/C1",
            "units": "-",
            "normalization": "ratio to C1",
            "description": "Chord ratio at section C3 relative to the root chord.",
        },
        "c4_c1_ratio": {
            "display_name": "Section C4 chord ratio",
            "symbol": "C4/C1",
            "units": "-",
            "normalization": "ratio to C1",
            "description": "Chord ratio at section C4 relative to the root chord.",
        },
        "s1_deg": {
            "display_name": "Sweep segment S1",
            "symbol": "S1",
            "units": "deg",
            "normalization": "angle",
            "description": "Leading-edge sweep of segment S1.",
        },
        "s2_deg": {
            "display_name": "Sweep segment S2",
            "symbol": "S2",
            "units": "deg",
            "normalization": "angle",
            "description": "Leading-edge sweep of segment S2.",
        },
        "s3_deg": {
            "display_name": "Sweep segment S3",
            "symbol": "S3",
            "units": "deg",
            "normalization": "angle",
            "description": "Leading-edge sweep of segment S3.",
        },
        "nose_blend_y": {
            "display_name": "Root nose blend length",
            "symbol": "y_blend,nose",
            "units": "m",
            "normalization": "absolute",
            "description": "Spanwise blending length used to round the root nose.",
        },
        "cst_n1": {
            "display_name": "CST class exponent N1",
            "symbol": "N1",
            "units": "-",
            "normalization": "class-function exponent",
            "description": "Leading-edge exponent of the Kulfan class function.",
        },
        "cst_n2": {
            "display_name": "CST class exponent N2",
            "symbol": "N2",
            "units": "-",
            "normalization": "class-function exponent",
            "description": "Trailing-edge exponent of the Kulfan class function.",
        },
        "dihedral_deg": {
            "display_name": "Global dihedral",
            "symbol": "Gamma",
            "units": "deg",
            "normalization": "angle",
            "description": "Global dihedral angle of the wing.",
        },
        "twist_c1_deg": {
            "display_name": "Twist at C1",
            "symbol": "twist_C1",
            "units": "deg",
            "normalization": "angle",
            "description": "Geometric twist at section C1.",
        },
        "twist_c2_deg": {
            "display_name": "Twist at C2",
            "symbol": "twist_C2",
            "units": "deg",
            "normalization": "angle",
            "description": "Geometric twist at section C2.",
        },
        "twist_c3_deg": {
            "display_name": "Twist at C3",
            "symbol": "twist_C3",
            "units": "deg",
            "normalization": "angle",
            "description": "Geometric twist at section C3.",
        },
        "twist_c4_deg": {
            "display_name": "Twist at C4",
            "symbol": "twist_C4",
            "units": "deg",
            "normalization": "angle",
            "description": "Geometric twist at section C4.",
        },
        "camber_c1": {
            "display_name": "Camber delta at C1",
            "symbol": "dcamber_C1",
            "units": "-",
            "normalization": "mode amplitude",
            "description": "Additional camber mode amplitude at section C1.",
        },
        "camber_c2": {
            "display_name": "Camber delta at C2",
            "symbol": "dcamber_C2",
            "units": "-",
            "normalization": "mode amplitude",
            "description": "Additional camber mode amplitude at section C2.",
        },
        "camber_c3": {
            "display_name": "Camber delta at C3",
            "symbol": "dcamber_C3",
            "units": "-",
            "normalization": "mode amplitude",
            "description": "Additional camber mode amplitude at section C3.",
        },
        "camber_c4": {
            "display_name": "Camber delta at C4",
            "symbol": "dcamber_C4",
            "units": "-",
            "normalization": "mode amplitude",
            "description": "Additional camber mode amplitude at section C4.",
        },
        "c1_tc_max": {
            "display_name": "Maximum thickness at C1",
            "symbol": "(t/c)_max,C1",
            "units": "-",
            "normalization": "t/c",
            "description": "Maximum thickness ratio target at section C1.",
        },
        "c2_tc_max": {
            "display_name": "Maximum thickness at C2",
            "symbol": "(t/c)_max,C2",
            "units": "-",
            "normalization": "t/c",
            "description": "Maximum thickness ratio target at section C2.",
        },
        "c3_tc_max": {
            "display_name": "Maximum thickness at C3",
            "symbol": "(t/c)_max,C3",
            "units": "-",
            "normalization": "t/c",
            "description": "Maximum thickness ratio target at section C3.",
        },
        "c4_tc_max": {
            "display_name": "Maximum thickness at C4",
            "symbol": "(t/c)_max,C4",
            "units": "-",
            "normalization": "t/c",
            "description": "Maximum thickness ratio target at section C4.",
        },
        "c1_x_tmax": {
            "display_name": "Thickness peak location at C1",
            "symbol": "x_tmax,C1/c",
            "units": "-",
            "normalization": "x/c",
            "description": "Chordwise position of the maximum thickness at section C1.",
        },
        "c2_x_tmax": {
            "display_name": "Thickness peak location at C2",
            "symbol": "x_tmax,C2/c",
            "units": "-",
            "normalization": "x/c",
            "description": "Chordwise position of the maximum thickness at section C2.",
        },
        "c3_x_tmax": {
            "display_name": "Thickness peak location at C3",
            "symbol": "x_tmax,C3/c",
            "units": "-",
            "normalization": "x/c",
            "description": "Chordwise position of the maximum thickness at section C3.",
        },
        "c4_x_tmax": {
            "display_name": "Thickness peak location at C4",
            "symbol": "x_tmax,C4/c",
            "units": "-",
            "normalization": "x/c",
            "description": "Chordwise position of the maximum thickness at section C4.",
        },
        "c1_te_thickness": {
            "display_name": "Trailing-edge thickness at C1",
            "symbol": "t_TE,C1/c",
            "units": "-",
            "normalization": "t_TE/c_local",
            "description": "Trailing-edge thickness ratio at section C1.",
        },
        "c2_te_thickness": {
            "display_name": "Trailing-edge thickness at C2",
            "symbol": "t_TE,C2/c",
            "units": "-",
            "normalization": "t_TE/c_local",
            "description": "Trailing-edge thickness ratio at section C2.",
        },
        "c3_te_thickness": {
            "display_name": "Trailing-edge thickness at C3",
            "symbol": "t_TE,C3/c",
            "units": "-",
            "normalization": "t_TE/c_local",
            "description": "Trailing-edge thickness ratio at section C3.",
        },
        "c4_te_thickness": {
            "display_name": "Trailing-edge thickness at C4",
            "symbol": "t_TE,C4/c",
            "units": "-",
            "normalization": "t_TE/c_local",
            "description": "Trailing-edge thickness ratio at section C4.",
        },
    }


def _cst_parameter_info(parameter_name: str) -> Optional[Dict[str, str]]:
    match = re.match(r"c([1-4])_(upper|lower)_cst_(\d+)$", parameter_name)
    if not match:
        return None

    section_idx = int(match.group(1))
    surface = match.group(2).upper()
    coeff_idx = int(match.group(3))
    return {
        "display_name": "C%d %s CST coefficient %d" % (section_idx, surface, coeff_idx),
        "symbol": "A_%s,%d^C%d" % (surface.lower(), coeff_idx, section_idx),
        "units": "-",
        "normalization": "Bernstein coefficient",
        "description": "Kulfan CST Bernstein coefficient %d on the %s surface of section C%d."
        % (coeff_idx, surface.lower(), section_idx),
    }


def parameter_info(parameter_name: str) -> Dict[str, str]:
    scalar_catalog = _scalar_parameter_catalog()
    if parameter_name in scalar_catalog:
        return scalar_catalog[parameter_name]
    cst_info = _cst_parameter_info(parameter_name)
    if cst_info is not None:
        return cst_info
    return {
        "display_name": parameter_name,
        "symbol": parameter_name,
        "units": "-",
        "normalization": "-",
        "description": "No description available.",
    }


def parameter_metadata(
    reference_design: Optional[SectionedBWBDesignVariables] = None,
) -> List[Dict[str, object]]:
    if reference_design is None:
        reference_design = SectionedBWBDesignVariables.reference_seed()
    flat_reference = flatten_design(reference_design)
    bounds = SectionedBWBDesignVariables.default_bounds()
    groups = parameter_groups()
    reverse_group = {}
    for group_name, names in groups.items():
        for name in names:
            reverse_group[name] = group_name

    rows: List[Dict[str, object]] = []
    for name in SectionedBWBDesignVariables.variable_names():
        lower, upper = bounds[name]
        info = parameter_info(name)
        rows.append(
            {
                "parameter": name,
                "group": reverse_group.get(name, "ungrouped"),
                "display_name": info["display_name"],
                "symbol": info["symbol"],
                "units": info["units"],
                "normalization": info["normalization"],
                "description": info["description"],
                "reference": float(flat_reference[name]),
                "lower_bound": float(lower),
                "upper_bound": float(upper),
            }
        )
    return rows


@dataclass
class DesignSpace:
    preset_name: str
    active_groups: Tuple[str, ...]
    active_variables: Tuple[str, ...]
    reference_design: SectionedBWBDesignVariables
    bounds: Dict[str, Tuple[float, float]]

    def reference_flat(self) -> Dict[str, float]:
        return flatten_design(self.reference_design)

    def active_metadata(self) -> List[Dict[str, object]]:
        metadata = {row["parameter"]: row for row in parameter_metadata(self.reference_design)}
        return [metadata[name] for name in self.active_variables]

    def sample_designs(
        self,
        count: int,
        seed: int = 7,
        variation_scale: float = 0.25,
    ) -> List[SectionedBWBDesignVariables]:
        rng = np.random.default_rng(int(seed))
        reference_flat = self.reference_flat()
        index_map = {
            name: idx
            for idx, name in enumerate(SectionedBWBDesignVariables.variable_names())
        }
        reference_vector = self.reference_design.as_vector()
        sampled: List[SectionedBWBDesignVariables] = []

        for _ in range(int(count)):
            vector = reference_vector.copy()
            active_set = set(self.active_variables)
            if {"b1_span_ratio", "b2_span_ratio", "b3_span_ratio"} & active_set:
                b1, b2, b3 = self._sample_b_ratios(rng, reference_flat, variation_scale)
                vector[index_map["b1_span_ratio"]] = b1
                vector[index_map["b2_span_ratio"]] = b2
                vector[index_map["b3_span_ratio"]] = b3

            for name in self.active_variables:
                if name in {"b1_span_ratio", "b2_span_ratio", "b3_span_ratio"}:
                    continue
                lower, upper = self._local_bounds(name, reference_flat[name], variation_scale)
                vector[index_map[name]] = rng.uniform(lower, upper)

            sampled.append(SectionedBWBDesignVariables.from_vector(vector))
        return sampled

    def _local_bounds(self, name: str, reference_value: float, variation_scale: float) -> Tuple[float, float]:
        global_lower, global_upper = self.bounds[name]
        span = global_upper - global_lower
        local_lower = max(global_lower, reference_value - variation_scale * span)
        local_upper = min(global_upper, reference_value + variation_scale * span)
        if local_upper <= local_lower:
            return global_lower, global_upper
        return float(local_lower), float(local_upper)

    def _sample_b_ratios(
        self,
        rng: np.random.Generator,
        reference_flat: Dict[str, float],
        variation_scale: float,
    ) -> Tuple[float, float, float]:
        for _ in range(200):
            b1_lower, b1_upper = self._local_bounds("b1_span_ratio", reference_flat["b1_span_ratio"], variation_scale)
            b2_lower, b2_upper = self._local_bounds("b2_span_ratio", reference_flat["b2_span_ratio"], variation_scale)
            b1 = float(rng.uniform(b1_lower, b1_upper))
            b2 = float(rng.uniform(b2_lower, b2_upper))
            b3 = 1.0 - b1 - b2
            lower_b3, upper_b3 = self.bounds["b3_span_ratio"]
            if lower_b3 <= b3 <= upper_b3:
                return b1, b2, float(b3)
        return (
            float(reference_flat["b1_span_ratio"]),
            float(reference_flat["b2_span_ratio"]),
            float(reference_flat["b3_span_ratio"]),
        )


def build_design_space(
    preset_name: str,
    reference_design: Optional[SectionedBWBDesignVariables] = None,
) -> DesignSpace:
    presets = preset_groups()
    if preset_name not in presets:
        raise ValueError(
            "Unknown design-space preset %r. Available presets: %s"
            % (preset_name, ", ".join(sorted(presets)))
        )

    if reference_design is None:
        reference_design = SectionedBWBDesignVariables.reference_seed()
    groups = parameter_groups()
    active_groups = presets[preset_name]
    active_variables: List[str] = []
    for group_name in active_groups:
        active_variables.extend(groups[group_name])

    return DesignSpace(
        preset_name=preset_name,
        active_groups=active_groups,
        active_variables=tuple(active_variables),
        reference_design=reference_design,
        bounds=SectionedBWBDesignVariables.default_bounds(),
    )


def recommended_design_space() -> DesignSpace:
    return build_design_space("presentation_core")


def available_presets() -> Tuple[str, ...]:
    return tuple(sorted(preset_groups().keys()))


def group_for_parameter(parameter_name: str) -> str:
    for group_name, names in parameter_groups().items():
        if parameter_name in names:
            return group_name
    return "ungrouped"


def preset_parameter_metadata(
    preset_name: str,
    reference_design: Optional[SectionedBWBDesignVariables] = None,
) -> List[Dict[str, object]]:
    design_space = build_design_space(preset_name, reference_design=reference_design)
    return design_space.active_metadata()
