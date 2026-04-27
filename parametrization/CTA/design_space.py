"""CTA design-space definition for AI sampling and optimization."""

from dataclasses import asdict, dataclass
import re
from typing import Dict, List, Optional, Tuple

import numpy as np

from .case import (
    CTA_C3_TRANSITION_CHORD_M,
    SWEEP_NAME_TO_VARIABLE,
    VARIABLE_TO_SWEEP_NAME,
    _chord_sweep_from_leading_edge_sweep,
    apply_cta_fixed_parameters,
    build_cta_design,
    cta_fixed_values,
)
from parametrization.bwb.design_space import DesignSpace, flatten_design, parameter_info
from parametrization.bwb.design_variables import SectionedBWBDesignVariables

CTA_FIXED_PARAMETERS: Tuple[str, ...] = (
    "b1_fixed_m",
    "s1_deg",
)


def _cst_names(section_index: int) -> Tuple[str, ...]:
    return tuple(
        [f"c{section_index}_upper_cst_{idx}" for idx in range(6)]
        + [f"c{section_index}_lower_cst_{idx}" for idx in range(6)]
    )


CTA_ACTIVE_VARIABLES: Tuple[str, ...] = (
    "span",
    "c1_root_chord",
    "c2_c1_ratio",
    "c4_c3_ratio",
    "b2_span_ratio",
    "c4_c1_ratio",
    "s2_deg",  # CTA label S1
    "s3_deg",  # CTA label S2
    "med_3_te_sweep_deg",
    "twist_c1_deg",
    "twist_c3_deg",
    "twist_c4_deg",
    *_cst_names(1),
    *_cst_names(2),
    *_cst_names(3),
    *_cst_names(4),
)

CTA_PUBLIC_ACTIVE_BOUNDS: Dict[str, Tuple[float, float]] = {
    "span": (28.0, 35.0),  # wing span = B2 + B3
    "c1_root_chord": (39.0, 43.0),  # C0
    "c2_c1_ratio": (13.0, 16.0),  # C3 in absolute meters at the CTA layer
    "c4_c3_ratio": (0.45, 0.60),  # C4/C3 transition taper ratio
    "b2_span_ratio": (0.13, 0.21),  # B2 / (B2 + B3)
    "c4_c1_ratio": (0.80, 1.80),  # C5 in absolute meters at the CTA layer
    "s2_deg": (15.0, 45.0),  # S1
    "s3_deg": (22.0, 33.0),  # S2
    "med_3_te_sweep_deg": (-10.0, 25.0),  # med_3_TEswp
    "twist_c1_deg": (0.2, 2.0),
    "twist_c3_deg": (0.1, 1.5),
    "twist_c4_deg": (0.1, 1.0),
}

_CTA_CST_BOUND_HALF_WIDTHS: Tuple[float, ...] = (0.06, 0.06, 0.05, 0.05, 0.04, 0.03)

_TWIST_PUBLIC_LABELS: Dict[str, str] = {
    "twist_c1_deg": "C0/C3",
    "twist_c3_deg": "C4",
    "twist_c4_deg": "C5",
}

_SECTION_PUBLIC_LABELS: Dict[str, str] = {
    "c1": "C0",
    "c2": "C3",
    "c3": "C4",
    "c4": "C5",
}

CTA_THICKNESS_CONSTRAINTS: Tuple[Dict[str, object], ...] = (
    {
        "parameter": "RThickness_CentreBody",
        "station": "C0",
        "display_name": "Relative thickness centre body",
        "symbol": "RThickness_C0",
        "units": "-",
        "normalization": "t/c",
        "reference": 0.141,
        "lower_bound": 0.12,
        "upper_bound": 0.18,
        "description": "Thickness/chord ratio at the centre body root chord.",
    },
    {
        "parameter": "RThickness_MidWing",
        "station": "C3",
        "display_name": "Relative thickness mid wing",
        "symbol": "RThickness_C3",
        "units": "-",
        "normalization": "t/c",
        "reference": 0.164,
        "lower_bound": 0.14,
        "upper_bound": 0.20,
        "description": "Thickness/chord ratio at the mid wing root chord.",
    },
    {
        "parameter": "RThickness_OutWing",
        "station": "C4",
        "display_name": "Relative thickness outer wing",
        "symbol": "RThickness_C4",
        "units": "-",
        "normalization": "t/c",
        "reference": 0.1027,
        "lower_bound": 0.08,
        "upper_bound": 0.16,
        "description": "Thickness/chord ratio at the outer wing root chord.",
    },
    {
        "parameter": "RThickness_OutWing_Tip",
        "station": "C5",
        "display_name": "Relative thickness outer wing tip",
        "symbol": "RThickness_C5",
        "units": "-",
        "normalization": "t/c",
        "reference": 0.095,
        "lower_bound": 0.08,
        "upper_bound": 0.13,
        "description": "Thickness/chord ratio at the outer wing tip chord.",
    },
)


def _wing_span_m(span: float, b1_fixed_m: float) -> float:
    wing_span = float(span - b1_fixed_m)
    if wing_span <= 0.0:
        raise ValueError(
            f"Wing span must remain positive after removing fixed B1={b1_fixed_m:.6f} m from span={span:.6f} m"
        )
    return wing_span


def _b2_wing_span_ratio(span: float, b2_span_ratio: float, b1_fixed_m: float) -> float:
    return float((float(span) * float(b2_span_ratio)) / _wing_span_m(float(span), float(b1_fixed_m)))


def _internal_b2_span_ratio(span: float, b2_wing_ratio: float, b1_fixed_m: float) -> float:
    return float(float(b2_wing_ratio) * _wing_span_m(float(span), float(b1_fixed_m)) / float(span))


def _total_span_from_wing_span(wing_span: float, b1_fixed_m: float) -> float:
    return float(wing_span + b1_fixed_m)


def _c5_absolute_chord(c0_body_chord: float, c5_ratio_to_c0: float) -> float:
    return float(float(c0_body_chord) * float(c5_ratio_to_c0))


def _internal_c5_ratio(c0_body_chord: float, c5_absolute_chord: float) -> float:
    if c0_body_chord <= 0.0:
        raise ValueError(f"C0/body chord must be positive, got {c0_body_chord:.6f}")
    return float(float(c5_absolute_chord) / float(c0_body_chord))


def _c4_absolute_chord(c0_body_chord: float, c4_ratio_to_c0: float) -> float:
    return float(float(c0_body_chord) * float(c4_ratio_to_c0))


def _internal_c4_ratio(c0_body_chord: float, c4_absolute_chord: float) -> float:
    if c0_body_chord <= 0.0:
        raise ValueError(f"C0/body chord must be positive, got {c0_body_chord:.6f}")
    return float(float(c4_absolute_chord) / float(c0_body_chord))


def _transition_taper_ratio(c3_transition_chord: float, c4_outer_chord: float) -> float:
    c3 = float(c3_transition_chord)
    if c3 <= 0.0:
        raise ValueError(f"C3/transition chord must be positive, got {c3:.6f}")
    return float(float(c4_outer_chord) / c3)


def _c4_from_transition_taper(c3_transition_chord: float, taper_ratio: float) -> float:
    return float(float(c3_transition_chord) * float(taper_ratio))


def _absolute_chord(c0_body_chord: float, chord_ratio_to_c0: float) -> float:
    return float(float(c0_body_chord) * float(chord_ratio_to_c0))


def _internal_chord_ratio(c0_body_chord: float, absolute_chord: float) -> float:
    if c0_body_chord <= 0.0:
        raise ValueError(f"C0/body chord must be positive, got {c0_body_chord:.6f}")
    return float(float(absolute_chord) / float(c0_body_chord))


def _cta_public_flat_from_design(
    design: SectionedBWBDesignVariables,
    cta_design: Optional[SectionedBWBDesignVariables] = None,
) -> Dict[str, float]:
    flat = flatten_design(design)
    b1_fixed_m = float(cta_fixed_values(cta_design=cta_design or design)["b1_fixed_m"])
    total_span = float(flat["span"])
    internal_b2_ratio = float(flat["b2_span_ratio"])
    internal_b3_ratio = float(flat["b3_span_ratio"])
    c0_body_chord = float(flat["c1_root_chord"])
    c3_transition_chord = _absolute_chord(c0_body_chord, float(flat["c2_c1_ratio"]))
    c4_outer_chord = _c4_absolute_chord(c0_body_chord, float(flat["c3_c1_ratio"]))
    c5_wing_tip = _c5_absolute_chord(c0_body_chord, float(flat["c4_c1_ratio"]))
    b2_length_m = float(total_span * internal_b2_ratio)
    b3_length_m = float(total_span * internal_b3_ratio)

    flat["span"] = _wing_span_m(total_span, b1_fixed_m)
    flat["b2_span_ratio"] = _b2_wing_span_ratio(total_span, internal_b2_ratio, b1_fixed_m)
    flat["c4_c1_ratio"] = c5_wing_tip
    flat["c2_c1_ratio"] = c3_transition_chord
    flat["c3_c1_ratio"] = c4_outer_chord
    flat["c4_c3_ratio"] = _transition_taper_ratio(c3_transition_chord, c4_outer_chord)
    flat["s2_deg"] = _chord_sweep_from_leading_edge_sweep(
        le_sweep_deg=float(flat["s2_deg"]),
        chord_fraction=0.50,
        dy_m=b2_length_m,
        chord_in_m=c3_transition_chord,
        chord_out_m=c4_outer_chord,
    )
    flat["s3_deg"] = _chord_sweep_from_leading_edge_sweep(
        le_sweep_deg=float(flat["s3_deg"]),
        chord_fraction=0.25,
        dy_m=b3_length_m,
        chord_in_m=c4_outer_chord,
        chord_out_m=c5_wing_tip,
    )
    return flat


def _cta_design_space_seed(
    cta_design: Optional[SectionedBWBDesignVariables] = None,
) -> SectionedBWBDesignVariables:
    base = build_cta_design() if cta_design is None else cta_design
    return apply_cta_fixed_parameters(base, cta_design=base)


def _cta_cst_bounds(cta_flat: Dict[str, float]) -> Dict[str, Tuple[float, float]]:
    bounds: Dict[str, Tuple[float, float]] = {}
    for name, value in cta_flat.items():
        match = re.fullmatch(r"c[1-4]_(upper|lower)_cst_(\d+)", name)
        if match is None:
            continue
        coeff_idx = int(match.group(2))
        if coeff_idx >= len(_CTA_CST_BOUND_HALF_WIDTHS):
            continue
        half_width = float(_CTA_CST_BOUND_HALF_WIDTHS[coeff_idx])
        cta_value = float(value)
        bounds[name] = (
            max(0.0, cta_value - half_width),
            cta_value + half_width,
        )
    return bounds


def _cta_public_bounds(cta_design: SectionedBWBDesignVariables) -> Dict[str, Tuple[float, float]]:
    bounds = dict(SectionedBWBDesignVariables.default_bounds())
    fixed = cta_fixed_values(cta_design=cta_design)
    b1_fixed_m = float(fixed["b1_fixed_m"])
    span_lower, span_upper = bounds["span"]
    b2_lower_internal, b2_upper_internal = bounds["b2_span_ratio"]
    c0_lower, c0_upper = bounds["c1_root_chord"]
    c5_ratio_lower, c5_ratio_upper = bounds["c4_c1_ratio"]

    bounds["span"] = (
        _wing_span_m(span_lower, b1_fixed_m),
        _wing_span_m(span_upper, b1_fixed_m),
    )
    bounds["b2_span_ratio"] = (
        _b2_wing_span_ratio(span_upper, b2_lower_internal, b1_fixed_m),
        _b2_wing_span_ratio(span_lower, b2_upper_internal, b1_fixed_m),
    )
    bounds["c4_c1_ratio"] = (
        _c5_absolute_chord(c0_lower, c5_ratio_lower),
        _c5_absolute_chord(c0_upper, c5_ratio_upper),
    )
    c3_ratio_lower, c3_ratio_upper = bounds["c2_c1_ratio"]
    bounds["c2_c1_ratio"] = (
        _absolute_chord(c0_lower, c3_ratio_lower),
        _absolute_chord(c0_upper, c3_ratio_upper),
    )
    c4_ratio_lower, c4_ratio_upper = bounds["c3_c1_ratio"]
    bounds["c3_c1_ratio"] = (
        _c4_absolute_chord(c0_lower, c4_ratio_lower),
        _c4_absolute_chord(c0_upper, c4_ratio_upper),
    )
    cta_flat = _cta_public_flat_from_design(cta_design, cta_design=cta_design)
    bounds.update(CTA_PUBLIC_ACTIVE_BOUNDS)
    bounds.update(_cta_cst_bounds(cta_flat))
    for name in CTA_ACTIVE_VARIABLES:
        lower, upper = bounds[name]
        cta_value = float(cta_flat[name])
        bounds[name] = (min(lower, cta_value), max(upper, cta_value))
    return bounds


def _rename_section_tokens(text: str) -> str:
    mapping = {internal.upper(): public for internal, public in _SECTION_PUBLIC_LABELS.items()}
    return re.sub(r"\bC[1-4]\b", lambda match: mapping.get(match.group(0), match.group(0)), str(text))


@dataclass
class CTADesignSpace(DesignSpace):
    def cta_flat(self) -> Dict[str, float]:
        return _cta_public_flat_from_design(self.seed_design, cta_design=self.seed_design)

    def seed_flat(self) -> Dict[str, float]:
        return self.cta_flat()

    def to_design(
        self,
        sample: Dict[str, float],
    ) -> SectionedBWBDesignVariables:
        cta_flat = self.cta_flat()
        index_map = {
            name: idx
            for idx, name in enumerate(SectionedBWBDesignVariables.variable_names())
        }
        vector = self.seed_design.as_vector().copy()
        b1_fixed_m = float(cta_fixed_values(cta_design=self.seed_design)["b1_fixed_m"])

        for name in self.active_variables:
            if name not in sample:
                continue
            value = float(sample[name])
            if name == "span":
                vector[index_map["span"]] = _total_span_from_wing_span(value, b1_fixed_m)
                continue
            if name == "b2_span_ratio":
                span_here = float(vector[index_map["span"]])
                vector[index_map["b2_span_ratio"]] = _internal_b2_span_ratio(span_here, value, b1_fixed_m)
                continue
            if name == "c2_c1_ratio":
                c0_here = float(vector[index_map["c1_root_chord"]])
                vector[index_map["c2_c1_ratio"]] = _internal_chord_ratio(c0_here, value)
                continue
            if name == "c4_c3_ratio":
                c0_here = float(vector[index_map["c1_root_chord"]])
                c3_here = _absolute_chord(c0_here, float(vector[index_map["c2_c1_ratio"]]))
                c4_here = _c4_from_transition_taper(c3_here, value)
                vector[index_map["c3_c1_ratio"]] = _internal_c4_ratio(c0_here, c4_here)
                continue
            if name == "c4_c1_ratio":
                c0_here = float(vector[index_map["c1_root_chord"]])
                vector[index_map["c4_c1_ratio"]] = _internal_c5_ratio(c0_here, value)
                continue
            vector[index_map[name]] = value

        design = SectionedBWBDesignVariables.from_vector(vector)
        return apply_cta_fixed_parameters(design, cta_design=self.seed_design)

    def sample_designs(
        self,
        count: int,
        seed: int = 7,
        variation_scale: float = 0.25,
    ) -> List[SectionedBWBDesignVariables]:
        rng = np.random.default_rng(int(seed))
        cta_flat = self.cta_flat()
        index_map = {
            name: idx
            for idx, name in enumerate(SectionedBWBDesignVariables.variable_names())
        }
        cta_vector = self.seed_design.as_vector()
        internal_bounds = SectionedBWBDesignVariables.default_bounds()
        b1_fixed_m = float(cta_fixed_values(cta_design=self.seed_design)["b1_fixed_m"])
        sampled: List[SectionedBWBDesignVariables] = []

        for _ in range(int(count)):
            last_error: Optional[ValueError] = None
            for _attempt in range(200):
                vector = cta_vector.copy()
                for name in self.active_variables:
                    if name == "span":
                        wing_span_lower, wing_span_upper = self._local_bounds(name, cta_flat[name], variation_scale)
                        sampled_wing_span = rng.uniform(wing_span_lower, wing_span_upper)
                        vector[index_map["span"]] = _total_span_from_wing_span(sampled_wing_span, b1_fixed_m)
                        continue
                    if name == "b2_span_ratio":
                        span_here = float(vector[index_map["span"]])
                        global_lower, global_upper = internal_bounds["b2_span_ratio"]
                        ref_internal = float(self.seed_design.b2_span_ratio)
                        span_internal = global_upper - global_lower
                        local_lower_internal = max(global_lower, ref_internal - variation_scale * span_internal)
                        local_upper_internal = min(global_upper, ref_internal + variation_scale * span_internal)
                        sampled_public = rng.uniform(
                            _b2_wing_span_ratio(span_here, local_lower_internal, b1_fixed_m),
                            _b2_wing_span_ratio(span_here, local_upper_internal, b1_fixed_m),
                        )
                        vector[index_map[name]] = _internal_b2_span_ratio(span_here, sampled_public, b1_fixed_m)
                        continue
                    if name == "c2_c1_ratio":
                        c3_lower, c3_upper = self._local_bounds(name, cta_flat[name], variation_scale)
                        sampled_c3 = rng.uniform(c3_lower, c3_upper)
                        c0_here = float(vector[index_map["c1_root_chord"]])
                        vector[index_map["c2_c1_ratio"]] = _internal_chord_ratio(c0_here, sampled_c3)
                        continue
                    if name == "c4_c3_ratio":
                        taper_lower, taper_upper = self._local_bounds(name, cta_flat[name], variation_scale)
                        sampled_taper = rng.uniform(taper_lower, taper_upper)
                        c0_here = float(vector[index_map["c1_root_chord"]])
                        c3_here = _absolute_chord(c0_here, float(vector[index_map["c2_c1_ratio"]]))
                        c4_here = _c4_from_transition_taper(c3_here, sampled_taper)
                        vector[index_map["c3_c1_ratio"]] = _internal_c4_ratio(c0_here, c4_here)
                        continue
                    if name == "c4_c1_ratio":
                        c5_lower, c5_upper = self._local_bounds(name, cta_flat[name], variation_scale)
                        sampled_c5 = rng.uniform(c5_lower, c5_upper)
                        c0_here = float(vector[index_map["c1_root_chord"]])
                        vector[index_map["c4_c1_ratio"]] = _internal_c5_ratio(c0_here, sampled_c5)
                        continue

                    lower, upper = self._local_bounds(name, cta_flat[name], variation_scale)
                    vector[index_map[name]] = rng.uniform(lower, upper)

                try:
                    sampled_design = SectionedBWBDesignVariables.from_vector(vector)
                    sampled.append(
                        apply_cta_fixed_parameters(sampled_design, cta_design=self.seed_design)
                    )
                    break
                except ValueError as exc:
                    last_error = exc
            else:
                raise ValueError(
                    "Could not sample a feasible CTA design within the requested public bounds. "
                    "Check the C0/C3/C4/C5, transition taper ratio, and wing-span combinations "
                    "because all active chords must remain positive."
                ) from last_error
        return sampled


def build_cta_design_space(
    cta_design: Optional[SectionedBWBDesignVariables] = None,
) -> DesignSpace:
    """
    Build CTA AI design space with requested fixed/variable split:
    - fixed: B1 and S
    - derived: C1 and B3
    - variable: wing span, B2 fraction, C0, C3, transition taper ratio C4/C3, C5,
      public chord sweeps S1/S2, twists, and CST.
    """
    design = _cta_design_space_seed(cta_design=cta_design)
    bounds = _cta_public_bounds(design)
    return CTADesignSpace(
        preset_name="cta_ai_core",
        active_groups=("cta_core",),
        active_variables=CTA_ACTIVE_VARIABLES,
        seed_design=design,
        bounds=bounds,
    )


def cta_thickness_constraints() -> List[Dict[str, object]]:
    """Thickness constraints tracked for CTA designs but not sampled as active variables."""
    return [dict(item) for item in CTA_THICKNESS_CONSTRAINTS]


def cta_fixed_parameters(
    cta_design: Optional[SectionedBWBDesignVariables] = None,
) -> Dict[str, object]:
    """Fixed CTA values plus naming aliases for presentation."""
    design = build_cta_design() if cta_design is None else cta_design
    fixed = cta_fixed_values(cta_design=design)
    cta_view = _cta_public_flat_from_design(design, cta_design=design)
    fixed.update(
        {
            "sweep_naming": dict(SWEEP_NAME_TO_VARIABLE),
            "legacy_cta_view": {
                "S(legacy S1)": fixed["s_deg"],
                "S1(50%-chord)": cta_view["s2_deg"],
                "S2(25%-chord)": cta_view["s3_deg"],
            },
            "notes": (
                "CTA keeps B1 fixed in absolute meters while wing span can vary. "
                "B2 is defined as a fraction of wing span, with wing span = B2 + B3. "
                "C0/body chord, C3/transition-wing chord, transition taper ratio C4/C3, and C5/wing tip chord remain active. "
                "C1 is derived from fixed sweep S and straight TE(C0->C1). "
                "There is no public C2 parameter; the inboard TE blend from C1 to C3 uses a hidden helper point. "
                "Transition-wing sweep S1 is defined on the 50% chord line and outer-wing sweep S2 on the 25% chord line. "
                "TE(C3->C4) is no longer forced to stay straight. "
                "Twist stays constant from root to C3."
            ),
            "nomenclature": {
                "C0": "Body chord",
                "C1": "Derived from C0 with fixed S and straight TE(C0->C1)",
                "C3": "Active transition-wing chord",
                "C4": "Derived from transition taper ratio C4/C3",
                "C5": "Wing tip (active)",
                "B2": "Transition wing fraction of wing span",
                "wing_span": "B2 + B3",
                "med_3_TEswp": "Transition-wing trailing-edge sweep shaping TE(C3->C4)",
            },
        }
    )
    return fixed


def cta_parameter_metadata(
    cta_design: Optional[SectionedBWBDesignVariables] = None,
) -> List[Dict[str, object]]:
    """Return active CTA variables with metadata and renamed sweep labels."""
    design_space = build_cta_design_space(cta_design=cta_design)
    cta_flat = design_space.cta_flat()
    rows: List[Dict[str, object]] = []
    for name in design_space.active_variables:
        info = dict(parameter_info(name))
        if name in VARIABLE_TO_SWEEP_NAME:
            cta_name = VARIABLE_TO_SWEEP_NAME[name]
            info["display_name"] = f"Sweep {cta_name}"
            info["symbol"] = cta_name
            info["description"] = f"CTA sweep {cta_name} (mapped from {name})."
        if name == "span":
            info["display_name"] = "Wing span"
            info["symbol"] = "B2+B3"
            info["description"] = (
                "Outer semispan of the CTA planform. B1 remains fixed in meters, "
                "so the total semispan is B1 + wing span and wing span = B2 + B3."
            )
        if name == "c1_root_chord":
            info["display_name"] = "Body chord"
            info["symbol"] = "C0"
            info["description"] = (
                "CTA body chord. C1 is derived from C0 with fixed sweep S and straight TE(C0->C1)."
            )
        if name == "c2_c1_ratio":
            info["display_name"] = "Transition-wing chord"
            info["symbol"] = "C3"
            info["units"] = "m"
            info["normalization"] = "absolute"
            info["description"] = (
                "CTA transition-wing chord C3. This section remains anchored at y=8.041 m."
            )
        if name == "c4_c3_ratio":
            info["display_name"] = "Transition taper ratio"
            info["symbol"] = "C4/C3"
            info["units"] = "-"
            info["normalization"] = "ratio"
            info["description"] = "CTA transition-wing taper ratio driving C4 from C3."
        if name == "b2_span_ratio":
            info["display_name"] = "Transition wing fraction B2"
            info["symbol"] = "B2/(B2+B3)"
            info["normalization"] = "fraction of wing span"
            info["description"] = (
                "CTA transition-wing fraction measured against wing span. "
                "Wing span is the outer semispan excluding the fixed body segment B1, "
                "so wing span = B2 + B3."
            )
        if name == "c4_c1_ratio":
            info["display_name"] = "Wing tip chord"
            info["symbol"] = "C5"
            info["units"] = "m"
            info["normalization"] = "absolute"
            info["description"] = "CTA wing-tip chord C5."
        if name == "s2_deg":
            info["display_name"] = "Transition-wing 50% chord sweep"
            info["symbol"] = "S1"
            info["description"] = "CTA transition-wing sweep measured on the 50% chord line."
        if name == "s3_deg":
            info["display_name"] = "Outer-wing 25% chord sweep"
            info["symbol"] = "S2"
            info["description"] = "CTA outer-wing sweep measured on the 25% chord line."
        if name == "med_3_te_sweep_deg":
            info["display_name"] = "Transition-wing trailing-edge sweep"
            info["symbol"] = "med_3_TEswp"
            info["description"] = (
                "CTA transition-wing trailing-edge sweep shaping the TE between C3 and C4."
            )
        if name in _TWIST_PUBLIC_LABELS:
            public_label = _TWIST_PUBLIC_LABELS[name]
            info["display_name"] = f"Twist at {public_label}"
            info["symbol"] = f"twist_{public_label}"
            info["description"] = f"Geometric twist at CTA section {public_label}."
        if name.startswith("c") and "_cst_" in name:
            info["display_name"] = _rename_section_tokens(info["display_name"])
            info["symbol"] = _rename_section_tokens(info["symbol"])
            info["description"] = _rename_section_tokens(info["description"])
        lower, upper = design_space.bounds[name]
        rows.append(
            {
                "parameter": name,
                "display_name": info["display_name"],
                "symbol": info["symbol"],
                "units": info["units"],
                "normalization": info["normalization"],
                "description": info["description"],
                "cta_value": float(cta_flat[name]),
                "lower_bound": float(lower),
                "upper_bound": float(upper),
            }
        )
    return rows


def sample_cta_designs(
    count: int,
    seed: int = 7,
    variation_scale: float = 0.25,
    cta_design: Optional[SectionedBWBDesignVariables] = None,
) -> List[SectionedBWBDesignVariables]:
    """
    Sample CTA designs from the active AI design-space variables only.
    Fixed parameters remain locked by construction.
    """
    design_space = build_cta_design_space(cta_design=cta_design)
    return design_space.sample_designs(count=count, seed=seed, variation_scale=variation_scale)


def cta_design_space_summary(
    cta_design: Optional[SectionedBWBDesignVariables] = None,
) -> Dict[str, object]:
    """Compact serializable summary used by docs/scripts."""
    ds = build_cta_design_space(cta_design=cta_design)
    return {
        "preset_name": ds.preset_name,
        "active_variable_count": len(ds.active_variables),
        "active_variables": list(ds.active_variables),
        "thickness_constraints": cta_thickness_constraints(),
        "fixed_parameters": cta_fixed_parameters(cta_design=ds.seed_design),
        "cta_active_view": ds.cta_flat(),
        "cta_design": asdict(ds.seed_design),
    }
