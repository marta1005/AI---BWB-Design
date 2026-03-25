"""CTA design-space definition for AI sampling and optimization."""

from dataclasses import asdict, dataclass
import re
from typing import Dict, List, Optional, Tuple

import numpy as np

from parametrization.bwb.design_space import DesignSpace, flatten_design, parameter_info
from parametrization.bwb.design_variables import SectionedBWBDesignVariables

from .reference import (
    SWEEP_NAME_TO_VARIABLE,
    VARIABLE_TO_SWEEP_NAME,
    apply_cta_fixed_parameters,
    build_reference_design,
    cta_fixed_values,
)

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
    "c3_c1_ratio",
    "b2_span_ratio",
    "c4_c1_ratio",
    "s2_deg",  # CTA label S1
    "s3_deg",  # CTA label S2
    "twist_c1_deg",
    "twist_c3_deg",
    "twist_c4_deg",
    *_cst_names(1),
    *_cst_names(2),
    *_cst_names(3),
    *_cst_names(4),
)

CTA_PUBLIC_ACTIVE_BOUNDS: Dict[str, Tuple[float, float]] = {
    "span": (30.0, 35.0),  # wing span = B2 + B3
    "c1_root_chord": (37.0, 45.0),  # C0
    "c2_c1_ratio": (13.0, 16.0),  # C3 in absolute meters at the CTA layer
    "c3_c1_ratio": (6.8, 9.8),  # C4 in absolute meters at the CTA layer
    "b2_span_ratio": (0.14, 0.20),  # B2 / (B2 + B3)
    "c4_c1_ratio": (0.80, 1.80),  # C5 in absolute meters at the CTA layer
    "s2_deg": (45.0, 66.0),  # S1
    "s3_deg": (27.0, 40.0),  # S2
    "twist_c1_deg": (0.2, 2.0),
    "twist_c3_deg": (0.1, 1.5),
    "twist_c4_deg": (0.1, 1.0),
}

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


def _absolute_chord(c0_body_chord: float, chord_ratio_to_c0: float) -> float:
    return float(float(c0_body_chord) * float(chord_ratio_to_c0))


def _internal_chord_ratio(c0_body_chord: float, absolute_chord: float) -> float:
    if c0_body_chord <= 0.0:
        raise ValueError(f"C0/body chord must be positive, got {c0_body_chord:.6f}")
    return float(float(absolute_chord) / float(c0_body_chord))


def _cta_public_bounds(reference_design: SectionedBWBDesignVariables) -> Dict[str, Tuple[float, float]]:
    bounds = dict(SectionedBWBDesignVariables.default_bounds())
    fixed = cta_fixed_values(reference_design=reference_design)
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
    reference_flat = flatten_design(reference_design)
    reference_flat["span"] = _wing_span_m(
        float(reference_flat["span"]),
        b1_fixed_m,
    )
    reference_flat["b2_span_ratio"] = _b2_wing_span_ratio(
        float(flatten_design(reference_design)["span"]),
        float(reference_flat["b2_span_ratio"]),
        b1_fixed_m,
    )
    reference_flat["c4_c1_ratio"] = _c5_absolute_chord(
        float(reference_flat["c1_root_chord"]),
        float(reference_flat["c4_c1_ratio"]),
    )
    reference_flat["c2_c1_ratio"] = _absolute_chord(
        float(reference_flat["c1_root_chord"]),
        float(reference_flat["c2_c1_ratio"]),
    )
    reference_flat["c3_c1_ratio"] = _c4_absolute_chord(
        float(reference_flat["c1_root_chord"]),
        float(reference_flat["c3_c1_ratio"]),
    )
    bounds.update(CTA_PUBLIC_ACTIVE_BOUNDS)
    for name in CTA_ACTIVE_VARIABLES:
        lower, upper = bounds[name]
        reference_value = float(reference_flat[name])
        bounds[name] = (min(lower, reference_value), max(upper, reference_value))
    return bounds


def _rename_section_tokens(text: str) -> str:
    mapping = {internal.upper(): public for internal, public in _SECTION_PUBLIC_LABELS.items()}
    return re.sub(r"\bC[1-4]\b", lambda match: mapping.get(match.group(0), match.group(0)), str(text))


@dataclass
class CTADesignSpace(DesignSpace):
    def reference_flat(self) -> Dict[str, float]:
        flat = flatten_design(self.reference_design)
        b1_fixed_m = float(cta_fixed_values(reference_design=self.reference_design)["b1_fixed_m"])
        total_span = float(flat["span"])
        flat["span"] = _wing_span_m(total_span, b1_fixed_m)
        flat["b2_span_ratio"] = _b2_wing_span_ratio(
            total_span,
            float(flat["b2_span_ratio"]),
            b1_fixed_m,
        )
        flat["c4_c1_ratio"] = _c5_absolute_chord(
            float(flat["c1_root_chord"]),
            float(flat["c4_c1_ratio"]),
        )
        flat["c2_c1_ratio"] = _absolute_chord(
            float(flat["c1_root_chord"]),
            float(flat["c2_c1_ratio"]),
        )
        flat["c3_c1_ratio"] = _c4_absolute_chord(
            float(flat["c1_root_chord"]),
            float(flat["c3_c1_ratio"]),
        )
        return flat

    def to_design(
        self,
        sample: Dict[str, float],
    ) -> SectionedBWBDesignVariables:
        reference_flat = self.reference_flat()
        index_map = {
            name: idx
            for idx, name in enumerate(SectionedBWBDesignVariables.variable_names())
        }
        vector = self.reference_design.as_vector().copy()
        b1_fixed_m = float(cta_fixed_values(reference_design=self.reference_design)["b1_fixed_m"])

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
            if name == "c3_c1_ratio":
                c0_here = float(vector[index_map["c1_root_chord"]])
                vector[index_map["c3_c1_ratio"]] = _internal_c4_ratio(c0_here, value)
                continue
            if name == "c4_c1_ratio":
                c0_here = float(vector[index_map["c1_root_chord"]])
                vector[index_map["c4_c1_ratio"]] = _internal_c5_ratio(c0_here, value)
                continue
            vector[index_map[name]] = value

        design = SectionedBWBDesignVariables.from_vector(vector)
        return apply_cta_fixed_parameters(design, reference_design=self.reference_design)

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
        internal_bounds = SectionedBWBDesignVariables.default_bounds()
        b1_fixed_m = float(cta_fixed_values(reference_design=self.reference_design)["b1_fixed_m"])
        sampled: List[SectionedBWBDesignVariables] = []

        for _ in range(int(count)):
            last_error: Optional[ValueError] = None
            for _attempt in range(200):
                vector = reference_vector.copy()
                for name in self.active_variables:
                    if name == "span":
                        wing_span_lower, wing_span_upper = self._local_bounds(name, reference_flat[name], variation_scale)
                        sampled_wing_span = rng.uniform(wing_span_lower, wing_span_upper)
                        vector[index_map["span"]] = _total_span_from_wing_span(sampled_wing_span, b1_fixed_m)
                        continue
                    if name == "b2_span_ratio":
                        span_here = float(vector[index_map["span"]])
                        global_lower, global_upper = internal_bounds["b2_span_ratio"]
                        ref_internal = float(self.reference_design.b2_span_ratio)
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
                        c3_lower, c3_upper = self._local_bounds(name, reference_flat[name], variation_scale)
                        sampled_c3 = rng.uniform(c3_lower, c3_upper)
                        c0_here = float(vector[index_map["c1_root_chord"]])
                        vector[index_map["c2_c1_ratio"]] = _internal_chord_ratio(c0_here, sampled_c3)
                        continue
                    if name == "c3_c1_ratio":
                        c4_lower, c4_upper = self._local_bounds(name, reference_flat[name], variation_scale)
                        sampled_c4 = rng.uniform(c4_lower, c4_upper)
                        c0_here = float(vector[index_map["c1_root_chord"]])
                        vector[index_map["c3_c1_ratio"]] = _internal_c4_ratio(c0_here, sampled_c4)
                        continue
                    if name == "c4_c1_ratio":
                        c5_lower, c5_upper = self._local_bounds(name, reference_flat[name], variation_scale)
                        sampled_c5 = rng.uniform(c5_lower, c5_upper)
                        c0_here = float(vector[index_map["c1_root_chord"]])
                        vector[index_map["c4_c1_ratio"]] = _internal_c5_ratio(c0_here, sampled_c5)
                        continue

                    lower, upper = self._local_bounds(name, reference_flat[name], variation_scale)
                    vector[index_map[name]] = rng.uniform(lower, upper)

                try:
                    sampled_design = SectionedBWBDesignVariables.from_vector(vector)
                    sampled.append(
                        apply_cta_fixed_parameters(sampled_design, reference_design=self.reference_design)
                    )
                    break
                except ValueError as exc:
                    last_error = exc
            else:
                raise ValueError(
                    "Could not sample a feasible CTA design within the requested public bounds. "
                    "Check the C0/C3/C4/C5 and wing-span combinations because all active chords must remain positive."
                ) from last_error
        return sampled


def build_cta_design_space(
    reference_design: Optional[SectionedBWBDesignVariables] = None,
) -> DesignSpace:
    """
    Build CTA AI design space with requested fixed/variable split:
    - fixed: B1 and S
    - derived: C1 and B3
    - variable: wing span, B2 fraction, C0, C3, C4, C5, following sweeps, twists, and CST.
    """
    reference = build_reference_design() if reference_design is None else reference_design
    reference = apply_cta_fixed_parameters(reference, reference_design=reference)
    bounds = _cta_public_bounds(reference)
    return CTADesignSpace(
        preset_name="cta_reference_ai_core",
        active_groups=("cta_core",),
        active_variables=CTA_ACTIVE_VARIABLES,
        reference_design=reference,
        bounds=bounds,
    )


def cta_fixed_parameters(
    reference_design: Optional[SectionedBWBDesignVariables] = None,
) -> Dict[str, object]:
    """Fixed CTA values plus naming aliases for presentation."""
    reference = build_reference_design() if reference_design is None else reference_design
    fixed = cta_fixed_values(reference_design=reference)
    fixed.update(
        {
            "sweep_naming": dict(SWEEP_NAME_TO_VARIABLE),
            "legacy_reference": {
                "S(legacy S1)": fixed["s_deg"],
                "S1(legacy S2)": float(reference.s2_deg),
                "S2(legacy S3)": float(reference.s3_deg),
            },
            "notes": (
                "CTA keeps B1 fixed in absolute meters while wing span can vary. "
                "B2 is defined as a fraction of wing span, with wing span = B2 + B3. "
                "C0/body chord, C3/transition-wing chord, C4/outer-wing chord, and C5/wing tip chord remain active. "
                "C1 is derived from fixed sweep S and straight TE(C0->C1). "
                "There is no public C2 parameter; the inboard TE blend from C1 to C3 uses a hidden helper point. "
                "TE(C3->C4) is no longer forced to stay straight. "
                "Twist stays constant from root to C3."
            ),
            "nomenclature": {
                "C0": "Body chord",
                "C1": "Derived from C0 with fixed S and straight TE(C0->C1)",
                "C3": "Active transition-wing chord",
                "C4": "Active outer-wing chord",
                "C5": "Wing tip (active)",
                "B2": "Transition wing fraction of wing span",
                "wing_span": "B2 + B3",
            },
        }
    )
    return fixed


def cta_parameter_metadata(
    reference_design: Optional[SectionedBWBDesignVariables] = None,
) -> List[Dict[str, object]]:
    """Return active CTA variables with metadata and renamed sweep labels."""
    design_space = build_cta_design_space(reference_design=reference_design)
    reference_flat = design_space.reference_flat()
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
        if name == "c3_c1_ratio":
            info["display_name"] = "Outer-wing chord"
            info["symbol"] = "C4"
            info["units"] = "m"
            info["normalization"] = "absolute"
            info["description"] = "CTA outer-wing chord C4."
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
                "reference": float(reference_flat[name]),
                "lower_bound": float(lower),
                "upper_bound": float(upper),
            }
        )
    return rows


def sample_cta_designs(
    count: int,
    seed: int = 7,
    variation_scale: float = 0.25,
    reference_design: Optional[SectionedBWBDesignVariables] = None,
) -> List[SectionedBWBDesignVariables]:
    """
    Sample CTA designs from the active AI design-space variables only.
    Fixed parameters remain locked by construction.
    """
    design_space = build_cta_design_space(reference_design=reference_design)
    return design_space.sample_designs(count=count, seed=seed, variation_scale=variation_scale)


def cta_design_space_summary(
    reference_design: Optional[SectionedBWBDesignVariables] = None,
) -> Dict[str, object]:
    """Compact serializable summary used by docs/scripts."""
    ds = build_cta_design_space(reference_design=reference_design)
    return {
        "preset_name": ds.preset_name,
        "active_variable_count": len(ds.active_variables),
        "active_variables": list(ds.active_variables),
        "fixed_parameters": cta_fixed_parameters(reference_design=ds.reference_design),
        "reference_active_view": ds.reference_flat(),
        "reference_design": asdict(ds.reference_design),
    }
