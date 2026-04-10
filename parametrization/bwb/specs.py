from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Optional, Tuple

import numpy as np

from .topology import SectionedBWBTopologySpec


@dataclass
class SectionCSTSpec:
    upper_coeffs: Tuple[float, ...]
    lower_coeffs: Tuple[float, ...]
    tc_max: float
    x_tmax: float
    te_thickness: float

    def validate(self, prefix: str, expected_size: int) -> None:
        if len(self.upper_coeffs) != expected_size:
            raise ValueError(
                f"{prefix}.upper_coeffs must have size {expected_size}, "
                f"got {len(self.upper_coeffs)}"
            )
        if len(self.lower_coeffs) != expected_size:
            raise ValueError(
                f"{prefix}.lower_coeffs must have size {expected_size}, "
                f"got {len(self.lower_coeffs)}"
            )
        if not all(np.isfinite(value) for value in self.upper_coeffs):
            raise ValueError(f"{prefix}.upper_coeffs must be finite")
        if not all(np.isfinite(value) for value in self.lower_coeffs):
            raise ValueError(f"{prefix}.lower_coeffs must be finite")
        if self.tc_max <= 0.0:
            raise ValueError(f"{prefix}.tc_max must be positive, got {self.tc_max}")
        if not (0.0 < self.x_tmax < 1.0):
            raise ValueError(f"{prefix}.x_tmax must lie in (0, 1), got {self.x_tmax}")
        if self.te_thickness < 0.0:
            raise ValueError(f"{prefix}.te_thickness must be non-negative, got {self.te_thickness}")


@dataclass(frozen=True)
class SectionProfileRelationSpec:
    shape_source_index: Optional[int] = None
    upper_source_index: Optional[int] = None
    lower_source_index: Optional[int] = None
    te_source_index: Optional[int] = None

    def validate(self, prefix: str, section_index: int, section_count: int) -> None:
        for label, source_index in (
            ("shape_source_index", self.shape_source_index),
            ("upper_source_index", self.upper_source_index),
            ("lower_source_index", self.lower_source_index),
            ("te_source_index", self.te_source_index),
        ):
            if source_index is None:
                continue
            if not (0 <= source_index < section_count):
                raise ValueError(
                    f"{prefix}.{label} must lie inside [0, {section_count - 1}], got {source_index}"
                )
            if source_index == section_index:
                raise ValueError(f"{prefix}.{label} cannot reference itself (section {section_index})")
            if source_index > section_index:
                raise ValueError(
                    f"{prefix}.{label} must reference an earlier section to avoid cyclic dependencies, "
                    f"got source={source_index}, section={section_index}"
                )


@dataclass
class PlanformSpec:
    le_root_x: float = 0.0
    c1_root_chord: float = 40.0
    c2_c1_ratio: float = 0.70
    c3_c1_ratio: float = 0.23
    c4_c1_ratio: float = 0.08
    s1_deg: float = 55.0
    s2_deg: float = 50.0
    s3_deg: float = 32.0
    med_3_te_sweep_deg: float = 0.0
    med_3_te_helper_fraction: float = 0.72
    te_c1_span_fraction: float = 5.0 / 8.0
    te_inboard_blend_fraction: float = 0.96
    te_inboard_blend_dx: float = 0.0
    te_inboard_radius_factor: float = 1.0
    te_outer_blend_fraction: float = 0.10
    body_le_fixed_points: Tuple[Tuple[float, float], ...] = ()
    continuity_order: int = 2
    blend_fraction: float = 0.18
    min_linear_core_fraction: float = 0.75
    te_blend_fraction: Optional[float] = None
    te_min_linear_core_fraction: Optional[float] = None
    le_linear_start_index: Optional[int] = None
    te_spline_bridge: Optional[Tuple[int, int]] = None
    symmetry_blend_y: float = 2.50
    te_exact_segments: Tuple[int, ...] = (0, 3)

    def validate(self) -> None:
        if self.c1_root_chord <= 0.0:
            raise ValueError(f"c1_root_chord must be positive, got {self.c1_root_chord}")
        ratios = (self.c2_c1_ratio, self.c3_c1_ratio, self.c4_c1_ratio)
        if any(value <= 0.0 for value in ratios):
            raise ValueError(f"Chord ratios must be positive, got {ratios}")
        if self.continuity_order not in {1, 2}:
            raise ValueError(
                f"continuity_order must be 1 or 2 for smooth planforms, got {self.continuity_order}"
            )
        if not (0.0 <= self.blend_fraction <= 0.45):
            raise ValueError(
                f"blend_fraction must lie in [0, 0.45], got {self.blend_fraction}"
            )
        if not (0.0 <= self.min_linear_core_fraction < 1.0):
            raise ValueError(
                "min_linear_core_fraction must lie in [0, 1), "
                f"got {self.min_linear_core_fraction}"
            )
        if self.te_blend_fraction is not None and not (0.0 <= self.te_blend_fraction <= 0.45):
            raise ValueError(
                f"te_blend_fraction must lie in [0, 0.45], got {self.te_blend_fraction}"
            )
        if self.te_min_linear_core_fraction is not None and not (
            0.0 <= self.te_min_linear_core_fraction < 1.0
        ):
            raise ValueError(
                "te_min_linear_core_fraction must lie in [0, 1), "
                f"got {self.te_min_linear_core_fraction}"
            )
        if self.le_linear_start_index is not None and self.le_linear_start_index < 1:
            raise ValueError(
                f"le_linear_start_index must be >= 1 when provided, got {self.le_linear_start_index}"
            )
        if self.te_spline_bridge is not None:
            if len(self.te_spline_bridge) != 2:
                raise ValueError(
                    f"te_spline_bridge must contain exactly 2 indices, got {self.te_spline_bridge}"
                )
            start_idx, end_idx = self.te_spline_bridge
            if not (0 < start_idx < end_idx):
                raise ValueError(
                    "te_spline_bridge must satisfy 0 < start_idx < end_idx, "
                    f"got {self.te_spline_bridge}"
                )
        if self.symmetry_blend_y < 0.0:
            raise ValueError(
                f"symmetry_blend_y must be non-negative, got {self.symmetry_blend_y}"
            )
        if any(index < 0 for index in self.te_exact_segments):
            raise ValueError(
                f"te_exact_segments must contain non-negative segment indices, got {self.te_exact_segments}"
            )
        if not (0.0 < self.te_c1_span_fraction < self.te_inboard_blend_fraction < 1.0):
            raise ValueError(
                "CTA TE fractions must satisfy 0 < te_c1_span_fraction < "
                "te_inboard_blend_fraction < 1, got "
                f"{(self.te_c1_span_fraction, self.te_inboard_blend_fraction)}"
            )
        if self.te_inboard_radius_factor <= 0.0:
            raise ValueError(
                "te_inboard_radius_factor must be positive, "
                f"got {self.te_inboard_radius_factor}"
            )
        if not (0.0 <= self.te_outer_blend_fraction < 1.0):
            raise ValueError(
                f"te_outer_blend_fraction must lie in [0, 1), got {self.te_outer_blend_fraction}"
            )
        if self.body_le_fixed_points:
            helper_y = [float(point[1]) for point in self.body_le_fixed_points]
            if any(value <= 0.0 for value in helper_y):
                raise ValueError(
                    "body_le_fixed_points must have positive spanwise y values, "
                    f"got {self.body_le_fixed_points}"
                )
            if any(left >= right for left, right in zip(helper_y[:-1], helper_y[1:])):
                raise ValueError(
                    "body_le_fixed_points must be strictly increasing in spanwise y, "
                    f"got {self.body_le_fixed_points}"
                )
        for label, value in (("s1_deg", self.s1_deg), ("s2_deg", self.s2_deg), ("s3_deg", self.s3_deg)):
            if not (0.0 < value < 85.0):
                raise ValueError(f"{label} must lie in (0, 85), got {value}")
        if not (-85.0 < self.med_3_te_sweep_deg < 85.0):
            raise ValueError(
                f"med_3_te_sweep_deg must lie in (-85, 85), got {self.med_3_te_sweep_deg}"
            )
        if not (0.0 < self.med_3_te_helper_fraction < 1.0):
            raise ValueError(
                "med_3_te_helper_fraction must lie in (0, 1), "
                f"got {self.med_3_te_helper_fraction}"
            )

    def section_chords(self) -> np.ndarray:
        return np.asarray(
            (
                self.c1_root_chord,
                self.c2_c1_ratio * self.c1_root_chord,
                self.c3_c1_ratio * self.c1_root_chord,
                self.c4_c1_ratio * self.c1_root_chord,
            ),
            dtype=float,
        )

    def leading_edge_x_sections(self, topology: SectionedBWBTopologySpec) -> np.ndarray:
        y_sections = topology.y_sections_array
        segment_dy = np.diff(y_sections)
        segment_sweeps = np.deg2rad([self.s1_deg, self.s2_deg, self.s3_deg])
        x_sections = np.zeros_like(y_sections)
        x_sections[0] = float(self.le_root_x)
        for idx, (dy, sweep_rad) in enumerate(zip(segment_dy, segment_sweeps), start=1):
            x_sections[idx] = x_sections[idx - 1] + np.tan(float(sweep_rad)) * float(dy)
        return x_sections.astype(float)

    def trailing_edge_x_sections(self, topology: SectionedBWBTopologySpec) -> np.ndarray:
        return self.leading_edge_x_sections(topology) + self.section_chords()

    def leading_edge_points(self, topology: SectionedBWBTopologySpec) -> np.ndarray:
        y_sections = topology.y_sections_array
        le_sections = self.leading_edge_x_sections(topology)
        points = [[float(le_sections[0]), float(y_sections[0]), 0.0]]
        for x_helper, y_helper in self.body_le_fixed_points:
            points.append([float(x_helper), float(y_helper), 0.0])
        for x_section, y_section in zip(le_sections[1:], y_sections[1:]):
            points.append([float(x_section), float(y_section), 0.0])
        return np.asarray(points, dtype=float)

    def trailing_edge_points(self, topology: SectionedBWBTopologySpec) -> np.ndarray:
        y_sections = topology.y_sections_array
        te_sections = self.trailing_edge_x_sections(topology)

        y_c3 = float(y_sections[1])
        if y_c3 <= 1e-12:
            raise ValueError("CTA-style TE auxiliary points require the first main section to lie away from the root")

        y_c1 = float(self.te_c1_span_fraction * y_c3)
        y_inboard_blend = float(self.te_inboard_blend_fraction * y_c3)

        te_root = float(te_sections[0])
        te_c3 = float(te_sections[1])
        te_c4 = float(te_sections[2])
        te_tip = float(te_sections[3])
        y_c4 = float(y_sections[2])
        y_tip = float(y_sections[3])

        # CTA reference:
        # - C0 -> C1 is a straight vertical TE segment.
        # - C1 -> C3 is a marked blend.
        # - C3 is a real geometric break.
        #
        # We keep one hidden helper point inside the C1 -> C3 blend so the
        # planform can stay smooth there without exposing a public "C2".
        te_c1 = te_root
        # The inboard hidden helper shapes the marked C1->C3 blend. By default
        # it lies on the public C1->C3 secant, but CTA can override it with a
        # fixed aft offset measured from C3 to reproduce the Airbus sketch.
        if abs(float(self.te_inboard_blend_dx)) > 1e-12:
            te_inboard_blend = te_c3 + float(self.te_inboard_blend_dx)
            te_inboard_blend = float(min(max(te_c3, te_inboard_blend), te_c1))
        else:
            blend_ratio = (y_inboard_blend - y_c1) / max(y_c3 - y_c1, 1e-12)
            te_inboard_blend = te_c1 + blend_ratio * (te_c3 - te_c1)
        points = [
            [te_root, 0.0, 0.0],
            [te_c1, y_c1, 0.0],
            [te_inboard_blend, y_inboard_blend, 0.0],
            [te_c3, float(y_sections[1]), 0.0],
        ]

        if abs(float(self.med_3_te_sweep_deg)) > 1e-12:
            y_med3 = y_c3 + float(self.med_3_te_helper_fraction) * (y_c4 - y_c3)
            te_med3 = te_c3 + np.tan(np.deg2rad(float(self.med_3_te_sweep_deg))) * (y_med3 - y_c3)
            points.append([te_med3, y_med3, 0.0])

        points.append([te_c4, float(y_sections[2]), 0.0])

        if self.te_outer_blend_fraction > 1e-12:
            # Optional short outer-wing helper point for a smooth C4->C5 TE.
            y_outer_blend = y_c4 + self.te_outer_blend_fraction * (y_tip - y_c4)
            outer_ratio = (y_outer_blend - y_c4) / max(y_tip - y_c4, 1e-12)
            te_outer_blend = te_c4 + outer_ratio * (te_tip - te_c4)
            points.append([te_outer_blend, y_outer_blend, 0.0])

        points.append([te_tip, float(y_sections[3]), 0.0])
        return np.asarray(points, dtype=float)

    def validate_with_topology(self, topology: SectionedBWBTopologySpec) -> None:
        le_x = self.leading_edge_x_sections(topology)
        te_x = self.trailing_edge_x_sections(topology)
        chord = te_x - le_x
        if np.any(chord <= 0.0):
            raise ValueError(
                "non-positive control-section chord detected; "
                f"LE={tuple(le_x)}, TE={tuple(te_x)}, chord={tuple(chord)}"
            )
        segment_count = self.trailing_edge_points(topology).shape[0] - 1
        if any(index >= segment_count for index in self.te_exact_segments):
            raise ValueError(
                f"te_exact_segments must lie inside [0, {segment_count - 1}], "
                f"got {self.te_exact_segments}"
            )
        if self.body_le_fixed_points:
            first_main_section_y = float(topology.y_sections_array[1])
            if float(self.body_le_fixed_points[-1][1]) >= first_main_section_y:
                raise ValueError(
                    "body_le_fixed_points must remain inside the center-body span before C3, "
                    f"got {self.body_le_fixed_points} with first main section at y={first_main_section_y:.6f}"
                )
        if self.le_linear_start_index is not None:
            point_count = self.leading_edge_points(topology).shape[0]
            if self.le_linear_start_index >= point_count - 1:
                raise ValueError(
                    f"le_linear_start_index must be < {point_count - 1}, got {self.le_linear_start_index}"
                )
        if self.te_spline_bridge is not None:
            point_count = self.trailing_edge_points(topology).shape[0]
            start_idx, end_idx = self.te_spline_bridge
            if end_idx >= point_count - 1:
                raise ValueError(
                    f"te_spline_bridge end index must be < {point_count - 1}, got {self.te_spline_bridge}"
                )


@dataclass
class AnchoredSpanwiseLaw:
    section_indices: Tuple[int, ...]
    values: Tuple[float, ...]
    interpolation: str = "pyspline"

    def validate(self, topology: SectionedBWBTopologySpec, label: str) -> None:
        if len(self.section_indices) != len(self.values):
            raise ValueError(
                f"{label} must have the same number of section_indices and values, "
                f"got {len(self.section_indices)} and {len(self.values)}"
            )
        if len(self.section_indices) < 2:
            raise ValueError(f"{label} must have at least 2 anchors")
        if self.interpolation not in {"linear", "pchip", "cubic", "pyspline"}:
            raise ValueError(
                f"{label}.interpolation must be 'linear', 'pchip', 'cubic' or 'pyspline', "
                f"got {self.interpolation!r}"
            )
        section_count = topology.anchor_y_array.size
        if any(index < 0 or index >= section_count for index in self.section_indices):
            raise ValueError(
                f"{label}.section_indices must lie inside [0, {section_count - 1}], "
                f"got {self.section_indices}"
            )
        if any(left >= right for left, right in zip(self.section_indices[:-1], self.section_indices[1:])):
            raise ValueError(f"{label}.section_indices must be strictly increasing, got {self.section_indices}")


@dataclass
class SectionFamilySpec:
    x_tc_window: Tuple[float, float] = (0.15, 0.65)
    x_valid_window: Tuple[float, float] = (0.02, 0.98)
    cst_degree: int = 5
    n1: float = 0.50
    n2: float = 1.0
    shared_leading_edge: bool = False
    profile_generation_mode: str = "cst_only"
    camber_mode_center: float = 3.5
    camber_mode_width: float = 2.0
    c1_spec: SectionCSTSpec = field(
        default_factory=lambda: SectionCSTSpec(
            upper_coeffs=(0.30, 0.36, 0.32, 0.24, 0.14, 0.05),
            lower_coeffs=(0.13, 0.10, 0.07, 0.045, 0.025, 0.010),
            tc_max=0.22,
            x_tmax=0.33,
            te_thickness=0.002,
        )
    )
    c2_spec: SectionCSTSpec = field(
        default_factory=lambda: SectionCSTSpec(
            upper_coeffs=(0.24, 0.30, 0.27, 0.20, 0.12, 0.04),
            lower_coeffs=(0.11, 0.085, 0.060, 0.040, 0.020, 0.008),
            tc_max=0.18,
            x_tmax=0.31,
            te_thickness=0.002,
        )
    )
    c3_spec: SectionCSTSpec = field(
        default_factory=lambda: SectionCSTSpec(
            upper_coeffs=(0.16, 0.21, 0.19, 0.14, 0.085, 0.03),
            lower_coeffs=(0.075, 0.055, 0.040, 0.025, 0.012, 0.005),
            tc_max=0.11,
            x_tmax=0.28,
            te_thickness=0.0015,
        )
    )
    c4_spec: SectionCSTSpec = field(
        default_factory=lambda: SectionCSTSpec(
            upper_coeffs=(0.12, 0.15, 0.13, 0.095, 0.060, 0.02),
            lower_coeffs=(0.055, 0.040, 0.028, 0.018, 0.009, 0.003),
            tc_max=0.08,
            x_tmax=0.24,
            te_thickness=0.001,
        )
    )
    section_specs_override: Optional[Tuple[SectionCSTSpec, ...]] = None
    profile_relations: Tuple[SectionProfileRelationSpec, ...] = field(
        default_factory=lambda: (
            SectionProfileRelationSpec(),
            SectionProfileRelationSpec(),
            SectionProfileRelationSpec(),
            SectionProfileRelationSpec(),
        )
    )

    @property
    def ncoeff(self) -> int:
        return self.cst_degree + 1

    @property
    def total_coeff_count(self) -> int:
        return 2 * self.ncoeff

    @property
    def base_section_specs(self) -> Tuple[SectionCSTSpec, ...]:
        if self.section_specs_override is not None:
            return tuple(self.section_specs_override)
        return (self.c1_spec, self.c2_spec, self.c3_spec, self.c4_spec)

    @property
    def section_specs(self) -> Tuple[SectionCSTSpec, ...]:
        resolved = [replace(spec) for spec in self.base_section_specs]
        if len(self.profile_relations) != len(resolved):
            raise ValueError(
                "profile_relations must match the number of control sections, "
                f"got {len(self.profile_relations)} relations for {len(resolved)} sections"
            )

        for idx, relation in enumerate(self.profile_relations):
            resolved_spec = resolved[idx]
            if relation.shape_source_index is not None:
                source = resolved[relation.shape_source_index]
                resolved_spec = replace(
                    resolved_spec,
                    upper_coeffs=source.upper_coeffs,
                    lower_coeffs=source.lower_coeffs,
                    tc_max=source.tc_max,
                    x_tmax=source.x_tmax,
                    te_thickness=source.te_thickness,
                )
            if relation.upper_source_index is not None:
                resolved_spec = replace(
                    resolved_spec,
                    upper_coeffs=resolved[relation.upper_source_index].upper_coeffs,
                )
            if relation.lower_source_index is not None:
                resolved_spec = replace(
                    resolved_spec,
                    lower_coeffs=resolved[relation.lower_source_index].lower_coeffs,
                )
            if relation.te_source_index is not None:
                resolved_spec = replace(
                    resolved_spec,
                    te_thickness=resolved[relation.te_source_index].te_thickness,
                )
            resolved[idx] = resolved_spec

        return tuple(resolved)

    def validate(self) -> None:
        if self.cst_degree < 2:
            raise ValueError(f"cst_degree must be at least 2, got {self.cst_degree}")
        if self.n1 <= 0.0 or self.n2 <= 0.0:
            raise ValueError(f"n1 and n2 must be positive, got {(self.n1, self.n2)}")
        if self.profile_generation_mode not in {"cst_only", "enforce_targets"}:
            raise ValueError(
                "profile_generation_mode must be 'cst_only' or 'enforce_targets', "
                f"got {self.profile_generation_mode!r}"
            )
        if self.camber_mode_width <= 0.0:
            raise ValueError(
                f"camber_mode_width must be positive, got {self.camber_mode_width}"
            )
        x0, x1 = self.x_tc_window
        if not (0.0 <= x0 < x1 <= 1.0):
            raise ValueError(f"x_tc_window must satisfy 0 <= x0 < x1 <= 1, got {self.x_tc_window}")
        x0, x1 = self.x_valid_window
        if not (0.0 <= x0 < x1 <= 1.0):
            raise ValueError(
                f"x_valid_window must satisfy 0 <= x0 < x1 <= 1, got {self.x_valid_window}"
            )
        expected_size = self.ncoeff
        base_specs = self.base_section_specs
        if len(self.profile_relations) != len(base_specs):
            raise ValueError(
                "profile_relations must match the number of control sections, "
                f"got {len(self.profile_relations)} relations for {len(base_specs)} sections"
            )
        for idx, relation in enumerate(self.profile_relations):
            relation.validate(f"profile_relations[{idx}]", idx, len(base_specs))
        for idx, spec in enumerate(base_specs, start=1):
            spec.validate(f"section_spec[{idx - 1}]", expected_size)
            if not (self.x_tc_window[0] <= spec.x_tmax <= self.x_tc_window[1]):
                raise ValueError(
                    f"section_spec[{idx - 1}].x_tmax={spec.x_tmax} must lie inside x_tc_window={self.x_tc_window}"
                )
        for idx, spec in enumerate(self.section_specs, start=1):
            spec.validate(f"resolved_section_spec[{idx - 1}]", expected_size)


@dataclass
class SpanwiseLawSpec:
    dihedral_deg: float = 0.0
    vertical_offset_m: float = 0.03
    vertical_offset_z: Optional[AnchoredSpanwiseLaw] = None
    twist_deg: AnchoredSpanwiseLaw = field(
        default_factory=lambda: AnchoredSpanwiseLaw(
            section_indices=(0, 1, 2, 3),
            values=(1.0, 1.0, 0.8, 0.6),
            interpolation="pyspline",
        )
    )
    camber_delta: AnchoredSpanwiseLaw = field(
        default_factory=lambda: AnchoredSpanwiseLaw(
            section_indices=(0, 1, 2, 3),
            values=(0.0, 0.0, 0.0, 0.0),
            interpolation="pyspline",
        )
    )

    def validate(self, topology: SectionedBWBTopologySpec) -> None:
        if not (-30.0 <= self.dihedral_deg <= 30.0):
            raise ValueError(f"dihedral_deg must lie in [-30, 30], got {self.dihedral_deg}")
        if not np.isfinite(self.vertical_offset_m):
            raise ValueError(f"vertical_offset_m must be finite, got {self.vertical_offset_m}")
        if self.vertical_offset_z is not None:
            self.vertical_offset_z.validate(topology, "vertical_offset_z")
        self.twist_deg.validate(topology, "twist_deg")
        self.camber_delta.validate(topology, "camber_delta")


@dataclass
class SamplingSpec:
    num_airfoil_points: int = 241
    num_base_stations: int = 41
    section_curve_n_ctl: int = 41
    section_interpolation: str = "pyspline"
    airfoil_distribution_mode: str = "all"
    k_span: int = 4

    def validate(self) -> None:
        if self.num_airfoil_points < 5:
            raise ValueError("num_airfoil_points must be at least 5")
        if self.num_base_stations < 2:
            raise ValueError("num_base_stations must be at least 2")
        if self.section_curve_n_ctl < 4:
            raise ValueError("section_curve_n_ctl must be at least 4")
        if self.section_interpolation not in {"cubic", "pchip", "linear", "pyspline"}:
            raise ValueError(
                "section_interpolation must be 'linear', 'cubic', 'pchip' or 'pyspline', "
                f"got {self.section_interpolation!r}"
            )
        if self.airfoil_distribution_mode not in {"all", "anchors"}:
            raise ValueError(
                "airfoil_distribution_mode must be 'all' or 'anchors', "
                f"got {self.airfoil_distribution_mode!r}"
            )
        if self.k_span not in {2, 3, 4}:
            raise ValueError(f"k_span must be one of {{2, 3, 4}}, got {self.k_span}")


@dataclass
class ExportSpec:
    out_dir: Path = Path("profiles_out")
    iges_path: Path = Path("bwb_from_cst_pygeo.igs")
    tip_style: str = "rounded"
    blunt_te: bool = True
    symmetric: bool = False

    def validate(self) -> None:
        if self.tip_style not in {"rounded", "pinched"}:
            raise ValueError(f"tip_style must be 'rounded' or 'pinched', got {self.tip_style!r}")


@dataclass
class VolumeConstraintSpec:
    enabled: bool = False
    required_volume_m3: float = 0.0
    span_samples: int = 161
    enforce_hard: bool = False

    def validate(self, topology: SectionedBWBTopologySpec) -> None:
        if not self.enabled:
            return
        if self.required_volume_m3 <= 0.0:
            raise ValueError(
                f"required_volume_m3 must be positive when volume constraint is enabled, got {self.required_volume_m3}"
            )
        if self.span_samples < 11:
            raise ValueError(f"span_samples must be >= 11, got {self.span_samples}")
        if topology.span <= 0.0:
            raise ValueError(f"topology span must be positive, got {topology.span}")


@dataclass
class SectionedBWBModelConfig:
    topology: SectionedBWBTopologySpec = field(default_factory=SectionedBWBTopologySpec)
    planform: PlanformSpec = field(default_factory=PlanformSpec)
    sections: SectionFamilySpec = field(default_factory=SectionFamilySpec)
    spanwise: SpanwiseLawSpec = field(default_factory=SpanwiseLawSpec)
    sampling: SamplingSpec = field(default_factory=SamplingSpec)
    export: ExportSpec = field(default_factory=ExportSpec)
    volume: VolumeConstraintSpec = field(default_factory=VolumeConstraintSpec)

    def validate(self) -> None:
        self.topology.validate()
        self.planform.validate()
        self.planform.validate_with_topology(self.topology)
        self.sections.validate()
        self.spanwise.validate(self.topology)
        self.sampling.validate()
        self.export.validate()
        self.volume.validate(self.topology)
        section_count = self.topology.anchor_y_array.size
        if len(self.sections.section_specs) != section_count:
            raise ValueError(
                "section CST specs must match topology.anchor_y, "
                f"got {len(self.sections.section_specs)} specs for {section_count} sections"
            )
