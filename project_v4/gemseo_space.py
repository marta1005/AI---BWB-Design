from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from .design_space import DesignSpace as LocalDesignSpace
from .design_space import build_design_space, flatten_design
from .design_variables import SectionedBWBDesignVariables


@dataclass(frozen=True)
class GemseoVariableSpec:
    name: str
    project_parameters: Tuple[str, ...]
    description: str
    units: str = "-"
    normalization: str = "-"

    @property
    def size(self) -> int:
        return len(self.project_parameters)

    def lower_bounds(self, bounds: Dict[str, Tuple[float, float]]) -> np.ndarray:
        return np.asarray([bounds[name][0] for name in self.project_parameters], dtype=float)

    def upper_bounds(self, bounds: Dict[str, Tuple[float, float]]) -> np.ndarray:
        return np.asarray([bounds[name][1] for name in self.project_parameters], dtype=float)

    def reference_values(self, reference_flat: Dict[str, float]) -> np.ndarray:
        return np.asarray([reference_flat[name] for name in self.project_parameters], dtype=float)


@dataclass
class GemseoDesignSpaceAdapter:
    preset_name: str
    project_space: LocalDesignSpace
    variable_specs: Tuple[GemseoVariableSpec, ...]
    gemseo_space: object

    def reference_sample(self) -> Dict[str, np.ndarray]:
        reference_flat = flatten_design(self.project_space.reference_design)
        sample = {}
        for spec in self.variable_specs:
            sample[spec.name] = spec.reference_values(reference_flat)
        return sample

    def summary_rows(self) -> List[Dict[str, object]]:
        reference_flat = flatten_design(self.project_space.reference_design)
        rows = []
        for spec in self.variable_specs:
            rows.append(
                {
                    "gemseo_variable": spec.name,
                    "size": spec.size,
                    "project_parameters": list(spec.project_parameters),
                    "units": spec.units,
                    "normalization": spec.normalization,
                    "description": spec.description,
                    "lower_bounds": spec.lower_bounds(self.project_space.bounds).tolist(),
                    "upper_bounds": spec.upper_bounds(self.project_space.bounds).tolist(),
                    "reference": spec.reference_values(reference_flat).tolist(),
                }
            )
        return rows

    def to_project_design(
        self,
        sample: Mapping[str, Sequence[float]],
    ) -> SectionedBWBDesignVariables:
        reference_flat = flatten_design(self.project_space.reference_design)
        values = dict(reference_flat)
        for spec in self.variable_specs:
            if spec.name not in sample:
                continue
            spec_values = np.asarray(sample[spec.name], dtype=float).ravel()
            if spec_values.size != spec.size:
                raise ValueError(
                    "GEMSEO variable %r expects size %d, received %d"
                    % (spec.name, spec.size, spec_values.size)
                )
            if spec.name == "b_span_weights":
                ratios = _project_to_bounded_simplex(
                    spec_values,
                    spec.lower_bounds(self.project_space.bounds),
                    spec.upper_bounds(self.project_space.bounds),
                )
                for parameter_name, value in zip(spec.project_parameters, ratios):
                    values[parameter_name] = float(value)
                continue
            for parameter_name, value in zip(spec.project_parameters, spec_values):
                values[parameter_name] = float(value)

        vector = np.asarray(
            [values[name] for name in SectionedBWBDesignVariables.variable_names()],
            dtype=float,
        )
        design = SectionedBWBDesignVariables.from_vector(vector)
        design.validate()
        return design


def _project_to_bounded_simplex(
    values: Sequence[float],
    lower: np.ndarray,
    upper: np.ndarray,
    total: float = 1.0,
) -> np.ndarray:
    projected = np.asarray(values, dtype=float).ravel()
    if projected.size != lower.size:
        raise ValueError("Simplex projection received incompatible sizes")

    projected = np.maximum(projected, 1e-12)
    projected = projected / np.sum(projected) * float(total)

    free = np.ones(projected.size, dtype=bool)
    for _ in range(20):
        below = free & (projected < lower)
        above = free & (projected > upper)
        if not (np.any(below) or np.any(above)):
            break
        projected[below] = lower[below]
        projected[above] = upper[above]
        free[below | above] = False
        if not np.any(free):
            break
        remaining = float(total - np.sum(projected[~free]))
        weights = np.maximum(projected[free], 1e-12)
        projected[free] = weights / np.sum(weights) * remaining

    projected = np.clip(projected, lower, upper)
    for _ in range(20):
        delta = float(total - np.sum(projected))
        if abs(delta) <= 1e-12:
            break
        if delta > 0.0:
            capacity = upper - projected
        else:
            capacity = projected - lower
        movable = capacity > 1e-12
        if not np.any(movable):
            break
        step = delta / float(np.sum(capacity[movable]))
        projected[movable] += step * capacity[movable]
        projected = np.clip(projected, lower, upper)
    return projected


def _gemseo_variable_specs(active_groups: Tuple[str, ...]) -> Tuple[GemseoVariableSpec, ...]:
    specs: List[GemseoVariableSpec] = []
    for group_name in active_groups:
        if group_name == "topology":
            specs.extend(
                [
                    GemseoVariableSpec(
                        name="span",
                        project_parameters=("span",),
                        units="m",
                        normalization="absolute",
                        description="Half-span of the aircraft.",
                    ),
                    GemseoVariableSpec(
                        name="b_span_weights",
                        project_parameters=("b1_span_ratio", "b2_span_ratio", "b3_span_ratio"),
                        units="-",
                        normalization="normalized to sum 1.0 in adapter",
                        description="Positive topology weights mapped to B1/B2/B3 semi-span ratios.",
                    ),
                ]
            )
        elif group_name == "planform":
            specs.extend(
                [
                    GemseoVariableSpec(
                        name="le_root_x",
                        project_parameters=("le_root_x",),
                        units="m",
                        normalization="absolute",
                        description="Root leading-edge x-position.",
                    ),
                    GemseoVariableSpec(
                        name="c1_root_chord",
                        project_parameters=("c1_root_chord",),
                        units="m",
                        normalization="absolute",
                        description="Root chord.",
                    ),
                    GemseoVariableSpec(
                        name="chord_ratios",
                        project_parameters=("c2_c1_ratio", "c3_c1_ratio", "c4_c1_ratio"),
                        units="-",
                        normalization="ratio to C1",
                        description="Chord ratios of sections C2, C3 and C4.",
                    ),
                    GemseoVariableSpec(
                        name="sweeps_deg",
                        project_parameters=("s1_deg", "s2_deg", "s3_deg"),
                        units="deg",
                        normalization="angle",
                        description="Leading-edge sweep angles S1, S2 and S3.",
                    ),
                    GemseoVariableSpec(
                        name="nose_blend_y",
                        project_parameters=("nose_blend_y",),
                        units="m",
                        normalization="absolute",
                        description="Root nose blending length.",
                    ),
                ]
            )
        elif group_name == "class_function":
            specs.append(
                GemseoVariableSpec(
                    name="cst_class_exponents",
                    project_parameters=("cst_n1", "cst_n2"),
                    units="-",
                    normalization="class-function exponent",
                    description="Kulfan class-function exponents N1 and N2.",
                )
            )
        elif group_name == "attitude":
            specs.append(
                GemseoVariableSpec(
                    name="dihedral_deg",
                    project_parameters=("dihedral_deg",),
                    units="deg",
                    normalization="angle",
                    description="Global dihedral angle.",
                )
            )
        elif group_name == "twist":
            specs.append(
                GemseoVariableSpec(
                    name="twist_deg",
                    project_parameters=("twist_c1_deg", "twist_c2_deg", "twist_c3_deg", "twist_c4_deg"),
                    units="deg",
                    normalization="angle",
                    description="Twist angles at sections C1 to C4.",
                )
            )
        elif group_name == "camber_mode":
            specs.append(
                GemseoVariableSpec(
                    name="camber_delta",
                    project_parameters=("camber_c1", "camber_c2", "camber_c3", "camber_c4"),
                    units="-",
                    normalization="mode amplitude",
                    description="Additional camber mode at sections C1 to C4.",
                )
            )
        elif group_name == "thickness_targets":
            specs.extend(
                [
                    GemseoVariableSpec(
                        name="tc_max",
                        project_parameters=("c1_tc_max", "c2_tc_max", "c3_tc_max", "c4_tc_max"),
                        units="-",
                        normalization="t/c",
                        description="Maximum thickness ratios at sections C1 to C4.",
                    ),
                    GemseoVariableSpec(
                        name="x_tmax",
                        project_parameters=("c1_x_tmax", "c2_x_tmax", "c3_x_tmax", "c4_x_tmax"),
                        units="-",
                        normalization="x/c",
                        description="Location of maximum thickness at sections C1 to C4.",
                    ),
                    GemseoVariableSpec(
                        name="te_thickness",
                        project_parameters=(
                            "c1_te_thickness",
                            "c2_te_thickness",
                            "c3_te_thickness",
                            "c4_te_thickness",
                        ),
                        units="-",
                        normalization="t_TE/c_local",
                        description="Trailing-edge thickness ratios at sections C1 to C4.",
                    ),
                ]
            )
        elif group_name.startswith("cst_c"):
            section_name = group_name[-2:]
            specs.extend(
                [
                    GemseoVariableSpec(
                        name="%s_upper_cst" % section_name,
                        project_parameters=tuple("%s_upper_cst_%d" % (section_name, idx) for idx in range(6)),
                        units="-",
                        normalization="Bernstein coefficient",
                        description="Upper CST coefficients of section %s." % section_name.upper(),
                    ),
                    GemseoVariableSpec(
                        name="%s_lower_cst" % section_name,
                        project_parameters=tuple("%s_lower_cst_%d" % (section_name, idx) for idx in range(6)),
                        units="-",
                        normalization="Bernstein coefficient",
                        description="Lower CST coefficients of section %s." % section_name.upper(),
                    ),
                ]
            )
    return tuple(specs)


def _new_gemseo_design_space() -> object:
    try:
        from gemseo.algos.design_space import DesignSpace

        return DesignSpace()
    except ModuleNotFoundError:
        try:
            from gemseo import create_design_space

            return create_design_space()
        except ModuleNotFoundError as exc:
            raise ModuleNotFoundError(
                "GEMSEO is not installed. Use Python 3.11 and install it with "
                "`python -m pip install 'gemseo>=6,<7'`."
            ) from exc


def _add_variable(design_space: object, spec: GemseoVariableSpec, bounds: Dict[str, Tuple[float, float]], reference_flat: Dict[str, float]) -> None:
    lower = spec.lower_bounds(bounds)
    upper = spec.upper_bounds(bounds)
    value = spec.reference_values(reference_flat)
    try:
        design_space.add_variable(
            name=spec.name,
            size=spec.size,
            l_b=lower,
            u_b=upper,
            value=value,
            type_="float",
        )
    except TypeError:
        design_space.add_variable(
            spec.name,
            spec.size,
            lower,
            upper,
            value,
        )


def build_gemseo_design_space(
    preset_name: str = "presentation_core",
    reference_design: Optional[SectionedBWBDesignVariables] = None,
) -> GemseoDesignSpaceAdapter:
    project_space = build_design_space(preset_name, reference_design=reference_design)
    variable_specs = _gemseo_variable_specs(project_space.active_groups)
    reference_flat = flatten_design(project_space.reference_design)
    gemseo_space = _new_gemseo_design_space()
    for spec in variable_specs:
        _add_variable(gemseo_space, spec, project_space.bounds, reference_flat)
    return GemseoDesignSpaceAdapter(
        preset_name=preset_name,
        project_space=project_space,
        variable_specs=variable_specs,
        gemseo_space=gemseo_space,
    )
