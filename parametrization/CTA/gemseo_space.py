from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from parametrization.bwb.builder import prepare_geometry
from parametrization.bwb.design_variables import SectionedBWBDesignVariables

from .design_space import CTADesignSpace, build_cta_design_space
from .reference import to_cta_model_config


@dataclass(frozen=True)
class CTAGemseoVariableSpec:
    name: str
    cta_parameters: Tuple[str, ...]
    description: str
    units: str = "-"
    normalization: str = "-"

    @property
    def size(self) -> int:
        return len(self.cta_parameters)

    def lower_bounds(self, bounds: Dict[str, Tuple[float, float]]) -> np.ndarray:
        return np.asarray([bounds[name][0] for name in self.cta_parameters], dtype=float)

    def upper_bounds(self, bounds: Dict[str, Tuple[float, float]]) -> np.ndarray:
        return np.asarray([bounds[name][1] for name in self.cta_parameters], dtype=float)

    def reference_values(self, reference_flat: Dict[str, float]) -> np.ndarray:
        return np.asarray([reference_flat[name] for name in self.cta_parameters], dtype=float)


@dataclass
class CTAGemseoDesignSpaceAdapter:
    cta_space: CTADesignSpace
    variable_specs: Tuple[CTAGemseoVariableSpec, ...]
    gemseo_space: Optional[object]

    def reference_sample(self) -> Dict[str, np.ndarray]:
        reference_flat = self.cta_space.reference_flat()
        return {
            spec.name: spec.reference_values(reference_flat)
            for spec in self.variable_specs
        }

    def summary_rows(self) -> List[Dict[str, object]]:
        reference_flat = self.cta_space.reference_flat()
        rows: List[Dict[str, object]] = []
        for spec in self.variable_specs:
            rows.append(
                {
                    "gemseo_variable": spec.name,
                    "size": spec.size,
                    "cta_parameters": list(spec.cta_parameters),
                    "units": spec.units,
                    "normalization": spec.normalization,
                    "description": spec.description,
                    "lower_bounds": spec.lower_bounds(self.cta_space.bounds).tolist(),
                    "upper_bounds": spec.upper_bounds(self.cta_space.bounds).tolist(),
                    "reference": spec.reference_values(reference_flat).tolist(),
                }
            )
        return rows

    def flat_gemseo_variable_names(self) -> List[str]:
        names: List[str] = []
        for spec in self.variable_specs:
            if spec.size == 1:
                names.append(spec.name)
            else:
                names.extend(f"{spec.name}[{idx}]" for idx in range(spec.size))
        return names

    def to_cta_public_sample(
        self,
        sample: Mapping[str, Sequence[float]],
    ) -> Dict[str, float]:
        public_values = dict(self.cta_space.reference_flat())
        for spec in self.variable_specs:
            if spec.name not in sample:
                continue
            values = np.asarray(sample[spec.name], dtype=float).ravel()
            if values.size != spec.size:
                raise ValueError(
                    "GEMSEO variable %r expects size %d, received %d"
                    % (spec.name, spec.size, values.size)
                )
            for parameter_name, value in zip(spec.cta_parameters, values):
                public_values[parameter_name] = float(value)
        return public_values

    def to_project_design(
        self,
        sample: Mapping[str, Sequence[float]],
    ) -> SectionedBWBDesignVariables:
        return self.cta_space.to_design(self.to_cta_public_sample(sample))

    def flat_vector_to_gemseo_sample(
        self,
        vector: Sequence[float],
    ) -> Dict[str, np.ndarray]:
        values = np.asarray(vector, dtype=float).ravel()
        expected_size = sum(spec.size for spec in self.variable_specs)
        if values.size != expected_size:
            raise ValueError(
                "Expected flattened GEMSEO vector of size %d, received %d"
                % (expected_size, values.size)
            )
        sample: Dict[str, np.ndarray] = {}
        offset = 0
        for spec in self.variable_specs:
            sample[spec.name] = values[offset : offset + spec.size]
            offset += spec.size
        return sample

    def flat_vector_to_project_design(
        self,
        vector: Sequence[float],
    ) -> SectionedBWBDesignVariables:
        return self.to_project_design(self.flat_vector_to_gemseo_sample(vector))


@dataclass(frozen=True)
class CTASampleGeometryEvaluation:
    geometry_valid: bool
    error_message: str
    profile_generation_mode: str
    section_interpolation: str
    volume_enabled: bool
    volume_satisfied: bool
    enclosed_volume_m3: float
    required_volume_m3: float
    volume_margin_m3: float
    volume_ratio: float
    mean_cross_section_area_m2: float
    max_cross_section_area_m2: float
    max_cross_section_area_y: float
    min_inner_tc: float
    min_inner_tc_y: float
    min_inner_tc_xc: float


def evaluate_cta_gemseo_sample_geometry(
    adapter: CTAGemseoDesignSpaceAdapter,
    sample: Mapping[str, Sequence[float]],
    profile_generation_mode: str = "enforce_targets",
    required_volume_m3: Optional[float] = None,
    volume_span_samples: int = 161,
    interpolation_override: Optional[str] = None,
) -> CTASampleGeometryEvaluation:
    design = adapter.to_project_design(sample)
    config = to_cta_model_config(design)
    config.sections.profile_generation_mode = profile_generation_mode
    if interpolation_override is not None:
        config.sampling.section_interpolation = interpolation_override
        config.spanwise.twist_deg.interpolation = interpolation_override
        config.spanwise.camber_delta.interpolation = interpolation_override
    if required_volume_m3 is not None:
        config.volume.enabled = True
        config.volume.required_volume_m3 = float(required_volume_m3)
        config.volume.span_samples = int(volume_span_samples)
        config.volume.enforce_hard = False

    try:
        prepared = prepare_geometry(config)
    except Exception as exc:
        volume_enabled = required_volume_m3 is not None
        nan = float("nan")
        return CTASampleGeometryEvaluation(
            geometry_valid=False,
            error_message=str(exc),
            profile_generation_mode=profile_generation_mode,
            section_interpolation=config.sampling.section_interpolation,
            volume_enabled=volume_enabled,
            volume_satisfied=False,
            enclosed_volume_m3=nan,
            required_volume_m3=float(required_volume_m3) if required_volume_m3 is not None else nan,
            volume_margin_m3=nan,
            volume_ratio=nan,
            mean_cross_section_area_m2=nan,
            max_cross_section_area_m2=nan,
            max_cross_section_area_y=nan,
            min_inner_tc=nan,
            min_inner_tc_y=nan,
            min_inner_tc_xc=nan,
        )

    return CTASampleGeometryEvaluation(
        geometry_valid=True,
        error_message="",
        profile_generation_mode=profile_generation_mode,
        section_interpolation=config.sampling.section_interpolation,
        volume_enabled=prepared.volume.enabled,
        volume_satisfied=prepared.volume.satisfied,
        enclosed_volume_m3=prepared.volume.enclosed_volume_m3,
        required_volume_m3=prepared.volume.required_volume_m3,
        volume_margin_m3=prepared.volume.volume_margin_m3,
        volume_ratio=prepared.volume.volume_ratio,
        mean_cross_section_area_m2=prepared.volume.mean_cross_section_area_m2,
        max_cross_section_area_m2=prepared.volume.max_cross_section_area_m2,
        max_cross_section_area_y=prepared.volume.max_cross_section_area_y,
        min_inner_tc=prepared.validation.min_inner_tc,
        min_inner_tc_y=prepared.validation.min_inner_tc_y,
        min_inner_tc_xc=prepared.validation.min_inner_tc_xc,
    )


def _cta_gemseo_variable_specs() -> Tuple[CTAGemseoVariableSpec, ...]:
    return (
        CTAGemseoVariableSpec(
            name="wing_span",
            cta_parameters=("span",),
            units="m",
            normalization="absolute",
            description="CTA outer semispan, i.e. wing span = B2 + B3.",
        ),
        CTAGemseoVariableSpec(
            name="c0_body_chord",
            cta_parameters=("c1_root_chord",),
            units="m",
            normalization="absolute",
            description="CTA body chord C0.",
        ),
        CTAGemseoVariableSpec(
            name="c3_transition_chord",
            cta_parameters=("c2_c1_ratio",),
            units="m",
            normalization="absolute",
            description="CTA transition-wing chord C3.",
        ),
        CTAGemseoVariableSpec(
            name="transition_taper_ratio",
            cta_parameters=("c4_c3_ratio",),
            units="-",
            normalization="ratio",
            description="CTA transition-wing taper ratio C4/C3 driving the C4 chord.",
        ),
        CTAGemseoVariableSpec(
            name="b2_wing_fraction",
            cta_parameters=("b2_span_ratio",),
            units="-",
            normalization="fraction of wing span",
            description="CTA transition-wing fraction B2/(B2+B3).",
        ),
        CTAGemseoVariableSpec(
            name="c5_wing_tip_chord",
            cta_parameters=("c4_c1_ratio",),
            units="m",
            normalization="absolute",
            description="CTA wing-tip chord C5.",
        ),
        CTAGemseoVariableSpec(
            name="sweeps_deg",
            cta_parameters=("s2_deg", "s3_deg"),
            units="deg",
            normalization="angle",
            description="CTA public sweeps S1 (50% chord) and S2 (25% chord).",
        ),
        CTAGemseoVariableSpec(
            name="transition_te_sweep_deg",
            cta_parameters=("med_3_te_sweep_deg",),
            units="deg",
            normalization="angle",
            description="CTA transition-wing trailing-edge sweep (med_3_TEswp).",
        ),
        CTAGemseoVariableSpec(
            name="twist_deg",
            cta_parameters=("twist_c1_deg", "twist_c3_deg", "twist_c4_deg"),
            units="deg",
            normalization="angle",
            description="Twist values at CTA sections C0/C3, C4 and C5.",
        ),
        CTAGemseoVariableSpec(
            name="c0_upper_cst",
            cta_parameters=tuple(f"c1_upper_cst_{idx}" for idx in range(6)),
            units="-",
            normalization="Bernstein coefficient",
            description="Upper CST coefficients of CTA section C0.",
        ),
        CTAGemseoVariableSpec(
            name="c0_lower_cst",
            cta_parameters=tuple(f"c1_lower_cst_{idx}" for idx in range(6)),
            units="-",
            normalization="Bernstein coefficient",
            description="Lower CST coefficients of CTA section C0.",
        ),
        CTAGemseoVariableSpec(
            name="c3_upper_cst",
            cta_parameters=tuple(f"c2_upper_cst_{idx}" for idx in range(6)),
            units="-",
            normalization="Bernstein coefficient",
            description="Upper CST coefficients of CTA section C3.",
        ),
        CTAGemseoVariableSpec(
            name="c3_lower_cst",
            cta_parameters=tuple(f"c2_lower_cst_{idx}" for idx in range(6)),
            units="-",
            normalization="Bernstein coefficient",
            description="Lower CST coefficients of CTA section C3.",
        ),
        CTAGemseoVariableSpec(
            name="c4_upper_cst",
            cta_parameters=tuple(f"c3_upper_cst_{idx}" for idx in range(6)),
            units="-",
            normalization="Bernstein coefficient",
            description="Upper CST coefficients of CTA section C4.",
        ),
        CTAGemseoVariableSpec(
            name="c4_lower_cst",
            cta_parameters=tuple(f"c3_lower_cst_{idx}" for idx in range(6)),
            units="-",
            normalization="Bernstein coefficient",
            description="Lower CST coefficients of CTA section C4.",
        ),
        CTAGemseoVariableSpec(
            name="c5_upper_cst",
            cta_parameters=tuple(f"c4_upper_cst_{idx}" for idx in range(6)),
            units="-",
            normalization="Bernstein coefficient",
            description="Upper CST coefficients of CTA section C5.",
        ),
        CTAGemseoVariableSpec(
            name="c5_lower_cst",
            cta_parameters=tuple(f"c4_lower_cst_{idx}" for idx in range(6)),
            units="-",
            normalization="Bernstein coefficient",
            description="Lower CST coefficients of CTA section C5.",
        ),
    )


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


def _add_variable(
    design_space: object,
    spec: CTAGemseoVariableSpec,
    bounds: Dict[str, Tuple[float, float]],
    reference_flat: Dict[str, float],
) -> None:
    lower = spec.lower_bounds(bounds)
    upper = spec.upper_bounds(bounds)
    value = spec.reference_values(reference_flat)
    try:
        design_space.add_variable(
            name=spec.name,
            size=spec.size,
            type_="float",
            lower_bound=lower,
            upper_bound=upper,
            value=value,
        )
    except TypeError:
        design_space.add_variable(
            spec.name,
            spec.size,
            "float",
            lower,
            upper,
            value,
        )


def _build_cta_gemseo_adapter(
    reference_design: Optional[SectionedBWBDesignVariables] = None,
    instantiate_backend: bool = True,
) -> CTAGemseoDesignSpaceAdapter:
    cta_space = build_cta_design_space(reference_design=reference_design)
    if not isinstance(cta_space, CTADesignSpace):
        raise TypeError("build_cta_design_space() must return CTADesignSpace for CTA GEMSEO integration.")
    variable_specs = _cta_gemseo_variable_specs()
    reference_flat = cta_space.reference_flat()
    gemseo_space = None
    if instantiate_backend:
        gemseo_space = _new_gemseo_design_space()
        for spec in variable_specs:
            _add_variable(gemseo_space, spec, cta_space.bounds, reference_flat)
    return CTAGemseoDesignSpaceAdapter(
        cta_space=cta_space,
        variable_specs=variable_specs,
        gemseo_space=gemseo_space,
    )


def build_cta_gemseo_design_space(
    reference_design: Optional[SectionedBWBDesignVariables] = None,
) -> CTAGemseoDesignSpaceAdapter:
    return _build_cta_gemseo_adapter(reference_design=reference_design, instantiate_backend=True)


def build_cta_gemseo_design_space_definition(
    reference_design: Optional[SectionedBWBDesignVariables] = None,
) -> CTAGemseoDesignSpaceAdapter:
    return _build_cta_gemseo_adapter(reference_design=reference_design, instantiate_backend=False)


def available_cta_gemseo_doe_algorithms() -> Tuple[str, ...]:
    try:
        from gemseo import get_available_doe_algorithms
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "GEMSEO is not installed. Use Python 3.11 and install it with "
            "`python -m pip install 'gemseo>=6,<7'`."
        ) from exc

    return tuple(sorted(get_available_doe_algorithms()))


def sample_cta_gemseo_doe(
    adapter: CTAGemseoDesignSpaceAdapter,
    n_samples: int,
    algo_name: str = "LHS",
    seed: int = 7,
) -> List[np.ndarray]:
    if adapter.gemseo_space is None:
        raise ValueError(
            "The GEMSEO backend is not instantiated. "
            "Use build_cta_gemseo_design_space(...), not build_cta_gemseo_design_space_definition(...)."
        )
    if int(n_samples) <= 0:
        raise ValueError("n_samples must be a strictly positive integer.")

    try:
        from gemseo import create_scenario
        try:
            from gemseo.core.discipline import Discipline
        except ImportError:
            from gemseo import Discipline
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "GEMSEO is not installed. Use Python 3.11 and install it with "
            "`python -m pip install 'gemseo>=6,<7'`."
        ) from exc

    available_algorithms = available_cta_gemseo_doe_algorithms()
    if algo_name not in available_algorithms:
        raise ValueError(
            "DOE algorithm %r is not available in this GEMSEO installation. Available algorithms: %s"
            % (algo_name, ", ".join(available_algorithms))
        )

    class _CTAGemseoSamplingDiscipline(Discipline):
        def __init__(self, variable_specs: Tuple[CTAGemseoVariableSpec, ...]):
            super().__init__()
            self._variable_specs = variable_specs
            default_inputs = {
                spec.name: np.zeros(spec.size, dtype=float)
                for spec in variable_specs
            }
            self.input_grammar.update_from_data(default_inputs)
            self.output_grammar.update_from_data(
                {"sampling_objective": np.zeros(1, dtype=float)}
            )
            self.default_input_data = default_inputs

        def _run(self, input_data):
            objective = 0.0
            for spec in self._variable_specs:
                values = np.asarray(input_data[spec.name], dtype=float).ravel()
                objective += float(np.dot(values, values))
            return {"sampling_objective": np.asarray([objective], dtype=float)}

    scenario = create_scenario(
        disciplines=[_CTAGemseoSamplingDiscipline(adapter.variable_specs)],
        objective_name="sampling_objective",
        design_space=adapter.gemseo_space,
        scenario_type="DOE",
        name="cta_design_space_sampling",
        formulation_name="DisciplinaryOpt",
    )
    scenario.execute(
        algo_name=algo_name,
        n_samples=int(n_samples),
        seed=int(seed),
    )

    optimization_problem = getattr(scenario, "optimization_problem", None)
    if optimization_problem is None:
        optimization_problem = scenario.formulation.optimization_problem
    database = optimization_problem.database
    return [np.asarray(vector, dtype=float).ravel() for vector in database.get_x_vect_history()]
