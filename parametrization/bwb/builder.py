from dataclasses import dataclass

import numpy as np

from .exporters import build_pygeo_surface, build_te_height_scaled, write_station_airfoils
from .planform import SectionedPlanform, build_sectioned_bwb_planform
from .sections import SectionModel, build_section_model
from .specs import SectionedBWBModelConfig
from .spanwise_laws import ResolvedSpanwiseLaws, resolve_spanwise_laws, vertical_offsets
from .validation import (
    LoftDefinition,
    ValidationSummary,
    build_validation_stations,
    validate_loft_definition,
    validate_section_geometry,
)
from .volume import VolumeConstraintSummary, evaluate_volume_constraint


@dataclass
class PreparedGeometry:
    planform: SectionedPlanform
    section_model: SectionModel
    spanwise_laws: ResolvedSpanwiseLaws
    loft: LoftDefinition
    validation: ValidationSummary
    volume: VolumeConstraintSummary


def build_span_stations(config: SectionedBWBModelConfig) -> np.ndarray:
    base_stations = np.linspace(0.0, config.topology.span, config.sampling.num_base_stations)
    planform_helper_stations = np.unique(
        np.concatenate(
            [
                config.planform.leading_edge_points(config.topology)[:, 1],
                config.planform.trailing_edge_points(config.topology)[:, 1],
            ]
        )
    )
    span_stations = np.unique(
        np.round(
            np.concatenate([base_stations, config.topology.anchor_y_array, planform_helper_stations]),
            decimals=12,
        )
    )
    return span_stations.astype(float)


def build_loft_definition(
    config: SectionedBWBModelConfig,
    planform: SectionedPlanform,
    laws: ResolvedSpanwiseLaws,
) -> LoftDefinition:
    span_stations = build_span_stations(config)
    span_z = span_stations.copy()
    leading_edge_x = np.array([planform.le_x(yy) for yy in span_stations], dtype=float)
    trailing_edge_x = np.array([planform.te_x(yy) for yy in span_stations], dtype=float)
    chord = trailing_edge_x - leading_edge_x
    vertical_y = vertical_offsets(config.topology, config.spanwise, span_z)
    twist_deg = np.array([laws.twist_deg(yy) for yy in span_stations], dtype=float)
    offset = np.zeros((span_stations.size, 2), dtype=float)
    return LoftDefinition(
        span_stations=span_stations,
        leading_edge_x=leading_edge_x,
        chord=chord,
        vertical_y=vertical_y,
        span_z=span_z,
        twist_deg=twist_deg,
        offset=offset,
    )


def prepare_geometry(config: SectionedBWBModelConfig) -> PreparedGeometry:
    config.validate()
    planform = build_sectioned_bwb_planform(config.topology, config.planform)
    laws = resolve_spanwise_laws(config)
    section_model = build_section_model(config, laws)
    loft = build_loft_definition(config, planform, laws)
    validate_loft_definition(loft)
    validation_stations = build_validation_stations(config.topology, loft)
    validation = validate_section_geometry(section_model, validation_stations)
    volume = evaluate_volume_constraint(config, planform, section_model)
    if config.volume.enforce_hard and not volume.satisfied:
        raise ValueError(
            "Volume constraint violated: "
            f"volume={volume.enclosed_volume_m3:.3f} m^3, "
            f"required={volume.required_volume_m3:.3f} m^3, "
            f"margin={volume.volume_margin_m3:.3f} m^3"
        )
    return PreparedGeometry(
        planform=planform,
        section_model=section_model,
        spanwise_laws=laws,
        loft=loft,
        validation=validation,
        volume=volume,
    )


def build_surface(config: SectionedBWBModelConfig):
    prepared = prepare_geometry(config)
    airfoil_list = write_station_airfoils(config, prepared.section_model, prepared.loft)
    te_height_scaled = build_te_height_scaled(prepared.section_model, prepared.loft)
    surface = build_pygeo_surface(config, prepared.loft, airfoil_list, te_height_scaled=te_height_scaled)
    return surface, prepared


def export_iges(config: SectionedBWBModelConfig) -> PreparedGeometry:
    surface, prepared = build_surface(config)
    surface.writeIGES(str(config.export.iges_path))
    print(f"IGES exported: {config.export.iges_path}")
    print(f"Airfoil .dat written to: {config.export.out_dir}/")
    print(
        "Spanwise interpolation: "
        f"{prepared.section_model.interpolation_name} | "
        f"stations={prepared.loft.span_stations.size} | "
        f"anchors={config.topology.anchor_y_array.size}"
    )
    print(
        "Geometry validation: "
        f"min inner t/c={prepared.validation.min_inner_tc:.6f} "
        f"at y={prepared.validation.min_inner_tc_y:.6f}, "
        f"x/c={prepared.validation.min_inner_tc_xc:.6f} "
        f"over {prepared.validation.num_samples} samples"
    )
    if prepared.volume.enabled:
        print(
            "Volume constraint: "
            f"satisfied={prepared.volume.satisfied} | "
            f"V={prepared.volume.enclosed_volume_m3:.2f}/{prepared.volume.required_volume_m3:.2f} m^3"
        )
    return prepared
