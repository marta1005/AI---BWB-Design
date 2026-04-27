from __future__ import annotations

from pathlib import Path
from typing import Optional

from parametrization.bwb.builder import PreparedGeometry, prepare_geometry
from parametrization.bwb.internal_volume_constraints import (
    CadReferenceFrame,
    InternalVolumeConstraintResult,
    InternalVolumeConstraintSet,
    evaluate_internal_volume_constraints,
    load_internal_volume_constraint_set,
)
from parametrization.bwb.specs import SectionedBWBModelConfig

from .case import build_cta_design, to_cta_model_config


CTA_INTERNAL_VOLUME_CONSTRAINTS_PATH = (
    Path(__file__).resolve().parent.parent.parent
    / "BWB B359-V0 - Internal Volume constraint - Set1 - BWB B359-V0 - Internal Volume constraint - Set1.csv"
)

# CAD / .geo frame offset relative to the internal CTA model frame.
# Current working alignment keeps the CAD and CTA model X origins coincident.
CTA_CAD_REFERENCE_FRAME = CadReferenceFrame(
    name="CTA CAD reference",
    offset_x_m=0.0,
    offset_y_m=0.0,
    offset_z_m=0.0,
    mirror_about_symmetry_plane=True,
)


def load_cta_internal_volume_constraint_set(
    path: str | Path | None = None,
    *,
    reference_frame: CadReferenceFrame = CTA_CAD_REFERENCE_FRAME,
) -> InternalVolumeConstraintSet:
    constraint_path = CTA_INTERNAL_VOLUME_CONSTRAINTS_PATH if path is None else Path(path)
    return load_internal_volume_constraint_set(
        constraint_path,
        reference_frame=reference_frame,
        name="CTA B359-V0 internal volume constraints",
    )


def evaluate_cta_internal_volume_constraints(
    *,
    prepared: Optional[PreparedGeometry] = None,
    config: Optional[SectionedBWBModelConfig] = None,
    path: str | Path | None = None,
    reference_frame: CadReferenceFrame = CTA_CAD_REFERENCE_FRAME,
    triangle_resolution: int = 18,
) -> InternalVolumeConstraintResult:
    if prepared is None:
        if config is None:
            config = to_cta_model_config(build_cta_design(), use_cta_anchor_twist=True)
        prepared = prepare_geometry(config)

    constraint_set = load_cta_internal_volume_constraint_set(
        path=path,
        reference_frame=reference_frame,
    )
    return evaluate_internal_volume_constraints(
        prepared,
        constraint_set,
        triangle_resolution=triangle_resolution,
    )
