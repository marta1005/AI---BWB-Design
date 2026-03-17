"""Separated workspace for future FFD and mesh-deformation development."""

from .ffd_box import (
    FFDBoxSpec,
    FFDBoxSummary,
    build_fullwing_ffd_box_points,
    build_reference_ffd_box,
    build_reference_surface_grids,
    build_reference_surface_pointset,
    build_semispan_ffd_box_points,
    load_dvgeometry,
    read_plot3d_ffd,
    reshape_surface_pointset,
)
from .dvgeo_tools import (
    DVGeoReferenceCase,
    apply_demo_local_deformation,
    create_reference_dvgeo_case,
    pointset_to_surface_grids,
    write_dvgeo_baseline_outputs,
)
