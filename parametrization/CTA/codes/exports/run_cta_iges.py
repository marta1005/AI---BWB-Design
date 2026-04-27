from pathlib import Path
import sys

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
CTA_DIR = SCRIPT_DIR.parent.parent
REPO_ROOT = CTA_DIR.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from parametrization.CTA.case import build_cta_design, to_cta_model_config
from parametrization.bwb.builder import prepare_geometry
from parametrization.bwb.exporters import (
    append_xy_symmetry_frame_surfaces,
    build_xy_symmetry_frame_surfaces,
    build_pygeo_surface,
    build_te_height_scaled,
    write_station_airfoils,
)


def main() -> None:
    output_dir = CTA_DIR / "outputs" / "wing"
    profiles_dir = output_dir / "station_airfoils"
    iges_path = output_dir / "cta.igs"
    meshing_iges_path = output_dir / "cta_meshing_xy_frame.igs"
    frame_only_iges_path = output_dir / "cta_xy_frame_only.igs"
    output_dir.mkdir(parents=True, exist_ok=True)

    design = build_cta_design()
    config = to_cta_model_config(design, use_cta_anchor_twist=True)
    config.export.out_dir = profiles_dir
    config.export.iges_path = iges_path
    config.export.symmetric = False
    # For the CTA IGES, pass the fully resolved spanwise geometry to pyGeo
    # instead of a reduced anchor subset.
    config.sampling.airfoil_distribution_mode = "all"

    prepared = prepare_geometry(config)
    airfoil_list = write_station_airfoils(config, prepared.section_model, prepared.loft)
    te_height_scaled = build_te_height_scaled(prepared.section_model, prepared.loft)
    surface = build_pygeo_surface(config, prepared.loft, airfoil_list, te_height_scaled=te_height_scaled)
    surface.writeIGES(str(config.export.iges_path))
    meshing_surface = build_pygeo_surface(config, prepared.loft, airfoil_list, te_height_scaled=te_height_scaled)
    append_xy_symmetry_frame_surfaces(
        meshing_surface,
        prepared,
        x_min_m=-400.0,
        x_max_m=400.0,
        y_min_m=-400.0,
        y_max_m=400.0,
    )
    meshing_surface.writeIGES(str(meshing_iges_path))

    frame_only_surface = build_pygeo_surface(config, prepared.loft, airfoil_list, te_height_scaled=te_height_scaled)
    frame_only_surface.surfs = build_xy_symmetry_frame_surfaces(
        prepared,
        x_min_m=-400.0,
        x_max_m=400.0,
        y_min_m=-400.0,
        y_max_m=400.0,
    )
    frame_only_surface.nSurf = len(frame_only_surface.surfs)
    frame_only_surface.writeIGES(str(frame_only_iges_path))
    print(f"IGES exported: {config.export.iges_path}")
    print(f"Meshing IGES exported: {meshing_iges_path}")
    print(f"Frame-only IGES exported: {frame_only_iges_path}")
    print(f"Airfoil .dat written to: {config.export.out_dir}/")
    print(
        "Spanwise interpolation: "
        f"{prepared.section_model.interpolation_name} | "
        f"stations={prepared.loft.span_stations.size} | "
        f"anchors={config.topology.anchor_y_array.size} | "
        "pyGeo receives all resolved stations"
    )
    print(
        "Geometry validation: "
        f"min inner t/c={prepared.validation.min_inner_tc:.6f} "
        f"at y={prepared.validation.min_inner_tc_y:.6f}, "
        f"x/c={prepared.validation.min_inner_tc_xc:.6f} "
        f"over {prepared.validation.num_samples} samples"
    )
    print(f"CTA IGES written to: {iges_path}")
    print(f"CTA meshing IGES written to: {meshing_iges_path}")
    print(f"CTA frame-only IGES written to: {frame_only_iges_path}")
    print(f"CTA station profiles written to: {profiles_dir}")


if __name__ == "__main__":
    main()
