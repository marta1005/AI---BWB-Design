import argparse
import json
from pathlib import Path
import sys

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from project_v4.ffd.dvgeo_tools import (
    apply_demo_local_deformation,
    evaluate_demo_local_deformation,
    create_reference_dvgeo_case,
    pointset_to_surface_grids,
    write_dvgeo_baseline_outputs,
)
from project_v4.ffd.ffd_box import FFDBoxSpec, read_plot3d_ffd
from project_v4.ffd.visualization import plot_ffd_3d, plot_ffd_with_surface, write_ffd_deformation_gif


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Deform the reference FFD box with DVGeometry and export visualization.")
    parser.add_argument("--full-wing", action="store_true", help="Use the full-wing FFD box instead of semispan.")
    parser.add_argument("--surface-spanwise", type=int, default=19, help="Surface spanwise samples per semispan.")
    parser.add_argument("--surface-chordwise", type=int, default=81, help="Surface chordwise samples.")
    parser.add_argument("--dv-scale", type=float, default=0.68, help="Magnitude of the sample local DV perturbation.")
    parser.add_argument("--n-frames", type=int, default=28, help="Number of frames in the deformation GIF.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = SCRIPT_DIR.parent / "example_outputs" / "reference_ffd_dvgeo"
    output_dir.mkdir(parents=True, exist_ok=True)

    spec = FFDBoxSpec()

    ffd_path = output_dir / ("reference_fullwing_ffd.xyz" if args.full_wing else "reference_semispan_ffd.xyz")
    case, dvgeo = create_reference_dvgeo_case(
        output_ffd_path=ffd_path,
        ffd_spec=spec,
        full_wing=bool(args.full_wing),
        surface_n_spanwise=int(args.surface_spanwise),
        surface_n_chordwise=int(args.surface_chordwise),
    )
    write_dvgeo_baseline_outputs(
        dvgeo,
        case.pointset_name,
        output_dir,
        stem="baseline",
    )
    stem = "deformed_fullwing" if args.full_wing else "deformed_semispan"
    deformed_pointset, deformed_ffd_points, target_shape_values = apply_demo_local_deformation(
        dvgeo,
        case.pointset_name,
        output_dir,
        stem=stem,
        dv_scale=float(args.dv_scale),
    )
    deformed_upper, deformed_lower = pointset_to_surface_grids(deformed_pointset, case.upper_grid)

    figure_path = output_dir / ("reference_fullwing_ffd_deformation.png" if args.full_wing else "reference_semispan_ffd_deformation.png")
    figure_3d_path = output_dir / ("reference_fullwing_ffd_deformation_3d.png" if args.full_wing else "reference_semispan_ffd_deformation_3d.png")
    gif_path = output_dir / ("reference_fullwing_ffd_deformation.gif" if args.full_wing else "reference_semispan_ffd_deformation.gif")
    baseline_ffd_points = read_plot3d_ffd(ffd_path)

    plot_ffd_with_surface(
        ffd_points=baseline_ffd_points,
        upper_grid=case.upper_grid,
        lower_grid=case.lower_grid,
        out_path=figure_path,
        title="Reference BWB DVGeometry deformation preview",
        baseline_pointset=case.pointset,
        deformed_ffd_points=deformed_ffd_points,
        deformed_pointset=deformed_pointset,
        deformed_upper_grid=deformed_upper,
        deformed_lower_grid=deformed_lower,
    )
    plot_ffd_3d(
        ffd_points=baseline_ffd_points,
        upper_grid=case.upper_grid,
        lower_grid=case.lower_grid,
        out_path=figure_3d_path,
        title="Reference BWB DVGeometry deformation in 3D",
        baseline_pointset=case.pointset,
        deformed_ffd_points=deformed_ffd_points,
        deformed_pointset=deformed_pointset,
        deformed_upper_grid=deformed_upper,
        deformed_lower_grid=deformed_lower,
    )

    animation_frames = []
    for frame_idx, fraction in enumerate(np.linspace(0.0, 1.0, int(args.n_frames), dtype=float)):
        frame_pointset, frame_ffd_points = evaluate_demo_local_deformation(
            dvgeo,
            pointset_name=case.pointset_name,
            target_shape_values=target_shape_values,
            fraction=float(fraction),
        )
        frame_upper, frame_lower = pointset_to_surface_grids(frame_pointset, case.upper_grid)
        animation_frames.append(
            {
                "baseline_ffd_points": baseline_ffd_points,
                "baseline_upper_grid": case.upper_grid,
                "baseline_lower_grid": case.lower_grid,
                "baseline_pointset": case.pointset,
                "deformed_ffd_points": frame_ffd_points,
                "deformed_pointset": frame_pointset,
                "deformed_upper_grid": frame_upper,
                "deformed_lower_grid": frame_lower,
                "frame_label": f"deformation = {fraction:.2f}",
            }
        )
    write_ffd_deformation_gif(
        out_path=gif_path,
        title="Reference BWB DVGeometry deformation",
        frames=animation_frames,
        duration_ms=135,
    )

    metadata = {
        "baseline_ffd": str(case.ffd_summary.xyz_path),
        "deformed_ffd": str(output_dir / f"{stem}_ffd.xyz"),
        "figure_path": str(figure_path),
        "figure_3d_path": str(figure_3d_path),
        "gif_path": str(gif_path),
        "full_wing": bool(args.full_wing),
        "surface_spanwise": int(args.surface_spanwise),
        "surface_chordwise": int(args.surface_chordwise),
        "dv_scale": float(args.dv_scale),
        "n_frames": int(args.n_frames),
    }
    metadata_path = output_dir / "dvgeo_deformation_summary.json"
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    print(f"Baseline FFD written to: {case.ffd_summary.xyz_path}")
    print(f"Deformed FFD written to: {output_dir / f'{stem}_ffd.xyz'}")
    print(f"Visualization written to: {figure_path}")
    print(f"3D visualization written to: {figure_3d_path}")
    print(f"GIF written to: {gif_path}")
    print(f"Metadata written to: {metadata_path}")


if __name__ == "__main__":
    main()
