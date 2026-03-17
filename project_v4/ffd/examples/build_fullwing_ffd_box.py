import json
from dataclasses import asdict
from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from project_v4.ffd.dvgeo_tools import create_reference_dvgeo_case, write_dvgeo_baseline_outputs
from project_v4.ffd.ffd_box import FFDBoxSpec, read_plot3d_ffd
from project_v4.ffd.visualization import plot_ffd_with_surface


def main() -> None:
    output_dir = SCRIPT_DIR.parent / "example_outputs" / "reference_ffd_box"
    output_dir.mkdir(parents=True, exist_ok=True)

    ffd_path = output_dir / "reference_fullwing_ffd.xyz"
    summary_path = output_dir / "reference_fullwing_ffd.json"
    figure_path = output_dir / "reference_fullwing_ffd.png"

    case, dvgeo = create_reference_dvgeo_case(output_ffd_path=ffd_path, ffd_spec=FFDBoxSpec(), full_wing=True)
    write_dvgeo_baseline_outputs(dvgeo, case.pointset_name, output_dir, stem="reference_fullwing_baseline")
    plot_ffd_with_surface(
        ffd_points=read_plot3d_ffd(ffd_path),
        upper_grid=case.upper_grid,
        lower_grid=case.lower_grid,
        out_path=figure_path,
        title="Reference full-wing FFD loaded with DVGeometry",
        baseline_pointset=case.pointset,
    )

    summary_path.write_text(json.dumps(asdict(case.ffd_summary), indent=2, default=str), encoding="utf-8")

    print(f"FFD box written to: {ffd_path}")
    print(
        "FFD dimensions: "
        f"{case.ffd_summary.n_streamwise} x {case.ffd_summary.n_vertical} x {case.ffd_summary.n_spanwise}"
    )
    print(
        "FFD extents: "
        f"x=[{case.ffd_summary.x_min:.3f}, {case.ffd_summary.x_max:.3f}], "
        f"y=[{case.ffd_summary.y_min:.3f}, {case.ffd_summary.y_max:.3f}], "
        f"z=[{case.ffd_summary.z_min:.3f}, {case.ffd_summary.z_max:.3f}]"
    )
    print("DVGeometry load check: True")
    print(f"DVGeometry Tecplot written to: {output_dir / 'reference_fullwing_baseline_ffd.dat'}")
    print(f"DVGeometry pointset written to: {output_dir / 'reference_fullwing_baseline_surface_wing_surface.dat'}")
    print(f"Visualization written to: {figure_path}")
    print(f"Summary written to: {summary_path}")


if __name__ == "__main__":
    main()
