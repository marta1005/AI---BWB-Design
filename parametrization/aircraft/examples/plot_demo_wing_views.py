from __future__ import annotations

from pathlib import Path
import os
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR.parent / ".mplconfig"))

import matplotlib

matplotlib.use("Agg")

from common_demo import build_demo_profiles, build_demo_wing, default_build_options, default_output_dir

from parametrization.aircraft import prepare_lifting_surface
from parametrization.aircraft.plotting import save_lifting_surface_3d, save_lifting_surface_overview


def main() -> None:
    profiles = build_demo_profiles()
    wing = build_demo_wing()
    component = wing.to_component_spec()
    output_dir = default_output_dir()
    profile_etas = tuple(float(station.eta) for station in wing.stations)

    prepared = prepare_lifting_surface(component, profiles, options=default_build_options())
    views_png = save_lifting_surface_overview(
        prepared,
        output_dir,
        stem=component.component_id,
        title="Demo main wing",
        profile_station_etas=profile_etas,
    )
    view3d_png = save_lifting_surface_3d(
        prepared,
        output_dir,
        stem=component.component_id,
        title="Demo main wing | 3D view",
    )

    print(f"Views PNG written to: {views_png}")
    print(f"3D PNG written to: {view3d_png}")


if __name__ == "__main__":
    main()
