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

from parametrization.aircraft import ContinuityOrder, prepare_lifting_surface
from parametrization.aircraft.plotting import save_lifting_surface_overview


def main() -> None:
    profiles = build_demo_profiles()
    output_dir = default_output_dir()

    for continuity in (ContinuityOrder.C1, ContinuityOrder.C2):
        wing = build_demo_wing(transition_continuity=continuity)
        component = wing.to_component_spec()
        prepared = prepare_lifting_surface(component, profiles, options=default_build_options())
        stem = f"{component.component_id}_transition_c{int(continuity)}"
        views_png = save_lifting_surface_overview(
            prepared,
            output_dir,
            stem=stem,
            title=f"Demo main wing | segmented C{int(continuity)} transition",
            profile_station_etas=tuple(float(station.eta) for station in wing.stations),
        )
        print(f"C{int(continuity)} PNG written to: {views_png}")


if __name__ == "__main__":
    main()
