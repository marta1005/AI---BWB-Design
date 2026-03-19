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

from common_demo import build_demo_aircraft, default_aircraft_output_dir

from parametrization.aircraft import prepare_aircraft_geometry
from parametrization.aircraft.plotting import save_aircraft_3d, save_aircraft_overview


def main() -> None:
    aircraft = build_demo_aircraft()
    output_dir = default_aircraft_output_dir()
    prepared = prepare_aircraft_geometry(aircraft)

    views_png, views_svg = save_aircraft_overview(
        prepared,
        output_dir,
        stem=aircraft.aircraft_id,
        title="Demo aircraft",
    )
    view3d_png = save_aircraft_3d(
        prepared,
        output_dir,
        stem=aircraft.aircraft_id,
        title="Demo aircraft | 3D view",
    )

    print(f"Views PNG written to: {views_png}")
    print(f"Views SVG written to: {views_svg}")
    print(f"3D PNG written to: {view3d_png}")


if __name__ == "__main__":
    main()
