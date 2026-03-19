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

from common_demo import build_demo_fuselage, default_fuselage_build_options, default_fuselage_output_dir

from parametrization.aircraft import prepare_fuselage
from parametrization.aircraft.plotting import save_fuselage_3d, save_fuselage_overview


def main() -> None:
    fuselage = build_demo_fuselage()
    output_dir = default_fuselage_output_dir()
    prepared = prepare_fuselage(fuselage, options=default_fuselage_build_options())

    views_png, views_svg = save_fuselage_overview(
        prepared,
        output_dir,
        stem=fuselage.fuselage_id,
        title="Demo fuselage",
    )
    view3d_png = save_fuselage_3d(
        prepared,
        output_dir,
        stem=fuselage.fuselage_id,
        title="Demo fuselage | 3D view",
    )

    print(f"Views PNG written to: {views_png}")
    print(f"Views SVG written to: {views_svg}")
    print(f"3D PNG written to: {view3d_png}")


if __name__ == "__main__":
    main()
