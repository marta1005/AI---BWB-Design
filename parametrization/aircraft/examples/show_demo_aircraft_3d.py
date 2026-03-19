from __future__ import annotations

from pathlib import Path
import os
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR.parent / ".mplconfig"))

import matplotlib.pyplot as plt

from common_demo import build_demo_aircraft

from parametrization.aircraft import prepare_aircraft_geometry
from parametrization.aircraft.plotting import create_aircraft_3d_figure


def main() -> None:
    aircraft = build_demo_aircraft()
    prepared = prepare_aircraft_geometry(aircraft)
    create_aircraft_3d_figure(prepared, title="Demo aircraft | interactive 3D")
    plt.show()


if __name__ == "__main__":
    main()
