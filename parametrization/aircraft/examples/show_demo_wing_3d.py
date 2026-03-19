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

from common_demo import build_demo_profiles, build_demo_wing, default_build_options

from parametrization.aircraft import prepare_lifting_surface
from parametrization.aircraft.plotting import create_lifting_surface_3d_figure


def main() -> None:
    profiles = build_demo_profiles()
    wing = build_demo_wing()
    component = wing.to_component_spec()
    prepared = prepare_lifting_surface(component, profiles, options=default_build_options())
    create_lifting_surface_3d_figure(prepared, title="Demo main wing | interactive 3D")
    plt.show()


if __name__ == "__main__":
    main()
