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

from cta_planform_match import (
    build_cta_planform_profiles,
    build_cta_planform_reference,
    build_cta_planform_wing,
    default_build_options,
)

from parametrization.aircraft import prepare_lifting_surface
from parametrization.aircraft.plotting import create_lifting_surface_3d_figure


def main() -> None:
    reference = build_cta_planform_reference()
    profiles = build_cta_planform_profiles()
    wing = build_cta_planform_wing(reference=reference)
    prepared = prepare_lifting_surface(
        wing.to_component_spec(),
        profiles,
        options=default_build_options(reference=reference),
    )
    create_lifting_surface_3d_figure(prepared, title="CTA planform match | interactive 3D")
    plt.show()


if __name__ == "__main__":
    main()
