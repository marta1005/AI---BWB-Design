from __future__ import annotations

from pathlib import Path
import os
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR.parent / ".mplconfig"))

from common_demo import build_demo_profiles, build_demo_wing, default_build_options, default_output_dir

from parametrization.aircraft import export_lifting_surface_iges


def main() -> None:
    profiles = build_demo_profiles()
    wing = build_demo_wing()
    component = wing.to_component_spec()
    output_dir = default_output_dir()

    result = export_lifting_surface_iges(
        component,
        profiles,
        out_dir=output_dir,
        options=default_build_options(),
    )

    print(f"IGES written to: {result.iges_path}")
    for path in result.airfoil_paths:
        print(f"Profile written to: {path}")


if __name__ == "__main__":
    main()
