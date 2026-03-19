from __future__ import annotations

from pathlib import Path
import os
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR.parent / ".mplconfig"))

from cta_planform_match import (
    build_cta_planform_profiles,
    build_cta_planform_reference,
    build_cta_planform_wing,
    default_output_dir,
    export_build_options,
)

from parametrization.aircraft import export_lifting_surface_iges


def main() -> None:
    reference = build_cta_planform_reference()
    profiles = build_cta_planform_profiles()
    wing = build_cta_planform_wing(reference=reference)
    output_dir = default_output_dir()

    result = export_lifting_surface_iges(
        wing.to_component_spec(),
        profiles,
        out_dir=output_dir,
        options=export_build_options(reference=reference),
        file_name="cta_planform_match.igs",
    )

    print(f"IGES written to: {result.iges_path}")
    print(f"Profile directory written to: {output_dir / 'profiles'}")


if __name__ == "__main__":
    main()
