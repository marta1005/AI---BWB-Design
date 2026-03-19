from __future__ import annotations

from pathlib import Path
import os
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR.parent / ".mplconfig"))

from common_demo import build_demo_fuselage, default_fuselage_build_options, default_fuselage_output_dir

from parametrization.aircraft import export_fuselage_iges


def main() -> None:
    fuselage = build_demo_fuselage()
    output_dir = default_fuselage_output_dir()

    result = export_fuselage_iges(
        fuselage,
        out_dir=output_dir,
        options=default_fuselage_build_options(),
    )

    print(f"IGES written to: {result.iges_path}")
    for path in result.section_paths:
        print(f"Section written to: {path}")


if __name__ == "__main__":
    main()
