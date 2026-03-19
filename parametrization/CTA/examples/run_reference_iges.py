from pathlib import Path
import os
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR.parent.parent / ".mplconfig"))

from parametrization.CTA.reference import build_reference_design, to_cta_model_config
from parametrization.bwb.builder import export_iges


def main() -> None:
    output_dir = SCRIPT_DIR.parent / "example_outputs" / "reference"
    profiles_dir = output_dir / "profiles_reference"
    iges_path = output_dir / "cta_reference.igs"
    output_dir.mkdir(parents=True, exist_ok=True)

    design = build_reference_design()
    config = to_cta_model_config(design)
    config.export.out_dir = profiles_dir
    config.export.iges_path = iges_path

    export_iges(config)
    print(f"CTA reference IGES written to: {iges_path}")
    print(f"CTA reference profiles written to: {profiles_dir}")


if __name__ == "__main__":
    main()
