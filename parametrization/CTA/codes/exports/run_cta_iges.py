from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
CTA_DIR = SCRIPT_DIR.parent.parent
REPO_ROOT = CTA_DIR.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from parametrization.CTA.case import build_cta_design, to_cta_model_config
from parametrization.bwb.builder import export_iges


def main() -> None:
    output_dir = CTA_DIR / "outputs" / "wing"
    profiles_dir = output_dir / "station_airfoils"
    iges_path = output_dir / "cta.igs"
    output_dir.mkdir(parents=True, exist_ok=True)

    design = build_cta_design()
    config = to_cta_model_config(design, use_cta_anchor_twist=True)
    config.export.out_dir = profiles_dir
    config.export.iges_path = iges_path
    config.export.symmetric = True

    export_iges(config)
    print(f"CTA IGES written to: {iges_path}")
    print(f"CTA station profiles written to: {profiles_dir}")


if __name__ == "__main__":
    main()
