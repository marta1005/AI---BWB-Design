import json
import os
from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
MPLCONFIG_DIR = REPO_ROOT / ".mplconfig"
MPLCONFIG_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG_DIR))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from parametrization.CTA.gemseo_space import (
    build_cta_gemseo_design_space,
    build_cta_gemseo_design_space_definition,
)


def main() -> None:
    try:
        adapter = build_cta_gemseo_design_space()
    except ModuleNotFoundError as exc:
        print(str(exc))
        print("Falling back to the CTA GEMSEO design-space definition only (without backend instantiation).")
        adapter = build_cta_gemseo_design_space_definition()

    output_dir = SCRIPT_DIR.parent / "example_outputs" / "cta_gemseo_design_space"
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "cta_gemseo_design_space.json"
    json_path.write_text(json.dumps(adapter.summary_rows(), indent=2), encoding="utf-8")

    print("CTA GEMSEO variables: %s" % ", ".join(spec.name for spec in adapter.variable_specs))
    print("Summary JSON written to: %s" % json_path)


if __name__ == "__main__":
    main()
