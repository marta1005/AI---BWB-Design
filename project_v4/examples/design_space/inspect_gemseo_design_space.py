import argparse
import json
from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from project_v4.gemseo_space import build_gemseo_design_space, build_gemseo_design_space_definition


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Inspect the GEMSEO design-space mapping for project_v4.")
    parser.add_argument("--preset", default="ai_geometry_core")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        adapter = build_gemseo_design_space(args.preset)
    except ModuleNotFoundError as exc:
        print(str(exc))
        print("Falling back to the GEMSEO design-space definition only (without backend instantiation).")
        adapter = build_gemseo_design_space_definition(args.preset)
    output_dir = SCRIPT_DIR.parent.parent / "example_outputs" / ("gemseo_design_space_%s" % args.preset)
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / ("gemseo_design_space_%s.json" % args.preset)
    json_path.write_text(json.dumps(adapter.summary_rows(), indent=2), encoding="utf-8")

    print("GEMSEO preset: %s" % args.preset)
    print("Variables: %s" % ", ".join(spec.name for spec in adapter.variable_specs))
    print("Summary JSON written to: %s" % json_path)


if __name__ == "__main__":
    main()
