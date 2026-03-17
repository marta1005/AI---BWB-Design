import csv
import argparse
from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from project_v4.gemseo_space import build_gemseo_design_space, build_gemseo_design_space_definition


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export GEMSEO bounds tables for a project_v4 design-space preset."
    )
    parser.add_argument("--preset", default="ai_geometry_core")
    return parser.parse_args()


def _format_scalar_or_vector(values) -> str:
    if len(values) == 1:
        return "%.6g" % float(values[0])
    return "[" + ", ".join("%.6g" % float(value) for value in values) + "]"


def _write_markdown(path: Path, rows) -> None:
    headers = [
        "GEMSEO Variable",
        "Project Parameters",
        "Units",
        "Normalization",
        "Lower Bounds",
        "Reference",
        "Upper Bounds",
        "Description",
    ]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        values = [
            row["gemseo_variable"],
            ", ".join(row["project_parameters"]),
            row["units"],
            row["normalization"],
            _format_scalar_or_vector(row["lower_bounds"]),
            _format_scalar_or_vector(row["reference"]),
            _format_scalar_or_vector(row["upper_bounds"]),
            row["description"],
        ]
        safe_values = [str(value).replace("|", "\\|") for value in values]
        lines.append("| " + " | ".join(safe_values) + " |")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    preset = args.preset
    try:
        adapter = build_gemseo_design_space(preset)
    except ModuleNotFoundError as exc:
        print(str(exc))
        print("Falling back to the GEMSEO design-space definition only (without backend instantiation).")
        adapter = build_gemseo_design_space_definition(preset)
    rows = adapter.summary_rows()

    docs_dir = SCRIPT_DIR.parent.parent / "docs"
    docs_dir.mkdir(parents=True, exist_ok=True)
    csv_path = docs_dir / f"{preset}_bounds.csv"
    md_path = docs_dir / f"{preset}_bounds.md"

    fieldnames = [
        "gemseo_variable",
        "size",
        "project_parameters",
        "units",
        "normalization",
        "description",
        "lower_bounds",
        "reference",
        "upper_bounds",
    ]

    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    _write_markdown(md_path, rows)

    print("GEMSEO bounds CSV written to: %s" % csv_path)
    print("GEMSEO bounds Markdown written to: %s" % md_path)


if __name__ == "__main__":
    main()
