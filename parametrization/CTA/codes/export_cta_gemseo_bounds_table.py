import csv
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


def _format_scalar_or_vector(values) -> str:
    if len(values) == 1:
        return "%.6g" % float(values[0])
    return "[" + ", ".join("%.6g" % float(value) for value in values) + "]"


def _write_markdown(path: Path, rows) -> None:
    headers = [
        "GEMSEO Variable",
        "CTA Parameters",
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
            ", ".join(row["cta_parameters"]),
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
    try:
        adapter = build_cta_gemseo_design_space()
    except ModuleNotFoundError as exc:
        print(str(exc))
        print("Falling back to the CTA GEMSEO design-space definition only (without backend instantiation).")
        adapter = build_cta_gemseo_design_space_definition()

    rows = adapter.summary_rows()
    docs_dir = SCRIPT_DIR.parent / "docs"
    docs_dir.mkdir(parents=True, exist_ok=True)
    csv_path = docs_dir / "cta_gemseo_bounds.csv"
    md_path = docs_dir / "cta_gemseo_bounds.md"

    fieldnames = [
        "gemseo_variable",
        "size",
        "cta_parameters",
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

    print("CTA GEMSEO bounds CSV written to: %s" % csv_path)
    print("CTA GEMSEO bounds Markdown written to: %s" % md_path)


if __name__ == "__main__":
    main()
