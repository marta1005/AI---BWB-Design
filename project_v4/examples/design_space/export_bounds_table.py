import csv
from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from project_v4.design_space import parameter_metadata


def _format_value(value: object) -> str:
    if isinstance(value, float):
        return "%.6g" % value
    return str(value)


def _write_markdown(path: Path, rows) -> None:
    headers = [
        "Group",
        "Parameter",
        "Display Name",
        "Symbol",
        "Units",
        "Normalization",
        "Lower Bound",
        "Reference",
        "Upper Bound",
        "Description",
    ]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        values = [
            row["group"],
            row["parameter"],
            row["display_name"],
            row["symbol"],
            row["units"],
            row["normalization"],
            _format_value(row["lower_bound"]),
            _format_value(row["reference"]),
            _format_value(row["upper_bound"]),
            row["description"],
        ]
        safe_values = [str(value).replace("|", "\\|") for value in values]
        lines.append("| " + " | ".join(safe_values) + " |")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    docs_dir = SCRIPT_DIR.parent.parent / "docs"
    docs_dir.mkdir(parents=True, exist_ok=True)
    csv_path = docs_dir / "design_variable_bounds.csv"
    md_path = docs_dir / "design_variable_bounds.md"

    rows = parameter_metadata()
    fieldnames = [
        "group",
        "parameter",
        "display_name",
        "symbol",
        "units",
        "normalization",
        "lower_bound",
        "reference",
        "upper_bound",
        "description",
    ]

    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    _write_markdown(md_path, rows)

    print("Bounds CSV written to: %s" % csv_path)
    print("Bounds Markdown written to: %s" % md_path)


if __name__ == "__main__":
    main()

