from pathlib import Path
import csv
import json
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from parametrization.CTA.design_space import cta_design_space_summary, cta_parameter_metadata


def write_csv(path: Path, rows):
    header = [
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
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=header)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_markdown(path: Path, rows, summary):
    fixed = summary["fixed_parameters"]
    lines = [
        "# CTA AI Design-Space Bounds",
        "",
        "## Fixed Parameters",
        "",
        f"- `S` (`s1_deg`): `{fixed['s_deg']}` deg",
        f"- `B1` (`b1_fixed_m`): `{fixed['b1_fixed_m']}` m",
        "- `C1` conditioned by fixed `S` and straight `TE(C0->C1)`",
        "- no public `C2`: the inboard `TE(C1->C3)` blend is built with a hidden helper point",
        "- `C4` derived from active `C3`, `B2`, `S1`, and straight `TE(C3->C4)`",
        f"- `te_exact_segments`: `{fixed['te_exact_segments']}`",
        "",
        "## Active Variables",
        "",
        "| Parameter | Display Name | Symbol | Units | Lower | Reference | Upper |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in rows:
        lines.append(
            "| {parameter} | {display_name} | {symbol} | {units} | {lower_bound:.6g} | "
            "{reference:.6g} | {upper_bound:.6g} |".format(**row)
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    docs_dir = SCRIPT_DIR.parent / "docs"
    docs_dir.mkdir(parents=True, exist_ok=True)
    csv_path = docs_dir / "cta_ai_bounds.csv"
    md_path = docs_dir / "cta_ai_bounds.md"
    json_path = docs_dir / "cta_ai_summary.json"

    rows = cta_parameter_metadata()
    summary = cta_design_space_summary()

    write_csv(csv_path, rows)
    write_markdown(md_path, rows, summary)
    json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(f"CTA bounds CSV written to: {csv_path}")
    print(f"CTA bounds Markdown written to: {md_path}")
    print(f"CTA summary JSON written to: {json_path}")


if __name__ == "__main__":
    main()
