from __future__ import annotations

import csv
import os
from pathlib import Path
import sys
from typing import Dict, List

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
CTA_DIR = SCRIPT_DIR.parent.parent
REPO_ROOT = CTA_DIR.parent.parent
MPLCONFIG_DIR = REPO_ROOT / ".mplconfig"
MPLCONFIG_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG_DIR))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from parametrization.CTA.design_space import build_cta_design_space
from parametrization.CTA.codes.plotting.animate_cta_planform_parameter_sweep import build_planform_sweep_samples


OUTPUT_PNG = CTA_DIR / "outputs" / "wing" / "cta_planform_parameter_frame_traces.png"
OUTPUT_CSV = CTA_DIR / "outputs" / "wing" / "cta_planform_parameter_frame_traces.csv"


def _c4_absolute(sample: Dict[str, float]) -> float:
    return float(sample["c2_c1_ratio"]) * float(sample["c4_c3_ratio"])


def _parameter_specs(bounds: Dict[str, tuple[float, float]]):
    return [
        ("c1_root_chord", "C0", "m", float(bounds["c1_root_chord"][0]), float(bounds["c1_root_chord"][1])),
        ("c2_c1_ratio", "C3", "m", float(bounds["c2_c1_ratio"][0]), float(bounds["c2_c1_ratio"][1])),
        ("c4_absolute", "C4", "m", 6.8, 9.8),
        ("c4_c1_ratio", "C5", "m", float(bounds["c4_c1_ratio"][0]), float(bounds["c4_c1_ratio"][1])),
        ("b2_span_ratio", "B2/Bw", "-", float(bounds["b2_span_ratio"][0]), float(bounds["b2_span_ratio"][1])),
        ("span", "Bw", "m", float(bounds["span"][0]), float(bounds["span"][1])),
    ]


def _sample_value(sample: Dict[str, float], name: str) -> float | None:
    if name == "c4_absolute":
        return _c4_absolute(sample)
    return float(sample[name])


def _write_csv(samples: List[Dict[str, float]]) -> None:
    fieldnames = [
        "frame",
        "C0_m",
        "C3_m",
        "C4_m",
        "C5_m",
        "B2_over_Bw",
        "Bw_m",
    ]
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for idx, sample in enumerate(samples):
            writer.writerow(
                {
                    "frame": idx,
                    "C0_m": float(sample["c1_root_chord"]),
                    "C3_m": float(sample["c2_c1_ratio"]),
                    "C4_m": _c4_absolute(sample),
                    "C5_m": float(sample["c4_c1_ratio"]),
                    "B2_over_Bw": float(sample["b2_span_ratio"]),
                    "Bw_m": float(sample["span"]),
                }
            )


def main() -> None:
    OUTPUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    samples = build_planform_sweep_samples(frame_count=51, seed=11)
    _write_csv(samples)

    space = build_cta_design_space()
    bounds = space.bounds
    cta_view = space.cta_flat()
    specs = _parameter_specs(bounds)
    x = np.arange(len(samples), dtype=int)

    fig, axes = plt.subplots(2, 3, figsize=(16.5, 7.9), constrained_layout=True, sharex=True)
    axes_flat = axes.ravel()

    for ax, (name, label, units, lower, upper) in zip(axes_flat, specs):
        y = np.array([float(_sample_value(sample, name)) for sample in samples], dtype=float)
        ref_key = {
            "c4_absolute": None,
        }.get(name, name)
        cta_value = _c4_absolute(cta_view) if name == "c4_absolute" else float(cta_view[ref_key])

        ax.axhline(lower, color="#94a3b8", linewidth=1.1, linestyle="--", label="Lower/Upper bound")
        ax.axhline(upper, color="#94a3b8", linewidth=1.1, linestyle="--")
        ax.axhline(cta_value, color="#475569", linewidth=1.0, linestyle=":", label="CTA")
        ax.plot(x, y, color="#0f4c5c", linewidth=1.8, marker="o", markersize=3.2)
        ax.fill_between(x, lower, upper, color="#bfdbfe", alpha=0.18)
        ax.set_title(label)
        ax.set_ylabel(units)
        ax.grid(True, linewidth=0.35, alpha=0.25)

    for ax in axes_flat:
        ax.set_xlabel("Frame")

    fig.suptitle("CTA planform animation parameter traces", fontsize=17)
    fig.savefig(OUTPUT_PNG, dpi=200)
    plt.close(fig)

    print(f"CTA parameter frame-trace PNG written to: {OUTPUT_PNG}")
    print(f"CTA parameter frame-trace CSV written to: {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
