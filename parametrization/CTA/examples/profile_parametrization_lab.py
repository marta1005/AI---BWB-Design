from __future__ import annotations

import argparse
import csv
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR.parent.parent / ".mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from parametrization.shared.cst import KulfanCSTAirfoil, cosine_spacing
from parametrization.bwb.design_variables import SectionedBWBDesignVariables


SECTION_KEYS: Dict[str, Tuple[str, str, str, str, str]] = {
    "c1": ("c1_upper_cst", "c1_lower_cst", "c1_tc_max", "c1_x_tmax", "c1_te_thickness"),
    "c2": ("c2_upper_cst", "c2_lower_cst", "c2_tc_max", "c2_x_tmax", "c2_te_thickness"),
    "c3": ("c3_upper_cst", "c3_lower_cst", "c3_tc_max", "c3_x_tmax", "c3_te_thickness"),
    "c4": ("c4_upper_cst", "c4_lower_cst", "c4_tc_max", "c4_x_tmax", "c4_te_thickness"),
}
THICKNESS_CHECK_POINTS = (0.20, 0.30, 0.40, 0.50, 0.60, 0.70)


@dataclass
class ProfileCase:
    case_id: str
    family: str
    section: str
    n1: float
    n2: float
    upper: np.ndarray
    lower: np.ndarray
    tc_target: float
    x_tmax_target: float
    te_thickness: float
    x: np.ndarray
    yu: np.ndarray
    yl: np.ndarray
    metrics: Dict[str, float]


def _section_data(
    design: SectionedBWBDesignVariables,
    section: str,
) -> Tuple[np.ndarray, np.ndarray, float, float, float]:
    upper_key, lower_key, tc_key, xt_key, te_key = SECTION_KEYS[section]
    upper = np.asarray(getattr(design, upper_key), dtype=float)
    lower = np.asarray(getattr(design, lower_key), dtype=float)
    tc_target = float(getattr(design, tc_key))
    x_tmax_target = float(getattr(design, xt_key))
    te_thickness = float(getattr(design, te_key))
    return upper, lower, tc_target, x_tmax_target, te_thickness


def _evaluate_profile(
    section: str,
    n1: float,
    n2: float,
    upper: np.ndarray,
    lower: np.ndarray,
    tc_target: float,
    x_tmax_target: float,
    te_thickness: float,
    enforce_targets: bool,
    sample_count: int = 501,
) -> ProfileCase:
    x = cosine_spacing(sample_count)
    airfoil = KulfanCSTAirfoil(degree=5, n1=float(n1), n2=float(n2), x_tc_window=(0.15, 0.65))
    coeffs = np.concatenate([upper, lower], dtype=float)
    kwargs = {"te_thickness": float(te_thickness)}
    if enforce_targets:
        kwargs["tc_target"] = float(tc_target)
        kwargs["x_tmax"] = float(x_tmax_target)
    yu, yl = airfoil.evaluate(x, coeffs, **kwargs)
    metrics = _compute_metrics(x, yu, yl)
    return ProfileCase(
        case_id="",
        family="",
        section=section,
        n1=float(n1),
        n2=float(n2),
        upper=upper.copy(),
        lower=lower.copy(),
        tc_target=float(tc_target),
        x_tmax_target=float(x_tmax_target),
        te_thickness=float(te_thickness),
        x=x,
        yu=yu,
        yl=yl,
        metrics=metrics,
    )


def _compute_metrics(x: np.ndarray, yu: np.ndarray, yl: np.ndarray) -> Dict[str, float]:
    thickness = yu - yl
    camber = 0.5 * (yu + yl)
    mask = (x >= 0.15) & (x <= 0.65)
    x_window = x[mask]
    t_window = thickness[mask]
    idx = int(np.argmax(t_window))
    metrics: Dict[str, float] = {
        "tc_max": float(t_window[idx]),
        "x_tmax": float(x_window[idx]),
        "te_thickness": float(thickness[-1]),
        "max_camber": float(np.max(np.abs(camber))),
        "area": float(np.trapezoid(thickness, x)),
    }
    for xc in THICKNESS_CHECK_POINTS:
        key = f"t_xc_{xc:.2f}"
        metrics[key] = float(np.interp(float(xc), x, thickness))
    return metrics


def _safe_coeff_variation(value: float, rel_delta: float, sign: int, min_floor: float = 1.0e-6) -> float:
    if value > 1.0e-8:
        candidate = value * (1.0 + sign * rel_delta)
    else:
        candidate = value + sign * rel_delta * 0.03
    return float(max(min_floor, candidate))


def _to_row(case: ProfileCase) -> Dict[str, float | str]:
    row: Dict[str, float | str] = {
        "case_id": case.case_id,
        "family": case.family,
        "section": case.section,
        "n1": case.n1,
        "n2": case.n2,
        "tc_target": case.tc_target,
        "x_tmax_target": case.x_tmax_target,
        "te_target": case.te_thickness,
    }
    row.update(case.metrics)
    return row


def _plot_profile(ax, case: ProfileCase, color: str, alpha: float, lw: float, ls: str = "-") -> None:
    ax.plot(case.x, case.yu, color=color, alpha=alpha, linewidth=lw, linestyle=ls)
    ax.plot(case.x, case.yl, color=color, alpha=alpha, linewidth=lw, linestyle=ls)


def _make_nclass_sweep(
    section: str,
    base_case: ProfileCase,
    n1_min: float,
    n1_max: float,
    n2_min: float,
    n2_max: float,
    n_grid: int,
    enforce_targets: bool,
) -> List[ProfileCase]:
    n1_values = np.linspace(float(n1_min), float(n1_max), int(n_grid))
    n2_values = np.linspace(float(n2_min), float(n2_max), int(n_grid))
    cases: List[ProfileCase] = []
    idx = 0
    for n1 in n1_values:
        for n2 in n2_values:
            case = _evaluate_profile(
                section=section,
                n1=float(n1),
                n2=float(n2),
                upper=base_case.upper,
                lower=base_case.lower,
                tc_target=base_case.tc_target,
                x_tmax_target=base_case.x_tmax_target,
                te_thickness=base_case.te_thickness,
                enforce_targets=enforce_targets,
            )
            case.case_id = f"nclass_{idx:03d}"
            case.family = "nclass"
            cases.append(case)
            idx += 1
    return cases


def _make_coeff_sweep(
    section: str,
    base_case: ProfileCase,
    coeff_rel_delta: float,
    enforce_targets: bool,
) -> List[ProfileCase]:
    cases: List[ProfileCase] = []
    idx = 0
    for surface_name in ("upper", "lower"):
        for coeff_idx in range(6):
            for sign in (-1, 1):
                upper = base_case.upper.copy()
                lower = base_case.lower.copy()
                if surface_name == "upper":
                    upper[coeff_idx] = _safe_coeff_variation(upper[coeff_idx], coeff_rel_delta, sign)
                else:
                    lower[coeff_idx] = _safe_coeff_variation(lower[coeff_idx], coeff_rel_delta, sign)
                case = _evaluate_profile(
                    section=section,
                    n1=base_case.n1,
                    n2=base_case.n2,
                    upper=upper,
                    lower=lower,
                    tc_target=base_case.tc_target,
                    x_tmax_target=base_case.x_tmax_target,
                    te_thickness=base_case.te_thickness,
                    enforce_targets=enforce_targets,
                )
                sign_label = "plus" if sign > 0 else "minus"
                case.case_id = f"coeff_{surface_name}_{coeff_idx}_{sign_label}"
                case.family = "coeff"
                cases.append(case)
                idx += 1
    return cases


def _plot_nclass_figure(
    output_path: Path,
    section: str,
    base_case: ProfileCase,
    nclass_cases: List[ProfileCase],
) -> None:
    fig, ax = plt.subplots(figsize=(11.0, 6.2), constrained_layout=True)
    cmap = plt.get_cmap("viridis")
    total = max(1, len(nclass_cases) - 1)
    for idx, case in enumerate(nclass_cases):
        color = cmap(idx / total)
        _plot_profile(ax, case, color=color, alpha=0.40, lw=0.9)
    _plot_profile(ax, base_case, color="#0f172a", alpha=1.0, lw=2.2)

    ax.set_title(f"{section.upper()} class-function sweep (N1/N2)")
    ax.set_xlabel("x/c")
    ax.set_ylabel("z/c")
    ax.grid(True, linewidth=0.35, alpha=0.30)
    ax.set_aspect("equal", adjustable="box")
    ax.text(
        0.02,
        0.98,
        (
            f"Baseline: N1={base_case.n1:.3f}, N2={base_case.n2:.3f}\n"
            f"tc={base_case.metrics['tc_max']:.4f}, x_tmax={base_case.metrics['x_tmax']:.4f}"
        ),
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9.0,
        bbox={"boxstyle": "round,pad=0.2", "facecolor": "white", "edgecolor": "#d1d5db", "alpha": 0.90},
    )
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _plot_coeff_sensitivity_figure(
    output_path: Path,
    section: str,
    base_case: ProfileCase,
    coeff_cases: List[ProfileCase],
    surface: str,
) -> None:
    fig, ax = plt.subplots(figsize=(11.0, 6.2), constrained_layout=True)
    colors = plt.get_cmap("tab10")
    for case in coeff_cases:
        if f"coeff_{surface}_" not in case.case_id:
            continue
        coeff_idx = int(case.case_id.split("_")[2])
        sign = case.case_id.split("_")[-1]
        ls = "-" if sign == "plus" else "--"
        _plot_profile(ax, case, color=colors(coeff_idx % 10), alpha=0.78, lw=1.1, ls=ls)
    _plot_profile(ax, base_case, color="#0f172a", alpha=1.0, lw=2.3)

    legend_items = [f"a{idx} (+ solid / - dashed)" for idx in range(6)]
    ax.text(
        0.98,
        0.98,
        "\n".join(legend_items),
        transform=ax.transAxes,
        va="top",
        ha="right",
        fontsize=8.8,
        bbox={"boxstyle": "round,pad=0.22", "facecolor": "white", "edgecolor": "#d1d5db", "alpha": 0.92},
    )
    ax.set_title(f"{section.upper()} CST sensitivity ({surface})")
    ax.set_xlabel("x/c")
    ax.set_ylabel("z/c")
    ax.grid(True, linewidth=0.35, alpha=0.30)
    ax.set_aspect("equal", adjustable="box")
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _build_sensitivity_scores(
    base_case: ProfileCase,
    coeff_cases: List[ProfileCase],
) -> Dict[str, np.ndarray]:
    weights = {
        "tc_max": 1.0,
        "x_tmax": 0.9,
        "max_camber": 0.8,
        "area": 0.6,
        "t_xc_0.40": 0.8,
        "t_xc_0.60": 0.8,
    }
    scores = {"upper": np.zeros(6, dtype=float), "lower": np.zeros(6, dtype=float)}
    counts = {"upper": np.zeros(6, dtype=float), "lower": np.zeros(6, dtype=float)}

    for case in coeff_cases:
        parts = case.case_id.split("_")
        if len(parts) < 4:
            continue
        surface = parts[1]
        idx = int(parts[2])
        score = 0.0
        for metric_name, weight in weights.items():
            score += float(weight) * abs(float(case.metrics[metric_name]) - float(base_case.metrics[metric_name]))
        scores[surface][idx] += score
        counts[surface][idx] += 1.0

    for surface in ("upper", "lower"):
        valid = counts[surface] > 0.0
        scores[surface][valid] /= counts[surface][valid]
    return scores


def _plot_score_bars(output_path: Path, section: str, scores: Dict[str, np.ndarray]) -> None:
    idx = np.arange(6)
    fig, ax = plt.subplots(figsize=(9.4, 5.4), constrained_layout=True)
    width = 0.37
    ax.bar(idx - 0.5 * width, scores["upper"], width=width, color="#1d4ed8", label="Upper coefficients")
    ax.bar(idx + 0.5 * width, scores["lower"], width=width, color="#9333ea", label="Lower coefficients")
    ax.set_xticks(idx)
    ax.set_xticklabels([f"a{i}" for i in idx])
    ax.set_ylabel("Sensitivity score")
    ax.set_title(f"{section.upper()} CST coefficient influence ranking")
    ax.grid(True, axis="y", linewidth=0.35, alpha=0.35)
    ax.legend(loc="upper right")
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _write_cases_csv(path: Path, rows: List[Dict[str, float | str]]) -> None:
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _write_report(path: Path, section: str, base_case: ProfileCase, scores: Dict[str, np.ndarray]) -> None:
    lines = []
    lines.append(f"Section: {section.upper()}")
    lines.append(f"Baseline N1={base_case.n1:.4f}, N2={base_case.n2:.4f}")
    lines.append(
        "Baseline metrics: "
        f"tc_max={base_case.metrics['tc_max']:.5f}, "
        f"x_tmax={base_case.metrics['x_tmax']:.5f}, "
        f"te={base_case.metrics['te_thickness']:.5f}, "
        f"max_camber={base_case.metrics['max_camber']:.5f}"
    )
    lines.append("")
    lines.append("Most influential UPPER coefficients:")
    upper_rank = np.argsort(scores["upper"])[::-1]
    for idx in upper_rank[:3]:
        lines.append(f"  a{idx}: score={scores['upper'][idx]:.6f}")
    lines.append("")
    lines.append("Most influential LOWER coefficients:")
    lower_rank = np.argsort(scores["lower"])[::-1]
    for idx in lower_rank[:3]:
        lines.append(f"  a{idx}: score={scores['lower'][idx]:.6f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Standalone profile parametrization laboratory for CTA CST sections."
    )
    parser.add_argument("--section", choices=tuple(SECTION_KEYS.keys()), default="c1")
    parser.add_argument("--n1-min", type=float, default=0.42)
    parser.add_argument("--n1-max", type=float, default=0.62)
    parser.add_argument("--n2-min", type=float, default=0.85)
    parser.add_argument("--n2-max", type=float, default=1.20)
    parser.add_argument("--n-grid", type=int, default=5)
    parser.add_argument("--coeff-rel-delta", type=float, default=0.15)
    parser.add_argument("--enforce-targets", action="store_true", default=False)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=SCRIPT_DIR.parent / "example_outputs" / "profile_lab",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    section = str(args.section).lower()
    output_dir = Path(args.output_dir).resolve() / section
    output_dir.mkdir(parents=True, exist_ok=True)

    design = SectionedBWBDesignVariables.reference_seed()
    upper, lower, tc_target, x_tmax_target, te_thickness = _section_data(design, section)
    baseline = _evaluate_profile(
        section=section,
        n1=float(design.cst_n1),
        n2=float(design.cst_n2),
        upper=upper,
        lower=lower,
        tc_target=tc_target,
        x_tmax_target=x_tmax_target,
        te_thickness=te_thickness,
        enforce_targets=bool(args.enforce_targets),
    )
    baseline.case_id = "baseline"
    baseline.family = "baseline"

    nclass_cases = _make_nclass_sweep(
        section=section,
        base_case=baseline,
        n1_min=float(args.n1_min),
        n1_max=float(args.n1_max),
        n2_min=float(args.n2_min),
        n2_max=float(args.n2_max),
        n_grid=int(args.n_grid),
        enforce_targets=bool(args.enforce_targets),
    )
    coeff_cases = _make_coeff_sweep(
        section=section,
        base_case=baseline,
        coeff_rel_delta=float(args.coeff_rel_delta),
        enforce_targets=bool(args.enforce_targets),
    )

    _plot_nclass_figure(output_dir / f"{section}_nclass_sweep.png", section, baseline, nclass_cases)
    _plot_coeff_sensitivity_figure(
        output_dir / f"{section}_coeff_sensitivity_upper.png", section, baseline, coeff_cases, surface="upper"
    )
    _plot_coeff_sensitivity_figure(
        output_dir / f"{section}_coeff_sensitivity_lower.png", section, baseline, coeff_cases, surface="lower"
    )

    scores = _build_sensitivity_scores(baseline, coeff_cases)
    _plot_score_bars(output_dir / f"{section}_coeff_scores.png", section, scores)

    rows = [_to_row(baseline)] + [_to_row(case) for case in nclass_cases] + [_to_row(case) for case in coeff_cases]
    _write_cases_csv(output_dir / f"{section}_profile_cases.csv", rows)
    _write_report(output_dir / f"{section}_summary.txt", section, baseline, scores)

    print(f"Profile lab outputs written to: {output_dir}")
    print(f"Section: {section.upper()}")
    print(f"Baseline: N1={baseline.n1:.3f}, N2={baseline.n2:.3f}")
    print(
        "Baseline metrics: "
        f"tc_max={baseline.metrics['tc_max']:.4f}, "
        f"x_tmax={baseline.metrics['x_tmax']:.4f}, "
        f"max_camber={baseline.metrics['max_camber']:.4f}"
    )
    for surface in ("upper", "lower"):
        rank = np.argsort(scores[surface])[::-1]
        top = ", ".join([f"a{idx}({scores[surface][idx]:.5f})" for idx in rank[:3]])
        print(f"Top {surface} coefficients to tune first: {top}")


if __name__ == "__main__":
    main()
