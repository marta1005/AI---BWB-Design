from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

from .builder import PreparedGeometry, prepare_geometry
from .specs import SectionedBWBModelConfig


@dataclass(frozen=True)
class CadReferenceFrame:
    name: str = "model"
    offset_x_m: float = 0.0
    offset_y_m: float = 0.0
    offset_z_m: float = 0.0
    mirror_about_symmetry_plane: bool = True

    def xy_to_model(self, x_m, y_m) -> tuple[np.ndarray, np.ndarray]:
        x_model = np.asarray(x_m, dtype=float) - float(self.offset_x_m)
        y_model = np.asarray(y_m, dtype=float) - float(self.offset_y_m)
        if self.mirror_about_symmetry_plane:
            y_model = np.abs(y_model)
        return x_model, y_model

    def z_to_cad(self, z_model_m) -> np.ndarray:
        return np.asarray(z_model_m, dtype=float) + float(self.offset_z_m)

    def xyz_to_cad(self, xyz_model_m: np.ndarray) -> np.ndarray:
        xyz = np.asarray(xyz_model_m, dtype=float).copy()
        xyz[..., 0] += float(self.offset_x_m)
        xyz[..., 1] += float(self.offset_y_m)
        xyz[..., 2] += float(self.offset_z_m)
        return xyz


@dataclass(frozen=True)
class IndicatorSurfaceSpec:
    category: str
    sub_category: str
    sense: str
    vertices_xyz_m: np.ndarray
    minimum_clearance_m: np.ndarray

    def __post_init__(self) -> None:
        vertices = np.asarray(self.vertices_xyz_m, dtype=float)
        clearance = np.asarray(self.minimum_clearance_m, dtype=float)
        object.__setattr__(self, "vertices_xyz_m", vertices)
        object.__setattr__(self, "minimum_clearance_m", clearance)
        if vertices.ndim != 2 or vertices.shape[1] != 3:
            raise ValueError(
                f"vertices_xyz_m must have shape (N, 3), got {vertices.shape}"
            )
        if vertices.shape[0] < 3:
            raise ValueError("indicator surfaces must define at least 3 vertices")
        if clearance.shape != (vertices.shape[0],):
            raise ValueError(
                "minimum_clearance_m must have one value per vertex, "
                f"got shape {clearance.shape} for {vertices.shape[0]} vertices"
            )
        if self.sense not in {"upper", "lower"}:
            raise ValueError(f"sense must be 'upper' or 'lower', got {self.sense!r}")

    @property
    def label(self) -> str:
        return f"{self.category}:{self.sub_category}"

    @property
    def polygon_xy_m(self) -> np.ndarray:
        return np.asarray(self.vertices_xyz_m[:, :2], dtype=float)

    @property
    def z_vertices_m(self) -> np.ndarray:
        return np.asarray(self.vertices_xyz_m[:, 2], dtype=float)


@dataclass(frozen=True)
class InternalVolumeConstraintSet:
    name: str
    source_path: Optional[Path]
    reference_frame: CadReferenceFrame
    surfaces: tuple[IndicatorSurfaceSpec, ...]


@dataclass(frozen=True)
class IndicatorSurfaceResult:
    label: str
    category: str
    sub_category: str
    sense: str
    satisfied: bool
    sample_count: int
    invalid_sample_count: int
    footprint_area_m2: float
    worst_margin_m: float
    mean_margin_m: float
    critical_x_m: float
    critical_y_m: float
    critical_target_z_m: float
    critical_geometry_z_m: float
    critical_clearance_m: float


@dataclass(frozen=True)
class InternalVolumeConstraintResult:
    name: str
    source_path: Optional[Path]
    reference_frame: CadReferenceFrame
    satisfied: bool
    minimum_margin_m: float
    total_enclosed_volume_m3: float
    surface_results: tuple[IndicatorSurfaceResult, ...]


def _integrate(values: np.ndarray, coordinates: np.ndarray) -> float:
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(values, coordinates))
    return float(np.trapz(values, coordinates))


def _polygon_area_xy(vertices_xy_m: np.ndarray) -> float:
    xy = np.asarray(vertices_xy_m, dtype=float)
    x = xy[:, 0]
    y = xy[:, 1]
    return float(
        0.5
        * abs(
            np.dot(x, np.roll(y, -1))
            - np.dot(y, np.roll(x, -1))
        )
    )


def _triangle_sample_points(
    vertices_xyz_m: np.ndarray,
    clearance_m: np.ndarray,
    resolution: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    tri_xyz = np.asarray(vertices_xyz_m, dtype=float)
    tri_clearance = np.asarray(clearance_m, dtype=float)
    if tri_xyz.shape != (3, 3):
        raise ValueError(f"triangle vertices must have shape (3, 3), got {tri_xyz.shape}")
    if tri_clearance.shape != (3,):
        raise ValueError(
            f"triangle clearance must have shape (3,), got {tri_clearance.shape}"
        )

    resolution = max(int(resolution), 1)
    xyz_samples = []
    clearance_samples = []
    weights = []

    for i in range(resolution + 1):
        for j in range(resolution + 1 - i):
            w0 = float(i) / float(resolution)
            w1 = float(j) / float(resolution)
            w2 = 1.0 - w0 - w1
            bary = np.asarray([w0, w1, w2], dtype=float)
            xyz_samples.append(bary @ tri_xyz)
            clearance_samples.append(float(bary @ tri_clearance))
            weights.append(1.0)

    return (
        np.asarray(xyz_samples, dtype=float),
        np.asarray(clearance_samples, dtype=float),
        np.asarray(weights, dtype=float),
    )


def sample_indicator_surface(
    surface: IndicatorSurfaceSpec,
    triangle_resolution: int = 18,
) -> tuple[np.ndarray, np.ndarray]:
    vertices = np.asarray(surface.vertices_xyz_m, dtype=float)
    clearance = np.asarray(surface.minimum_clearance_m, dtype=float)
    sample_xyz = []
    sample_weights = []

    for idx in range(1, vertices.shape[0] - 1):
        tri_vertices = np.vstack([vertices[0], vertices[idx], vertices[idx + 1]])
        tri_clearance = np.asarray(
            [clearance[0], clearance[idx], clearance[idx + 1]],
            dtype=float,
        )
        tri_area = _polygon_area_xy(tri_vertices[:, :2])
        tri_xyz_samples, tri_clearance_samples, tri_weights = _triangle_sample_points(
            tri_vertices,
            tri_clearance,
            triangle_resolution,
        )
        tri_xyz_samples = tri_xyz_samples.copy()
        tri_xyz_samples[:, 2] = tri_xyz_samples[:, 2]
        sample_xyz.append(
            np.column_stack(
                [
                    tri_xyz_samples[:, 0],
                    tri_xyz_samples[:, 1],
                    tri_xyz_samples[:, 2],
                    tri_clearance_samples,
                ]
            )
        )
        sample_weights.append(
            np.full(tri_weights.shape, tri_area / max(float(tri_weights.size), 1.0), dtype=float)
        )

    if not sample_xyz:
        raise ValueError(f"Could not triangulate indicator surface {surface.label}")

    return np.vstack(sample_xyz), np.concatenate(sample_weights)


def _triangle_barycentric_coordinates_xy(
    x_m: float,
    y_m: float,
    triangle_xy_m: np.ndarray,
) -> tuple[float, float, float] | None:
    tri = np.asarray(triangle_xy_m, dtype=float)
    x0, y0 = tri[0]
    x1, y1 = tri[1]
    x2, y2 = tri[2]
    det = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2)
    if abs(float(det)) <= 1.0e-14:
        return None
    l0 = ((y1 - y2) * (float(x_m) - x2) + (x2 - x1) * (float(y_m) - y2)) / det
    l1 = ((y2 - y0) * (float(x_m) - x2) + (x0 - x2) * (float(y_m) - y2)) / det
    l2 = 1.0 - l0 - l1
    tol = 1.0e-10
    if l0 < -tol or l1 < -tol or l2 < -tol:
        return None
    return float(l0), float(l1), float(l2)


def _surface_required_z_at_points(
    surface: IndicatorSurfaceSpec,
    x_m: np.ndarray,
    y_m: np.ndarray,
) -> np.ndarray:
    x_arr = np.asarray(x_m, dtype=float)
    y_arr = np.asarray(y_m, dtype=float)
    required = np.full(x_arr.shape, np.nan, dtype=float)
    vertices = np.asarray(surface.vertices_xyz_m, dtype=float)
    clearance = np.asarray(surface.minimum_clearance_m, dtype=float)

    for tri_idx in range(1, vertices.shape[0] - 1):
        tri_vertices = np.vstack([vertices[0], vertices[tri_idx], vertices[tri_idx + 1]])
        tri_xy = tri_vertices[:, :2]
        tri_z = tri_vertices[:, 2]
        tri_clearance = np.asarray(
            [clearance[0], clearance[tri_idx], clearance[tri_idx + 1]],
            dtype=float,
        )
        for idx, (xx, yy) in enumerate(zip(np.ravel(x_arr), np.ravel(y_arr))):
            flat_idx = np.unravel_index(idx, x_arr.shape)
            if np.isfinite(required[flat_idx]):
                continue
            bary = _triangle_barycentric_coordinates_xy(float(xx), float(yy), tri_xy)
            if bary is None:
                continue
            weights = np.asarray(bary, dtype=float)
            z_target = float(weights @ tri_z)
            local_clearance = float(weights @ tri_clearance)
            if surface.sense == "upper":
                required[flat_idx] = z_target + local_clearance
            else:
                required[flat_idx] = z_target - local_clearance
    return required


def required_constraint_bounds_at_points(
    surfaces: tuple[IndicatorSurfaceSpec, ...] | list[IndicatorSurfaceSpec],
    x_m: np.ndarray,
    y_m: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    x_arr = np.asarray(x_m, dtype=float)
    y_arr = np.asarray(y_m, dtype=float)
    upper_required = np.full(x_arr.shape, -np.inf, dtype=float)
    lower_required = np.full(x_arr.shape, np.inf, dtype=float)
    has_upper = np.zeros(x_arr.shape, dtype=bool)
    has_lower = np.zeros(x_arr.shape, dtype=bool)

    for surface in tuple(surfaces):
        required = _surface_required_z_at_points(surface, x_arr, y_arr)
        valid = np.isfinite(required)
        if not np.any(valid):
            continue
        if surface.sense == "upper":
            upper_required[valid] = np.maximum(upper_required[valid], required[valid])
            has_upper |= valid
        else:
            lower_required[valid] = np.minimum(lower_required[valid], required[valid])
            has_lower |= valid

    upper_required[~has_upper] = np.nan
    lower_required[~has_lower] = np.nan
    return upper_required, lower_required


class GeometryEnvelopeEvaluator:
    def __init__(
        self,
        prepared: PreparedGeometry,
        reference_frame: CadReferenceFrame = CadReferenceFrame(),
    ) -> None:
        self.prepared = prepared
        self.reference_frame = reference_frame
        self._curve_cache: dict[float, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}

    @staticmethod
    def _sorted_curve(x_coords: np.ndarray, z_coords: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        order = np.argsort(np.asarray(x_coords, dtype=float))
        x_sorted = np.asarray(x_coords, dtype=float)[order]
        z_sorted = np.asarray(z_coords, dtype=float)[order]
        unique_x, inverse = np.unique(x_sorted, return_inverse=True)
        if unique_x.size == x_sorted.size:
            return x_sorted, z_sorted
        z_unique = np.zeros_like(unique_x, dtype=float)
        counts = np.zeros_like(unique_x, dtype=float)
        np.add.at(z_unique, inverse, z_sorted)
        np.add.at(counts, inverse, 1.0)
        z_unique /= np.maximum(counts, 1.0)
        return unique_x, z_unique

    def _surface_curves_at_y_model(
        self,
        y_model_m: float,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        cache_key = round(float(y_model_m), 9)
        if cache_key in self._curve_cache:
            return self._curve_cache[cache_key]

        y_model = float(y_model_m)
        if y_model < -1.0e-12 or y_model > float(self.prepared.loft.span_stations[-1]) + 1.0e-12:
            raise ValueError(
                f"y_model={y_model:.6f} lies outside the prepared semispan "
                f"[0, {float(self.prepared.loft.span_stations[-1]):.6f}]"
            )

        x_air = np.asarray(self.prepared.section_model.x_air, dtype=float)
        chord = float(np.interp(y_model, self.prepared.loft.span_stations, self.prepared.loft.chord))
        le_x = float(np.interp(y_model, self.prepared.loft.span_stations, self.prepared.loft.leading_edge_x))
        vertical_z = float(np.interp(y_model, self.prepared.loft.span_stations, self.prepared.loft.vertical_y))
        twist_deg = float(self.prepared.spanwise_laws.twist_deg(y_model))

        upper, lower, _ = self.prepared.section_model.coordinates_at_y(y_model)
        x_local = x_air * chord
        z_upper_local = np.asarray(upper, dtype=float) * chord
        z_lower_local = np.asarray(lower, dtype=float) * chord
        twist_rad = np.deg2rad(twist_deg)
        cos_twist = float(np.cos(twist_rad))
        sin_twist = float(np.sin(twist_rad))

        x_upper = le_x + x_local * cos_twist - z_upper_local * sin_twist
        x_lower = le_x + x_local * cos_twist - z_lower_local * sin_twist
        z_upper = vertical_z + x_local * sin_twist + z_upper_local * cos_twist
        z_lower = vertical_z + x_local * sin_twist + z_lower_local * cos_twist

        curves = (
            *self._sorted_curve(x_upper, z_upper),
            *self._sorted_curve(x_lower, z_lower),
        )
        self._curve_cache[cache_key] = curves
        return curves

    def _surface_z_at_xy_model(
        self,
        x_model_m: float,
        y_model_m: float,
        sense: str,
    ) -> float:
        x_upper, z_upper, x_lower, z_lower = self._surface_curves_at_y_model(y_model_m)
        if sense == "upper":
            return float(np.interp(float(x_model_m), x_upper, z_upper, left=np.nan, right=np.nan))
        if sense == "lower":
            return float(np.interp(float(x_model_m), x_lower, z_lower, left=np.nan, right=np.nan))
        raise ValueError(f"Unknown sense {sense!r}")

    def evaluate_points(
        self,
        x_cad_m: np.ndarray,
        y_cad_m: np.ndarray,
        sense: str,
    ) -> np.ndarray:
        x_model, y_model = self.reference_frame.xy_to_model(x_cad_m, y_cad_m)
        result = np.empty_like(x_model, dtype=float)
        flat_x = np.ravel(x_model)
        flat_y = np.ravel(y_model)
        flat_out = np.empty_like(flat_x, dtype=float)
        for idx, (xx, yy) in enumerate(zip(flat_x, flat_y)):
            try:
                z_model = self._surface_z_at_xy_model(float(xx), float(yy), sense)
            except ValueError:
                z_model = np.nan
            flat_out[idx] = float(self.reference_frame.z_to_cad(z_model))
        result[:] = flat_out.reshape(result.shape)
        return result


def load_internal_volume_constraint_set(
    csv_path: str | Path,
    *,
    reference_frame: CadReferenceFrame = CadReferenceFrame(),
    name: Optional[str] = None,
) -> InternalVolumeConstraintSet:
    path = Path(csv_path)
    rows = list(csv.reader(path.open(encoding="utf-8")))
    surfaces = []
    idx = 1

    while idx < len(rows):
        row = rows[idx]
        if len(row) >= 3 and row[0] and row[1] and row[2] in {"upper", "lower"}:
            category = str(row[0]).strip()
            sub_category = str(row[1]).strip()
            sense = str(row[2]).strip()
            idx += 1
            if idx >= len(rows):
                break
            header = rows[idx]
            if len(header) < 4 or header[0] != "X" or header[1] != "Y" or header[2] != "Z":
                raise ValueError(
                    f"Unexpected indicator-surface header at row {idx + 1}: {header}"
                )
            idx += 1
            vertices = []
            clearance = []
            while idx < len(rows):
                point_row = rows[idx]
                if not any(str(value).strip() for value in point_row):
                    break
                x_value = float(point_row[0])
                y_value = float(point_row[1])
                z_value = float(point_row[2])
                clearance_value = 0.0 if len(point_row) < 4 or point_row[3] == "" else float(point_row[3])
                vertices.append((x_value, y_value, z_value))
                clearance.append(clearance_value)
                idx += 1

            surfaces.append(
                IndicatorSurfaceSpec(
                    category=category,
                    sub_category=sub_category,
                    sense=sense,
                    vertices_xyz_m=np.asarray(vertices, dtype=float),
                    minimum_clearance_m=np.asarray(clearance, dtype=float),
                )
            )
        idx += 1

    return InternalVolumeConstraintSet(
        name=path.stem if name is None else str(name),
        source_path=path,
        reference_frame=reference_frame,
        surfaces=tuple(surfaces),
    )


def evaluate_indicator_surface(
    surface: IndicatorSurfaceSpec,
    evaluator: GeometryEnvelopeEvaluator,
    *,
    triangle_resolution: int = 18,
) -> IndicatorSurfaceResult:
    samples_xyzc, sample_weights = sample_indicator_surface(
        surface,
        triangle_resolution=triangle_resolution,
    )
    x_samples = samples_xyzc[:, 0]
    y_samples = samples_xyzc[:, 1]
    z_target = samples_xyzc[:, 2]
    clearance = samples_xyzc[:, 3]
    z_geometry = evaluator.evaluate_points(x_samples, y_samples, surface.sense)

    if surface.sense == "upper":
        required_z = z_target + clearance
        margin = z_geometry - required_z
    else:
        required_z = z_target - clearance
        margin = required_z - z_geometry

    invalid_mask = ~np.isfinite(z_geometry)
    if np.any(invalid_mask):
        margin = margin.copy()
        margin[invalid_mask] = -np.inf

    critical_idx = int(np.argmin(margin))
    valid_mask = np.isfinite(margin)
    if np.any(valid_mask):
        mean_margin = float(np.sum(margin[valid_mask] * sample_weights[valid_mask]) / np.sum(sample_weights[valid_mask]))
        worst_margin = float(np.min(margin[valid_mask]))
        satisfied = bool(np.all(margin[valid_mask] >= 0.0) and not np.any(invalid_mask))
    else:
        mean_margin = float("nan")
        worst_margin = float("-inf")
        satisfied = False

    return IndicatorSurfaceResult(
        label=surface.label,
        category=surface.category,
        sub_category=surface.sub_category,
        sense=surface.sense,
        satisfied=satisfied,
        sample_count=int(margin.size),
        invalid_sample_count=int(np.count_nonzero(invalid_mask)),
        footprint_area_m2=_polygon_area_xy(surface.polygon_xy_m),
        worst_margin_m=worst_margin,
        mean_margin_m=mean_margin,
        critical_x_m=float(x_samples[critical_idx]),
        critical_y_m=float(y_samples[critical_idx]),
        critical_target_z_m=float(required_z[critical_idx]),
        critical_geometry_z_m=float(z_geometry[critical_idx]),
        critical_clearance_m=float(clearance[critical_idx]),
    )


def evaluate_internal_volume_constraints(
    prepared: PreparedGeometry,
    constraint_set: InternalVolumeConstraintSet,
    *,
    triangle_resolution: int = 18,
) -> InternalVolumeConstraintResult:
    evaluator = GeometryEnvelopeEvaluator(
        prepared,
        reference_frame=constraint_set.reference_frame,
    )
    surface_results = tuple(
        evaluate_indicator_surface(
            surface,
            evaluator,
            triangle_resolution=triangle_resolution,
        )
        for surface in constraint_set.surfaces
    )
    minimum_margin = min(result.worst_margin_m for result in surface_results)
    satisfied = all(result.satisfied for result in surface_results)
    return InternalVolumeConstraintResult(
        name=constraint_set.name,
        source_path=constraint_set.source_path,
        reference_frame=constraint_set.reference_frame,
        satisfied=satisfied,
        minimum_margin_m=float(minimum_margin),
        total_enclosed_volume_m3=float(prepared.volume.enclosed_volume_m3),
        surface_results=surface_results,
    )


def evaluate_internal_volume_constraints_from_config(
    config: SectionedBWBModelConfig,
    constraint_set: InternalVolumeConstraintSet,
    *,
    triangle_resolution: int = 18,
) -> InternalVolumeConstraintResult:
    prepared = prepare_geometry(config)
    return evaluate_internal_volume_constraints(
        prepared,
        constraint_set,
        triangle_resolution=triangle_resolution,
    )
