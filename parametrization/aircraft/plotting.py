from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil
import tempfile

import matplotlib.pyplot as plt
import numpy as np

from .aircraft import PreparedAircraftGeometry
from .fuselage import PreparedFuselage
from .lifting_surface import PreparedLiftingSurface


@dataclass(frozen=True)
class LiftingSurfaceMesh:
    x_upper: np.ndarray
    y_upper: np.ndarray
    z_upper: np.ndarray
    x_lower: np.ndarray
    y_lower: np.ndarray
    z_lower: np.ndarray

    @property
    def x_all(self) -> np.ndarray:
        return np.concatenate([self.x_upper.ravel(), self.x_lower.ravel()])

    @property
    def y_all(self) -> np.ndarray:
        return np.concatenate([self.y_upper.ravel(), self.y_lower.ravel()])

    @property
    def z_all(self) -> np.ndarray:
        return np.concatenate([self.z_upper.ravel(), self.z_lower.ravel()])


@dataclass(frozen=True)
class FuselageMesh:
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray

    @property
    def x_all(self) -> np.ndarray:
        return self.x.ravel()

    @property
    def y_all(self) -> np.ndarray:
        return self.y.ravel()

    @property
    def z_all(self) -> np.ndarray:
        return self.z.ravel()


def _save_figure(fig, path: Path, **kwargs) -> None:
    try:
        fig.savefig(path, **kwargs)
    except OSError as exc:
        if "Resource deadlock avoided" not in str(exc):
            raise
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir) / path.name
            fig.savefig(tmp_path, **kwargs)
            shutil.copy2(tmp_path, path)


def _rotation_matrix_z(angle_rad: float) -> np.ndarray:
    c = float(np.cos(angle_rad))
    s = float(np.sin(angle_rad))
    return np.asarray([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]], dtype=float)


def _rotation_matrix_x(angle_rad: float) -> np.ndarray:
    c = float(np.cos(angle_rad))
    s = float(np.sin(angle_rad))
    return np.asarray([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]], dtype=float)


def _rotation_matrix_y(angle_rad: float) -> np.ndarray:
    c = float(np.cos(angle_rad))
    s = float(np.sin(angle_rad))
    return np.asarray([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]], dtype=float)


def _transform_section_points(
    x_local: np.ndarray,
    y_local: np.ndarray,
    station,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    local = np.column_stack(
        [
            x_local,
            y_local,
            np.zeros_like(x_local, dtype=float),
        ]
    )
    rotation = (
        _rotation_matrix_y(np.deg2rad(float(station.pitch_deg)))
        @ _rotation_matrix_x(np.deg2rad(float(station.roll_deg)))
        @ _rotation_matrix_z(np.deg2rad(float(station.twist_deg)))
    )
    world = local @ rotation.T
    world[:, 0] += float(station.x)
    world[:, 1] += float(station.y)
    world[:, 2] += float(station.z)
    return world[:, 0], world[:, 1], world[:, 2]


def build_lifting_surface_mesh(prepared: PreparedLiftingSurface) -> LiftingSurfaceMesh:
    stations = prepared.stations
    if not stations:
        raise ValueError("prepared lifting surface contains no stations")

    n_span = len(stations)
    n_airfoil = stations[0].x_air.size
    x_upper = np.zeros((n_span, n_airfoil), dtype=float)
    y_upper = np.zeros((n_span, n_airfoil), dtype=float)
    z_upper = np.zeros((n_span, n_airfoil), dtype=float)
    x_lower = np.zeros((n_span, n_airfoil), dtype=float)
    y_lower = np.zeros((n_span, n_airfoil), dtype=float)
    z_lower = np.zeros((n_span, n_airfoil), dtype=float)

    for index, station in enumerate(stations):
        x_local = float(station.chord) * np.asarray(station.x_air, dtype=float)
        yu_local = float(station.chord) * np.asarray(station.yu, dtype=float)
        yl_local = float(station.chord) * np.asarray(station.yl, dtype=float)

        xu, yu, zu = _transform_section_points(x_local, yu_local, station)
        xl, yl, zl = _transform_section_points(x_local, yl_local, station)

        x_upper[index, :] = xu
        y_upper[index, :] = yu
        z_upper[index, :] = zu
        x_lower[index, :] = xl
        y_lower[index, :] = yl
        z_lower[index, :] = zl

    return LiftingSurfaceMesh(
        x_upper=x_upper,
        y_upper=y_upper,
        z_upper=z_upper,
        x_lower=x_lower,
        y_lower=y_lower,
        z_lower=z_lower,
    )


def build_fuselage_mesh(prepared: PreparedFuselage) -> FuselageMesh:
    return FuselageMesh(
        x=prepared.x_grid,
        y=prepared.y_grid,
        z=prepared.z_grid,
    )


def _set_equal_2d(ax, x_values: np.ndarray, y_values: np.ndarray, pad_ratio: float = 0.06) -> None:
    x0 = float(np.min(x_values))
    x1 = float(np.max(x_values))
    y0 = float(np.min(y_values))
    y1 = float(np.max(y_values))
    dx = max(x1 - x0, 1e-9)
    dy = max(y1 - y0, 1e-9)
    x_pad = pad_ratio * dx
    y_pad = pad_ratio * dy
    ax.set_xlim(x0 - x_pad, x1 + x_pad)
    ax.set_ylim(y0 - y_pad, y1 + y_pad)
    ax.set_aspect("equal", adjustable="box")


def plot_planform(ax, mesh: LiftingSurfaceMesh, title: str = "Top view") -> None:
    ax.plot(mesh.x_upper[:, 0], mesh.z_upper[:, 0], color="#0f4c5c", linewidth=2.2, label="Leading edge")
    ax.plot(mesh.x_upper[:, -1], mesh.z_upper[:, -1], color="#c44536", linewidth=2.2, label="Trailing edge")
    ax.fill_between(
        mesh.x_upper[:, 0],
        mesh.z_upper[:, 0],
        mesh.z_upper[:, -1],
        color="#dbeafe",
        alpha=0.0,
    )
    ax.plot(mesh.x_upper[0, :], mesh.z_upper[0, :], color="#64748b", linewidth=1.2, alpha=0.9)
    ax.plot(mesh.x_upper[-1, :], mesh.z_upper[-1, :], color="#64748b", linewidth=1.2, alpha=0.9)
    for idx in range(mesh.x_upper.shape[0]):
        ax.plot(mesh.x_upper[idx, :], mesh.z_upper[idx, :], color="#93c5fd", linewidth=0.55, alpha=0.35)
        ax.plot(mesh.x_lower[idx, :], mesh.z_lower[idx, :], color="#93c5fd", linewidth=0.55, alpha=0.35)
    ax.set_title(title)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("z spanwise [m]")
    ax.grid(True, linewidth=0.35, alpha=0.25)
    _set_equal_2d(ax, mesh.x_all, mesh.z_all)


def plot_side_view(ax, mesh: LiftingSurfaceMesh, title: str = "Side view") -> None:
    stride = max(1, mesh.x_upper.shape[0] // 9)
    for idx in range(0, mesh.x_upper.shape[0], stride):
        ax.plot(mesh.x_upper[idx, :], mesh.y_upper[idx, :], color="#0f766e", linewidth=1.0, alpha=0.55)
        ax.plot(mesh.x_lower[idx, :], mesh.y_lower[idx, :], color="#0f766e", linewidth=1.0, alpha=0.55)
    ax.plot(mesh.x_upper[:, 0], mesh.y_upper[:, 0], color="#0f4c5c", linewidth=2.0)
    ax.plot(mesh.x_upper[:, -1], mesh.y_upper[:, -1], color="#c44536", linewidth=2.0)
    ax.plot(mesh.x_lower[:, 0], mesh.y_lower[:, 0], color="#0f4c5c", linewidth=1.4, alpha=0.65)
    ax.plot(mesh.x_lower[:, -1], mesh.y_lower[:, -1], color="#c44536", linewidth=1.4, alpha=0.65)
    ax.set_title(title)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y vertical [m]")
    ax.grid(True, linewidth=0.35, alpha=0.25)
    _set_equal_2d(ax, mesh.x_all, mesh.y_all)


def plot_front_view(ax, mesh: LiftingSurfaceMesh, title: str = "Front view") -> None:
    stride = max(1, mesh.x_upper.shape[0] // 9)
    for idx in range(0, mesh.x_upper.shape[0], stride):
        ax.plot(mesh.z_upper[idx, :], mesh.y_upper[idx, :], color="#0f766e", linewidth=1.0, alpha=0.55)
        ax.plot(mesh.z_lower[idx, :], mesh.y_lower[idx, :], color="#0f766e", linewidth=1.0, alpha=0.55)
    ax.plot(mesh.z_upper[:, 0], mesh.y_upper[:, 0], color="#0f4c5c", linewidth=2.0)
    ax.plot(mesh.z_upper[:, -1], mesh.y_upper[:, -1], color="#c44536", linewidth=2.0)
    ax.plot(mesh.z_lower[:, 0], mesh.y_lower[:, 0], color="#0f4c5c", linewidth=1.4, alpha=0.65)
    ax.plot(mesh.z_lower[:, -1], mesh.y_lower[:, -1], color="#c44536", linewidth=1.4, alpha=0.65)
    ax.set_title(title)
    ax.set_xlabel("z spanwise [m]")
    ax.set_ylabel("y vertical [m]")
    ax.grid(True, linewidth=0.35, alpha=0.25)
    _set_equal_2d(ax, mesh.z_all, mesh.y_all)


def plot_section_profiles(
    ax,
    prepared: PreparedLiftingSurface,
    title: str = "Section profiles (true scale)",
    selected_etas: tuple[float, ...] | None = None,
    max_profiles: int = 4,
) -> None:
    stations = prepared.stations
    if not stations:
        raise ValueError("prepared lifting surface contains no section profiles")

    def select_indices() -> list[int]:
        n = len(stations)
        if selected_etas:
            chosen = []
            for eta in selected_etas:
                target = float(eta)
                nearest = min(range(n), key=lambda idx: abs(float(stations[idx].eta) - target))
                if nearest not in chosen:
                    chosen.append(nearest)
            return chosen

        count = min(max(1, int(max_profiles)), n)
        if count >= n:
            return list(range(n))
        return sorted(set(np.linspace(0, n - 1, count, dtype=int).tolist()))

    indices = select_indices()
    raw_profiles = [
        (
            float(station.chord) * np.asarray(station.x_air, dtype=float),
            float(station.chord) * np.asarray(station.yu, dtype=float),
            float(station.chord) * np.asarray(station.yl, dtype=float),
            station,
        )
        for idx, station in enumerate(stations)
        if idx in indices
    ]

    height_spans = [float(np.max(yu) - np.min(yl)) for _, yu, yl, _ in raw_profiles]
    base_gap = max(max(height_spans, default=0.0), 1.0) * 0.55
    offsets = []
    current = 0.0
    for _, yu, yl, _ in raw_profiles:
        offsets.append(current)
        current += float(np.max(yu) - np.min(yl)) + base_gap

    for offset, (x_vals, yu, yl, station) in zip(offsets, raw_profiles):
        color = "#0f4c5c"
        ax.plot(x_vals, yu + offset, color=color, linewidth=1.6)
        ax.plot(x_vals, yl + offset, color=color, linewidth=1.6)
        ax.fill_between(x_vals, yl + offset, yu + offset, color="#dbeafe", alpha=0.65)
        ax.text(
            float(x_vals[-1]) + 0.2,
            offset,
            f"eta={station.eta:.2f}  c={station.chord:.2f} m",
            fontsize=8.6,
            va="center",
            color="#334155",
        )

    ax.set_title(title)
    ax.set_xlabel("local chordwise x [m]")
    ax.set_ylabel("stacked local y [m]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, linewidth=0.35, alpha=0.25)
    ax.set_yticks([])


def plot_lifting_surface_3d(ax, mesh: LiftingSurfaceMesh, title: str = "3D view") -> None:
    stride_span = max(1, mesh.x_upper.shape[0] // 45)
    stride_airfoil = max(1, mesh.x_upper.shape[1] // 45)

    xu = mesh.x_upper[::stride_span, ::stride_airfoil]
    yu = mesh.y_upper[::stride_span, ::stride_airfoil]
    zu = mesh.z_upper[::stride_span, ::stride_airfoil]
    xl = mesh.x_lower[::stride_span, ::stride_airfoil]
    yl = mesh.y_lower[::stride_span, ::stride_airfoil]
    zl = mesh.z_lower[::stride_span, ::stride_airfoil]

    ax.plot_surface(xu, zu, yu, rstride=1, cstride=1, linewidth=0.25, edgecolor="#2563eb", alpha=0.82, color="#dbeafe")
    ax.plot_surface(xl, zl, yl, rstride=1, cstride=1, linewidth=0.25, edgecolor="#2563eb", alpha=0.82, color="#dbeafe")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("z spanwise [m]")
    ax.set_zlabel("y vertical [m]")
    ax.set_title(title)
    ax.view_init(elev=22, azim=-126)

    x_min = float(np.min(mesh.x_all))
    x_max = float(np.max(mesh.x_all))
    z_min = float(np.min(mesh.z_all))
    z_max = float(np.max(mesh.z_all))
    y_min = float(np.min(mesh.y_all))
    y_max = float(np.max(mesh.y_all))
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(z_min, z_max)
    ax.set_zlim(y_min, y_max)
    try:
        ax.set_box_aspect((max(x_max - x_min, 1e-9), max(z_max - z_min, 1e-9), max(y_max - y_min, 1e-9)))
    except Exception:
        pass


def create_lifting_surface_overview_figure(
    prepared: PreparedLiftingSurface,
    title: str | None = None,
    profile_station_etas: tuple[float, ...] | None = None,
    max_section_profiles: int = 4,
) -> tuple[plt.Figure, np.ndarray]:
    mesh = build_lifting_surface_mesh(prepared)
    figure_title = title or prepared.component_id

    fig, axes = plt.subplots(2, 2, figsize=(16.0, 11.0), constrained_layout=True)
    plot_planform(axes[0, 0], mesh, title=f"{figure_title} | top view")
    plot_side_view(axes[0, 1], mesh, title=f"{figure_title} | side view")
    plot_front_view(axes[1, 0], mesh, title=f"{figure_title} | front view")
    plot_section_profiles(
        axes[1, 1],
        prepared,
        title=f"{figure_title} | section profiles",
        selected_etas=profile_station_etas,
        max_profiles=max_section_profiles,
    )
    return fig, axes


def create_lifting_surface_3d_figure(
    prepared: PreparedLiftingSurface,
    title: str | None = None,
) -> tuple[plt.Figure, plt.Axes]:
    mesh = build_lifting_surface_mesh(prepared)
    fig = plt.figure(figsize=(12.0, 8.0), constrained_layout=True)
    ax = fig.add_subplot(111, projection="3d")
    plot_lifting_surface_3d(ax, mesh, title=title or f"{prepared.component_id} | 3D view")
    return fig, ax


def save_lifting_surface_overview(
    prepared: PreparedLiftingSurface,
    out_dir: str | Path,
    stem: str | None = None,
    title: str | None = None,
    profile_station_etas: tuple[float, ...] | None = None,
    max_section_profiles: int = 4,
) -> Path:
    output_dir = Path(out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    file_stem = stem or prepared.component_id
    fig, _ = create_lifting_surface_overview_figure(
        prepared,
        title=title,
        profile_station_etas=profile_station_etas,
        max_section_profiles=max_section_profiles,
    )
    png_path = output_dir / f"{file_stem}_views.png"
    _save_figure(fig, png_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return png_path


def save_lifting_surface_3d(
    prepared: PreparedLiftingSurface,
    out_dir: str | Path,
    stem: str | None = None,
    title: str | None = None,
) -> Path:
    output_dir = Path(out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    file_stem = stem or prepared.component_id
    fig, _ = create_lifting_surface_3d_figure(prepared, title=title)
    png_path = output_dir / f"{file_stem}_3d.png"
    _save_figure(fig, png_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return png_path


def plot_fuselage_top_view(ax, mesh: FuselageMesh, title: str = "Top view") -> None:
    right = np.max(mesh.z, axis=1)
    left = np.min(mesh.z, axis=1)
    x_axis = mesh.x[:, 0]
    ax.fill_between(x_axis, left, right, color="#dbeafe", alpha=0.70)
    ax.plot(x_axis, right, color="#0f4c5c", linewidth=2.2)
    ax.plot(x_axis, left, color="#0f4c5c", linewidth=2.2)
    stride = max(1, mesh.x.shape[0] // 8)
    for idx in range(0, mesh.x.shape[0], stride):
        ax.plot(mesh.x[idx, :], mesh.z[idx, :], color="#93c5fd", linewidth=0.65, alpha=0.45)
    ax.set_title(title)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("z lateral [m]")
    ax.grid(True, linewidth=0.35, alpha=0.25)
    _set_equal_2d(ax, mesh.x_all, mesh.z_all)


def plot_fuselage_side_view(ax, mesh: FuselageMesh, title: str = "Side view") -> None:
    upper = np.max(mesh.y, axis=1)
    lower = np.min(mesh.y, axis=1)
    x_axis = mesh.x[:, 0]
    ax.fill_between(x_axis, lower, upper, color="#dbeafe", alpha=0.70)
    ax.plot(x_axis, upper, color="#0f4c5c", linewidth=2.2)
    ax.plot(x_axis, lower, color="#0f4c5c", linewidth=2.2)
    stride = max(1, mesh.x.shape[0] // 8)
    for idx in range(0, mesh.x.shape[0], stride):
        ax.plot(mesh.x[idx, :], mesh.y[idx, :], color="#93c5fd", linewidth=0.65, alpha=0.45)
    ax.set_title(title)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y vertical [m]")
    ax.grid(True, linewidth=0.35, alpha=0.25)
    _set_equal_2d(ax, mesh.x_all, mesh.y_all)


def plot_fuselage_front_sections(ax, mesh: FuselageMesh, title: str = "Cross sections") -> None:
    n_sections = mesh.x.shape[0]
    selected = sorted(set(np.linspace(0, n_sections - 1, min(6, n_sections), dtype=int).tolist()))
    areas = (np.max(mesh.y, axis=1) - np.min(mesh.y, axis=1)) * (np.max(mesh.z, axis=1) - np.min(mesh.z, axis=1))
    main_idx = int(np.argmax(areas))

    for idx in selected:
        alpha = 0.8 if idx == main_idx else 0.45
        width = 2.1 if idx == main_idx else 1.0
        color = "#0f4c5c" if idx == main_idx else "#93c5fd"
        ax.plot(mesh.z[idx, :], mesh.y[idx, :], color=color, linewidth=width, alpha=alpha)

    ax.set_title(title)
    ax.set_xlabel("z lateral [m]")
    ax.set_ylabel("y vertical [m]")
    ax.grid(True, linewidth=0.35, alpha=0.25)
    _set_equal_2d(ax, mesh.z_all, mesh.y_all)


def plot_fuselage_section_family(
    ax,
    prepared: PreparedFuselage,
    title: str = "Section family",
) -> None:
    stations = prepared.stations
    selected = sorted(set(np.linspace(0, len(stations) - 1, min(6, len(stations)), dtype=int).tolist()))
    height_spans = [float(np.max(station.loop_y) - np.min(station.loop_y)) for station in stations]
    base_gap = max(max(height_spans, default=0.0), 1.0) * 0.55
    offsets = []
    current = 0.0
    for idx in selected:
        station = stations[idx]
        offsets.append((idx, current))
        current += float(np.max(station.loop_y) - np.min(station.loop_y)) + base_gap

    for idx, offset in offsets:
        station = stations[idx]
        ax.plot(station.loop_z, station.loop_y + offset, color="#0f4c5c", linewidth=1.5)
        ax.fill(station.loop_z, station.loop_y + offset, color="#dbeafe", alpha=0.65)
        ax.text(
            float(np.max(station.loop_z)) + 0.15,
            offset,
            f"{station.section_id}  eta={station.eta:.2f}  w={station.width:.2f} m  h={station.height:.2f} m",
            fontsize=8.0,
            va="center",
            color="#334155",
        )

    ax.set_title(title)
    ax.set_xlabel("local z [m]")
    ax.set_ylabel("stacked local y [m]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, linewidth=0.35, alpha=0.25)
    ax.set_yticks([])


def plot_fuselage_3d(ax, mesh: FuselageMesh, title: str = "3D view") -> None:
    stride_u = max(1, mesh.x.shape[0] // 50)
    stride_v = max(1, mesh.x.shape[1] // 60)
    xs = mesh.x[::stride_u, ::stride_v]
    ys = mesh.y[::stride_u, ::stride_v]
    zs = mesh.z[::stride_u, ::stride_v]

    ax.plot_surface(xs, zs, ys, rstride=1, cstride=1, linewidth=0.2, edgecolor="#2563eb", alpha=0.84, color="#dbeafe")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("z lateral [m]")
    ax.set_zlabel("y vertical [m]")
    ax.set_title(title)
    ax.view_init(elev=18, azim=-122)

    x_min = float(np.min(mesh.x_all))
    x_max = float(np.max(mesh.x_all))
    z_min = float(np.min(mesh.z_all))
    z_max = float(np.max(mesh.z_all))
    y_min = float(np.min(mesh.y_all))
    y_max = float(np.max(mesh.y_all))
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(z_min, z_max)
    ax.set_zlim(y_min, y_max)
    try:
        ax.set_box_aspect((max(x_max - x_min, 1e-9), max(z_max - z_min, 1e-9), max(y_max - y_min, 1e-9)))
    except Exception:
        pass


def create_fuselage_overview_figure(
    prepared: PreparedFuselage,
    title: str | None = None,
) -> tuple[plt.Figure, np.ndarray]:
    mesh = build_fuselage_mesh(prepared)
    figure_title = title or prepared.component_id

    fig, axes = plt.subplots(2, 2, figsize=(16.0, 11.0), constrained_layout=True)
    plot_fuselage_top_view(axes[0, 0], mesh, title=f"{figure_title} | top view")
    plot_fuselage_side_view(axes[0, 1], mesh, title=f"{figure_title} | side view")
    plot_fuselage_front_sections(axes[1, 0], mesh, title=f"{figure_title} | cross sections")
    plot_fuselage_section_family(axes[1, 1], prepared, title=f"{figure_title} | section family")
    return fig, axes


def create_fuselage_3d_figure(
    prepared: PreparedFuselage,
    title: str | None = None,
) -> tuple[plt.Figure, plt.Axes]:
    mesh = build_fuselage_mesh(prepared)
    fig = plt.figure(figsize=(12.0, 8.0), constrained_layout=True)
    ax = fig.add_subplot(111, projection="3d")
    plot_fuselage_3d(ax, mesh, title=title or f"{prepared.component_id} | 3D view")
    return fig, ax


def save_fuselage_overview(
    prepared: PreparedFuselage,
    out_dir: str | Path,
    stem: str | None = None,
    title: str | None = None,
) -> tuple[Path, Path]:
    output_dir = Path(out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    file_stem = stem or prepared.component_id
    fig, _ = create_fuselage_overview_figure(prepared, title=title)
    png_path = output_dir / f"{file_stem}_views.png"
    svg_path = output_dir / f"{file_stem}_views.svg"
    _save_figure(fig, png_path, dpi=220, bbox_inches="tight")
    _save_figure(fig, svg_path, bbox_inches="tight")
    plt.close(fig)
    return png_path, svg_path


def save_fuselage_3d(
    prepared: PreparedFuselage,
    out_dir: str | Path,
    stem: str | None = None,
    title: str | None = None,
) -> Path:
    output_dir = Path(out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    file_stem = stem or prepared.component_id
    fig, _ = create_fuselage_3d_figure(prepared, title=title)
    png_path = output_dir / f"{file_stem}_3d.png"
    _save_figure(fig, png_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return png_path


def _mirrored_z(values: np.ndarray) -> np.ndarray:
    return -np.asarray(values, dtype=float)


def _aircraft_bounds(prepared: PreparedAircraftGeometry) -> tuple[np.ndarray, np.ndarray]:
    xyz_min = np.asarray([np.inf, np.inf, np.inf], dtype=float)
    xyz_max = np.asarray([-np.inf, -np.inf, -np.inf], dtype=float)

    for fuselage_entry in prepared.fuselages:
        mesh = build_fuselage_mesh(fuselage_entry.prepared)
        xyz_min = np.minimum(xyz_min, np.asarray([np.min(mesh.x_all), np.min(mesh.y_all), np.min(mesh.z_all)]))
        xyz_max = np.maximum(xyz_max, np.asarray([np.max(mesh.x_all), np.max(mesh.y_all), np.max(mesh.z_all)]))

    for wing_entry in prepared.wings:
        mesh = build_lifting_surface_mesh(wing_entry.prepared)
        xyz_min = np.minimum(xyz_min, np.asarray([np.min(mesh.x_all), np.min(mesh.y_all), np.min(mesh.z_all)]))
        xyz_max = np.maximum(xyz_max, np.asarray([np.max(mesh.x_all), np.max(mesh.y_all), np.max(mesh.z_all)]))

    for tail_entry in prepared.vertical_tails:
        mesh = build_lifting_surface_mesh(tail_entry.prepared)
        xyz_min = np.minimum(xyz_min, np.asarray([np.min(mesh.x_all), np.min(mesh.y_all), np.min(mesh.z_all)]))
        xyz_max = np.maximum(xyz_max, np.asarray([np.max(mesh.x_all), np.max(mesh.y_all), np.max(mesh.z_all)]))

    return xyz_min, xyz_max


def plot_aircraft_top_view(ax, prepared: PreparedAircraftGeometry, title: str = "Top view") -> None:
    for fuselage_entry in prepared.fuselages:
        mesh = build_fuselage_mesh(fuselage_entry.prepared)
        right = np.max(mesh.z, axis=1)
        left = np.min(mesh.z, axis=1)
        x_axis = mesh.x[:, 0]
        ax.fill_between(x_axis, left, right, color="#dbeafe", alpha=0.75)
        ax.plot(x_axis, right, color="#0f4c5c", linewidth=2.0)
        ax.plot(x_axis, left, color="#0f4c5c", linewidth=2.0)

    for wing_entry in prepared.wings:
        mesh = build_lifting_surface_mesh(wing_entry.prepared)
        ax.plot(mesh.x_upper[:, 0], mesh.z_upper[:, 0], color="#1d4ed8", linewidth=1.8)
        ax.plot(mesh.x_upper[:, -1], mesh.z_upper[:, -1], color="#dc2626", linewidth=1.8)
        for idx in range(mesh.x_upper.shape[0]):
            ax.plot(mesh.x_upper[idx, :], mesh.z_upper[idx, :], color="#93c5fd", linewidth=0.5, alpha=0.30)
            ax.plot(mesh.x_lower[idx, :], mesh.z_lower[idx, :], color="#93c5fd", linewidth=0.5, alpha=0.30)

    for tail_entry in prepared.vertical_tails:
        mesh = build_lifting_surface_mesh(tail_entry.prepared)
        for idx in range(mesh.x_upper.shape[0]):
            ax.plot(mesh.x_upper[idx, :], mesh.z_upper[idx, :], color="#60a5fa", linewidth=0.55, alpha=0.35)
            ax.plot(mesh.x_lower[idx, :], mesh.z_lower[idx, :], color="#60a5fa", linewidth=0.55, alpha=0.35)

    xyz_min, xyz_max = _aircraft_bounds(prepared)
    ax.set_title(title)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("z lateral [m]")
    ax.grid(True, linewidth=0.35, alpha=0.25)
    _set_equal_2d(ax, np.asarray([xyz_min[0], xyz_max[0]]), np.asarray([xyz_min[2], xyz_max[2]]))


def plot_aircraft_side_view(ax, prepared: PreparedAircraftGeometry, title: str = "Side view") -> None:
    for fuselage_entry in prepared.fuselages:
        mesh = build_fuselage_mesh(fuselage_entry.prepared)
        upper = np.max(mesh.y, axis=1)
        lower = np.min(mesh.y, axis=1)
        x_axis = mesh.x[:, 0]
        ax.fill_between(x_axis, lower, upper, color="#dbeafe", alpha=0.75)
        ax.plot(x_axis, upper, color="#0f4c5c", linewidth=2.0)
        ax.plot(x_axis, lower, color="#0f4c5c", linewidth=2.0)

    for wing_entry in prepared.wings:
        mesh = build_lifting_surface_mesh(wing_entry.prepared)
        stride = max(1, mesh.x_upper.shape[0] // 8)
        for idx in range(0, mesh.x_upper.shape[0], stride):
            ax.plot(mesh.x_upper[idx, :], mesh.y_upper[idx, :], color="#1d4ed8", linewidth=0.8, alpha=0.55)
            ax.plot(mesh.x_lower[idx, :], mesh.y_lower[idx, :], color="#1d4ed8", linewidth=0.8, alpha=0.55)

    for tail_entry in prepared.vertical_tails:
        mesh = build_lifting_surface_mesh(tail_entry.prepared)
        stride = max(1, mesh.x_upper.shape[0] // 8)
        for idx in range(0, mesh.x_upper.shape[0], stride):
            ax.plot(mesh.x_upper[idx, :], mesh.y_upper[idx, :], color="#2563eb", linewidth=0.8, alpha=0.55)
            ax.plot(mesh.x_lower[idx, :], mesh.y_lower[idx, :], color="#2563eb", linewidth=0.8, alpha=0.55)

    xyz_min, xyz_max = _aircraft_bounds(prepared)
    ax.set_title(title)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y vertical [m]")
    ax.grid(True, linewidth=0.35, alpha=0.25)
    _set_equal_2d(ax, np.asarray([xyz_min[0], xyz_max[0]]), np.asarray([xyz_min[1], xyz_max[1]]))


def plot_aircraft_front_view(ax, prepared: PreparedAircraftGeometry, title: str = "Front view") -> None:
    for fuselage_entry in prepared.fuselages:
        mesh = build_fuselage_mesh(fuselage_entry.prepared)
        areas = (np.max(mesh.y, axis=1) - np.min(mesh.y, axis=1)) * (np.max(mesh.z, axis=1) - np.min(mesh.z, axis=1))
        idx = int(np.argmax(areas))
        ax.plot(mesh.z[idx, :], mesh.y[idx, :], color="#0f4c5c", linewidth=2.0)
        ax.fill(mesh.z[idx, :], mesh.y[idx, :], color="#dbeafe", alpha=0.55)

    for wing_entry in prepared.wings:
        mesh = build_lifting_surface_mesh(wing_entry.prepared)
        stride = max(1, mesh.x_upper.shape[0] // 10)
        for idx in range(0, mesh.x_upper.shape[0], stride):
            ax.plot(mesh.z_upper[idx, :], mesh.y_upper[idx, :], color="#1d4ed8", linewidth=0.8, alpha=0.45)
            ax.plot(mesh.z_lower[idx, :], mesh.y_lower[idx, :], color="#1d4ed8", linewidth=0.8, alpha=0.45)

    for tail_entry in prepared.vertical_tails:
        mesh = build_lifting_surface_mesh(tail_entry.prepared)
        stride = max(1, mesh.x_upper.shape[0] // 10)
        for idx in range(0, mesh.x_upper.shape[0], stride):
            ax.plot(mesh.z_upper[idx, :], mesh.y_upper[idx, :], color="#2563eb", linewidth=0.8, alpha=0.45)
            ax.plot(mesh.z_lower[idx, :], mesh.y_lower[idx, :], color="#2563eb", linewidth=0.8, alpha=0.45)

    xyz_min, xyz_max = _aircraft_bounds(prepared)
    ax.set_title(title)
    ax.set_xlabel("z lateral [m]")
    ax.set_ylabel("y vertical [m]")
    ax.grid(True, linewidth=0.35, alpha=0.25)
    _set_equal_2d(ax, np.asarray([xyz_min[2], xyz_max[2]]), np.asarray([xyz_min[1], xyz_max[1]]))


def plot_aircraft_3d(ax, prepared: PreparedAircraftGeometry, title: str = "3D view") -> None:
    for fuselage_entry in prepared.fuselages:
        mesh = build_fuselage_mesh(fuselage_entry.prepared)
        stride_u = max(1, mesh.x.shape[0] // 50)
        stride_v = max(1, mesh.x.shape[1] // 60)
        xs = mesh.x[::stride_u, ::stride_v]
        ys = mesh.y[::stride_u, ::stride_v]
        zs = mesh.z[::stride_u, ::stride_v]
        ax.plot_surface(xs, zs, ys, rstride=1, cstride=1, linewidth=0.2, edgecolor="#2563eb", alpha=0.84, color="#dbeafe")

    for wing_entry in prepared.wings:
        mesh = build_lifting_surface_mesh(wing_entry.prepared)
        stride_span = max(1, mesh.x_upper.shape[0] // 45)
        stride_airfoil = max(1, mesh.x_upper.shape[1] // 45)
        xu = mesh.x_upper[::stride_span, ::stride_airfoil]
        yu = mesh.y_upper[::stride_span, ::stride_airfoil]
        zu = mesh.z_upper[::stride_span, ::stride_airfoil]
        xl = mesh.x_lower[::stride_span, ::stride_airfoil]
        yl = mesh.y_lower[::stride_span, ::stride_airfoil]
        zl = mesh.z_lower[::stride_span, ::stride_airfoil]
        ax.plot_surface(xu, zu, yu, rstride=1, cstride=1, linewidth=0.2, edgecolor="#1d4ed8", alpha=0.80, color="#e0f2fe")
        ax.plot_surface(xl, zl, yl, rstride=1, cstride=1, linewidth=0.2, edgecolor="#1d4ed8", alpha=0.80, color="#e0f2fe")

    for tail_entry in prepared.vertical_tails:
        mesh = build_lifting_surface_mesh(tail_entry.prepared)
        stride_span = max(1, mesh.x_upper.shape[0] // 45)
        stride_airfoil = max(1, mesh.x_upper.shape[1] // 45)
        xu = mesh.x_upper[::stride_span, ::stride_airfoil]
        yu = mesh.y_upper[::stride_span, ::stride_airfoil]
        zu = mesh.z_upper[::stride_span, ::stride_airfoil]
        xl = mesh.x_lower[::stride_span, ::stride_airfoil]
        yl = mesh.y_lower[::stride_span, ::stride_airfoil]
        zl = mesh.z_lower[::stride_span, ::stride_airfoil]
        ax.plot_surface(xu, zu, yu, rstride=1, cstride=1, linewidth=0.2, edgecolor="#2563eb", alpha=0.82, color="#bfdbfe")
        ax.plot_surface(xl, zl, yl, rstride=1, cstride=1, linewidth=0.2, edgecolor="#2563eb", alpha=0.82, color="#bfdbfe")

    xyz_min, xyz_max = _aircraft_bounds(prepared)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("z lateral [m]")
    ax.set_zlabel("y vertical [m]")
    ax.set_title(title)
    ax.view_init(elev=19, azim=-124)
    ax.set_xlim(float(xyz_min[0]), float(xyz_max[0]))
    ax.set_ylim(float(xyz_min[2]), float(xyz_max[2]))
    ax.set_zlim(float(xyz_min[1]), float(xyz_max[1]))
    try:
        ax.set_box_aspect(
            (
                max(float(xyz_max[0] - xyz_min[0]), 1e-9),
                max(float(xyz_max[2] - xyz_min[2]), 1e-9),
                max(float(xyz_max[1] - xyz_min[1]), 1e-9),
            )
        )
    except Exception:
        pass


def create_aircraft_overview_figure(
    prepared: PreparedAircraftGeometry,
    title: str | None = None,
) -> tuple[plt.Figure, np.ndarray]:
    figure_title = title or prepared.aircraft_id
    fig = plt.figure(figsize=(16.0, 11.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 2)
    ax_top = fig.add_subplot(gs[0, 0])
    ax_side = fig.add_subplot(gs[0, 1])
    ax_front = fig.add_subplot(gs[1, 0])
    ax_3d = fig.add_subplot(gs[1, 1], projection="3d")
    plot_aircraft_top_view(ax_top, prepared, title=f"{figure_title} | top view")
    plot_aircraft_side_view(ax_side, prepared, title=f"{figure_title} | side view")
    plot_aircraft_front_view(ax_front, prepared, title=f"{figure_title} | front view")
    plot_aircraft_3d(ax_3d, prepared, title=f"{figure_title} | 3D view")
    return fig, np.asarray([[ax_top, ax_side], [ax_front, ax_3d]], dtype=object)


def create_aircraft_3d_figure(
    prepared: PreparedAircraftGeometry,
    title: str | None = None,
) -> tuple[plt.Figure, plt.Axes]:
    fig = plt.figure(figsize=(12.0, 8.0), constrained_layout=True)
    ax = fig.add_subplot(111, projection="3d")
    plot_aircraft_3d(ax, prepared, title=title or f"{prepared.aircraft_id} | 3D view")
    return fig, ax


def save_aircraft_overview(
    prepared: PreparedAircraftGeometry,
    out_dir: str | Path,
    stem: str | None = None,
    title: str | None = None,
) -> tuple[Path, Path]:
    output_dir = Path(out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    file_stem = stem or prepared.aircraft_id
    fig, _ = create_aircraft_overview_figure(prepared, title=title)
    png_path = output_dir / f"{file_stem}_views.png"
    svg_path = output_dir / f"{file_stem}_views.svg"
    _save_figure(fig, png_path, dpi=220, bbox_inches="tight")
    _save_figure(fig, svg_path, bbox_inches="tight")
    plt.close(fig)
    return png_path, svg_path


def save_aircraft_3d(
    prepared: PreparedAircraftGeometry,
    out_dir: str | Path,
    stem: str | None = None,
    title: str | None = None,
) -> Path:
    output_dir = Path(out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    file_stem = stem or prepared.aircraft_id
    fig, _ = create_aircraft_3d_figure(prepared, title=title)
    png_path = output_dir / f"{file_stem}_3d.png"
    _save_figure(fig, png_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return png_path
