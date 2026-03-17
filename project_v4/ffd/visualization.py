from pathlib import Path
from typing import Optional
import os
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib-bwb"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


def _to_plot_coords(points: np.ndarray) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    return points[..., [0, 2, 1]]


def _plot_ffd_lattice_3d(ax, points: np.ndarray, color: str, alpha: float, linewidth: float) -> None:
    points = _to_plot_coords(points)
    ni, nj, nk, _ = points.shape
    for j in range(nj):
        for k in range(nk):
            ax.plot(points[:, j, k, 0], points[:, j, k, 1], points[:, j, k, 2], color=color, alpha=alpha, lw=linewidth)
    for i in range(ni):
        for k in range(nk):
            ax.plot(points[i, :, k, 0], points[i, :, k, 1], points[i, :, k, 2], color=color, alpha=alpha, lw=linewidth)
    for i in range(ni):
        for j in range(nj):
            ax.plot(points[i, j, :, 0], points[i, j, :, 1], points[i, j, :, 2], color=color, alpha=alpha, lw=linewidth)


def _plot_ffd_lattice_top(ax, points: np.ndarray, color: str, alpha: float, linewidth: float) -> None:
    ni, nj, nk, _ = points.shape
    for j in range(nj):
        for k in range(nk):
            ax.plot(points[:, j, k, 2], points[:, j, k, 0], color=color, alpha=alpha, lw=linewidth)
    for i in range(ni):
        for k in range(nk):
            ax.plot(points[i, :, k, 2], points[i, :, k, 0], color=color, alpha=alpha, lw=linewidth)
    for i in range(ni):
        for j in range(nj):
            ax.plot(points[i, j, :, 2], points[i, j, :, 0], color=color, alpha=alpha, lw=linewidth)


def _plot_surface_wireframe_3d(ax, upper: np.ndarray, lower: np.ndarray, color: str, alpha: float) -> None:
    upper = _to_plot_coords(upper)
    lower = _to_plot_coords(lower)
    for grid in (upper, lower):
        for idx in range(grid.shape[0]):
            ax.plot(grid[idx, :, 0], grid[idx, :, 1], grid[idx, :, 2], color=color, alpha=alpha, lw=0.45)
        for idx in range(grid.shape[1]):
            ax.plot(grid[:, idx, 0], grid[:, idx, 1], grid[:, idx, 2], color=color, alpha=alpha * 0.75, lw=0.35)


def _plot_surface_shell_3d(ax, upper: np.ndarray, lower: np.ndarray, color: str, alpha: float) -> None:
    upper = _to_plot_coords(upper)
    lower = _to_plot_coords(lower)
    for grid in (upper, lower):
        ax.plot_surface(
            grid[:, :, 0],
            grid[:, :, 1],
            grid[:, :, 2],
            color=color,
            alpha=alpha,
            linewidth=0.0,
            shade=False,
            antialiased=True,
        )


def _plot_surface_top(ax, upper: np.ndarray, lower: np.ndarray, color: str, alpha: float) -> None:
    ax.plot(upper[:, 0, 2], upper[:, 0, 0], color=color, alpha=alpha, lw=1.4)
    ax.plot(upper[:, -1, 2], upper[:, -1, 0], color=color, alpha=alpha, lw=1.4)
    ax.plot(lower[:, -1, 2], lower[:, -1, 0], color=color, alpha=alpha, lw=1.0)
    for idx in range(0, upper.shape[0], max(1, upper.shape[0] // 12)):
        ax.plot(upper[idx, :, 2], upper[idx, :, 0], color=color, alpha=alpha * 0.45, lw=0.7)


def _scatter_points_3d(
    ax,
    points: np.ndarray,
    color: str,
    alpha: float,
    size: float,
    step: int = 1,
) -> None:
    sampled = _to_plot_coords(points).reshape(-1, 3)[:: max(1, int(step))]
    ax.scatter(sampled[:, 0], sampled[:, 1], sampled[:, 2], s=size, c=color, alpha=alpha, depthshade=False)


def _scatter_points_top(
    ax,
    points: np.ndarray,
    color: str,
    alpha: float,
    size: float,
    step: int = 1,
) -> None:
    sampled = points.reshape(-1, 3)[:: max(1, int(step))]
    ax.scatter(sampled[:, 2], sampled[:, 0], s=size, c=color, alpha=alpha)


def _set_equal_3d_box(ax, point_groups: list[np.ndarray]) -> None:
    stacked = np.vstack([_to_plot_coords(points).reshape(-1, 3) for points in point_groups])
    mins = np.min(stacked, axis=0)
    maxs = np.max(stacked, axis=0)
    centers = 0.5 * (mins + maxs)
    radius = 0.5 * np.max(maxs - mins)
    ax.set_xlim(centers[0] - radius, centers[0] + radius)
    ax.set_ylim(centers[1] - radius, centers[1] + radius)
    ax.set_zlim(centers[2] - radius, centers[2] + radius)
    ax.set_box_aspect((1.25, 1.35, 0.45))


def _draw_ffd_scene_3d(
    ax,
    ffd_points: np.ndarray,
    upper_grid: np.ndarray,
    lower_grid: np.ndarray,
    baseline_pointset: Optional[np.ndarray] = None,
    deformed_ffd_points: Optional[np.ndarray] = None,
    deformed_pointset: Optional[np.ndarray] = None,
    deformed_upper_grid: Optional[np.ndarray] = None,
    deformed_lower_grid: Optional[np.ndarray] = None,
    frame_label: Optional[str] = None,
) -> None:
    _plot_surface_shell_3d(ax, upper_grid, lower_grid, color="#D7DDD9", alpha=0.22)
    _plot_surface_wireframe_3d(ax, upper_grid, lower_grid, color="#6C7A80", alpha=0.55)
    _plot_ffd_lattice_3d(ax, ffd_points, color="#B51E3C", alpha=0.80, linewidth=0.55)
    _scatter_points_3d(ax, ffd_points, color="#B51E3C", alpha=0.90, size=6.0, step=1)

    point_groups = [ffd_points, upper_grid, lower_grid]
    if baseline_pointset is not None:
        _scatter_points_3d(ax, baseline_pointset, color="#6E7E83", alpha=0.10, size=1.2, step=4)
        point_groups.append(baseline_pointset)

    if (
        deformed_ffd_points is not None
        and deformed_pointset is not None
        and deformed_upper_grid is not None
        and deformed_lower_grid is not None
    ):
        _plot_surface_shell_3d(ax, deformed_upper_grid, deformed_lower_grid, color="#E2C8CD", alpha=0.18)
        _plot_surface_wireframe_3d(ax, deformed_upper_grid, deformed_lower_grid, color="#B56171", alpha=0.55)
        _plot_ffd_lattice_3d(ax, deformed_ffd_points, color="#D11945", alpha=0.92, linewidth=0.60)
        _scatter_points_3d(ax, deformed_ffd_points, color="#D11945", alpha=0.95, size=6.0, step=1)
        _scatter_points_3d(ax, deformed_pointset, color="#A15A68", alpha=0.10, size=1.2, step=4)
        point_groups.extend([deformed_ffd_points, deformed_pointset, deformed_upper_grid, deformed_lower_grid])

    _set_equal_3d_box(ax, point_groups)
    ax.view_init(elev=24, azim=-66)
    ax.set_xlabel("x")
    ax.set_ylabel("spanwise")
    ax.set_zlabel("vertical")
    ax.set_title("DVGeometry FFD deformation")
    if frame_label:
        ax.text2D(0.04, 0.96, frame_label, transform=ax.transAxes)


def plot_ffd_3d(
    ffd_points: np.ndarray,
    upper_grid: np.ndarray,
    lower_grid: np.ndarray,
    out_path: Path,
    title: str,
    baseline_pointset: Optional[np.ndarray] = None,
    deformed_ffd_points: Optional[np.ndarray] = None,
    deformed_pointset: Optional[np.ndarray] = None,
    deformed_upper_grid: Optional[np.ndarray] = None,
    deformed_lower_grid: Optional[np.ndarray] = None,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(8.8, 7.0))
    ax_3d = fig.add_subplot(1, 1, 1, projection="3d")
    _draw_ffd_scene_3d(
        ax_3d,
        ffd_points=ffd_points,
        upper_grid=upper_grid,
        lower_grid=lower_grid,
        baseline_pointset=baseline_pointset,
        deformed_ffd_points=deformed_ffd_points,
        deformed_pointset=deformed_pointset,
        deformed_upper_grid=deformed_upper_grid,
        deformed_lower_grid=deformed_lower_grid,
    )
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _render_ffd_scene_image(
    ffd_points: np.ndarray,
    upper_grid: np.ndarray,
    lower_grid: np.ndarray,
    title: str,
    baseline_pointset: Optional[np.ndarray] = None,
    deformed_ffd_points: Optional[np.ndarray] = None,
    deformed_pointset: Optional[np.ndarray] = None,
    deformed_upper_grid: Optional[np.ndarray] = None,
    deformed_lower_grid: Optional[np.ndarray] = None,
    frame_label: Optional[str] = None,
) -> Image.Image:
    fig = plt.figure(figsize=(8.8, 7.0))
    ax_3d = fig.add_subplot(1, 1, 1, projection="3d")
    _draw_ffd_scene_3d(
        ax_3d,
        ffd_points=ffd_points,
        upper_grid=upper_grid,
        lower_grid=lower_grid,
        baseline_pointset=baseline_pointset,
        deformed_ffd_points=deformed_ffd_points,
        deformed_pointset=deformed_pointset,
        deformed_upper_grid=deformed_upper_grid,
        deformed_lower_grid=deformed_lower_grid,
        frame_label=frame_label,
    )
    fig.suptitle(title)
    fig.tight_layout()
    fig.canvas.draw()
    width, height = fig.canvas.get_width_height()
    image = Image.frombuffer(
        "RGBA",
        (width, height),
        fig.canvas.buffer_rgba(),
        "raw",
        "RGBA",
        0,
        1,
    ).convert("P")
    plt.close(fig)
    return image


def write_ffd_deformation_gif(
    out_path: Path,
    title: str,
    frames: list[dict],
    duration_ms: int = 90,
) -> None:
    if not frames:
        raise ValueError("frames must contain at least one deformation state")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    images = [
        _render_ffd_scene_image(
            ffd_points=frame["baseline_ffd_points"],
            upper_grid=frame["baseline_upper_grid"],
            lower_grid=frame["baseline_lower_grid"],
            title=title,
            baseline_pointset=frame.get("baseline_pointset"),
            deformed_ffd_points=frame.get("deformed_ffd_points"),
            deformed_pointset=frame.get("deformed_pointset"),
            deformed_upper_grid=frame.get("deformed_upper_grid"),
            deformed_lower_grid=frame.get("deformed_lower_grid"),
            frame_label=frame.get("frame_label"),
        )
        for frame in frames
    ]
    images[0].save(
        out_path,
        save_all=True,
        append_images=images[1:],
        loop=0,
        duration=int(duration_ms),
        optimize=False,
    )


def plot_ffd_with_surface(
    ffd_points: np.ndarray,
    upper_grid: np.ndarray,
    lower_grid: np.ndarray,
    out_path: Path,
    title: str,
    baseline_pointset: Optional[np.ndarray] = None,
    deformed_ffd_points: Optional[np.ndarray] = None,
    deformed_pointset: Optional[np.ndarray] = None,
    deformed_upper_grid: Optional[np.ndarray] = None,
    deformed_lower_grid: Optional[np.ndarray] = None,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(14, 6.8))
    ax_3d = fig.add_subplot(1, 2, 1, projection="3d")
    ax_top = fig.add_subplot(1, 2, 2)

    _plot_surface_shell_3d(ax_3d, upper_grid, lower_grid, color="#D7DDD9", alpha=0.22)
    _plot_surface_wireframe_3d(ax_3d, upper_grid, lower_grid, color="#6C7A80", alpha=0.62)
    _plot_ffd_lattice_3d(ax_3d, ffd_points, color="#B51E3C", alpha=0.82, linewidth=0.55)
    _plot_surface_top(ax_top, upper_grid, lower_grid, color="#0F5D73", alpha=0.90)
    _plot_ffd_lattice_top(ax_top, ffd_points, color="#B51E3C", alpha=0.80, linewidth=0.55)
    _scatter_points_3d(ax_3d, ffd_points, color="#B51E3C", alpha=0.90, size=6.0, step=1)
    _scatter_points_top(ax_top, ffd_points, color="#B51E3C", alpha=0.90, size=4.0, step=1)

    if baseline_pointset is not None:
        _scatter_points_3d(ax_3d, baseline_pointset, color="#6E7E83", alpha=0.10, size=1.2, step=4)
        _scatter_points_top(ax_top, baseline_pointset, color="#6E7E83", alpha=0.10, size=1.0, step=4)

    point_groups = [ffd_points, upper_grid, lower_grid]
    if (
        deformed_ffd_points is not None
        and deformed_pointset is not None
        and deformed_upper_grid is not None
        and deformed_lower_grid is not None
    ):
        _plot_surface_shell_3d(ax_3d, deformed_upper_grid, deformed_lower_grid, color="#E2C8CD", alpha=0.18)
        _plot_surface_wireframe_3d(ax_3d, deformed_upper_grid, deformed_lower_grid, color="#B56171", alpha=0.55)
        _plot_ffd_lattice_3d(ax_3d, deformed_ffd_points, color="#D11945", alpha=0.90, linewidth=0.60)
        _plot_surface_top(ax_top, deformed_upper_grid, deformed_lower_grid, color="#D1495B", alpha=0.70)
        _plot_ffd_lattice_top(ax_top, deformed_ffd_points, color="#D11945", alpha=0.82, linewidth=0.60)
        _scatter_points_3d(ax_3d, deformed_ffd_points, color="#D11945", alpha=0.95, size=6.0, step=1)
        _scatter_points_top(ax_top, deformed_ffd_points, color="#D11945", alpha=0.95, size=4.0, step=1)
        _scatter_points_3d(ax_3d, deformed_pointset, color="#A15A68", alpha=0.10, size=1.2, step=4)
        _scatter_points_top(ax_top, deformed_pointset, color="#A15A68", alpha=0.10, size=1.0, step=4)
        point_groups.extend([deformed_ffd_points, deformed_upper_grid, deformed_lower_grid, deformed_pointset])
    if baseline_pointset is not None:
        point_groups.append(baseline_pointset)

    _set_equal_3d_box(ax_3d, point_groups)
    ax_3d.view_init(elev=24, azim=-66)
    ax_3d.set_xlabel("x")
    ax_3d.set_ylabel("spanwise")
    ax_3d.set_zlabel("vertical")
    ax_3d.set_title("FFD box and wing wireframe")

    ax_top.set_xlabel("spanwise z")
    ax_top.set_ylabel("streamwise x")
    ax_top.set_title("Top projection")
    ax_top.set_aspect("equal", adjustable="box")
    ax_top.grid(True, alpha=0.25)

    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
