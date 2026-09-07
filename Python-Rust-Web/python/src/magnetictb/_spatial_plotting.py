"""Matplotlib display of Rust-generated crystal/BZ geometry, without recompilation.

Only artist construction, sphere tessellation, and viewport layout live here.
Sites, cell edges, magnetic arrow endpoints, BZ facets, and folded k paths are
consumed in their returned order. No geometry or symmetry is reconstructed.
"""

from __future__ import annotations

import io
import math
from collections.abc import Sequence
from typing import Any


def _image_size(value: Any) -> tuple[float, float]:
    """Convert the display-size option to pixel dimensions, not geometry units."""
    if value is None:
        value = "Automatic"
    if isinstance(value, str):
        presets = {
            "Automatic": 600,
            "Tiny": 200,
            "Small": 320,
            "Medium": 480,
            "Large": 700,
        }
        if value not in presets:
            raise ValueError(
                "ImageSize must be Automatic/Tiny/Small/Medium/Large, a width, or (width, height)"
            )
        value = presets[value]
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        if len(value) != 2:
            raise ValueError("ImageSize must contain exactly width and height")
        dimensions = tuple(value)
    else:
        dimensions = (value, value)
    if any(isinstance(item, (bool, complex)) for item in dimensions):
        raise ValueError("ImageSize dimensions must be positive finite real numbers")
    try:
        width, height = (float(item) for item in dimensions)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError(
            "ImageSize dimensions must be positive finite real numbers"
        ) from error
    if any(not math.isfinite(item) or item <= 0 for item in (width, height)):
        raise ValueError("ImageSize dimensions must be positive finite real numbers")
    return width, height


def _equal_limits(ax: Any, points: Any) -> None:
    """Equal Cartesian scale, including arrowheads and physical atom radii."""
    lower = [min(point[axis] for point in points) for axis in range(3)]
    upper = [max(point[axis] for point in points) for axis in range(3)]
    span = max(high - low for low, high in zip(lower, upper))
    half = (span if span > 0 else 1.0) * 0.56
    for setter, low, high in zip((ax.set_xlim, ax.set_ylim, ax.set_zlim), lower, upper):
        center = (low + high) / 2
        setter(center - half, center + half)
    ax.set_box_aspect((1, 1, 1))


def _moment_arrow(ax: Any, start: Any, end: Any, index: int) -> Any:
    """Draw returned endpoints with a solid, view-independent arrowhead.

    Axes3D.quiver uses two planar head strokes. For an in-plane moment their
    plane becomes edge-on in a strict top view. A closed cone has nonzero
    projected area from every direction, including an end-on view (a disk).
    Only display tessellation is performed; neither endpoint is changed.
    """
    import numpy as np
    from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection

    start, end = np.asarray(start, dtype=float), np.asarray(end, dtype=float)
    displacement = end - start
    length = math.hypot(*displacement)
    if length == 0:
        return ()  # MomentScale=0 must not invent a visible arrow or direction.
    direction = displacement / length
    # Pick a nonparallel auxiliary axis even for x/y/z-aligned moments.
    auxiliary = np.eye(3)[np.argmin(np.abs(direction))]
    transverse = np.cross(direction, auxiliary)
    transverse /= np.linalg.norm(transverse)
    other = np.cross(direction, transverse)
    # A head narrower than the 2-point shaft still looks like a line in a
    # normal-size multi-cell figure. Widen the glyph, not the moment vector:
    # the cone occupies the final quarter of the unchanged shaft length.
    head_length = 0.25 * length
    base = end - head_length * direction
    radius = 0.125 * length
    angles = np.linspace(0, 2 * np.pi, 24, endpoint=False)
    ring = base + radius * (
        np.cos(angles)[:, None] * transverse + np.sin(angles)[:, None] * other
    )
    faces = []
    for vertex, successor in zip(ring, np.roll(ring, -1, axis=0)):
        faces.append((end, vertex, successor))
        faces.append((base, successor, vertex))

    shaft = Line3DCollection(((start, end),), colors="#d62728", linewidths=2)
    shaft.set_gid(f"crystal-moment-{index}")
    ax.add_collection3d(shaft)
    head = Poly3DCollection(faces, facecolors="#d62728", linewidths=0)
    head.set_gid(f"crystal-moment-head-{index}")
    ax.add_collection3d(head)
    return ring


class _SpatialPlot:
    """Shared display methods; no fields are added to the Web data record."""

    def _draw(self, ax: Any) -> None:
        raise NotImplementedError

    def _configure(self, ax: Any, elev: float, azim: float) -> None:
        if getattr(ax, "name", None) != "3d":
            raise TypeError("plot requires a Matplotlib 3D Axes (projection='3d')")
        if any(
            isinstance(value, bool) or not math.isfinite(float(value))
            for value in (elev, azim)
        ):
            raise ValueError("elev and azim must be finite angles in degrees")
        ax.set_proj_type("ortho")
        ax.view_init(elev=float(elev), azim=float(azim))
        ax.set_axis_off()

    def _figure(self, elev: float = 24, azim: float = -60) -> Any:
        from matplotlib.figure import Figure

        figure = Figure(
            figsize=tuple(size / 100 for size in _image_size(self.image_size)), dpi=100
        )
        ax = figure.add_subplot(projection="3d")
        self._configure(ax, elev, azim)
        self._draw(ax)
        figure.subplots_adjust(left=0, bottom=0, right=1, top=1)
        return figure

    def plot(self, ax: Any = None, *, elev: float = 24, azim: float = -60) -> Any:
        """Draw and return a 3D Axes with orthographic projection.

        Mouse rotation/zoom uses the active Matplotlib interactive backend.
        Jupyter's default inline SVG is a static view; an interactive backend
        may be selected separately without changing the computed result.
        """
        if ax is None:
            from matplotlib import pyplot as plt

            figure = plt.figure(
                figsize=tuple(size / 100 for size in _image_size(self.image_size)),
                dpi=100,
            )
            figure.subplots_adjust(left=0, bottom=0, right=1, top=1)
            ax = figure.add_subplot(projection="3d")
        self._configure(ax, elev, azim)
        self._draw(ax)
        return ax

    def show(
        self, *, block: bool | None = None, elev: float = 24, azim: float = -60
    ) -> None:
        """Open/display a rotatable view using the current Matplotlib backend."""
        from matplotlib import pyplot as plt

        self.plot(elev=elev, azim=azim)
        plt.show(block=block)

    def savefig(
        self, filename: Any, *, elev: float = 24, azim: float = -60, **options: Any
    ) -> None:
        """Save a chosen 3D view without opening a window or calling Rust again."""
        self._figure(elev=elev, azim=azim).savefig(filename, **options)

    def _repr_svg_(self) -> str:
        output = io.StringIO()
        self._figure().savefig(output, format="svg", metadata={"Date": None})
        return output.getvalue()


class _CrystalPlot(_SpatialPlot):
    def _draw(self, ax: Any) -> None:
        import numpy as np
        from matplotlib import colormaps
        from mpl_toolkits.mplot3d.art3d import Line3DCollection

        edges = Line3DCollection(self.cell_edges, colors="#595959", linewidths=1.2)
        edges.set_gid("crystal-cell-edges")
        ax.add_collection3d(edges)
        points = [point for edge in self.cell_edges for point in edge]
        # Tessellating a display sphere does not generate or alter atom sites.
        longitude, latitude = np.meshgrid(
            np.linspace(0, 2 * np.pi, 25), np.linspace(0, np.pi, 17)
        )
        sphere = (
            np.cos(longitude) * np.sin(latitude),
            np.sin(longitude) * np.sin(latitude),
            np.cos(latitude),
        )
        palette = colormaps["tab10"]
        radius = self.atom_radius
        for index, record in enumerate(self.atom_records):
            center = record.cartesian_position
            surface = ax.plot_surface(
                *(center[axis] + radius * sphere[axis] for axis in range(3)),
                color=palette((record.orbit_index - 1) % 10),
                linewidth=0,
                antialiased=True,
                shade=True,
            )
            surface.set_gid(f"crystal-atom-{index}")
            points.extend(
                tuple(value + sign * radius for value in center) for sign in (-1, 1)
            )
            if record.magnetic and record.arrow_end is not None:
                points.extend(_moment_arrow(ax, center, record.arrow_end, index))
                points.append(record.arrow_end)
            if self.show_atom_labels:
                label_position = (center[0], center[1], center[2] + 1.5 * radius)
                label = ax.text(
                    *label_position,
                    f"{record.orbit_index}.{record.equivalent_index}",
                    color="black",
                    fontsize=self.font_size,
                    ha="center",
                )
                label.set_gid(f"crystal-label-{index}")
                points.append(label_position)
        _equal_limits(ax, points)


class _BrillouinPlot(_SpatialPlot):
    def _draw(self, ax: Any) -> None:
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection

        faces = [[self.vertices[index] for index in facet] for facet in self.facets]
        mesh = Poly3DCollection(
            faces,
            facecolors=(0.68, 0.85, 0.90, 0.18),
            edgecolors="#404040",
            linewidths=1.2,
        )
        mesh.set_gid("brillouin-zone-facets")
        ax.add_collection3d(mesh)
        points = list(self.vertices)
        if self.show_k_path:
            labels_seen = set()
            # Cartesian endpoints already include stable first-BZ folding.
            # Use the paired labels only; do not transform fractional points here.
            for index, (segment, displayed) in enumerate(
                zip(self.cartesian_k_path, self.displayed_k_path)
            ):
                (line,) = ax.plot(*zip(*segment), color="#d62728", linewidth=2)
                line.set_gid(f"brillouin-path-{index}")
                points.extend(segment)
                for point, text in zip(segment, displayed[1]):
                    key = (str(text), tuple(point))
                    if key in labels_seen:
                        continue
                    labels_seen.add(key)
                    ax.scatter(*point, color="#d62728", s=25, depthshade=False)
                    label = ax.text(
                        *point,
                        str(text).replace("\\Gamma", "Γ"),
                        color="black",
                        fontsize=self.font_size,
                        fontweight="bold",
                        va="bottom",
                    )
                    label.set_gid(f"brillouin-label-{len(labels_seen) - 1}")
        _equal_limits(ax, points)
