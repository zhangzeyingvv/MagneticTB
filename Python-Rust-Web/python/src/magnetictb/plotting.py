"""Stable 2.0.10 band-path and interactive-band data backed by Rust."""

from __future__ import annotations

import csv
import io
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Any, Optional

from .api import ExactError, ModelError, PropertiesError, geometry, properties
from .io import _encode_symbolic_matrix
from .model import PreparedModel, current_model
from .tight_binding import Hamiltonian, HamiltonianExpression
from .properties import _complex_matrix


BandPoint = tuple[float, float, float]
BandSegment = tuple[tuple[BandPoint, BandPoint], tuple[str, str]]


def _finite_real(value: Any, field: str) -> float:
    if isinstance(value, bool):
        raise TypeError(f"{field} must be a finite real number")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise TypeError(f"{field} must be a finite real number") from error
    if not math.isfinite(result):
        raise ValueError(f"{field} must be finite")
    return result


def _path(value: Any) -> tuple[BandSegment, ...]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise PropertiesError(
            "InvalidBandPath",
            "the band path must contain one or more ((q1,q2),(label1,label2)) segments",
        )
    segments: list[BandSegment] = []
    for segment_index, segment in enumerate(value):
        if (
            isinstance(segment, (str, bytes))
            or not isinstance(segment, Sequence)
            or len(segment) != 2
        ):
            raise PropertiesError(
                "InvalidBandPath", f"path segment {segment_index + 1} is malformed"
            )
        endpoints, labels = segment
        if (
            isinstance(endpoints, (str, bytes))
            or not isinstance(endpoints, Sequence)
            or len(endpoints) != 2
            or any(
                isinstance(point, (str, bytes))
                or not isinstance(point, Sequence)
                or len(point) != 3
                for point in endpoints
            )
        ):
            raise PropertiesError(
                "InvalidBandPath",
                f"path segment {segment_index + 1} must contain two three-vectors",
            )
        if (
            isinstance(labels, (str, bytes))
            or not isinstance(labels, Sequence)
            or len(labels) != 2
            or any(not isinstance(label, str) for label in labels)
        ):
            raise PropertiesError(
                "InvalidBandPath",
                f"path segment {segment_index + 1} must contain two string labels",
            )
        points = tuple(
            tuple(
                _finite_real(coordinate, f"path[{segment_index}][0][{point_index}][{axis}]")
                for axis, coordinate in enumerate(point)
            )
            for point_index, point in enumerate(endpoints)
        )
        segments.append(((points[0], points[1]), (labels[0], labels[1])))
    return tuple(segments)


def _subdivisions(value: Any) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise PropertiesError(
            "InvalidBandSubdivisionCount",
            f"npoint must be a positive integer; received {value!r}",
        )
    return value


def _sample_path(
    path: tuple[BandSegment, ...], npoint: int
) -> tuple[tuple[BandPoint, ...], tuple[float, ...], tuple[BandPoint, ...]]:
    segment_points: list[tuple[BandPoint, ...]] = []
    x_coordinates: list[float] = []
    flat: list[BandPoint] = []
    for segment_index, (endpoints, _labels) in enumerate(path):
        start, end = endpoints
        points = tuple(
            tuple(
                start[axis] + (end[axis] - start[axis]) * index / npoint
                for axis in range(3)
            )
            for index in range(npoint + 1)
        )
        segment_points.append(points)
        flat.extend(points)
        x_coordinates.extend(
            float(npoint * segment_index + index) for index in range(npoint + 1)
        )
    return tuple(segment_points), tuple(x_coordinates), tuple(flat)


def _boundary_labels(path: tuple[BandSegment, ...]) -> tuple[str, ...]:
    labels = [path[0][1][0]]
    for index, segment in enumerate(path):
        label = segment[1][1]
        if index + 1 < len(path) and label != path[index + 1][1][0]:
            label = f"{label}|{path[index + 1][1][0]}"
        labels.append(label)
    return tuple(label.replace("\\Gamma", "Γ") for label in labels)


def _rules(value: Any) -> dict[str, Any]:
    if value is None:
        return {}
    if isinstance(value, Mapping):
        return {str(key): item for key, item in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        try:
            return {str(key): item for key, item in value}
        except (TypeError, ValueError) as error:
            raise TypeError("rules must be a mapping or sequence of (name, value) pairs") from error
    raise TypeError("rules must be a mapping or sequence of (name, value) pairs")


def _symbolic_matrix(value: Any) -> bool:
    if isinstance(value, Hamiltonian):
        return True
    if (
        not isinstance(value, Sequence)
        or isinstance(value, (str, bytes))
        or not value
        or not all(
            isinstance(row, Sequence) and not isinstance(row, (str, bytes))
            for row in value
        )
    ):
        return False
    return any(
        isinstance(item, HamiltonianExpression) for row in value for item in row
    )


def _matrices(
    hamiltonian: Any,
    rules: Mapping[str, Any],
    fractional_points: Sequence[BandPoint],
) -> list[list[list[dict[str, float]]]]:
    momenta = [[math.tau * coordinate for coordinate in point] for point in fractional_points]
    if _symbolic_matrix(hamiltonian):
        matrix, context, encoded_parameters = _encode_symbolic_matrix(
            hamiltonian, rules
        )
        try:
            result = geometry(
                "evaluate_symbolic_expression_matrices",
                matrix=matrix,
                context=context,
                parameters=encoded_parameters,
                momenta=momenta,
            )
        except ExactError as error:
            raise ModelError(error.tag, error.detail) from None
        return result["matrices"]
    if callable(hamiltonian):
        if rules:
            raise TypeError("rules apply only to a symbolic MagneticTB Hamiltonian")
        return [
            _complex_matrix(hamiltonian(list(momentum)), f"Hamiltonian at point {index + 1}")
            for index, momentum in enumerate(momenta)
        ]
    if rules:
        raise TypeError("rules apply only to a symbolic MagneticTB Hamiltonian")
    matrix = _complex_matrix(hamiltonian, "hamiltonian")
    return [matrix for _ in fractional_points]


@dataclass(frozen=True)
class _BandPlotStyle:
    plot_range: Optional[tuple[float, float]] = None
    y_ticks: Any = "Automatic"
    font_size: float = 24.0
    font_family: str = "Times"
    image_size: tuple[float, float] = (640.0, 420.0)


def _band_plot_style(options: Mapping[str, Any]) -> _BandPlotStyle:
    def invalid(detail: str) -> PropertiesError:
        return PropertiesError("InvalidBandPlotOption", detail)

    plot_range = options.get("plotRange", "All")
    if plot_range is None or isinstance(plot_range, str) and plot_range in ("All", "Automatic"):
        plot_range = None
    else:
        if not isinstance(plot_range, Sequence) or isinstance(plot_range, str) or len(plot_range) != 2:
            raise invalid("plotRange must be All, Automatic, None, or (minimum, maximum)")
        plot_range = tuple(_finite_real(value, "plotRange") for value in plot_range)
        if plot_range[0] >= plot_range[1]:
            raise invalid("plotRange minimum must be less than maximum")
    ticks = options.get("yTicks", "Automatic")
    if ticks is None or isinstance(ticks, str) and ticks == "None":
        ticks = ()
    elif isinstance(ticks, str):
        if ticks != "Automatic":
            raise invalid("yTicks must be Automatic, None, positions, or (position, label) pairs")
    elif isinstance(ticks, Sequence):
        parsed = []
        for tick in ticks:
            if isinstance(tick, Sequence) and not isinstance(tick, (str, bytes)):
                if len(tick) != 2:
                    raise invalid("a custom y tick must be a (position, label) pair")
                parsed.append((_finite_real(tick[0], "yTicks position"), str(tick[1])))
            else:
                position = _finite_real(tick, "yTicks position")
                parsed.append((position, f"{position:g}"))
        ticks = tuple(parsed)
    else:
        raise invalid("yTicks must be Automatic, None, positions, or (position, label) pairs")
    font_size = _finite_real(options.get("FontSize", 24), "FontSize")
    font_family = options.get("FontFamily", "Times")
    if font_size <= 0 or not isinstance(font_family, str) or not font_family:
        raise invalid("FontSize must be positive and FontFamily must be a nonempty string")
    size = options.get("ImageSize", "Automatic")
    if size is None or isinstance(size, str) and size == "Automatic":
        size = (640.0, 420.0)
    elif isinstance(size, Sequence) and not isinstance(size, (str, bytes)):
        if len(size) != 2:
            raise invalid("ImageSize must be Automatic, a positive width, or (width, height)")
        size = tuple(_finite_real(value, "ImageSize") for value in size)
    else:
        width = _finite_real(size, "ImageSize")
        size = (width, width * 420 / 640)
    if any(value <= 0 for value in size):
        raise invalid("ImageSize dimensions must be positive")
    return _BandPlotStyle(plot_range, ticks, font_size, font_family, size)


class _BandPlotDisplay:
    """Shared rendering lifecycle for cached bands and band comparisons."""

    _style: _BandPlotStyle

    def _draw(self, ax: Any) -> None:
        raise NotImplementedError

    def _figure(self) -> Any:
        from matplotlib.figure import Figure

        # Notebook display and export must not create a GUI figure manager.
        figure = Figure(
            figsize=tuple(value / 100 for value in self._style.image_size),
            dpi=100,
            layout="constrained",
        )
        self._draw(figure.subplots())
        return figure

    def plot(self, ax: Any = None) -> Any:
        """Draw on a Matplotlib Axes (or create one) and return that Axes."""
        if ax is None:
            from matplotlib import pyplot as plt

            _, ax = plt.subplots(
                figsize=tuple(value / 100 for value in self._style.image_size),
                dpi=100,
                layout="constrained",
            )
        self._draw(ax)
        return ax

    def show(self, *, block: Optional[bool] = None) -> None:
        """Show the figure using the active Matplotlib notebook/GUI backend."""
        from matplotlib import pyplot as plt

        self.plot()
        plt.show(block=block)

    def savefig(self, filename: Any, **options: Any) -> None:
        """Save PNG, SVG, PDF, etc. without opening a window or recomputing bands."""
        self._figure().savefig(filename, **options)

    def _repr_svg_(self) -> str:
        output = io.StringIO()
        self._figure().savefig(output, format="svg", metadata={"Date": None})
        return output.getvalue()


@dataclass(frozen=True)
class BandPlotResult(_BandPlotDisplay):
    """Rust-computed bands with automatic Jupyter graphics and explicit plotting.

    A last-expression ``bandplot(...)`` displays an SVG in Jupyter. In scripts,
    use ``result.show()`` or ``result.savefig('bands.png')``. Rendering never
    reevaluates the Hamiltonian or calls the eigensolver.
    """

    path: tuple[BandSegment, ...]
    npoint: int
    segment_points: tuple[tuple[BandPoint, ...], ...]
    x_coordinates: tuple[float, ...]
    boundary_labels: tuple[str, ...]
    eigenvalues: tuple[tuple[float, ...], ...]
    bands: tuple[tuple[float, ...], ...]
    _style: _BandPlotStyle = field(default_factory=_BandPlotStyle, repr=False)

    @property
    def dimension(self) -> int:
        return len(self.bands)

    @property
    def boundary_positions(self) -> tuple[int, ...]:
        return tuple(self.npoint * index for index in range(len(self.path) + 1))

    def _draw(self, ax: Any) -> None:
        # Stable ListLinePlot draws each band/segment independently. Do not
        # connect X and M across a disconnected X|M path boundary.
        for band in self.bands:
            for segment in range(len(self.path)):
                start = segment * (self.npoint + 1)
                stop = start + self.npoint + 1
                ax.plot(self.x_coordinates[start:stop], band[start:stop], color="black")
        # Gridlines and explicit tick positions must not enlarge the energy
        # range selected from the bands (or from the user's plotRange).
        y_limits = self._style.plot_range or ax.get_ylim()
        for position in self.boundary_positions:
            ax.axvline(position, color="black", linewidth=0.6)
        ax.axhline(0, color="black", linewidth=0.6)
        ax.set_xticks(self.boundary_positions, self.boundary_labels)
        ax.set_xlim(self.boundary_positions[0], self.boundary_positions[-1])
        if self._style.y_ticks != "Automatic":
            ax.set_yticks(
                [position for position, _ in self._style.y_ticks],
                [label for _, label in self._style.y_ticks],
            )
        ax.set_ylim(*y_limits)
        # Map Mathematica's Times/Helvetica families to portable Matplotlib
        # categories; other explicitly supplied font names are preserved.
        family = {"Times": "serif", "Helvetica": "sans-serif"}.get(
            self._style.font_family, self._style.font_family
        )
        ax.tick_params(axis="both", colors="black", labelsize=self._style.font_size)
        for label in (*ax.get_xticklabels(), *ax.get_yticklabels()):
            label.set_fontfamily(family)
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color("black")

    def to_dict(self) -> dict[str, Any]:
        return {
            "path": [
                [[list(start), list(end)], [left, right]]
                for ((start, end), (left, right)) in self.path
            ],
            "npoint": self.npoint,
            "segment_points": [
                [list(point) for point in segment] for segment in self.segment_points
            ],
            "x_coordinates": list(self.x_coordinates),
            "boundary_positions": list(self.boundary_positions),
            "boundary_labels": list(self.boundary_labels),
            "eigenvalues": [list(values) for values in self.eigenvalues],
            "bands": [list(values) for values in self.bands],
        }


def bandplot(
    pathstr: Any,
    npoint: int,
    ham: Any,
    rules: Any = (),
    *,
    HermitianTolerance: float = 1.0e-10,
    **plot_options: Any,
) -> BandPlotResult:
    """Sample the stable 2.0.10 path convention and diagonalize in Rust.

    The result displays automatically in Jupyter and supports show/savefig/plot.
    Plot-only options affect rendering, never the Rust-computed energy data.
    ``"Automatic"`` resolves the current model's standard path, as in stable MMA.
    """

    supported_plot_options = {
        "plotRange",
        "yTicks",
        "FontSize",
        "FontFamily",
        "ImageSize",
    }
    unknown = sorted(set(plot_options) - supported_plot_options)
    if unknown:
        raise PropertiesError(
            "UnknownBandPlotOption", f"unknown bandplot option(s): {', '.join(unknown)}"
        )
    style = _band_plot_style(plot_options)
    path = standardKPath() if isinstance(pathstr, str) and pathstr == "Automatic" else _path(pathstr)
    subdivisions = _subdivisions(npoint)
    segment_points, x_coordinates, points = _sample_path(path, subdivisions)
    tolerance = _finite_real(HermitianTolerance, "HermitianTolerance")
    if tolerance < 0:
        raise PropertiesError(
            "InvalidPropertiesOption", "HermitianTolerance must be nonnegative"
        )
    result = properties(
        "band_eigenvalues",
        matrices=_matrices(ham, _rules(rules), points),
        hermitian_tolerance=tolerance,
    )
    eigenvalues = tuple(
        tuple(float(value) for value in values) for values in result["Eigenvalues"]
    )
    dimension = len(eigenvalues[0])
    if any(len(values) != dimension for values in eigenvalues):
        raise PropertiesError(
            "InconsistentBandDimension", "Hamiltonian dimension changed along the path"
        )
    bands = tuple(
        tuple(eigenvalues[point][band] for point in range(len(eigenvalues)))
        for band in range(dimension)
    )
    return BandPlotResult(
        path=path,
        npoint=subdivisions,
        segment_points=segment_points,
        x_coordinates=x_coordinates,
        boundary_labels=_boundary_labels(path),
        eigenvalues=eigenvalues,
        bands=bands,
        _style=style,
    )


def standardKPath(
    *,
    BravaisType: str = "Automatic",
    Tolerance: float = 1.0e-6,
    model: Optional[PreparedModel] = None,
) -> tuple[BandSegment, ...]:
    """Return the stable 2.0.10 conventional path for the current model."""

    active = current_model() if model is None else model
    if not isinstance(active, PreparedModel):
        raise TypeError("model must be a PreparedModel")
    raw = active.to_canonical_dict(1)
    requested = BravaisType
    if not isinstance(requested, str):
        raise TypeError("BravaisType must be Automatic or a string")
    try:
        result = geometry(
            "standard_k_path",
            context=raw["field_context"],
            lattice=raw["model_lattice"],
            bravais_type=requested,
            tolerance=_finite_real(Tolerance, "Tolerance"),
        )
    except ExactError as error:
        raise ModelError(error.tag, error.detail) from None
    return _path(result["Path"])


def showband(
    npoint: int,
    ham: Any,
    rules: Any = (),
    **options: Any,
) -> BandPlotResult:
    """Stable ``showband`` using ``standardKPath`` from the current session."""

    path = standardKPath(
        BravaisType=options.pop("BravaisType", "Automatic"),
        Tolerance=options.pop("Tolerance", 1.0e-6),
    )
    return bandplot(path, npoint, ham, rules, **options)


def banddata(
    pathstr: Any,
    npoint: int,
    ham: Any,
    rules: Any,
    save: str | Path,
    *,
    HermitianTolerance: float = 1.0e-10,
) -> str:
    """Write the stable point-major band-energy CSV and return its path."""

    result = bandplot(
        pathstr,
        npoint,
        ham,
        rules,
        HermitianTolerance=HermitianTolerance,
    )
    target = Path(save)
    with target.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerows(result.eigenvalues)
    return str(target)


@dataclass(frozen=True)
class BandExplorer:
    """Python/Web counterpart of stable ``bandManipulate`` parameter controls."""

    path: tuple[BandSegment, ...]
    npoint: int
    hamiltonian: Any
    parameter_names: tuple[str, ...]
    parameter_range: tuple[float, float] = (-1.0, 1.0)
    hermitian_tolerance: float = 1.0e-10

    def defaults(self) -> Mapping[str, float]:
        return MappingProxyType({name: 0.0 for name in self.parameter_names})

    def evaluate(self, rules: Optional[Mapping[str, Any]] = None) -> BandPlotResult:
        assigned = dict(self.defaults())
        if rules is not None:
            unknown = sorted(set(rules) - set(self.parameter_names))
            if unknown:
                raise PropertiesError(
                    "UnknownBandParameter",
                    f"unknown band parameter(s): {', '.join(unknown)}",
                )
            assigned.update(rules)
        return bandplot(
            self.path,
            self.npoint,
            self.hamiltonian,
            assigned,
            HermitianTolerance=self.hermitian_tolerance,
        )


def _parameter_names(hamiltonian: Any) -> tuple[str, ...]:
    if isinstance(hamiltonian, Hamiltonian):
        return hamiltonian.parameter_names
    if (
        isinstance(hamiltonian, (str, bytes))
        or not isinstance(hamiltonian, Sequence)
        or not hamiltonian
        or any(
            isinstance(row, (str, bytes))
            or not isinstance(row, Sequence)
            or len(row) != len(hamiltonian)
            for row in hamiltonian
        )
    ):
        raise PropertiesError(
            "InvalidHamiltonian", "Hamiltonian must be a nonempty square matrix"
        )
    if not _symbolic_matrix(hamiltonian):
        _complex_matrix(hamiltonian, "hamiltonian")
        return ()
    names: list[str] = []
    for row in hamiltonian:
        for cell in row:
            if isinstance(cell, HamiltonianExpression):
                for contribution in cell.contributions:
                    if contribution.parameter_name not in names:
                        names.append(contribution.parameter_name)
    return tuple(names)


def bandManipulate(
    pathstr: Any,
    npoint: int,
    h: Any,
    *,
    HermitianTolerance: float = 1.0e-10,
) -> BandExplorer:
    """Return immutable controls for the stable ordinary-symbol band workflow."""

    path = _path(pathstr)
    subdivisions = _subdivisions(npoint)
    tolerance = _finite_real(HermitianTolerance, "HermitianTolerance")
    if tolerance < 0:
        raise PropertiesError(
            "InvalidPropertiesOption", "HermitianTolerance must be nonnegative"
        )
    return BandExplorer(
        path=path,
        npoint=subdivisions,
        hamiltonian=h,
        parameter_names=_parameter_names(h),
        hermitian_tolerance=tolerance,
    )


def standard_k_path(**kwargs: Any) -> tuple[BandSegment, ...]:
    return standardKPath(**kwargs)


def band_manipulate(*args: Any, **kwargs: Any) -> BandExplorer:
    return bandManipulate(*args, **kwargs)
