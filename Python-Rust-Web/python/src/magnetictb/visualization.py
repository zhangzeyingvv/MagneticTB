"""Stable 2.0.10 numerical plot-data APIs backed by the Rust Properties core."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
from typing import Any, Optional

from .api import DataCatalog, ExactError, ModelError, PropertiesError, geometry
from .linear_algebra import _decode_element, _thaw_json
from .model import PreparedModel, current_model
from .plotting import _path, standardKPath
from .properties import berryCurvature, surfaceSpectralFunction, wilsonLoop
from ._spatial_plotting import _BrillouinPlot, _CrystalPlot, _image_size


def _real(value: Any, name: str) -> float:
    if isinstance(value, (bool, complex)):
        raise TypeError(f"{name} must be a finite real number")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise TypeError(f"{name} must be a finite real number") from error
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _positive_integer(value: Any, name: str, minimum: int = 1) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise PropertiesError(
            "InvalidPropertiesOption", f"{name} must be an integer >= {minimum}"
        )
    return value


def _point(value: Any, dimension: int, name: str) -> tuple[float, ...]:
    if (
        isinstance(value, (str, bytes))
        or not isinstance(value, Sequence)
        or len(value) != dimension
    ):
        raise PropertiesError(
            "InvalidPlotPath", f"{name} must be a real {dimension}-vector"
        )
    return tuple(_real(item, f"{name}[{index}]") for index, item in enumerate(value))


def _linspace(start: float, stop: float, count: int) -> tuple[float, ...]:
    return tuple(start + (stop - start) * index / (count - 1) for index in range(count))


@dataclass(frozen=True)
class _SampledPath:
    points: tuple[tuple[float, ...], ...]
    coordinates: tuple[float, ...]
    boundary_coordinates: Optional[tuple[float, ...]]
    boundary_labels: Optional[tuple[str, ...]]


@dataclass(frozen=True)
class CrystalAtomRecord:
    orbit_index: int
    equivalent_index: int
    fractional_position: tuple[Any, Any, Any]
    fractional_moment: tuple[Any, Any, Any]
    cell_translation: tuple[int, int, int]
    cartesian_position: tuple[float, float, float]
    cartesian_moment: tuple[float, float, float]
    magnetic: bool
    arrow_end: Optional[tuple[float, float, float]]


@dataclass(frozen=True)
class CrystalStructureResult(_CrystalPlot):
    """Rust crystal data with Jupyter SVG display and plot/show/savefig methods."""
    lattice: tuple[tuple[float, float, float], ...]
    cell_range: tuple[tuple[int, int], tuple[int, int], tuple[int, int]]
    translations: tuple[tuple[int, int, int], ...]
    atom_radius: float
    moment_scale: float
    atom_records: tuple[CrystalAtomRecord, ...]
    atom_count: int
    magnetic_atom_count: int
    cell_edges: tuple[
        tuple[tuple[float, float, float], tuple[float, float, float]], ...
    ]
    show_atom_labels: bool
    font_size: float
    image_size: Any

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class BrillouinZoneResult(_BrillouinPlot):
    """Rust first-BZ data with Jupyter SVG display and plot/show/savefig methods."""
    reciprocal_lattice: tuple[tuple[float, float, float], ...]
    vertices: tuple[tuple[float, float, float], ...]
    facets: tuple[tuple[int, ...], ...]
    facet_count: int
    bravais_type: Optional[str]
    bz_type: Optional[str]
    k_path: tuple[Any, ...]
    displayed_k_path: tuple[Any, ...]
    cartesian_k_path: tuple[Any, ...]
    translation_range: int
    tolerance: float
    show_k_path: bool
    font_size: float
    image_size: Any

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _active_model(model: Optional[PreparedModel]) -> PreparedModel:
    active = current_model() if model is None else model
    if not isinstance(active, PreparedModel):
        raise TypeError("model must be a PreparedModel")
    return active


def _cell_range(value: Any) -> tuple[tuple[int, int], tuple[int, int], tuple[int, int]]:
    if isinstance(value, int) and not isinstance(value, bool) and value >= 0:
        return ((-value, value), (-value, value), (-value, value))
    if (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and len(value) == 3
    ):
        result = []
        for axis, bounds in enumerate(value):
            if (
                not isinstance(bounds, Sequence)
                or isinstance(bounds, (str, bytes))
                or len(bounds) != 2
                or any(isinstance(item, bool) or not isinstance(item, int) for item in bounds)
                or bounds[0] > bounds[1]
            ):
                raise PropertiesError(
                    "InvalidCrystalCellRange",
                    f"CellRange axis {axis + 1} must contain ordered integer bounds",
                )
            result.append((bounds[0], bounds[1]))
        return (result[0], result[1], result[2])
    raise PropertiesError(
        "InvalidCrystalCellRange",
        "CellRange must be a nonnegative integer or three ordered integer pairs",
    )


def showCrystalStructure(
    *,
    CellRange: Any = ((0, 0), (0, 0), (0, 0)),
    MomentRules: Any = (),
    MomentScale: Any = "Automatic",
    AtomRadius: Any = "Automatic",
    ShowAtomLabels: bool = False,
    FontSize: Any = 14,
    ImageSize: Any = "Large",
    model: Optional[PreparedModel] = None,
) -> CrystalStructureResult:
    """Show Rust-generated crystal geometry as an orthographic 3D plot.

    The returned object displays directly as a Jupyter last expression, retains
    its numerical fields/to_dict(), and supports plot(), show(), and savefig().
    """

    active = _active_model(model)
    recipe = active._compilation_recipe
    if recipe is None:
        raise ModelError(
            "ModelRecompilationUnavailable",
            "showCrystalStructure requires the original stable init recipe",
        )
    if MomentRules not in ((), [], {}, None):
        raise ModelError(
            "UnsupportedCrystalMomentRules",
            "the current exact Python input domain cannot retain unresolved moment symbols",
        )
    if not isinstance(ShowAtomLabels, bool):
        raise TypeError("ShowAtomLabels must be bool")
    font_size = _real(FontSize, "FontSize")
    if font_size <= 0.0:
        raise PropertiesError("InvalidCrystalStructureOption", "FontSize must be positive")
    _image_size(ImageSize)
    cell_range = _cell_range(CellRange)
    moment_scale = None if MomentScale == "Automatic" else _real(MomentScale, "MomentScale")
    atom_radius = None if AtomRadius == "Automatic" else _real(AtomRadius, "AtomRadius")
    model_result = active.to_canonical_dict(1)
    dispatcher = DataCatalog().geometry if model_result.get("data_loaded") else geometry
    try:
        raw = dispatcher(
            "crystal_structure_data",
            model_result=model_result,
            compiler_input=_thaw_json(recipe.compiler_input),
            cell_range=cell_range,
            moment_scale=moment_scale,
            atom_radius=atom_radius,
        )
    except ExactError as error:
        raise ModelError(error.tag, error.detail) from None
    records = tuple(
        CrystalAtomRecord(
            orbit_index=int(record["OrbitIndex"]),
            equivalent_index=int(record["EquivalentIndex"]),
            fractional_position=tuple(
                _decode_element(value, active._field)
                for value in record["FractionalPosition"]
            ),
            fractional_moment=tuple(
                _decode_element(value, active._field)
                for value in record["FractionalMoment"]
            ),
            cell_translation=tuple(int(value) for value in record["CellTranslation"]),
            cartesian_position=tuple(float(value) for value in record["CartesianPosition"]),
            cartesian_moment=tuple(float(value) for value in record["CartesianMoment"]),
            magnetic=bool(record["Magnetic"]),
            arrow_end=(
                None
                if record["ArrowEnd"] is None
                else tuple(float(value) for value in record["ArrowEnd"])
            ),
        )
        for record in raw["AtomRecords"]
    )
    return CrystalStructureResult(
        lattice=tuple(tuple(float(value) for value in row) for row in raw["Lattice"]),
        cell_range=cell_range,
        translations=tuple(tuple(int(value) for value in row) for row in raw["Translations"]),
        atom_radius=float(raw["AtomRadius"]),
        moment_scale=float(raw["MomentScale"]),
        atom_records=records,
        atom_count=int(raw["AtomCount"]),
        magnetic_atom_count=int(raw["MagneticAtomCount"]),
        cell_edges=tuple(
            tuple(tuple(float(value) for value in point) for point in edge)
            for edge in raw["CellEdges"]
        ),
        show_atom_labels=ShowAtomLabels,
        font_size=font_size,
        image_size=ImageSize,
    )


def showBrillouinZone(
    *,
    KPath: Any = "Automatic",
    ShowKPath: bool = True,
    TranslationRange: int = 2,
    Tolerance: Any = 1.0e-8,
    FontSize: Any = 14,
    ImageSize: Any = "Large",
    model: Optional[PreparedModel] = None,
) -> BrillouinZoneResult:
    """Show Rust-generated first-BZ facets and folded paths in parallel projection.

    The result displays directly in Jupyter; plot()/show() use Matplotlib's
    active backend and savefig() exports a view. Numerical/Web fields remain
    available through to_dict().
    """

    active = _active_model(model)
    if not isinstance(ShowKPath, bool):
        raise TypeError("ShowKPath must be bool")
    translation_range = _positive_integer(TranslationRange, "TranslationRange")
    tolerance = _real(Tolerance, "Tolerance")
    font_size = _real(FontSize, "FontSize")
    if tolerance <= 0.0 or font_size <= 0.0:
        raise PropertiesError(
            "InvalidBrillouinZoneOption", "Tolerance and FontSize must be positive"
        )
    _image_size(ImageSize)
    if KPath == "Automatic":
        path = standardKPath(Tolerance=max(tolerance, 1.0e-6), model=active)
        bravais_type: Optional[str] = None
        bz_type: Optional[str] = None
        path_raw = active.to_canonical_dict(1)
        standard = geometry(
            "standard_k_path",
            context=path_raw["field_context"],
            lattice=path_raw["model_lattice"],
            bravais_type="Automatic",
            tolerance=max(tolerance, 1.0e-6),
        )
        bravais_type = str(standard["BravaisType"])
        bz_type = str(standard["BZType"])
    elif KPath is None:
        path = ()
        bravais_type = None
        bz_type = None
    else:
        path = _path(KPath)
        bravais_type = "UserSupplied"
        bz_type = None
    canonical = active.to_canonical_dict(1)
    try:
        raw = geometry(
            "brillouin_zone_data",
            context=canonical["field_context"],
            lattice=canonical["model_lattice"],
            path=path,
            bravais_type=bravais_type,
            bz_type=bz_type,
            translation_range=translation_range,
            tolerance=tolerance,
        )
    except ExactError as error:
        raise PropertiesError(error.tag, error.detail) from None
    return BrillouinZoneResult(
        reciprocal_lattice=tuple(
            tuple(float(value) for value in row) for row in raw["ReciprocalLattice"]
        ),
        vertices=tuple(tuple(float(value) for value in row) for row in raw["Vertices"]),
        facets=tuple(tuple(int(value) for value in facet) for facet in raw["Facets"]),
        facet_count=int(raw["FacetCount"]),
        bravais_type=None if raw["BravaisType"] is None else str(raw["BravaisType"]),
        bz_type=None if raw["BZType"] is None else str(raw["BZType"]),
        k_path=tuple(raw["KPath"]),
        displayed_k_path=tuple(raw["DisplayedKPath"]),
        cartesian_k_path=tuple(raw["CartesianKPath"]),
        translation_range=translation_range,
        tolerance=tolerance,
        show_k_path=ShowKPath,
        font_size=font_size,
        image_size=ImageSize,
    )


def _sample_path(path: Any, subdivisions: int, dimension: int) -> _SampledPath:
    if isinstance(path, (str, bytes)) or not isinstance(path, Sequence) or not path:
        raise PropertiesError("InvalidPlotPath", "path must be nonempty")
    try:
        points = tuple(_point(item, dimension, "path point") for item in path)
    except (PropertiesError, TypeError, ValueError):
        points = ()
    if len(points) >= 2:
        distances = tuple(math.dist(left, right) for left, right in zip(points, points[1:]))
        if any(distance <= 0.0 for distance in distances):
            raise PropertiesError("InvalidPlotPath", "adjacent path points must differ")
        coordinate = 0.0
        coordinates = [coordinate]
        for distance in distances:
            coordinate += distance
            coordinates.append(coordinate)
        return _SampledPath(points, tuple(coordinates), None, None)

    segments: list[tuple[tuple[float, ...], tuple[float, ...], str, str]] = []
    for index, segment in enumerate(path):
        if (
            isinstance(segment, (str, bytes))
            or not isinstance(segment, Sequence)
            or len(segment) != 2
        ):
            raise PropertiesError("InvalidPlotPath", f"segment {index + 1} is malformed")
        endpoints, labels = segment
        if (
            isinstance(endpoints, (str, bytes))
            or not isinstance(endpoints, Sequence)
            or len(endpoints) != 2
            or isinstance(labels, (str, bytes))
            or not isinstance(labels, Sequence)
            or len(labels) != 2
            or not all(isinstance(label, str) for label in labels)
        ):
            raise PropertiesError("InvalidPlotPath", f"segment {index + 1} is malformed")
        start = _point(endpoints[0], dimension, f"segment {index + 1} start")
        end = _point(endpoints[1], dimension, f"segment {index + 1} end")
        if math.dist(start, end) <= 0.0:
            raise PropertiesError("InvalidPlotPath", "path segments must have positive length")
        if segments and math.dist(segments[-1][1], start) > 1.0e-10:
            raise PropertiesError("InvalidPlotPath", "labeled path segments must be connected")
        segments.append((start, end, labels[0], labels[1]))
    sampled: list[tuple[float, ...]] = []
    coordinates: list[float] = []
    labels = [segments[0][2]]
    for segment_index, (start, end, _left, right) in enumerate(segments):
        current = [
            tuple(
                start[axis] + (end[axis] - start[axis]) * index / subdivisions
                for axis in range(dimension)
            )
            for index in range(subdivisions + 1)
        ]
        x_values = [
            segment_index + index / subdivisions for index in range(subdivisions + 1)
        ]
        if segment_index:
            current = current[1:]
            x_values = x_values[1:]
        sampled.extend(current)
        coordinates.extend(x_values)
        next_left = segments[segment_index + 1][2] if segment_index + 1 < len(segments) else None
        labels.append(f"{right}|{next_left}" if next_left is not None and right != next_left else right)
    return _SampledPath(
        tuple(sampled),
        tuple(coordinates),
        tuple(float(index) for index in range(len(segments) + 1)),
        tuple(label.replace("\\Gamma", "Γ") for label in labels),
    )


@dataclass(frozen=True)
class WilsonLoopPlotResult:
    parameter_path: tuple[tuple[float, ...], ...]
    path_coordinate: tuple[float, ...]
    boundary_coordinates: Optional[tuple[float, ...]]
    boundary_labels: Optional[tuple[str, ...]]
    phase_convention: str
    values: tuple[tuple[float, ...], ...]

    @property
    def branch_count(self) -> int:
        return len(self.values[0]) if self.values else 0

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def plotWilsonLoop(
    hamiltonian: Any,
    centers: Sequence[Sequence[Any]],
    occupied: int,
    start: Sequence[Any],
    end: Sequence[Any],
    parameter_path: Any,
    *,
    parameters: Optional[Mapping[str, Any]] = None,
    ParameterSubdivisions: int = 60,
    LoopSubdivisions: int = 50,
    HermitianTolerance: float = 1.0e-10,
    GapTolerance: float = 1.0e-9,
    CovarianceTolerance: float = 1.0e-8,
    OverlapTolerance: float = 1.0e-10,
    PhaseConvention: str = "PhaseOverPi",
) -> WilsonLoopPlotResult:
    subdivisions = _positive_integer(ParameterSubdivisions, "ParameterSubdivisions")
    _positive_integer(LoopSubdivisions, "LoopSubdivisions")
    start_point = tuple(_real(value, "start") for value in start)
    end_point = tuple(_real(value, "end") for value in end)
    if not start_point or len(start_point) != len(end_point):
        raise PropertiesError("InvalidPlotPath", "start and end dimensions must agree")
    if PhaseConvention not in ("PhaseOverPi", "WannierCenters"):
        raise PropertiesError("InvalidPropertiesOption", "invalid PhaseConvention")
    sampled = _sample_path(parameter_path, subdivisions, len(start_point))
    values = tuple(
        tuple(
            float(value)
            for value in wilsonLoop(
                hamiltonian,
                centers,
                occupied,
                tuple(left + offset for left, offset in zip(start_point, point)),
                tuple(left + offset for left, offset in zip(end_point, point)),
                parameters=parameters,
                PathSubdivisions=LoopSubdivisions,
                HermitianTolerance=HermitianTolerance,
                GapTolerance=GapTolerance,
                CovarianceTolerance=CovarianceTolerance,
                OverlapTolerance=OverlapTolerance,
                Output=PhaseConvention,
            )
        )
        for point in sampled.points
    )
    return WilsonLoopPlotResult(
        sampled.points,
        sampled.coordinates,
        sampled.boundary_coordinates,
        sampled.boundary_labels,
        PhaseConvention,
        values,
    )


@dataclass(frozen=True)
class SurfaceSpectrumPlotResult:
    surface_momentum_path: tuple[tuple[float, float], ...]
    path_coordinate: tuple[float, ...]
    boundary_coordinates: Optional[tuple[float, ...]]
    boundary_labels: Optional[tuple[str, ...]]
    energies: tuple[float, ...]
    spectral_weight: tuple[tuple[float, ...], ...]
    broadening: float
    surface: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def plotSurfaceSpectrum(
    hopping_data: Mapping[str, Any],
    momentum_path: Any,
    energy_range: Sequence[Any],
    *,
    MomentumSubdivisions: int = 60,
    EnergyPoints: int = 201,
    CellMatrix: Any = "Automatic",
    Broadening: float = 1.0e-3,
    Tolerance: float = 1.0e-10,
    MaxIterations: int = 200,
    Surface: str = "Positive",
    HermiticityTolerance: float = 1.0e-9,
) -> SurfaceSpectrumPlotResult:
    subdivisions = _positive_integer(MomentumSubdivisions, "MomentumSubdivisions")
    energy_points = _positive_integer(EnergyPoints, "EnergyPoints", 2)
    if not isinstance(energy_range, Sequence) or len(energy_range) != 2:
        raise PropertiesError("InvalidEnergyRange", "energy range must have two endpoints")
    minimum, maximum = (_real(value, "energy range") for value in energy_range)
    if minimum >= maximum:
        raise PropertiesError("InvalidEnergyRange", "energy range must be increasing")
    sampled = _sample_path(momentum_path, subdivisions, 2)
    energies = _linspace(minimum, maximum, energy_points)
    weights = tuple(
        tuple(
            float(
                surfaceSpectralFunction(
                    hopping_data,
                    point,
                    energy,
                    CellMatrix=CellMatrix,
                    Broadening=Broadening,
                    Tolerance=Tolerance,
                    MaxIterations=MaxIterations,
                    Surface=Surface,
                    HermiticityTolerance=HermiticityTolerance,
                )
            )
            for energy in energies
        )
        for point in sampled.points
    )
    return SurfaceSpectrumPlotResult(
        tuple((point[0], point[1]) for point in sampled.points),
        sampled.coordinates,
        sampled.boundary_coordinates,
        sampled.boundary_labels,
        energies,
        weights,
        _real(Broadening, "Broadening"),
        Surface,
    )


def _ranges(value: Any, dimension: int) -> tuple[tuple[float, float], ...]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) != dimension:
        raise PropertiesError("InvalidPlotRange", f"expected {dimension} ordered ranges")
    ranges = tuple(tuple(_real(item, "plot range") for item in interval) for interval in value)
    if any(len(interval) != 2 or interval[0] >= interval[1] for interval in ranges):
        raise PropertiesError("InvalidPlotRange", "plot ranges must be finite and increasing")
    return ranges


def _grid_size(value: int | Sequence[int], dimension: int) -> tuple[int, ...]:
    if isinstance(value, int) and not isinstance(value, bool):
        return (_positive_integer(value, "GridSize", 2),) * dimension
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)) and len(value) == dimension:
        return tuple(_positive_integer(item, "GridSize", 2) for item in value)
    raise PropertiesError("InvalidPropertiesOption", "GridSize has the wrong shape")


@dataclass(frozen=True)
class BerryCurvature2DResult:
    directions: tuple[int, int]
    fixed_coordinates: tuple[float, ...]
    grid_values: tuple[tuple[float, ...], tuple[float, ...]]
    grid_size: tuple[int, int]
    step_size: tuple[float, float]
    curvature: tuple[tuple[float, ...], ...]
    maximum_absolute_curvature: float

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def plotBerryCurvature2D(
    hamiltonian: Any,
    centers: Sequence[Sequence[Any]],
    occupied: int,
    ranges: Any,
    *,
    parameters: Optional[Mapping[str, Any]] = None,
    Directions: Sequence[int] = (1, 2),
    FixedCoordinates: Any = "Automatic",
    GridSize: int | Sequence[int] = 31,
    StepSize: Any = "Automatic",
    HermitianTolerance: float = 1.0e-10,
    GapTolerance: float = 1.0e-9,
    CovarianceTolerance: float = 1.0e-8,
    OverlapTolerance: float = 1.0e-10,
) -> BerryCurvature2DResult:
    if not centers:
        raise PropertiesError("InvalidWannierCenters", "centers must be nonempty")
    dimension = len(centers[0])
    if dimension < 2 or any(len(center) != dimension for center in centers):
        raise PropertiesError("InvalidWannierCenters", "center dimensions must agree")
    plot_ranges = _ranges(ranges, 2)
    if (
        not isinstance(Directions, Sequence)
        or len(Directions) != 2
        or any(isinstance(index, bool) or not isinstance(index, int) for index in Directions)
        or len(set(Directions)) != 2
        or any(index < 1 or index > dimension for index in Directions)
    ):
        raise PropertiesError("InvalidBerryCurvatureDirections", "Directions are invalid")
    directions = (Directions[0], Directions[1])
    if FixedCoordinates == "Automatic":
        fixed = (0.0,) * dimension
    else:
        fixed = _point(FixedCoordinates, dimension, "FixedCoordinates")
    size = _grid_size(GridSize, 2)
    grids = tuple(_linspace(interval[0], interval[1], count) for interval, count in zip(plot_ranges, size))
    spacings = tuple((interval[1] - interval[0]) / (count - 1) for interval, count in zip(plot_ranges, size))
    if StepSize == "Automatic":
        steps = (spacings[0] / 2.0, spacings[1] / 2.0)
    elif isinstance(StepSize, Sequence) and not isinstance(StepSize, (str, bytes)):
        raw_steps = tuple(_real(value, "StepSize") for value in StepSize)
        if len(raw_steps) == dimension:
            steps = (raw_steps[directions[0] - 1], raw_steps[directions[1] - 1])
        elif len(raw_steps) == 2:
            steps = (raw_steps[0], raw_steps[1])
        else:
            raise PropertiesError("InvalidPropertiesOption", "StepSize has the wrong length")
    else:
        step = _real(StepSize, "StepSize")
        steps = (step, step)
    if any(step <= 0.0 for step in steps):
        raise PropertiesError("InvalidPropertiesOption", "StepSize must be positive")
    curvature_rows = []
    for first in grids[0]:
        row = []
        for second in grids[1]:
            point = list(fixed)
            point[directions[0] - 1] = first
            point[directions[1] - 1] = second
            row.append(
                float(
                    berryCurvature(
                        hamiltonian,
                        centers,
                        occupied,
                        point,
                        parameters=parameters,
                        Directions=directions,
                        StepSize=steps,
                        HermitianTolerance=HermitianTolerance,
                        GapTolerance=GapTolerance,
                        CovarianceTolerance=CovarianceTolerance,
                        OverlapTolerance=OverlapTolerance,
                    )
                )
            )
        curvature_rows.append(tuple(row))
    curvature = tuple(curvature_rows)
    maximum = max(abs(value) for row in curvature for value in row)
    return BerryCurvature2DResult(
        directions,
        fixed,
        (grids[0], grids[1]),
        (size[0], size[1]),
        steps,
        curvature,
        maximum,
    )


@dataclass(frozen=True)
class BerryCurvature3DResult:
    component_directions: tuple[tuple[int, int], ...]
    grid_values: tuple[tuple[float, ...], ...]
    grid_size: tuple[int, int, int]
    step_size: tuple[float, float, float]
    vectors: tuple[Any, ...]
    magnitudes: tuple[Any, ...]
    maximum_magnitude: float
    visible_vector_count: int
    vector_scale: float

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def plotBerryCurvature3D(
    hamiltonian: Any,
    centers: Sequence[Sequence[Any]],
    occupied: int,
    ranges: Any,
    *,
    parameters: Optional[Mapping[str, Any]] = None,
    GridSize: int | Sequence[int] = 7,
    StepSize: Any = "Automatic",
    VectorScale: Any = "Automatic",
    MagnitudeThreshold: float = 0.0,
    HermitianTolerance: float = 1.0e-10,
    GapTolerance: float = 1.0e-9,
    CovarianceTolerance: float = 1.0e-8,
    OverlapTolerance: float = 1.0e-10,
) -> BerryCurvature3DResult:
    if not centers or any(len(center) != 3 for center in centers):
        raise PropertiesError("InvalidWannierCenters", "centers must be real three-vectors")
    plot_ranges = _ranges(ranges, 3)
    size = _grid_size(GridSize, 3)
    grids = tuple(_linspace(interval[0], interval[1], count) for interval, count in zip(plot_ranges, size))
    spacings = tuple((interval[1] - interval[0]) / (count - 1) for interval, count in zip(plot_ranges, size))
    if StepSize == "Automatic":
        steps = tuple(spacing / 2.0 for spacing in spacings)
    elif isinstance(StepSize, Sequence) and not isinstance(StepSize, (str, bytes)):
        steps = tuple(_real(value, "StepSize") for value in StepSize)
        if len(steps) != 3:
            raise PropertiesError("InvalidPropertiesOption", "StepSize must have length three")
    else:
        step = _real(StepSize, "StepSize")
        steps = (step, step, step)
    if any(step <= 0.0 for step in steps):
        raise PropertiesError("InvalidPropertiesOption", "StepSize must be positive")
    component_directions = ((2, 3), (3, 1), (1, 2))
    vectors_nested = []
    magnitudes_nested = []
    for first in grids[0]:
        vector_plane = []
        magnitude_plane = []
        for second in grids[1]:
            vector_line = []
            magnitude_line = []
            for third in grids[2]:
                point = (first, second, third)
                vector = tuple(
                    float(
                        berryCurvature(
                            hamiltonian,
                            centers,
                            occupied,
                            point,
                            parameters=parameters,
                            Directions=directions,
                            StepSize=(steps[directions[0] - 1], steps[directions[1] - 1]),
                            HermitianTolerance=HermitianTolerance,
                            GapTolerance=GapTolerance,
                            CovarianceTolerance=CovarianceTolerance,
                            OverlapTolerance=OverlapTolerance,
                        )
                    )
                    for directions in component_directions
                )
                vector_line.append(vector)
                magnitude_line.append(math.sqrt(sum(value * value for value in vector)))
            vector_plane.append(tuple(vector_line))
            magnitude_plane.append(tuple(magnitude_line))
        vectors_nested.append(tuple(vector_plane))
        magnitudes_nested.append(tuple(magnitude_plane))
    vectors = tuple(vectors_nested)
    magnitudes = tuple(magnitudes_nested)
    flat_magnitudes = [value for plane in magnitudes for line in plane for value in line]
    maximum = max(flat_magnitudes)
    threshold = _real(MagnitudeThreshold, "MagnitudeThreshold")
    if threshold < 0.0:
        raise PropertiesError("InvalidPropertiesOption", "MagnitudeThreshold must be nonnegative")
    scale = 0.45 * min(spacings) if VectorScale == "Automatic" else _real(VectorScale, "VectorScale")
    if scale <= 0.0:
        raise PropertiesError("InvalidPropertiesOption", "VectorScale must be positive")
    return BerryCurvature3DResult(
        component_directions,
        grids,
        (size[0], size[1], size[2]),
        (steps[0], steps[1], steps[2]),
        vectors,
        magnitudes,
        maximum,
        sum(value > threshold for value in flat_magnitudes),
        scale,
    )


@dataclass(frozen=True)
class WavefunctionSiteRecord:
    site_index: int
    basis_indices: tuple[int, ...]
    cell_indices: tuple[int, ...]
    cells: tuple[tuple[int, int, int], ...]
    orbital_indices: tuple[int, ...]
    cartesian_position: tuple[float, float, float]
    weight: float
    representative_amplitude: complex
    phase: float


@dataclass(frozen=True)
class RealSpaceWavefunctionResult:
    aggregation: str
    normalized: bool
    input_norm: float
    amplitudes: tuple[complex, ...]
    site_records: tuple[WavefunctionSiteRecord, ...]
    visible_site_count: int
    maximum_weight: float
    phase_coloring: bool

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def plotRealSpaceWavefunction(
    real_space_data: Mapping[str, Any],
    state: Sequence[Any],
    *,
    Aggregation: str = "Atom",
    Normalize: bool = True,
    WeightThreshold: float = 1.0e-8,
    PositionTolerance: float = 1.0e-8,
    PhaseColoring: bool = False,
) -> RealSpaceWavefunctionResult:
    if (
        not isinstance(real_space_data, Mapping)
        or real_space_data.get("Schema") != "MagneticTBFiniteRealSpaceHamiltonian"
        or not real_space_data.get("PositionDataAvailable", False)
    ):
        raise PropertiesError(
            "MissingRealSpacePositionData",
            "buildRealSpaceHamiltonian(..., Output='Data') position metadata is required",
        )
    dimension = int(real_space_data.get("Dimension", 0))
    records = real_space_data.get("BasisRecords")
    if not isinstance(records, Sequence) or len(records) != dimension:
        raise PropertiesError("InvalidRealSpaceData", "BasisRecords are inconsistent")
    if isinstance(state, (str, bytes)) or not isinstance(state, Sequence) or len(state) != dimension:
        raise PropertiesError("InvalidWavefunction", f"state must have length {dimension}")
    amplitudes = tuple(complex(value) for value in state)
    if any(not math.isfinite(value.real) or not math.isfinite(value.imag) for value in amplitudes):
        raise PropertiesError("InvalidWavefunction", "state amplitudes must be finite")
    norm = math.sqrt(sum(abs(value) ** 2 for value in amplitudes))
    if norm == 0.0:
        raise PropertiesError("ZeroWavefunction", "the wavefunction has zero norm")
    if not isinstance(Normalize, bool) or not isinstance(PhaseColoring, bool):
        raise TypeError("Normalize and PhaseColoring must be bool")
    normalized = tuple(value / norm for value in amplitudes) if Normalize else amplitudes
    if Aggregation not in ("Atom", "Orbital"):
        raise PropertiesError("InvalidPropertiesOption", "Aggregation must be Atom or Orbital")
    tolerance = _real(PositionTolerance, "PositionTolerance")
    threshold = _real(WeightThreshold, "WeightThreshold")
    if tolerance <= 0.0 or threshold < 0.0:
        raise PropertiesError("InvalidPropertiesOption", "invalid wavefunction plot tolerance")
    positions = [tuple(_real(value, "CartesianPosition") for value in record["CartesianPosition"]) for record in records]
    groups: list[list[int]] = []
    for index, position in enumerate(positions):
        if Aggregation == "Orbital":
            groups.append([index])
            continue
        matched = next((group for group in groups if math.dist(positions[group[0]], position) <= tolerance), None)
        if matched is None:
            groups.append([index])
        else:
            matched.append(index)
    site_records = []
    for site_index, group in enumerate(groups, start=1):
        dominant = max(group, key=lambda index: abs(normalized[index]))
        position = tuple(sum(positions[index][axis] for index in group) / len(group) for axis in range(3))
        cells = tuple(dict.fromkeys(tuple(int(value) for value in records[index]["Cell"]) for index in group))
        site_records.append(
            WavefunctionSiteRecord(
                site_index,
                tuple(index + 1 for index in group),
                tuple(dict.fromkeys(int(records[index]["CellIndex"]) for index in group)),
                cells,
                tuple(int(records[index]["OrbitalIndex"]) for index in group),
                (position[0], position[1], position[2]),
                sum(abs(normalized[index]) ** 2 for index in group),
                normalized[dominant],
                math.atan2(normalized[dominant].imag, normalized[dominant].real),
            )
        )
    maximum = max(record.weight for record in site_records)
    return RealSpaceWavefunctionResult(
        Aggregation,
        Normalize,
        norm,
        normalized,
        tuple(site_records),
        sum(record.weight > threshold for record in site_records),
        maximum,
        PhaseColoring,
    )


plot_wilson_loop = plotWilsonLoop
plot_surface_spectrum = plotSurfaceSpectrum
plot_berry_curvature_2d = plotBerryCurvature2D
plot_berry_curvature_3d = plotBerryCurvature3D
plot_real_space_wavefunction = plotRealSpaceWavefunction
show_crystal_structure = showCrystalStructure
show_brillouin_zone = showBrillouinZone
