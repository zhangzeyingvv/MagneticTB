"""Stable 2.0.10 band fitting with all numerical algorithms in Rust."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from numbers import Number
from pathlib import Path
from types import MappingProxyType
from typing import Any, Optional

from .core_bindings import fitting as _rust_fitting
from .errors import FittingError
from .plotting import (
    BandPlotResult,
    _BandPlotDisplay,
    _BandPlotStyle,
    _matrices,
    _parameter_names,
    _path,
    _subdivisions,
    bandplot,
)
from .properties import _finite_real


def _real(value: Any, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Number):
        raise FittingError("InvalidNumericFittingInput", f"{field} must be real numeric")
    number = complex(value)
    if number.imag != 0.0 or not math.isfinite(number.real):
        raise FittingError("InvalidNumericFittingInput", f"{field} must be finite and real")
    return number.real


def _reference_data(value: Any) -> tuple[tuple[float, float, float], tuple[tuple[float, ...], ...]]:
    if hasattr(value, "tolist"):
        value = value.tolist()
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise FittingError(
            "InvalidReferenceBands",
            "reference data must be a nonempty sequence of (k-point, energies) records",
        )
    points: list[tuple[float, float, float]] = []
    bands: list[tuple[float, ...]] = []
    band_count: Optional[int] = None
    for record_index, record in enumerate(value):
        if hasattr(record, "tolist"):
            record = record.tolist()
        if (
            isinstance(record, (str, bytes))
            or not isinstance(record, Sequence)
            or len(record) != 2
        ):
            raise FittingError(
                "InvalidReferenceBands",
                f"reference record {record_index + 1} must contain a k-point and energies",
            )
        point, energies = record
        if hasattr(point, "tolist"):
            point = point.tolist()
        if (
            isinstance(point, (str, bytes))
            or not isinstance(point, Sequence)
            or len(point) != 3
        ):
            raise FittingError(
                "InvalidReferenceBands",
                f"reference k-point {record_index + 1} must have length three",
            )
        if hasattr(energies, "tolist"):
            energies = energies.tolist()
        if (
            isinstance(energies, (str, bytes))
            or not isinstance(energies, Sequence)
            or not energies
        ):
            raise FittingError(
                "InvalidReferenceBands",
                f"reference energies {record_index + 1} must be nonempty",
            )
        numeric_energies = tuple(
            _real(item, f"reference[{record_index + 1}].energies[{band + 1}]")
            for band, item in enumerate(energies)
        )
        if band_count is None:
            band_count = len(numeric_energies)
        elif len(numeric_energies) != band_count:
            raise FittingError(
                "InvalidReferenceBands", "all reference records must have one band count"
            )
        points.append(
            tuple(
                _real(item, f"reference[{record_index + 1}].k[{axis + 1}]")
                for axis, item in enumerate(point)
            )
        )
        bands.append(numeric_energies)
    return tuple(points), tuple(bands)


def _rules(value: Any, field: str) -> dict[str, Any]:
    if isinstance(value, Mapping):
        return {str(name): item for name, item in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        try:
            return {str(name): item for name, item in value}
        except (TypeError, ValueError) as error:
            raise FittingError(
                "InvalidInitialParameters",
                f"{field} must be a mapping or sequence of (name, value) pairs",
            ) from error
    raise FittingError(
        "InvalidInitialParameters",
        f"{field} must be a mapping or sequence of (name, value) pairs",
    )


def _frozen_matrix(value: Any) -> Any:
    if hasattr(value, "to_list"):
        value = value.to_list()
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return tuple(_frozen_matrix(item) for item in value)
    return value


def _band_selection(value: Any, band_count: int) -> tuple[int, ...]:
    if value is None or value == "All":
        return tuple(range(1, band_count + 1))
    if isinstance(value, bool):
        raise FittingError("InvalidBandSelection", "band indices are 1-based integers")
    if isinstance(value, int):
        indices = (value,)
    elif isinstance(value, BandSpan):
        indices = value.indices()
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        indices = tuple(value)
    else:
        raise FittingError("InvalidBandSelection", "invalid BandSelection")
    if (
        not indices
        or any(isinstance(item, bool) or not isinstance(item, int) for item in indices)
        or any(item < 1 or item > band_count for item in indices)
        or any(left >= right for left, right in zip(indices, indices[1:]))
    ):
        raise FittingError(
            "InvalidBandSelection",
            "BandSelection must be strictly increasing, unique, and within the 1-based band range",
        )
    return indices


def _energy_window(value: Any) -> Optional[list[float]]:
    if value is None or value == "All":
        return None
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) != 2:
        raise FittingError(
            "InvalidEnergyWindow", "EnergyWindow must be All/None or {minimum, maximum}"
        )
    result = [_real(item, f"EnergyWindow[{index + 1}]") for index, item in enumerate(value)]
    if result[0] > result[1]:
        raise FittingError("InvalidEnergyWindow", "EnergyWindow minimum exceeds maximum")
    return result


def _k_neighborhood(value: Any) -> dict[str, Any]:
    if value is None or value == "All":
        return {"mode": "All"}
    if not isinstance(value, Mapping):
        raise FittingError(
            "InvalidKPointNeighborhood", "KPointNeighborhood must be All/None or a mapping"
        )
    allowed = {"Center", "Radius", "Range", "Periodic"}
    if set(value) - allowed or "Center" not in value:
        raise FittingError("InvalidKPointNeighborhood", "invalid neighborhood fields")
    center = value["Center"]
    if isinstance(center, (str, bytes)) or not isinstance(center, Sequence) or len(center) != 3:
        raise FittingError("InvalidKPointNeighborhood", "Center must have length three")
    center_values = [_real(item, f"Center[{index + 1}]") for index, item in enumerate(center)]
    periodic = value.get("Periodic", True)
    if not isinstance(periodic, bool):
        raise FittingError("InvalidKPointNeighborhood", "Periodic must be Boolean")
    has_radius = "Radius" in value
    has_range = "Range" in value
    if has_radius == has_range:
        raise FittingError(
            "InvalidKPointNeighborhood", "specify exactly one of Radius or Range"
        )
    if has_radius:
        return {
            "mode": "Radius",
            "center": center_values,
            "radius": _real(value["Radius"], "Radius"),
            "periodic": periodic,
        }
    ranges = value["Range"]
    if isinstance(ranges, (str, bytes)) or not isinstance(ranges, Sequence) or len(ranges) != 3:
        raise FittingError("InvalidKPointNeighborhood", "Range must have length three")
    return {
        "mode": "Range",
        "center": center_values,
        "range": [_real(item, f"Range[{index + 1}]") for index, item in enumerate(ranges)],
        "periodic": periodic,
    }


def _weights(value: Any, records: int, bands: int) -> Optional[list[list[float]]]:
    if value is None or value == "Automatic":
        return None
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) != records:
        raise FittingError("InvalidResidualWeights", "ResidualWeights shape is invalid")
    result: list[list[float]] = []
    for row_index, row in enumerate(value):
        if isinstance(row, (str, bytes)) or not isinstance(row, Sequence) or len(row) != bands:
            raise FittingError("InvalidResidualWeights", "ResidualWeights shape is invalid")
        numeric = [_real(item, f"ResidualWeights[{row_index + 1}][{column + 1}]") for column, item in enumerate(row)]
        if any(item < 0.0 for item in numeric):
            raise FittingError("InvalidResidualWeights", "ResidualWeights must be nonnegative")
        result.append(numeric)
    return result


@dataclass(frozen=True)
class BandSpan:
    """Inclusive positive-step counterpart of a Mathematica fitting band Span."""

    start: int
    end: int
    step: int = 1

    def indices(self) -> tuple[int, ...]:
        if (
            isinstance(self.start, bool)
            or isinstance(self.end, bool)
            or isinstance(self.step, bool)
            or not all(isinstance(item, int) for item in (self.start, self.end, self.step))
            or self.step <= 0
            or self.end < self.start
        ):
            raise FittingError("InvalidBandSelection", "BandSpan is invalid")
        return tuple(range(self.start, self.end + 1, self.step))


@dataclass(frozen=True)
class BandFitComparison:
    k_point_indices: tuple[int, ...]
    reference_bands: tuple[tuple[float, ...], ...]
    initial_bands: tuple[tuple[float, ...], ...]
    fitted_bands: tuple[tuple[float, ...], ...]
    selection_mask: tuple[tuple[bool, ...], ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "k_point_indices": list(self.k_point_indices),
            "reference_bands": [list(row) for row in self.reference_bands],
            "initial_bands": [list(row) for row in self.initial_bands],
            "fitted_bands": [list(row) for row in self.fitted_bands],
            "selection_mask": [list(row) for row in self.selection_mask],
        }


@dataclass(frozen=True)
class BandFittingResult:
    fitted_parameters: Mapping[str, float]
    objective: str
    optimizer: str
    parameter_model: str
    converged: bool
    iterations: int
    final_damping: float
    initial_loss: float
    fitted_loss: float
    candidate_k_point_count: int
    used_k_point_indices: tuple[int, ...]
    unique_k_point_count: int
    residual_count: int
    selection_restricted: bool
    selected_band_indices_by_k_point: tuple[tuple[int, ...], ...]
    band_selection: tuple[int, ...]
    energy_window: Optional[tuple[float, float]]
    k_point_neighborhood: Mapping[str, Any]
    weighting: str
    band_pairing: str
    comparison: BandFitComparison

    def to_dict(self) -> dict[str, Any]:
        return {
            "FittedParams": dict(self.fitted_parameters),
            "Objective": self.objective,
            "Optimizer": self.optimizer,
            "ParameterModel": self.parameter_model,
            "Converged": self.converged,
            "Iterations": self.iterations,
            "FinalDamping": self.final_damping,
            "LSQForInitParams": self.initial_loss,
            "LSQForFittedParams": self.fitted_loss,
            "CandidateKPointCount": self.candidate_k_point_count,
            "UsedKPointCount": len(self.used_k_point_indices),
            "UniqueKPointCount": self.unique_k_point_count,
            "UsedKPointIndices": list(self.used_k_point_indices),
            "ResidualCount": self.residual_count,
            "SelectionRestricted": self.selection_restricted,
            "SelectedBandIndicesByKPoint": [
                list(indices) for indices in self.selected_band_indices_by_k_point
            ],
            "BandSelection": list(self.band_selection),
            "EnergyWindow": "All" if self.energy_window is None else list(self.energy_window),
            "KPointNeighborhood": dict(self.k_point_neighborhood),
            "Weighting": self.weighting,
            "BandPairing": self.band_pairing,
            "ComparisonData": self.comparison.to_dict(),
        }


def _fit_result(value: Mapping[str, Any]) -> BandFittingResult:
    fitted_parameters = MappingProxyType(
        {str(item["name"]): float(item["value"]) for item in value["FittedParams"]}
    )
    energy_window = value["EnergyWindow"]
    comparison = BandFitComparison(
        k_point_indices=tuple(int(item) for item in value["UsedKPointIndices"]),
        reference_bands=tuple(tuple(map(float, row)) for row in value["ReferenceBands"]),
        initial_bands=tuple(tuple(map(float, row)) for row in value["InitialBands"]),
        fitted_bands=tuple(tuple(map(float, row)) for row in value["FittedBands"]),
        selection_mask=tuple(tuple(map(bool, row)) for row in value["SelectionMask"]),
    )
    return BandFittingResult(
        fitted_parameters=fitted_parameters,
        objective=str(value["Objective"]),
        optimizer=str(value["Optimizer"]),
        parameter_model=str(value["ParameterModel"]),
        converged=bool(value["Converged"]),
        iterations=int(value["Iterations"]),
        final_damping=float(value["FinalDamping"]),
        initial_loss=float(value["LSQForInitParams"]),
        fitted_loss=float(value["LSQForFittedParams"]),
        candidate_k_point_count=int(value["CandidateKPointCount"]),
        used_k_point_indices=tuple(int(item) for item in value["UsedKPointIndices"]),
        unique_k_point_count=int(value["UniqueKPointCount"]),
        residual_count=int(value["ResidualCount"]),
        selection_restricted=bool(value["SelectionRestricted"]),
        selected_band_indices_by_k_point=tuple(
            tuple(int(item) for item in row)
            for row in value["SelectedBandIndicesByKPoint"]
        ),
        band_selection=tuple(int(item) for item in value["BandSelection"]),
        energy_window=None
        if energy_window is None or energy_window == "All"
        else (float(energy_window[0]), float(energy_window[1])),
        k_point_neighborhood=MappingProxyType(dict(value["KPointNeighborhood"])),
        weighting=str(value["Weighting"]),
        band_pairing=str(value["BandPairing"]),
        comparison=comparison,
    )


def vaspEig(
    filename: str | Path,
    efermi: Any,
    spin: int,
    startband: int,
    endband: int,
) -> list[list[list[float]]]:
    """Parse a VASP EIGENVAL with the stable 1-based spin and band convention."""

    if not isinstance(filename, (str, Path)):
        raise FittingError("VaspEigenvalFile", "filename must be a path")
    if any(isinstance(value, bool) or not isinstance(value, int) for value in (spin, startband, endband)):
        raise FittingError("InvalidVaspBandRange", "spin and band bounds must be integers")
    result = _rust_fitting(
        "parse_vasp_eigenval",
        path=str(filename),
        fermi_energy=_real(efermi, "efermi"),
        spin=spin,
        start_band=startband,
        end_band=endband,
    )
    return [[list(record[0]), list(record[1])] for record in result["Records"]]


def fittingTB(
    h: Any,
    eigdata: Any,
    krange: Sequence[int],
    initparms: Any,
    *,
    MaxIterations: int = 100,
    TimeConstraint: Optional[Any] = 60,
    EnergyWindow: Any = None,
    KPointNeighborhood: Any = None,
    BandSelection: Any = None,
    ResidualWeights: Any = None,
    FiniteDifferenceStep: Any = 1.0e-5,
    InitialDamping: Any = 1.0e-3,
    FitTolerance: Any = 1.0e-8,
    HermitianTolerance: Any = 1.0e-10,
) -> BandFittingResult:
    """Fit an affine MagneticTB Hamiltonian to reference bands entirely in Rust."""

    points, reference_bands = _reference_data(eigdata)
    parameter_names = _parameter_names(h)
    initial = _rules(initparms, "initparms")
    missing = [name for name in parameter_names if name not in initial]
    if missing:
        raise FittingError(
            "InvalidInitialParameters",
            f"missing initial values for: {', '.join(missing)}",
        )
    initial_values = [_real(initial[name], f"initial parameter {name}") for name in parameter_names]
    if isinstance(krange, (str, bytes)) or not isinstance(krange, Sequence) or not krange:
        raise FittingError("InvalidKRange", "krange must be a nonempty 1-based sequence")
    if any(isinstance(index, bool) or not isinstance(index, int) for index in krange):
        raise FittingError("InvalidKRange", "krange must contain 1-based integers")
    if any(index < 1 or index > len(points) for index in krange):
        raise FittingError("InvalidKRange", "krange contains an out-of-range index")
    if not isinstance(MaxIterations, int) or isinstance(MaxIterations, bool):
        raise FittingError("InvalidFittingOption", "MaxIterations must be an integer")
    if TimeConstraint is not None:
        TimeConstraint = _real(TimeConstraint, "TimeConstraint")
    band_count = len(reference_bands[0])
    selection = _band_selection(BandSelection, band_count)
    zero_rules = {name: 0.0 for name in parameter_names}
    constant_matrices = _matrices(h, zero_rules, points)
    unit_matrices = []
    for name in parameter_names:
        unit_rules = dict(zero_rules)
        unit_rules[name] = 1.0
        unit_matrices.append(_matrices(h, unit_rules, points))
    result = _rust_fitting(
        "fit_bands",
        parameter_names=list(parameter_names),
        k_points=[list(point) for point in points],
        reference_bands=[list(row) for row in reference_bands],
        k_range=[index - 1 for index in krange],
        initial_values=initial_values,
        constant_matrices=constant_matrices,
        unit_parameter_matrices=unit_matrices,
        options={
            "max_iterations": MaxIterations,
            "time_constraint": TimeConstraint,
            "energy_window": _energy_window(EnergyWindow),
            "k_point_neighborhood": _k_neighborhood(KPointNeighborhood),
            "band_selection": [index - 1 for index in selection],
            "residual_weights": _weights(ResidualWeights, len(points), band_count),
            "finite_difference_step": _real(FiniteDifferenceStep, "FiniteDifferenceStep"),
            "initial_damping": _real(InitialDamping, "InitialDamping"),
            "fit_tolerance": _real(FitTolerance, "FitTolerance"),
            "hermitian_tolerance": _real(HermitianTolerance, "HermitianTolerance"),
        },
    )
    return _fit_result(result)


@dataclass(frozen=True)
class BandFittingExplorer:
    """Immutable Python/Web counterpart of stable ``bandManipulateEig``."""

    hamiltonian: Any
    k_points: tuple[tuple[float, float, float], ...]
    reference_bands: tuple[tuple[float, ...], ...]
    parameter_names: tuple[str, ...]
    parameter_range: tuple[float, float] = (-1.0, 1.0)
    hermitian_tolerance: float = 1.0e-10

    def defaults(self) -> Mapping[str, float]:
        return MappingProxyType({name: 0.0 for name in self.parameter_names})

    def evaluate(self, rules: Optional[Mapping[str, Any]] = None) -> BandFitComparison:
        assigned = dict(self.defaults())
        if rules is not None:
            unknown = sorted(set(rules) - set(self.parameter_names))
            if unknown:
                raise FittingError(
                    "UnsupportedFittingParameter",
                    f"unknown fitting parameter(s): {', '.join(unknown)}",
                )
            assigned.update(
                {name: _real(value, f"parameter {name}") for name, value in rules.items()}
            )
        matrices = _matrices(self.hamiltonian, assigned, self.k_points)
        from .core_bindings import properties

        values = properties(
            "band_eigenvalues",
            matrices=matrices,
            hermitian_tolerance=self.hermitian_tolerance,
        )["Eigenvalues"]
        model = tuple(tuple(map(float, row)) for row in values)
        mask = tuple(tuple(True for _ in row) for row in self.reference_bands)
        indices = tuple(range(1, len(self.k_points) + 1))
        return BandFitComparison(indices, self.reference_bands, model, model, mask)


def bandManipulateEig(
    h: Any,
    eigdata: Any,
    *,
    HermitianTolerance: Any = 1.0e-10,
) -> BandFittingExplorer:
    """Create immutable fitting controls over stable reference-band records."""

    points, reference_bands = _reference_data(eigdata)
    return BandFittingExplorer(
        hamiltonian=_frozen_matrix(h),
        k_points=points,
        reference_bands=reference_bands,
        parameter_names=_parameter_names(h),
        hermitian_tolerance=_finite_real(HermitianTolerance, "HermitianTolerance"),
    )


@dataclass(frozen=True)
class BandPathComparison(_BandPlotDisplay):
    """Cached TB/reference bands with notebook display and show/plot/savefig."""

    model: BandPlotResult
    reference_k_points: tuple[tuple[float, float, float], ...]
    reference_bands: tuple[tuple[float, ...], ...]
    reference_x_coordinates: tuple[float, ...]
    plot_range: Optional[tuple[float, float]]

    @property
    def _style(self) -> _BandPlotStyle:
        return self.model._style

    def _draw(self, ax: Any) -> None:
        # Draw the reference first so BandPlotResult's automatic energy limits
        # include both datasets, before adding its zero/boundary gridlines.
        # Keep the existing reference x coordinates and energy order unchanged.
        reference_lines = []
        for band in range(len(self.reference_bands[0])):
            for segment in range(len(self.model.path)):
                start = segment * self.model.npoint
                stop = start + self.model.npoint
                reference_lines.extend(
                    ax.plot(
                        self.reference_x_coordinates[start:stop],
                        [row[band] for row in self.reference_bands[start:stop]],
                        color="#8080ff",
                        marker="." if self.model.npoint == 1 else None,
                    )
                )
        model_start = len(ax.lines)
        self.model._draw(ax)
        model_count = self.model.dimension * len(self.model.path)
        model_lines = ax.lines[model_start:model_start + model_count]
        for line in model_lines:
            line.set_color("purple")
        if self.plot_range is not None:
            ax.set_ylim(*self.plot_range)
        ax.legend(
            (model_lines[0], reference_lines[0]),
            ("TB", "Reference"),
            fontsize=12,
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "model": self.model.to_dict(),
            "reference_k_points": [list(point) for point in self.reference_k_points],
            "reference_bands": [list(row) for row in self.reference_bands],
            "reference_x_coordinates": list(self.reference_x_coordinates),
            "plot_range": None if self.plot_range is None else list(self.plot_range),
        }


def compareBand(
    pathstr: Any,
    npoint: int,
    ham: Any,
    rules: Any,
    vaspband: Any,
    *,
    plotRange: Any = None,
    HermitianTolerance: Any = 1.0e-10,
) -> BandPathComparison:
    """Compare cached TB/reference bands and display their overlay in Jupyter.

    Use show(), plot(), or savefig() for explicit graphics; to_dict() retains
    the existing data-only interface. Rendering does not refit or solve bands.
    """

    path = _path(pathstr)
    subdivisions = _subdivisions(npoint)
    points, reference_bands = _reference_data(vaspband)
    if len(points) != len(path) * subdivisions:
        raise FittingError(
            "InvalidReferenceBands",
            "compareBand requires exactly npoint reference records per path segment",
        )
    model = bandplot(
        path,
        subdivisions,
        ham,
        rules,
        HermitianTolerance=HermitianTolerance,
    )
    x_values: list[float] = []
    for segment in range(len(path)):
        start = subdivisions * segment
        if subdivisions == 1:
            x_values.append(float(start))
        else:
            x_values.extend(
                start + subdivisions * index / (subdivisions - 1)
                for index in range(subdivisions)
            )
    parsed_range = _energy_window(plotRange)
    return BandPathComparison(
        model=model,
        reference_k_points=points,
        reference_bands=reference_bands,
        reference_x_coordinates=tuple(x_values),
        plot_range=None
        if parsed_range is None
        else (parsed_range[0], parsed_range[1]),
    )


def vasp_eig(*args: Any, **kwargs: Any) -> list[list[list[float]]]:
    return vaspEig(*args, **kwargs)


def fitting_tb(*args: Any, **kwargs: Any) -> BandFittingResult:
    return fittingTB(*args, **kwargs)


def band_manipulate_eig(*args: Any, **kwargs: Any) -> BandFittingExplorer:
    return bandManipulateEig(*args, **kwargs)


def compare_band(*args: Any, **kwargs: Any) -> BandPathComparison:
    return compareBand(*args, **kwargs)
