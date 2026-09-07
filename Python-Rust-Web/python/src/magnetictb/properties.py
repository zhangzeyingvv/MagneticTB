"""Stable-aligned numerical properties backed exclusively by Rust algorithms."""

from __future__ import annotations

import math
import json
from collections.abc import Callable, Mapping, Sequence
from numbers import Number
from typing import Any, Optional

from .api import ExactError, PropertiesError, gapless_points_json, geometry, properties
from .linear_algebra import _encode_evaluation_scalar
from .model import current_model
from .tight_binding import Hamiltonian


def _finite_real(value: Any, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Number) or isinstance(value, complex):
        raise TypeError(f"{field} must be a finite real number")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{field} must be finite")
    return result


def _real_vector(value: Any, field: str) -> list[float]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise TypeError(f"{field} must be a nonempty real sequence")
    return [_finite_real(item, f"{field}[{index}]") for index, item in enumerate(value)]


def _complex_scalar(value: Any, field: str) -> dict[str, float]:
    if isinstance(value, bool) or not isinstance(value, Number):
        raise TypeError(f"{field} must be numeric")
    result = complex(value)
    if not (math.isfinite(result.real) and math.isfinite(result.imag)):
        raise ValueError(f"{field} must be finite")
    return {"real": result.real, "imaginary": result.imag}


def _complex_matrix(value: Any, field: str) -> list[list[dict[str, float]]]:
    if hasattr(value, "tolist"):
        value = value.tolist()
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise TypeError(f"{field} must be a nonempty numeric matrix")
    rows = []
    columns: Optional[int] = None
    for row_index, row in enumerate(value):
        if hasattr(row, "tolist"):
            row = row.tolist()
        if isinstance(row, (str, bytes)) or not isinstance(row, Sequence) or not row:
            raise TypeError(f"{field}[{row_index}] must be a nonempty row")
        if columns is None:
            columns = len(row)
        elif len(row) != columns:
            raise ValueError(f"{field} rows must have equal length")
        rows.append(
            [
                _complex_scalar(entry, f"{field}[{row_index}][{column_index}]")
                for column_index, entry in enumerate(row)
            ]
        )
    return rows


def _decode_numeric(value: Any) -> Any:
    if isinstance(value, dict):
        if set(value) == {"real", "imaginary"}:
            return complex(float(value["real"]), float(value["imaginary"]))
        return {key: _decode_numeric(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_decode_numeric(item) for item in value]
    return value


def _evaluator(
    hamiltonian: Any,
    parameters: Optional[Mapping[str, Any]],
) -> Callable[[Sequence[float]], Any]:
    if isinstance(hamiltonian, Hamiltonian):
        if parameters is None:
            raise TypeError("parameters are required for a symbolic Hamiltonian")
        return lambda point: hamiltonian.evaluate(parameters, point)
    if callable(hamiltonian):
        if parameters is not None:
            raise TypeError("parameters apply only to a MagneticTB Hamiltonian")
        return lambda point: hamiltonian(list(point))
    if parameters is not None:
        raise TypeError("parameters apply only to a MagneticTB Hamiltonian")
    _complex_matrix(hamiltonian, "hamiltonian")
    return lambda _point: hamiltonian


def _sample(
    evaluator: Callable[[Sequence[float]], Any],
    path: Sequence[Sequence[float]],
) -> list[list[list[dict[str, float]]]]:
    return [
        _complex_matrix(evaluator(point), f"Hamiltonian at path point {index + 1}")
        for index, point in enumerate(path)
    ]


def _wilson_options(
    hermitian_tolerance: float,
    gap_tolerance: float,
    covariance_tolerance: float,
    overlap_tolerance: float,
) -> dict[str, float]:
    return {
        "hermitian_tolerance": _finite_real(
            hermitian_tolerance, "HermitianTolerance"
        ),
        "gap_tolerance": _finite_real(gap_tolerance, "GapTolerance"),
        "covariance_tolerance": _finite_real(
            covariance_tolerance, "CovarianceTolerance"
        ),
        "overlap_tolerance": _finite_real(
            overlap_tolerance, "OverlapTolerance"
        ),
    }


def wilsonLoop(
    hamiltonian: Any,
    centers: Sequence[Sequence[Any]],
    occupied: int,
    start: Sequence[Any],
    end: Sequence[Any],
    *,
    parameters: Optional[Mapping[str, Any]] = None,
    PathSubdivisions: int = 50,
    HermitianTolerance: float = 1.0e-10,
    GapTolerance: float = 1.0e-9,
    CovarianceTolerance: float = 1.0e-8,
    OverlapTolerance: float = 1.0e-10,
    Output: str = "PhaseOverPi",
) -> Any:
    """Match stable ``wilsonLoop`` using Rust eigensystems and SVD links."""

    if isinstance(PathSubdivisions, bool) or not isinstance(PathSubdivisions, int):
        raise TypeError("PathSubdivisions must be a positive integer")
    if PathSubdivisions < 1:
        raise PropertiesError(
            "InvalidPropertiesOption", "PathSubdivisions must be a positive integer"
        )
    first = _real_vector(start, "start")
    last = _real_vector(end, "end")
    if len(first) != len(last):
        raise PropertiesError("InvalidWilsonPath", "start and end dimensions differ")
    path = [
        [
            left + (right - left) * index / PathSubdivisions
            for left, right in zip(first, last)
        ]
        for index in range(PathSubdivisions + 1)
    ]
    result = _decode_numeric(
        properties(
            "wilson_loop",
            matrices=_sample(_evaluator(hamiltonian, parameters), path),
            path=path,
            centers=[_real_vector(center, "centers") for center in centers],
            occupied=occupied,
            options=_wilson_options(
                HermitianTolerance,
                GapTolerance,
                CovarianceTolerance,
                OverlapTolerance,
            ),
        )
    )
    outputs = {
        "PhaseOverPi": "PhasesOverPi",
        "WannierCenters": "WannierCenters",
        "Eigenvalues": "Eigenvalues",
        "Data": None,
    }
    if Output not in outputs:
        raise PropertiesError("InvalidPropertiesOption", f"unknown Output {Output}")
    return result if outputs[Output] is None else result[outputs[Output]]


def berryPhase(
    hamiltonian: Any,
    centers: Sequence[Sequence[Any]],
    occupied: int,
    path: Sequence[Sequence[Any]],
    *,
    parameters: Optional[Mapping[str, Any]] = None,
    HermitianTolerance: float = 1.0e-10,
    GapTolerance: float = 1.0e-9,
    CovarianceTolerance: float = 1.0e-8,
    OverlapTolerance: float = 1.0e-10,
    Output: str = "Phase",
) -> Any:
    """Match stable ``berryPhase`` with Rust occupied-subspace links."""

    numeric_path = [_real_vector(point, "path") for point in path]
    result = _decode_numeric(
        properties(
            "berry_phase",
            matrices=_sample(_evaluator(hamiltonian, parameters), numeric_path),
            path=numeric_path,
            centers=[_real_vector(center, "centers") for center in centers],
            occupied=occupied,
            options=_wilson_options(
                HermitianTolerance,
                GapTolerance,
                CovarianceTolerance,
                OverlapTolerance,
            ),
        )
    )
    outputs = {
        "Phase": "Phase",
        "PhaseOverPi": "PhaseOverPi",
        "WilsonDeterminant": "WilsonDeterminant",
        "Data": None,
    }
    if Output not in outputs:
        raise PropertiesError("InvalidPropertiesOption", f"unknown Output {Output}")
    return result if outputs[Output] is None else result[outputs[Output]]


def z2path(
    hamiltonian: Any,
    occupied: int,
    path: Sequence[Sequence[Any]],
    *,
    parameters: Optional[Mapping[str, Any]] = None,
    HermitianTolerance: float = 1.0e-10,
) -> list[float]:
    """Reproduce the stable 2.0.10 compatibility ``z2path`` routine."""

    numeric_path = [_real_vector(point, "path") for point in path]
    if len(numeric_path) < 2:
        raise PropertiesError("InvalidWilsonPath", "z2path requires at least two path points")
    matrices = _sample(_evaluator(hamiltonian, parameters), numeric_path[:-1])
    result = properties(
        "legacy_z2_path",
        matrices=matrices,
        occupied=occupied,
        hermitian_tolerance=_finite_real(HermitianTolerance, "HermitianTolerance"),
    )
    return [float(value) for value in result["Values"]]


def berryph(
    hamiltonian: Any,
    occupied: int,
    path: Sequence[Sequence[Any]],
    *,
    parameters: Optional[Mapping[str, Any]] = None,
    HermitianTolerance: float = 1.0e-10,
) -> list[float]:
    """Reproduce the stable 2.0.10 compatibility ``berryph`` routine."""

    numeric_path = [_real_vector(point, "path") for point in path]
    if len(numeric_path) < 2:
        raise PropertiesError("InvalidWilsonPath", "berryph requires at least two path points")
    result = properties(
        "legacy_berryph",
        matrices=_sample(_evaluator(hamiltonian, parameters), numeric_path),
        occupied=occupied,
        hermitian_tolerance=_finite_real(HermitianTolerance, "HermitianTolerance"),
    )
    return [float(value) for value in result["Values"]]


def _legacy_piecewise_path(
    segments: Sequence[Sequence[Sequence[Any]]],
    subdivisions: int,
    field: str,
    initial: Optional[Sequence[float]] = None,
) -> list[list[float]]:
    if isinstance(segments, (str, bytes)) or not isinstance(segments, Sequence) or not segments:
        raise TypeError(f"{field} must be a nonempty sequence of start/end segments")
    parsed: list[tuple[list[float], list[float]]] = []
    for index, segment in enumerate(segments):
        if (
            isinstance(segment, (str, bytes))
            or not isinstance(segment, Sequence)
            or len(segment) != 2
        ):
            raise TypeError(f"{field}[{index}] must contain start and end points")
        start = _real_vector(segment[0], f"{field}[{index}][0]")
        end = _real_vector(segment[1], f"{field}[{index}][1]")
        if len(start) != len(end):
            raise PropertiesError("InvalidWilsonPath", f"{field}[{index}] dimensions differ")
        parsed.append((start, end))
    dimension = len(parsed[0][0])
    if any(len(start) != dimension for start, _end in parsed):
        raise PropertiesError("InvalidWilsonPath", f"{field} point dimensions differ")
    result = [list(parsed[0][0] if initial is None else initial)]
    for start, end in parsed:
        result.extend(
            [
                [
                    left + (right - left) * sample / subdivisions
                    for left, right in zip(start, end)
                ]
                for sample in range(1, subdivisions + 1)
            ]
        )
    return result


def wLoop(
    hamiltonian: Any,
    occupied: int,
    path1: Sequence[Sequence[Sequence[Any]]],
    path2: Sequence[Sequence[Sequence[Any]]],
    *,
    parameters: Optional[Mapping[str, Any]] = None,
    HermitianTolerance: float = 1.0e-10,
) -> list[list[float]]:
    """Reproduce stable 2.0.10 ``wLoop`` with its fixed 20/1020 sampling."""

    loop_path = _legacy_piecewise_path(path1, 20, "path1")
    offset_path = _legacy_piecewise_path(path2, 1020, "path2", initial=loop_path[0])
    if any(len(point) != len(loop_path[0]) for point in offset_path):
        raise PropertiesError("InvalidWilsonPath", "path1 and path2 dimensions differ")
    evaluator = _evaluator(hamiltonian, parameters)
    matrix_grid = [
        _sample(
            evaluator,
            [
                [left + right for left, right in zip(point, offset)]
                for point in loop_path
            ],
        )
        for offset in offset_path
    ]
    result = properties(
        "legacy_wloop",
        matrix_grid=matrix_grid,
        occupied=occupied,
        hermitian_tolerance=_finite_real(HermitianTolerance, "HermitianTolerance"),
    )
    return [[float(value) for value in row] for row in result["Values"]]


def berryCurvature(
    hamiltonian: Any,
    centers: Sequence[Sequence[Any]],
    occupied: int,
    point: Sequence[Any],
    *,
    parameters: Optional[Mapping[str, Any]] = None,
    Directions: Sequence[int] = (1, 2),
    StepSize: Any = 1.0e-3,
    HermitianTolerance: float = 1.0e-10,
    GapTolerance: float = 1.0e-9,
    CovarianceTolerance: float = 1.0e-8,
    OverlapTolerance: float = 1.0e-10,
    Output: str = "Curvature",
) -> Any:
    """Match stable ``berryCurvature`` on an oriented Rust-defined plaquette."""

    numeric_point = _real_vector(point, "point")
    if (
        isinstance(Directions, (str, bytes))
        or not isinstance(Directions, Sequence)
        or len(Directions) != 2
        or any(isinstance(index, bool) or not isinstance(index, int) for index in Directions)
    ):
        raise TypeError("Directions must contain two one-based integer indices")
    directions = [index - 1 for index in Directions]
    if isinstance(StepSize, Number) and not isinstance(StepSize, (bool, complex)):
        steps = [_finite_real(StepSize, "StepSize")] * 2
    else:
        steps = _real_vector(StepSize, "StepSize")
        if len(steps) != 2:
            raise TypeError("StepSize must be one number or a pair")
    options = _wilson_options(
        HermitianTolerance,
        GapTolerance,
        CovarianceTolerance,
        OverlapTolerance,
    )
    path = properties(
        "berry_plaquette_path",
        point=numeric_point,
        directions=directions,
        step_size=steps,
        options=options,
    )["PlaquettePath"]
    result = _decode_numeric(
        properties(
            "berry_curvature",
            matrices=_sample(_evaluator(hamiltonian, parameters), path),
            centers=[_real_vector(center, "centers") for center in centers],
            occupied=occupied,
            point=numeric_point,
            directions=directions,
            step_size=steps,
            options=options,
        )
    )
    outputs = {"Curvature": "Curvature", "Flux": "Flux", "Data": None}
    if Output not in outputs:
        raise PropertiesError("InvalidPropertiesOption", f"unknown Output {Output}")
    return result if outputs[Output] is None else result[outputs[Output]]


def _hopping_mapping(value: Any) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise TypeError(
            "hopping input must be a mapping returned by hoppingData, "
            "transformHoppings, or a compatible stable data record"
        )
    return _encode_numeric_tree(dict(value))


def _encode_numeric_tree(value: Any) -> Any:
    if isinstance(value, complex):
        if not (math.isfinite(value.real) and math.isfinite(value.imag)):
            raise ValueError("hopping data contains a non-finite complex value")
        return {"real": value.real, "imaginary": value.imag}
    if isinstance(value, Mapping):
        return {str(key): _encode_numeric_tree(item) for key, item in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return [_encode_numeric_tree(item) for item in value]
    return value


def _stable_output(result: dict[str, Any], output: str, *, scalar: str = "Hamiltonian") -> Any:
    if output == "Data":
        return result
    if output == "Matrix":
        return result[scalar]
    raise PropertiesError("InvalidPropertiesOption", f"unknown Output {output}")


def hoppingData(
    selection: int | Sequence[int],
    rules: Mapping[str, Any] | Sequence[tuple[str, Any]],
    *,
    Hermitian: bool = True,
    KernelMethod: str = "Iterative",
    ValidationLevel: str = "Basic",
) -> dict[str, Any]:
    """Return the stable numeric real-space cache for solved session shells."""

    model = current_model()
    if isinstance(selection, bool):
        raise PropertiesError("InvalidShellSelection", "shell selection must be positive")
    if isinstance(selection, int):
        shells = list(range(1, selection + 1)) if selection > 0 else []
    elif isinstance(selection, Sequence) and not isinstance(selection, (str, bytes)):
        shells = []
        for shell in selection:
            if isinstance(shell, bool) or not isinstance(shell, int) or shell < 1:
                shells = []
                break
            if shell not in shells:
                shells.append(shell)
    else:
        shells = []
    if not shells:
        raise PropertiesError(
            "InvalidShellSelection",
            "selection must be a positive integer or nonempty positive-integer sequence",
        )
    if any(shell > model.initial_bond_shells for shell in shells):
        raise PropertiesError(
            "ShellOutOfRange",
            f"only {model.initial_bond_shells} shells were prepared",
        )
    if not isinstance(Hermitian, bool) or not isinstance(KernelMethod, str) or not isinstance(ValidationLevel, str):
        raise PropertiesError("InvalidPropertiesOption", "invalid hoppingData options")
    if isinstance(rules, Mapping):
        pairs = list(rules.items())
    elif isinstance(rules, Sequence) and not isinstance(rules, (str, bytes)):
        pairs = list(rules)
        if any(not isinstance(pair, Sequence) or len(pair) != 2 for pair in pairs):
            raise TypeError("rules must be a mapping or sequence of (name, value) pairs")
    else:
        raise TypeError("rules must be a mapping or sequence of (name, value) pairs")
    parameters = {
        str(name): _encode_evaluation_scalar(value, f"parameter {name}")
        for name, value in pairs
    }
    try:
        shell_results = model._cached_hopping_shell_results(
            shells,
            hermitian=Hermitian,
            kernel_method=KernelMethod,
            validation_level=ValidationLevel,
        )
        result = geometry(
            "hopping_data_from_shell_results",
            shell_results=shell_results,
            parameters=parameters,
        )
    except ExactError as error:
        raise PropertiesError(error.tag, error.detail) from None
    return _decode_numeric(result)


def transformHoppings(
    hopping_data: Mapping[str, Any],
    cell_matrix: Sequence[Sequence[int]],
    *,
    WannierCenters: Any = "Automatic",
    Lattice: Any = "Automatic",
    HermiticityTolerance: float = 1.0e-9,
) -> dict[str, Any]:
    """Transform stable hopping data to an integer supercell in Rust."""

    data = _hopping_mapping(hopping_data)
    for key, setting in (("WannierCenters", WannierCenters), ("Lattice", Lattice)):
        if isinstance(setting, str) and setting == "Automatic":
            continue
        if setting is None:
            data.pop(key, None)
        else:
            data[key] = _encode_numeric_tree(setting)
    return _decode_numeric(
        properties(
            "transform_hoppings",
            data=data,
            cell_matrix=cell_matrix,
            hermiticity_tolerance=_finite_real(
                HermiticityTolerance, "HermiticityTolerance"
            ),
        )
    )


def buildBlochHamiltonian(
    hopping_data: Mapping[str, Any],
    momentum: Sequence[Any],
    *,
    CellMatrix: Any = "Automatic",
    HermiticityTolerance: float = 1.0e-9,
    Output: str = "Matrix",
) -> Any:
    """Build ``H(k)=Sum_R H(R) exp(2 pi i k.R)`` in Rust."""

    arguments: dict[str, Any] = {
        "data": _hopping_mapping(hopping_data),
        "momentum": _real_vector(momentum, "momentum"),
        "hermiticity_tolerance": _finite_real(
            HermiticityTolerance, "HermiticityTolerance"
        ),
    }
    if not (isinstance(CellMatrix, str) and CellMatrix == "Automatic"):
        arguments["cell_matrix"] = CellMatrix
    result = _decode_numeric(properties("build_bloch_hamiltonian", **arguments))
    return _stable_output(result, Output)


def buildRealSpaceHamiltonian(
    hopping_data: Mapping[str, Any],
    geometry: Sequence[Any],
    *,
    BoundaryConditions: Sequence[str] = ("Open", "Open", "Open"),
    HermiticityTolerance: float = 1.0e-9,
    Output: str = "Matrix",
) -> Any:
    """Assemble a rectangular or explicit-cell finite Hamiltonian in Rust."""

    result = _decode_numeric(
        properties(
            "build_real_space_hamiltonian",
            data=_hopping_mapping(hopping_data),
            geometry=geometry,
            boundary_conditions=list(BoundaryConditions),
            hermiticity_tolerance=_finite_real(
                HermiticityTolerance, "HermiticityTolerance"
            ),
        )
    )
    return _stable_output(result, Output)


def buildSlabHamiltonian(
    hopping_data: Mapping[str, Any],
    size: Sequence[int],
    momentum: Sequence[Any],
    *,
    PeriodicDirections: Sequence[int] = (1, 2),
    CellMatrix: Any = "Automatic",
    HermiticityTolerance: float = 1.0e-9,
    Output: str = "Matrix",
) -> Any:
    """Build the stable hybrid finite/periodic Hamiltonian in Rust."""

    arguments: dict[str, Any] = {
        "data": _hopping_mapping(hopping_data),
        "size": list(size),
        "momentum": _real_vector(momentum, "momentum"),
        "periodic_directions": list(PeriodicDirections),
        "hermiticity_tolerance": _finite_real(
            HermiticityTolerance, "HermiticityTolerance"
        ),
    }
    if not (isinstance(CellMatrix, str) and CellMatrix == "Automatic"):
        arguments["cell_matrix"] = CellMatrix
    result = _decode_numeric(properties("build_slab_hamiltonian", **arguments))
    return _stable_output(result, Output)


def surfaceGreenFunction(
    hopping_data: Mapping[str, Any],
    momentum: Sequence[Any],
    energy: Any,
    *,
    CellMatrix: Any = "Automatic",
    Broadening: float = 1.0e-3,
    Tolerance: float = 1.0e-10,
    MaxIterations: int = 200,
    Surface: str = "Positive",
    HermiticityTolerance: float = 1.0e-9,
    Output: str = "GreenFunction",
) -> Any:
    """Compute the stable iterative surface Green function in Rust."""

    if isinstance(MaxIterations, bool) or not isinstance(MaxIterations, int):
        raise TypeError("MaxIterations must be a positive integer")
    arguments: dict[str, Any] = {
        "data": _hopping_mapping(hopping_data),
        "momentum": _real_vector(momentum, "momentum"),
        "energy": _finite_real(energy, "energy"),
        "broadening": _finite_real(Broadening, "Broadening"),
        "tolerance": _finite_real(Tolerance, "Tolerance"),
        "max_iterations": MaxIterations,
        "surface": Surface,
        "hermiticity_tolerance": _finite_real(
            HermiticityTolerance, "HermiticityTolerance"
        ),
    }
    if not (isinstance(CellMatrix, str) and CellMatrix == "Automatic"):
        arguments["cell_matrix"] = CellMatrix
    result = _decode_numeric(properties("surface_green_function", **arguments))
    outputs = {
        "GreenFunction": "GreenFunction",
        "SpectralWeight": "SpectralWeight",
        "Data": None,
    }
    if Output not in outputs:
        raise PropertiesError("InvalidPropertiesOption", f"unknown Output {Output}")
    return result if outputs[Output] is None else result[outputs[Output]]


def surfaceSpectralFunction(*args: Any, **kwargs: Any) -> float:
    """Return only ``-Im Tr G_surface / pi``."""

    kwargs["Output"] = "SpectralWeight"
    return surfaceGreenFunction(*args, **kwargs)


def pointChernNumber(
    hamiltonian: Any,
    occupied: int,
    point: Sequence[Any],
    radius: Any,
    *,
    parameters: Optional[Mapping[str, Any]] = None,
    SurfaceSubdivisions: int = 8,
    HermitianTolerance: float = 1.0e-10,
    SurfaceGapTolerance: float = 1.0e-8,
    CenterGapTolerance: float = 1.0e-6,
    OverlapTolerance: float = 1.0e-10,
    IntegerTolerance: float = 5.0e-3,
    RequireGaplessCenter: bool = True,
    Output: str = "ChernNumber",
) -> Any:
    """Compute the stable oriented-cube occupied-subspace Chern charge."""

    if isinstance(SurfaceSubdivisions, bool) or not isinstance(SurfaceSubdivisions, int):
        raise TypeError("SurfaceSubdivisions must be an integer")
    if not isinstance(RequireGaplessCenter, bool):
        raise TypeError("RequireGaplessCenter must be bool")
    numeric_point = _real_vector(point, "point")
    if len(numeric_point) != 3:
        raise TypeError("point must be a three-vector")
    options = {
        "surface_subdivisions": SurfaceSubdivisions,
        "hermitian_tolerance": _finite_real(HermitianTolerance, "HermitianTolerance"),
        "surface_gap_tolerance": _finite_real(SurfaceGapTolerance, "SurfaceGapTolerance"),
        "center_gap_tolerance": _finite_real(CenterGapTolerance, "CenterGapTolerance"),
        "overlap_tolerance": _finite_real(OverlapTolerance, "OverlapTolerance"),
        "integer_tolerance": _finite_real(IntegerTolerance, "IntegerTolerance"),
        "require_gapless_center": RequireGaplessCenter,
    }
    numeric_radius = _finite_real(radius, "radius")
    mesh = properties(
        "point_chern_mesh",
        point=numeric_point,
        radius=numeric_radius,
        **options,
    )
    evaluator = _evaluator(hamiltonian, parameters)
    result = _decode_numeric(
        properties(
            "point_chern_number",
            point=numeric_point,
            radius=numeric_radius,
            occupied=occupied,
            points=mesh["Points"],
            triangles=mesh["Triangles"],
            center_matrix=_complex_matrix(evaluator(numeric_point), "Hamiltonian at center"),
            matrices=_sample(evaluator, mesh["Points"]),
            **options,
        )
    )
    if Output == "ChernNumber":
        return result["ChernNumber"]
    if Output == "Data":
        return result
    raise PropertiesError("InvalidPropertiesOption", f"unknown Output {Output}")


def findGaplessPoints(
    hamiltonian: Any,
    occupied: int,
    *,
    parameters: Optional[Mapping[str, Any]] = None,
    BrillouinZone: Sequence[Sequence[Any]] = (
        (-math.pi, math.pi),
        (-math.pi, math.pi),
        (-math.pi, math.pi),
    ),
    GridSize: int | Sequence[int] = 15,
    CandidateCount: int = 32,
    GapTolerance: float = 1.0e-7,
    MergeTolerance: float = 1.0e-4,
    HermitianTolerance: float = 1.0e-10,
    MaxIterations: int = 500,
    RefinementMethod: str = "PrincipalAxis",
    Output: str = "Points",
) -> Any:
    """Search periodic grid minima and refine their gap squared in Rust."""

    zone = [_real_vector(interval, "BrillouinZone") for interval in BrillouinZone]
    if not 1 <= len(zone) <= 3 or any(len(interval) != 2 for interval in zone):
        raise PropertiesError(
            "InvalidBrillouinZone",
            "BrillouinZone must contain one to three minimum/maximum intervals",
        )
    if isinstance(GridSize, bool):
        raise TypeError("GridSize must be an integer or integer sequence")
    if isinstance(GridSize, int):
        grid_size = [GridSize] * len(zone)
    elif isinstance(GridSize, Sequence) and not isinstance(GridSize, (str, bytes)):
        grid_size = list(GridSize)
    else:
        raise TypeError("GridSize must be an integer or integer sequence")
    for name, value in (
        ("CandidateCount", CandidateCount),
        ("MaxIterations", MaxIterations),
    ):
        if isinstance(value, bool) or not isinstance(value, int):
            raise TypeError(f"{name} must be an integer")
    evaluator = _evaluator(hamiltonian, parameters)
    request = {
        "occupied": occupied,
        "brillouin_zone": zone,
        "grid_size": grid_size,
        "candidate_count": CandidateCount,
        "gap_tolerance": _finite_real(GapTolerance, "GapTolerance"),
        "merge_tolerance": _finite_real(MergeTolerance, "MergeTolerance"),
        "hermitian_tolerance": _finite_real(
            HermitianTolerance, "HermitianTolerance"
        ),
        "max_iterations": MaxIterations,
        "refinement_method": RefinementMethod,
    }
    result = _decode_numeric(
        json.loads(
            gapless_points_json(
                lambda point: evaluator(point),
                json.dumps(
                    request,
                    ensure_ascii=True,
                    allow_nan=False,
                    separators=(",", ":"),
                ),
            )
        )
    )
    if Output == "Points":
        return result["Points"]
    if Output == "Data":
        return result
    raise PropertiesError("InvalidPropertiesOption", f"unknown Output {Output}")


def wilson_loop(*args: Any, **kwargs: Any) -> Any:
    return wilsonLoop(*args, **kwargs)


def berry_phase(*args: Any, **kwargs: Any) -> Any:
    return berryPhase(*args, **kwargs)


def berry_curvature(*args: Any, **kwargs: Any) -> Any:
    return berryCurvature(*args, **kwargs)


def transform_hoppings(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return transformHoppings(*args, **kwargs)


def hopping_data(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return hoppingData(*args, **kwargs)


def build_bloch_hamiltonian(*args: Any, **kwargs: Any) -> Any:
    return buildBlochHamiltonian(*args, **kwargs)


def build_real_space_hamiltonian(*args: Any, **kwargs: Any) -> Any:
    return buildRealSpaceHamiltonian(*args, **kwargs)


def build_slab_hamiltonian(*args: Any, **kwargs: Any) -> Any:
    return buildSlabHamiltonian(*args, **kwargs)


def surface_green_function(*args: Any, **kwargs: Any) -> Any:
    return surfaceGreenFunction(*args, **kwargs)


def surface_spectral_function(*args: Any, **kwargs: Any) -> float:
    return surfaceSpectralFunction(*args, **kwargs)


def point_chern_number(*args: Any, **kwargs: Any) -> Any:
    return pointChernNumber(*args, **kwargs)


def find_gapless_points(*args: Any, **kwargs: Any) -> Any:
    return findGaplessPoints(*args, **kwargs)
