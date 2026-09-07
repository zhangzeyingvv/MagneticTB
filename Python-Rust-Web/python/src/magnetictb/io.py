"""Stable 2.0.10 Wannier90 I/O with parsing and formatting in Rust."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from fractions import Fraction
from math import lcm
from numbers import Number
from pathlib import Path
from typing import Any

from .api import ExactError, PropertiesError, cyclotomic_context, geometry, properties
from .linear_algebra import CyclotomicElement, QuadraticNumber, _encode_evaluation_scalar
from .model import current_model
from .tight_binding import Hamiltonian, HamiltonianExpression
from .properties import _decode_numeric, _hopping_mapping


def _rational(value: Fraction) -> dict[str, str]:
    value = Fraction(value)
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def _element(value: CyclotomicElement) -> dict[str, Any]:
    """Encode an element in the matrix endpoint's selected context."""

    return {"coefficients": [_rational(item) for item in value.coefficients]}


def _encoded_element(value: CyclotomicElement) -> dict[str, Any]:
    """Encode an element as a context-independent exact expression for Rust."""

    terms: list[dict[str, Any]] = []
    for power, coefficient in enumerate(value.coefficients):
        if coefficient == 0:
            continue
        rational = _encode_evaluation_scalar(coefficient, "cyclotomic coefficient")
        if power == 0:
            terms.append(rational)
            continue
        root = {
            "kind": "root_of_unity",
            "order": value.conductor,
            "power": power,
        }
        if coefficient == 1:
            terms.append(root)
        else:
            terms.append(
                {
                    "kind": "call",
                    "head": "System`Times",
                    "arguments": [rational, root],
                }
            )
    if not terms:
        return {"kind": "integer", "value": "0"}
    if len(terms) == 1:
        return terms[0]
    return {"kind": "call", "head": "System`Plus", "arguments": terms}


def _cyclotomic_conductors(value: Any) -> set[int]:
    if isinstance(value, CyclotomicElement):
        return {value.conductor}
    if isinstance(value, Mapping):
        result: set[int] = set()
        for item in value.values():
            result.update(_cyclotomic_conductors(item))
        return result
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        result = set()
        for item in value:
            result.update(_cyclotomic_conductors(item))
        return result
    return set()


def _quadratic(value: Any) -> dict[str, Any]:
    if isinstance(value, QuadraticNumber):
        coefficients = value.coefficients
    elif isinstance(value, (int, Fraction)) and not isinstance(value, bool):
        coefficients = (Fraction(value), Fraction(0), Fraction(0), Fraction(0))
    else:
        raise PropertiesError(
            "InvalidHamiltonianTerm",
            "Wannier90 extraction requires a numerical quadratic displacement",
        )
    return {
        "basis": ["1", "sqrt(2)", "sqrt(3)", "sqrt(6)"],
        "coefficients": [_rational(item) for item in coefficients],
    }


def _matrix_rows(value: Any) -> list[list[Any]]:
    if isinstance(value, Hamiltonian):
        return value.to_list()
    if (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and value
    ):
        if any(
            not isinstance(row, Sequence)
            or isinstance(row, (str, bytes))
            or len(row) != len(value)
            for row in value
        ):
            raise PropertiesError(
                "InvalidHamiltonian", "the hand-written Hamiltonian must be square"
            )
        rows = [list(row) for row in value]
        return rows
    raise PropertiesError(
        "InvalidHamiltonian", "the hand-written Hamiltonian must be nonempty and square"
    )


def _centers(setting: Any, dimension: int) -> list[list[Any]]:
    if setting is None:
        return [[0, 0, 0] for _ in range(dimension)]
    if isinstance(setting, str) and setting == "Automatic":
        try:
            model = current_model()
        except Exception:
            return [[0, 0, 0] for _ in range(dimension)]
        if len(model.orbital_table) != dimension:
            return [[0, 0, 0] for _ in range(dimension)]
        return [list(record.fractional_position) for record in model.orbital_table]
    if (
        isinstance(setting, Sequence)
        and not isinstance(setting, (str, bytes))
        and len(setting) == dimension
        and all(
            isinstance(center, Sequence)
            and not isinstance(center, (str, bytes))
            and len(center) == 3
            for center in setting
        )
    ):
        return [list(center) for center in setting]
    raise PropertiesError(
        "InvalidWannierCenters",
        "wcc must be None, Automatic, or one three-vector per Hamiltonian row",
    )


def _encode_symbolic_matrix(
    value: Any,
    rules: Mapping[str, Any] | Sequence[tuple[str, Any]],
) -> tuple[list[list[Any]], dict[str, Any], dict[str, Any]]:
    rows = _matrix_rows(value)
    if isinstance(rules, Mapping):
        parameters = dict(rules)
    elif isinstance(rules, Sequence) and not isinstance(rules, (str, bytes)):
        parameters = dict(rules)
    else:
        raise TypeError("rules must be a mapping or sequence of (name, value) pairs")
    conductors = {
        contribution.coefficient.conductor
        for row in rows
        for cell in row
        if isinstance(cell, HamiltonianExpression)
        for contribution in cell.contributions
    }
    if len(conductors) > 1:
        raise PropertiesError(
            "IncompatibleHamiltonianField",
            "all Hamiltonian coefficients must belong to one cyclotomic field",
        )
    context = cyclotomic_context(next(iter(conductors), 1))
    one = {
        "coefficients": [
            _rational(Fraction(1 if index == 0 else 0))
            for index in range(int(context["degree"]))
        ]
    }

    def encoded_scalar(item: Any, name: str) -> Any:
        return (
            _element(item)
            if isinstance(item, CyclotomicElement)
            else _encode_evaluation_scalar(item, name)
        )

    encoded_parameters = {
        str(name): encoded_scalar(item, f"parameter {name}")
        for name, item in parameters.items()
    }
    encoded_rows: list[list[list[dict[str, Any]]]] = []
    for row_index, row in enumerate(rows):
        encoded_row = []
        for column_index, cell in enumerate(row):
            if isinstance(cell, HamiltonianExpression):
                encoded_row.append(
                    [
                        {
                            "parameter": contribution.parameter_name,
                            "coefficient": _element(contribution.coefficient),
                            "displacement": [
                                _quadratic(item) for item in contribution.displacement
                            ],
                        }
                        for contribution in cell.contributions
                    ]
                )
            elif isinstance(cell, Number) and not isinstance(cell, bool):
                name = f"__constant_{row_index}_{column_index}"
                encoded_parameters[name] = _encode_evaluation_scalar(cell, name)
                encoded_row.append(
                    [
                        {
                            "parameter": name,
                            "coefficient": one,
                            "displacement": [_quadratic(0)] * 3,
                        }
                    ]
                    if complex(cell) != 0
                    else []
                )
            else:
                raise PropertiesError(
                    "InvalidHamiltonianTerm",
                    "matrix entries must be HamiltonianExpression or numeric scalars",
                )
        encoded_rows.append(encoded_row)
    return encoded_rows, context, encoded_parameters


def _encode_hopping_symbolic_matrix(
    value: Any,
    rules: Mapping[str, Any] | Sequence[tuple[str, Any]],
    centers: Sequence[Sequence[Any]],
) -> tuple[list[list[Any]], dict[str, Any], dict[str, Any]]:
    """Encode hopping input without assuming all exact values share one context."""

    rows = _matrix_rows(value)
    if isinstance(rules, Mapping):
        parameters = dict(rules)
    elif isinstance(rules, Sequence) and not isinstance(rules, (str, bytes)):
        parameters = dict(rules)
    else:
        raise TypeError("rules must be a mapping or sequence of (name, value) pairs")

    conductors = {
        contribution.coefficient.conductor
        for row in rows
        for cell in row
        if isinstance(cell, HamiltonianExpression)
        for contribution in cell.contributions
    }
    conductors.update(_cyclotomic_conductors(parameters))
    conductors.update(_cyclotomic_conductors(centers))
    conductor = 1
    for required_conductor in conductors:
        conductor = lcm(conductor, required_conductor)
    context = cyclotomic_context(conductor)

    def encoded_scalar(item: Any, name: str) -> Any:
        return (
            _encoded_element(item)
            if isinstance(item, CyclotomicElement)
            else _encode_evaluation_scalar(item, name)
        )

    encoded_parameters = {
        str(name): encoded_scalar(item, f"parameter {name}")
        for name, item in parameters.items()
    }
    encoded_rows: list[list[list[dict[str, Any]]]] = []
    for row_index, row in enumerate(rows):
        encoded_row = []
        for column_index, cell in enumerate(row):
            if isinstance(cell, HamiltonianExpression):
                encoded_row.append(
                    [
                        {
                            "parameter": contribution.parameter_name,
                            "coefficient": _encoded_element(contribution.coefficient),
                            "displacement": [
                                _quadratic(item) for item in contribution.displacement
                            ],
                        }
                        for contribution in cell.contributions
                    ]
                )
            elif isinstance(cell, Number) and not isinstance(cell, bool):
                name = f"__constant_{row_index}_{column_index}"
                encoded_parameters[name] = _encode_evaluation_scalar(cell, name)
                encoded_row.append(
                    [
                        {
                            "parameter": name,
                            "coefficient": {"kind": "integer", "value": "1"},
                            "displacement": [_quadratic(0)] * 3,
                        }
                    ]
                    if complex(cell) != 0
                    else []
                )
            else:
                raise PropertiesError(
                    "InvalidHamiltonianTerm",
                    "matrix entries must be HamiltonianExpression or numeric scalars",
                )
        encoded_rows.append(encoded_row)
    return encoded_rows, context, encoded_parameters


def _symbolic_hopping_data(
    value: Any,
    rules: Mapping[str, Any] | Sequence[tuple[str, Any]],
    centers: Any,
    tolerance: float,
) -> dict[str, Any]:
    rows = _matrix_rows(value)
    resolved_centers = _centers(centers, len(rows))
    encoded_rows, context, encoded_parameters = _encode_hopping_symbolic_matrix(
        value, rules, resolved_centers
    )

    def encoded_scalar(item: Any) -> Any:
        return (
            _encoded_element(item)
            if isinstance(item, CyclotomicElement)
            else _encode_evaluation_scalar(item, "Wannier center")
        )

    try:
        return _decode_numeric(
            geometry(
                "hopping_data_from_symbolic_matrix",
                context=context,
                matrix=encoded_rows,
                parameters=encoded_parameters,
                wannier_centers=[
                    [encoded_scalar(item) for item in center]
                    for center in resolved_centers
                ],
                translation_tolerance=tolerance,
            )
        )
    except ExactError as error:
        raise PropertiesError(error.tag, error.detail) from None


def _emit_hr(text: str, target: Any) -> str:
    if target is None:
        print(text, end="")
        return text
    if not isinstance(target, (str, Path)):
        raise PropertiesError("InvalidWannier90Option", "hrExport must be a path or None")
    path = Path(target)
    if path.is_dir():
        path = path / "wannier90_hr.dat"
    if not path.parent.is_dir():
        raise PropertiesError(
            "InvalidExportTarget", "the target directory must already exist"
        )
    try:
        path.write_text(text, encoding="utf-8")
    except OSError as error:
        raise PropertiesError("ExportFailed", str(error)) from None
    return str(path)


def hop(
    selection: Any,
    rules: Mapping[str, Any] | Sequence[tuple[str, Any]],
    *,
    Hermitian: bool = True,
    KernelMethod: str = "Iterative",
    ValidationLevel: str = "Basic",
    hrExport: Any = None,
    RealDigits: int = 12,
    wcc: Any = "Automatic",
    TranslationTolerance: float = 1.0e-9,
) -> str:
    """Match stable ``hop`` and emit canonical wannier90_hr.dat text."""

    if isinstance(selection, Hamiltonian):
        data = _symbolic_hopping_data(selection, rules, wcc, TranslationTolerance)
    elif isinstance(selection, int) and not isinstance(selection, bool):
        from .properties import hoppingData

        data = hoppingData(
            selection,
            rules,
            Hermitian=Hermitian,
            KernelMethod=KernelMethod,
            ValidationLevel=ValidationLevel,
        )
    elif (
        isinstance(selection, Sequence)
        and not isinstance(selection, (str, bytes))
        and selection
        and all(isinstance(item, int) and not isinstance(item, bool) for item in selection)
    ):
        from .properties import hoppingData

        data = hoppingData(
            selection,
            rules,
            Hermitian=Hermitian,
            KernelMethod=KernelMethod,
            ValidationLevel=ValidationLevel,
        )
    else:
        data = _symbolic_hopping_data(selection, rules, wcc, TranslationTolerance)
    result = properties(
        "format_wannier90_hr", data=_hopping_mapping(data), real_digits=RealDigits
    )
    return _emit_hr(str(result["Text"]), hrExport)


def readHR(file: str | Path, *, prec: float = 1.0e-10, ncell: Any = "All") -> dict[str, Any]:
    """Strictly parse one complete wannier90_hr.dat file in Rust."""

    if not isinstance(file, (str, Path)):
        raise PropertiesError("UnreadableHRFile", "file must be a path")
    path = Path(file)
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as error:
        raise PropertiesError("UnreadableHRFile", str(error)) from None
    if ncell == "All" or ncell is None:
        cutoff = None
    elif isinstance(ncell, int) and not isinstance(ncell, bool) and ncell >= 0:
        cutoff = ncell
    else:
        raise PropertiesError(
            "InvalidWannier90Option", "ncell must be All or a nonnegative integer"
        )
    return _decode_numeric(
        properties(
            "parse_wannier90_hr",
            text=text,
            source=str(path.resolve()),
            precision_cutoff=prec,
            cell_cutoff=cutoff,
        )
    )


def read_hr(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return readHR(*args, **kwargs)
