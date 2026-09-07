"""Linear Algebra boundary for the Rust-backed MagneticTB Python API.

This module performs only immutable input/output adaptation; mathematical
algorithms remain in the Rust core.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from contextvars import ContextVar
from dataclasses import dataclass, field
from fractions import Fraction
from functools import lru_cache
from math import isfinite
from threading import RLock
from types import MappingProxyType
from typing import Any, Dict, Literal, Optional, Tuple, Union

ExactAtom = Union[int, float, complex, Fraction, "ExactExpr"]


MatrixInput = Sequence[Sequence[ExactAtom]]


def _freeze_json(value: Any) -> Any:
    if isinstance(value, MappingProxyType):
        return value
    if isinstance(value, Mapping):
        return MappingProxyType({key: _freeze_json(item) for key, item in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_freeze_json(item) for item in value)
    return value


def _thaw_json(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _thaw_json(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_thaw_json(item) for item in value]
    return value


def _hashable_json(value: Any) -> Any:
    if isinstance(value, Mapping):
        return tuple((key, _hashable_json(item)) for key, item in value.items())
    if isinstance(value, tuple):
        return tuple(_hashable_json(item) for item in value)
    return value


def _freeze_matrix_input(value: MatrixInput) -> Tuple[Tuple[ExactAtom, ...], ...]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or not value:
        raise TypeError("matrix must be a nonempty sequence of rows")
    rows = tuple(tuple(row) for row in value)
    if not rows[0] or any(len(row) != len(rows[0]) for row in rows):
        raise ValueError("matrix rows must be nonempty and rectangular")
    return rows


@dataclass(frozen=True)
class ExactExpr:
    """An immutable exact-expression syntax node; evaluation remains in Rust."""

    _value: Mapping[str, Any]

    def __post_init__(self) -> None:
        object.__setattr__(self, "_value", _freeze_json(self._value))

    def __hash__(self) -> int:
        return hash(_hashable_json(self._value))

    def __str__(self) -> str:
        return _format_exact_payload(self._value)

    def __repr__(self) -> str:
        return str(self)

    def __add__(self, other: ExactAtom) -> "ExactExpr":
        return _call("System`Plus", self, other)

    def __radd__(self, other: ExactAtom) -> "ExactExpr":
        return _call("System`Plus", other, self)

    def __sub__(self, other: ExactAtom) -> "ExactExpr":
        return _call("System`Plus", self, -_coerce_expr(other))

    def __rsub__(self, other: ExactAtom) -> "ExactExpr":
        return _call("System`Plus", other, -self)

    def __mul__(self, other: ExactAtom) -> "ExactExpr":
        return _call("System`Times", self, other)

    def __rmul__(self, other: ExactAtom) -> "ExactExpr":
        return _call("System`Times", other, self)

    def __truediv__(self, other: ExactAtom) -> "ExactExpr":
        if isinstance(other, (int, Fraction)) and not isinstance(other, bool):
            rational = Fraction(other)
            if rational == 0:
                raise ZeroDivisionError("exact expression division by zero")
            return _call("System`Times", self, 1 / rational)
        return _call("System`Times", self, _coerce_expr(other) ** -1)

    def __rtruediv__(self, other: ExactAtom) -> "ExactExpr":
        return _call("System`Times", other, self**-1)

    def __pow__(self, exponent: int | Fraction) -> "ExactExpr":
        if isinstance(exponent, bool) or not isinstance(exponent, (int, Fraction)):
            raise TypeError("exact exponent must be int or Fraction")
        return _call("System`Power", self, exponent)

    def __neg__(self) -> "ExactExpr":
        return _call("System`Times", -1, self)


def _format_exact_payload(value: Any) -> str:
    """Render the public exact AST without exposing its canonical JSON tags."""

    if not isinstance(value, Mapping):
        return str(value)
    kind = value.get("kind")
    if kind == "integer":
        return str(value["value"])
    if kind == "rational":
        return f"{value['numerator']}/{value['denominator']}"
    if kind == "symbol":
        return str(value["name"]).rsplit("`", 1)[-1]
    if kind == "root_of_unity":
        return _root_of_unity_text(int(value["order"]), int(value["power"]))
    if kind != "call":
        return f"ExactExpr({kind or 'unknown'})"

    head = str(value["head"]).rsplit("`", 1)[-1]
    arguments = tuple(_format_exact_payload(item) for item in value["arguments"])
    if head == "Plus":
        return "(" + " + ".join(arguments).replace("+ -", "- ") + ")"
    if head == "Times":
        return "(" + "*".join(arguments) + ")"
    if head == "Power" and len(arguments) == 2:
        return f"({arguments[0]})^({arguments[1]})"
    if head == "Complex" and len(arguments) == 2:
        return f"({arguments[0]} + I*{arguments[1]})"
    return f"{head}({', '.join(arguments)})"


def _root_of_unity_text(order: int, power: int) -> str:
    """Render an exact root of unity as a reduced exponential expression."""

    normalized_power = power % order
    if normalized_power == 0:
        return "1"
    phase = Fraction(2 * normalized_power, order)
    numerator = phase.numerator
    coefficient = "pi*i" if numerator == 1 else f"{numerator}*pi*i"
    if phase.denominator != 1:
        coefficient = f"{coefficient}/{phase.denominator}"
    return f"e^({coefficient})"


def _coerce_expr(value: ExactAtom) -> ExactExpr:
    if isinstance(value, ExactExpr):
        return value
    return ExactExpr(_freeze_json(_encode_exact(value)))


def _call(head: str, *arguments: ExactAtom) -> ExactExpr:
    return ExactExpr(
        _freeze_json(
            {
                "kind": "call",
                "head": head,
                "arguments": [_encode_exact(argument) for argument in arguments],
            }
        )
    )


def symbol(name: str) -> ExactExpr:
    """Create an exact lattice or continuous parameter symbol."""

    if not isinstance(name, str) or not name:
        raise TypeError("symbol name must be a nonempty string")
    qualified = name if "`" in name else f"MagneticTB`{name}"
    return ExactExpr(_freeze_json({"kind": "symbol", "name": qualified}))


def sqrt(radicand: int) -> ExactExpr:
    """Create an unevaluated exact square root for Rust compilation."""

    if isinstance(radicand, bool) or not isinstance(radicand, int) or radicand <= 0:
        raise TypeError("sqrt radicand must be a positive integer")
    return _call("System`Power", radicand, Fraction(1, 2))


def root_of_unity(order: int, power: int = 1) -> ExactExpr:
    if (
        isinstance(order, bool)
        or not isinstance(order, int)
        or order <= 0
        or isinstance(power, bool)
        or not isinstance(power, int)
    ):
        raise TypeError("root_of_unity requires a positive integer order and integer power")
    return ExactExpr(
        _freeze_json({"kind": "root_of_unity", "order": order, "power": power})
    )


def sin(argument: ExactAtom) -> ExactExpr:
    return _call("System`Sin", argument)


def cos(argument: ExactAtom) -> ExactExpr:
    return _call("System`Cos", argument)


def exact_complex(real: ExactAtom, imaginary: ExactAtom) -> ExactExpr:
    """Create an exact complex scalar syntax node for a Rust expression."""

    return _call("System`Complex", real, imaginary)


def exp(argument: ExactAtom) -> ExactExpr:
    """Create an exact exponential syntax node for continuous representations."""

    return _call("System`Power", symbol("System`E"), argument)


@dataclass(frozen=True)
class ExactField:
    """An explicit exact cyclotomic field; no automatic enlargement occurs."""

    conductor: int
    cyclotomic_polynomial: Tuple[int, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "cyclotomic_polynomial", tuple(self.cyclotomic_polynomial))
        if (
            isinstance(self.conductor, bool)
            or not isinstance(self.conductor, int)
            or self.conductor < 1
        ):
            raise TypeError("conductor must be a positive integer")
        if not self.cyclotomic_polynomial:
            raise ValueError("cyclotomic_polynomial must be nonempty")

    @classmethod
    def rational(cls) -> "ExactField":
        return cls(1, (-1, 1))

    @classmethod
    def gaussian(cls) -> "ExactField":
        return cls(4, (1, 0, 1))

    @classmethod
    def crystallographic(cls) -> "ExactField":
        return cls(24, (1, 0, 0, 0, -1, 0, 0, 0, 1))

    def _context(self) -> Dict[str, Any]:
        return {
            "conductor": self.conductor,
            "degree": len(self.cyclotomic_polynomial) - 1,
            "cyclotomic_polynomial": [
                {"numerator": str(coefficient), "denominator": "1"}
                for coefficient in self.cyclotomic_polynomial
            ],
        }


@dataclass(frozen=True)
class CyclotomicElement:
    coefficients: Tuple[Fraction, ...]
    conductor: int

    def __str__(self) -> str:
        if all(value == 0 for value in self.coefficients[1:]):
            return str(self.coefficients[0])
        terms = []
        for power, coefficient in enumerate(self.coefficients):
            if coefficient == 0:
                continue
            if power == 0:
                atom = "1"
            elif self.conductor % 4 == 0 and power == self.conductor // 4:
                atom = "I"
            else:
                atom = _root_of_unity_text(self.conductor, power)
            if atom == "1":
                term = str(coefficient)
            elif coefficient == 1:
                term = atom
            elif coefficient == -1:
                term = f"-{atom}"
            else:
                term = f"{coefficient}*{atom}"
            terms.append(term)
        return " + ".join(terms).replace("+ -", "- ")


@dataclass(frozen=True)
class ExactMatrix:
    rows: Tuple[Tuple[CyclotomicElement, ...], ...]

    @property
    def shape(self) -> Tuple[int, int]:
        return (len(self.rows), len(self.rows[0]) if self.rows else 0)

    def to_list(self) -> list[list[CyclotomicElement]]:
        return [list(row) for row in self.rows]


@dataclass(frozen=True)
class QuadraticNumber:
    coefficients: Tuple[Fraction, Fraction, Fraction, Fraction]

    def __str__(self) -> str:
        labels = ("", "sqrt(2)", "sqrt(3)", "sqrt(6)")
        terms = []
        for coefficient, label in zip(self.coefficients, labels):
            if coefficient == 0:
                continue
            if not label:
                terms.append(str(coefficient))
            elif coefficient == 1:
                terms.append(label)
            elif coefficient == -1:
                terms.append(f"-{label}")
            else:
                terms.append(f"{coefficient}*{label}")
        return " + ".join(terms).replace("+ -", "- ") if terms else "0"


PhaseCoordinate = Union[QuadraticNumber, int, Fraction, ExactExpr]


def _encode_lattpar(value: Any) -> Dict[str, Any]:
    items: Iterable[Tuple[Any, Any]]
    if isinstance(value, Mapping):
        items = value.items()
    else:
        items = value
    rules = []
    for key, item in items:
        parameter = key if isinstance(key, ExactExpr) else symbol(key)
        rules.append(
            {
                "kind": "call",
                "head": "System`Rule",
                "arguments": [_thaw_json(parameter._value), _encode_exact(item)],
            }
        )
    return _tagged_list(rules)


def _encode_exact(value: Any) -> Dict[str, Any]:
    if isinstance(value, bool):
        raise TypeError("bool is not an exact scalar")
    if isinstance(value, ExactExpr):
        return _thaw_json(value._value)
    if isinstance(value, int):
        return {"kind": "integer", "value": str(value)}
    if isinstance(value, Fraction):
        if value.denominator == 1:
            return {"kind": "integer", "value": str(value.numerator)}
        return {
            "kind": "rational",
            "numerator": str(value.numerator),
            "denominator": str(value.denominator),
        }
    if isinstance(value, float):
        return {"kind": "machine_real", "decimal": repr(value)}
    if isinstance(value, complex):
        return {
            "kind": "machine_complex",
            "real_decimal": repr(value.real),
            "imaginary_decimal": repr(value.imag),
        }
    raise TypeError(f"unsupported exact scalar {type(value).__name__}")


def _encode_evaluation_scalar(value: ExactAtom, name: str) -> Any:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be numeric, not bool")
    if isinstance(value, complex):
        if not isfinite(value.real) or not isfinite(value.imag):
            raise ValueError(f"{name} must be finite")
        return {
            "kind": "approx_complex",
            "real": value.real,
            "imaginary": value.imag,
        }
    if isinstance(value, float):
        if not isfinite(value):
            raise ValueError(f"{name} must be finite")
        return value
    if isinstance(value, (int, Fraction, ExactExpr)):
        return _encode_exact(value)
    raise TypeError(f"{name} must be int, Fraction, ExactExpr, float, or complex")


def _encode_tree(value: Any) -> Any:
    if isinstance(value, (int, float, complex, Fraction, ExactExpr)) and not isinstance(value, bool):
        return _encode_exact(value)
    if isinstance(value, Mapping):
        return _ordered([(str(key), _encode_tree(item)) for key, item in value.items()])
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return _tagged_list([_encode_tree(item) for item in value])
    return value


def _matrix(rows: MatrixInput) -> Dict[str, Any]:
    if not isinstance(rows, Sequence) or isinstance(rows, (str, bytes)) or not rows:
        raise TypeError("matrix must be a nonempty sequence of rows")
    columns = len(rows[0])
    if columns == 0 or any(not isinstance(row, Sequence) or len(row) != columns for row in rows):
        raise ValueError("matrix rows must be nonempty and rectangular")
    return {
        "kind": "matrix",
        "rows": len(rows),
        "columns": columns,
        "entries": [_encode_exact(item) for row in rows for item in row],
    }


def _tagged_list(items: Sequence[Any]) -> Dict[str, Any]:
    return {"kind": "list", "items": list(items)}


def _ordered(items: Iterable[Tuple[str, Any]]) -> Dict[str, Any]:
    return {
        "kind": "ordered_association",
        "entries": [{"key": key, "value": value} for key, value in items],
    }


def _strict_bool(value: Any, name: str) -> bool:
    if not isinstance(value, bool):
        raise TypeError(f"{name} must be bool")
    return value


def _decode_element(value: Mapping[str, Any], field_spec: ExactField) -> CyclotomicElement:
    return CyclotomicElement(
        tuple(
            Fraction(int(coefficient["numerator"]), int(coefficient["denominator"]))
            for coefficient in value["coefficients"]
        ),
        field_spec.conductor,
    )


def _decode_nested_elements(value: Any, field_spec: ExactField) -> Any:
    if isinstance(value, list):
        return tuple(_decode_nested_elements(item, field_spec) for item in value)
    if isinstance(value, Mapping) and "coefficients" in value:
        return _decode_element(value, field_spec)
    return value


def _decode_matrix(value: Mapping[str, Any], field_spec: ExactField) -> ExactMatrix:
    columns = int(value["columns"])
    entries = [_decode_element(item, field_spec) for item in value["entries"]]
    return ExactMatrix(
        tuple(
            tuple(entries[offset : offset + columns])
            for offset in range(0, len(entries), columns)
        )
    )


def _decode_quadratic(value: Mapping[str, Any]) -> QuadraticNumber:
    coefficients = tuple(
        Fraction(int(item["numerator"]), int(item["denominator"]))
        for item in value["coefficients"]
    )
    if len(coefficients) != 4:
        raise RuntimeError("Rust returned a malformed quadratic coefficient vector")
    return QuadraticNumber(coefficients)  # type: ignore[arg-type]


def _decode_phase_coordinate(value: Any) -> PhaseCoordinate:
    if isinstance(value, Mapping) and "basis" in value:
        return _decode_quadratic(value)
    if not isinstance(value, Mapping):
        if isinstance(value, bool) or not isinstance(value, int):
            raise RuntimeError("Rust returned a malformed exact phase coordinate")
        return value
    kind = value.get("kind")
    if kind == "integer":
        return int(value["value"])
    if kind == "rational":
        return Fraction(int(value["numerator"]), int(value["denominator"]))
    return ExactExpr(value)


def _encode_cyclotomic_element(value: CyclotomicElement) -> Dict[str, Any]:
    return {
        "coefficients": [
            {
                "numerator": str(coefficient.numerator),
                "denominator": str(coefficient.denominator),
            }
            for coefficient in value.coefficients
        ]
    }


def _encode_phase_coordinate(value: PhaseCoordinate) -> Dict[str, Any]:
    if isinstance(value, QuadraticNumber):
        return {
            "basis": ["1", "sqrt(2)", "sqrt(3)", "sqrt(6)"],
            "coefficients": [
                {
                    "numerator": str(coefficient.numerator),
                    "denominator": str(coefficient.denominator),
                }
                for coefficient in value.coefficients
            ],
        }
    return _encode_exact(value)


def _decode_input_exact(value: Any) -> Any:
    if not isinstance(value, Mapping):
        return value
    kind = value.get("kind")
    if kind == "integer":
        return int(value["value"])
    if kind == "rational":
        return Fraction(int(value["numerator"]), int(value["denominator"]))
    if kind in ("symbol", "root_of_unity", "call"):
        return ExactExpr(value)
    if kind == "list":
        return tuple(_decode_input_exact(item) for item in value["items"])
    if kind == "matrix":
        rows = int(value["rows"])
        columns = int(value["columns"])
        entries = tuple(_decode_input_exact(item) for item in value["entries"])
        if len(entries) != rows * columns:
            raise RuntimeError("Rust returned a malformed exact input matrix")
        return tuple(
            tuple(entries[row * columns : (row + 1) * columns])
            for row in range(rows)
        )
    if kind == "ordered_association":
        return tuple(
            (str(entry["key"]), _decode_input_exact(entry["value"]))
            for entry in value["entries"]
        )
    raise RuntimeError("Rust returned an unsupported exact input encoding")


def _decode_lattice_parameters(value: Mapping[str, Any]) -> Tuple[Tuple[Any, Any], ...]:
    if value.get("kind") != "list":
        raise RuntimeError("Rust returned malformed lattice parameters")
    output = []
    for rule in value["items"]:
        if rule.get("kind") != "call" or rule.get("head") != "System`Rule":
            raise RuntimeError("Rust returned a malformed lattice-parameter rule")
        arguments = rule.get("arguments", ())
        if len(arguments) != 2:
            raise RuntimeError("Rust returned a malformed lattice-parameter rule")
        parameter = arguments[0]
        if parameter.get("kind") != "symbol":
            raise RuntimeError("Rust returned a nonsymbol lattice parameter")
        output.append(
            (
                str(parameter["name"]).rsplit("`", 1)[-1],
                _decode_input_exact(arguments[1]),
            )
        )
    return tuple(output)


def _field_from_context(value: Mapping[str, Any]) -> ExactField:
    polynomial = []
    for coefficient in value["cyclotomic_polynomial"]:
        rational = Fraction(
            int(coefficient["numerator"]), int(coefficient["denominator"])
        )
        if rational.denominator != 1:
            raise RuntimeError("Rust returned a non-integral cyclotomic polynomial")
        polynomial.append(rational.numerator)
    return ExactField(int(value["conductor"]), tuple(polynomial))


def _cyclotomic_to_user(value: CyclotomicElement) -> Any:
    if all(coefficient == 0 for coefficient in value.coefficients[1:]):
        return value.coefficients[0]
    return {
        "conductor": value.conductor,
        "coefficients": list(value.coefficients),
    }


def _quadratic_to_user(value: QuadraticNumber) -> Any:
    if all(coefficient == 0 for coefficient in value.coefficients[1:]):
        return value.coefficients[0]
    return {
        "basis": ["1", "sqrt(2)", "sqrt(3)", "sqrt(6)"],
        "coefficients": list(value.coefficients),
    }


def _phase_coordinate_to_user(value: PhaseCoordinate) -> Any:
    return _quadratic_to_user(value) if isinstance(value, QuadraticNumber) else value


def _matrix_to_user(value: ExactMatrix) -> list[list[Any]]:
    return [
        [_cyclotomic_to_user(element) for element in row]
        for row in value.rows
    ]
