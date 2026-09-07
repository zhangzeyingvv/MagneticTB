"""Strict JSON boundary for the MagneticTB web API.

This module only reconstructs public immutable Python input objects and makes
public results JSON-safe.  It contains no group, representation, bond, kernel,
Hamiltonian, or band mathematics.
"""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
from fractions import Fraction
from math import isfinite
from typing import Any, Mapping, Sequence

from ..modeling import (
    CyclotomicElement,
    DirectProductRepresentation,
    ExactExpr,
    ExactField,
    ExactMatrix,
    ExplicitPolynomialBasis,
    InducedOrbit,
    InducedRepresentation,
    MSG,
    MagneticGroupSelector,
    QuadraticNumber,
    SpinSpaceOperation,
    SymmetryOperation,
    Wyckoff,
    exact_complex,
    root_of_unity,
    symbol,
)


def guard_json_tree(
    value: Any,
    *,
    maximum_nodes: int = 100_000,
    maximum_depth: int = 32,
    maximum_sequence: int = 8_192,
) -> None:
    """Reject structurally excessive input before it reaches Rust."""

    nodes = 0

    def visit(item: Any, depth: int) -> None:
        nonlocal nodes
        nodes += 1
        if nodes > maximum_nodes:
            raise ValueError("request contains too many JSON values")
        if depth > maximum_depth:
            raise ValueError("request nesting is too deep")
        if isinstance(item, Mapping):
            if len(item) > maximum_sequence:
                raise ValueError("request object contains too many fields")
            for key, nested in item.items():
                if not isinstance(key, str) or len(key) > 256:
                    raise ValueError("request object keys must be short strings")
                visit(nested, depth + 1)
        elif isinstance(item, Sequence) and not isinstance(item, (str, bytes)):
            if len(item) > maximum_sequence:
                raise ValueError("request array is too long")
            for nested in item:
                visit(nested, depth + 1)
        elif isinstance(item, str) and len(item) > 1_000_000:
            raise ValueError("request string is too long")
        elif isinstance(item, float) and not isfinite(item):
            raise ValueError("NaN and Infinity are not valid inputs")

    visit(value, 0)


def _required(value: Mapping[str, Any], name: str) -> Any:
    if name not in value:
        raise ValueError(f"tagged value is missing {name}")
    return value[name]


def decode_tagged(value: Any) -> Any:
    """Decode explicit ``$type`` web records into public Python API objects."""

    if isinstance(value, list):
        return [decode_tagged(item) for item in value]
    if not isinstance(value, Mapping):
        return value

    kind = value.get("$type")
    if kind is None:
        return {str(key): decode_tagged(item) for key, item in value.items()}
    if kind == "fraction":
        return Fraction(int(_required(value, "numerator")), int(_required(value, "denominator")))
    if kind == "symbol":
        return symbol(str(_required(value, "name")))
    if kind == "root_of_unity":
        return root_of_unity(int(_required(value, "order")), int(value.get("power", 1)))
    if kind == "exact_complex":
        return exact_complex(
            decode_tagged(_required(value, "real")),
            decode_tagged(_required(value, "imaginary")),
        )
    if kind == "complex_number":
        result = complex(float(_required(value, "real")), float(_required(value, "imaginary")))
        if not (isfinite(result.real) and isfinite(result.imag)):
            raise ValueError("complex_number must be finite")
        return result
    if kind == "exact_expr":
        payload = _required(value, "value")
        if not isinstance(payload, Mapping):
            raise TypeError("exact_expr value must be a tagged exact AST object")
        return ExactExpr(payload)
    if kind == "exact_field":
        preset = value.get("preset")
        if preset == "rational":
            return ExactField.rational()
        if preset == "gaussian":
            return ExactField.gaussian()
        if preset == "crystallographic":
            return ExactField.crystallographic()
        return ExactField(
            int(_required(value, "conductor")),
            tuple(int(item) for item in _required(value, "cyclotomic_polynomial")),
        )
    if kind == "msg":
        stable_id = value.get("stable_id")
        bns = value.get("bns")
        return MSG(
            stable_id=str(stable_id) if stable_id is not None else None,
            bns=tuple(int(item) for item in bns) if bns is not None else None,
        )
    if kind == "selector":
        return MagneticGroupSelector(
            str(_required(value, "classification")),
            tuple(int(item) for item in _required(value, "source_key")),
        )
    if kind == "wyckoff":
        return Wyckoff(
            letter=str(_required(value, "letter")),
            position=tuple(decode_tagged(value.get("position", [0, 0, 0]))),
            moment=tuple(decode_tagged(value.get("moment", [0, 0, 0]))),
            source_ordinal=value.get("source_ordinal"),
        )
    if kind == "symmetry_operation":
        return SymmetryOperation(
            label=str(_required(value, "label")),
            rotation=decode_tagged(_required(value, "rotation")),
            translation=decode_tagged(value.get("translation", [0, 0, 0])),
            antiunitary=bool(value.get("antiunitary", False)),
        )
    if kind == "spin_space_operation":
        parameter = value.get("continuous_parameter")
        return SpinSpaceOperation(
            label=str(_required(value, "label")),
            space_rotation=decode_tagged(_required(value, "space_rotation")),
            space_translation=decode_tagged(value.get("space_translation", [0, 0, 0])),
            spin_rotation=decode_tagged(_required(value, "spin_rotation")),
            antiunitary=bool(value.get("antiunitary", False)),
            continuous_parameter=(decode_tagged(parameter) if parameter is not None else None),
        )
    if kind == "explicit_polynomial_basis":
        labels = value.get("labels")
        return ExplicitPolynomialBasis(
            decode_tagged(_required(value, "functions")),
            tuple(str(item) for item in labels) if labels is not None else None,
        )
    if kind == "induced_orbit":
        return InducedOrbit(
            reference_site_index=int(_required(value, "reference_site_index")),
            site_symmetry_operation_indices=tuple(
                int(item) for item in _required(value, "site_symmetry_operation_indices")
            ),
            site_symmetry_matrices=decode_tagged(
                _required(value, "site_symmetry_matrices")
            ),
        )
    if kind == "induced_representation":
        continuous = value.get("continuous_site_matrices")
        return InducedRepresentation(
            tuple(decode_tagged(item) for item in _required(value, "orbits")),
            decode_tagged(continuous) if continuous is not None else None,
        )
    if kind == "direct_product_representation":
        continuous = value.get("continuous_site_matrices")
        return DirectProductRepresentation(
            decode_tagged(_required(value, "matrices_by_orbit")),
            decode_tagged(continuous) if continuous is not None else None,
        )
    raise ValueError(f"unsupported tagged web input type {kind!r}")


def _fraction(value: Fraction) -> Any:
    if value.denominator == 1:
        return value.numerator
    return {
        "type": "rational",
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
        "text": str(value),
    }


def jsonable(value: Any) -> Any:
    """Convert public immutable result objects to lossless JSON display data."""

    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not isfinite(value):
            raise ValueError("result contains NaN or Infinity")
        return value
    if isinstance(value, Fraction):
        return _fraction(value)
    if isinstance(value, complex):
        if not (isfinite(value.real) and isfinite(value.imag)):
            raise ValueError("result contains a non-finite complex value")
        return {"type": "complex", "real": value.real, "imaginary": value.imag}
    if isinstance(value, CyclotomicElement):
        return {
            "type": "cyclotomic",
            "conductor": value.conductor,
            "coefficients": [jsonable(item) for item in value.coefficients],
            "text": str(value),
        }
    if isinstance(value, QuadraticNumber):
        return {
            "type": "quadratic",
            "coefficients": [jsonable(item) for item in value.coefficients],
            "text": str(value),
        }
    if isinstance(value, ExactExpr):
        return {"type": "exact_expression", "text": str(value), "ast": jsonable(value._value)}
    if isinstance(value, ExactMatrix):
        return {
            "type": "exact_matrix",
            "shape": list(value.shape),
            "rows": jsonable(value.rows),
        }
    if isinstance(value, Mapping):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return [jsonable(item) for item in value]
    if hasattr(value, "to_dict") and callable(value.to_dict):
        return jsonable(value.to_dict())
    if is_dataclass(value):
        return jsonable(asdict(value))
    if hasattr(value, "_asdict"):
        return jsonable(value._asdict())
    raise TypeError(f"web result contains unsupported type {type(value).__name__}")
