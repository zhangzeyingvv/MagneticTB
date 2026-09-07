"""Non-evaluating exact form grammar for the browser init workflow."""

from __future__ import annotations

import ast
import math
import re
from fractions import Fraction
from typing import TYPE_CHECKING, Any, List, Mapping, Optional

from ..modeling import (
    ExactField,
    SpinSpaceOperation,
    SymmetryOperation,
    cos,
    exact_complex,
    root_of_unity,
    sin,
    sqrt,
    symbol,
)
from .codec import decode_tagged
if TYPE_CHECKING:
    from .model_schemas import InitFieldRequest, InitOperationRequest

_INTEGER_PATTERN = re.compile(r"^[+-]?\d+$")


def _parse_scalar(value: Any, field_name: str) -> Any:
    """Parse a small, non-evaluating exact-expression grammar for form fields."""

    if isinstance(value, bool):
        raise TypeError(f"{field_name} must be a scalar, not bool")
    if isinstance(value, (int, float, complex)):
        return value
    if isinstance(value, Mapping):
        return decode_tagged(value)
    if not isinstance(value, str) or not value.strip():
        raise TypeError(f"{field_name} must be a nonempty scalar expression")
    text = value.strip()
    if _INTEGER_PATTERN.fullmatch(text):
        return int(text)
    try:
        expression = ast.parse(text, mode="eval").body
    except SyntaxError as error:
        raise ValueError(f"{field_name} has invalid scalar syntax") from error

    def convert(node: ast.AST) -> Any:
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float, complex)):
            if isinstance(node.value, bool):
                raise ValueError
            return node.value
        if isinstance(node, ast.Name):
            if node.id == "I":
                return exact_complex(0, 1)
            if node.id in {"Pi", "pi"}:
                return symbol("System`Pi")
            return symbol(node.id)
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            operand = convert(node.operand)
            return operand if isinstance(node.op, ast.UAdd) else -operand
        if isinstance(node, ast.BinOp):
            left = convert(node.left)
            right = convert(node.right)
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
            if isinstance(node.op, ast.Div):
                if isinstance(left, (int, Fraction)) and isinstance(right, (int, Fraction)):
                    return Fraction(left) / Fraction(right)
                return left / right
            if isinstance(node.op, ast.Pow):
                if isinstance(right, bool) or not isinstance(right, int):
                    raise ValueError(f"{field_name} exponent must be an integer")
                return left**right
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            name = node.func.id
            arguments = [convert(argument) for argument in node.args]
            if node.keywords:
                raise ValueError(f"{field_name} does not accept keyword calls")
            if name == "sqrt" and len(arguments) == 1 and isinstance(arguments[0], int):
                return sqrt(arguments[0])
            if name == "sin" and len(arguments) == 1:
                return sin(arguments[0])
            if name == "cos" and len(arguments) == 1:
                return cos(arguments[0])
            if name in {"root", "root_of_unity"} and len(arguments) in {1, 2} and all(
                isinstance(item, int) for item in arguments
            ):
                return root_of_unity(arguments[0], arguments[1] if len(arguments) == 2 else 1)
            if name in {"complex", "exact_complex"} and len(arguments) == 2:
                return exact_complex(arguments[0], arguments[1])
        raise ValueError(
            f"{field_name} supports numbers, symbols, + - * / **, sqrt, sin, cos, root, and exact_complex"
        )

    return convert(expression)


def _parse_lattice_scalar(value: Any, field_name: str) -> Any:
    """Preserve the lattice AST for Rust's exact/7-significant-digit policy."""

    parsed = _parse_scalar(value, field_name)
    if isinstance(parsed, complex):
        if parsed.imag != 0.0:
            raise ValueError(f"{field_name} must be real")
        parsed = parsed.real
    if isinstance(parsed, float):
        if not math.isfinite(parsed):
            raise ValueError(f"{field_name} must be finite")
    return parsed


def _parse_vector(values: List[Any], field_name: str) -> tuple[Any, Any, Any]:
    if len(values) != 3:
        raise ValueError(f"{field_name} must contain exactly three values")
    return tuple(
        _parse_scalar(value, f"{field_name}[{index}]")
        for index, value in enumerate(values)
    )  # type: ignore[return-value]


def _parse_matrix(values: List[List[Any]], field_name: str) -> tuple[tuple[Any, ...], ...]:
    if len(values) != 3 or any(len(row) != 3 for row in values):
        raise ValueError(f"{field_name} must be a 3 by 3 matrix")
    return tuple(
        tuple(
            _parse_scalar(value, f"{field_name}[{row_index},{column_index}]")
            for column_index, value in enumerate(row)
        )
        for row_index, row in enumerate(values)
    )


def _parse_representation_matrix(
    values: List[List[Any]], field_name: str
) -> tuple[tuple[Any, ...], ...]:
    if not values or not values[0] or any(len(row) != len(values[0]) for row in values):
        raise ValueError(f"{field_name} must be a nonempty rectangular matrix")
    return tuple(
        tuple(
            _parse_scalar(value, f"{field_name}[{row_index},{column_index}]")
            for column_index, value in enumerate(row)
        )
        for row_index, row in enumerate(values)
    )


def _parse_lattice_matrix(
    values: List[List[Any]], field_name: str
) -> tuple[tuple[Any, ...], ...]:
    if len(values) != 3 or any(len(row) != 3 for row in values):
        raise ValueError(f"{field_name} must be a 3 by 3 matrix")
    return tuple(
        tuple(
            _parse_lattice_scalar(
                value, f"{field_name}[{row_index},{column_index}]"
            )
            for column_index, value in enumerate(row)
        )
        for row_index, row in enumerate(values)
    )


def _field_from_form(value: InitFieldRequest) -> ExactField:
    if value.preset == "crystallographic":
        return ExactField.crystallographic()
    if value.preset == "rational":
        return ExactField.rational()
    if value.preset == "gaussian":
        return ExactField.gaussian()
    if value.conductor is None or value.cyclotomic_polynomial is None:
        raise ValueError("custom exact field requires conductor and cyclotomic_polynomial")
    return ExactField(value.conductor, tuple(value.cyclotomic_polynomial))


def _operation_from_form(value: InitOperationRequest) -> Any:
    rotation = _parse_matrix(value.rotation, f"operation {value.label} rotation")
    translation = _parse_vector(value.translation, f"operation {value.label} translation")
    if value.kind == "spatial":
        if value.spin_rotation is not None or value.continuous_parameter is not None:
            raise ValueError("spatial operation cannot contain spin/continuous fields")
        return SymmetryOperation(value.label, rotation, translation, value.antiunitary)
    if value.spin_rotation is None:
        raise ValueError("spin-space operation requires spin_rotation")
    parameter = (
        symbol(value.continuous_parameter.strip())
        if value.continuous_parameter is not None and value.continuous_parameter.strip()
        else None
    )
    return SpinSpaceOperation(
        value.label,
        rotation,
        translation,
        _parse_matrix(value.spin_rotation, f"operation {value.label} spin rotation"),
        value.antiunitary,
        parameter,
    )
