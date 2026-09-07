"""Data-linked Bravais, MSG, and Wyckoff options for browser init."""

from __future__ import annotations

from fractions import Fraction
from typing import Any, Dict, List, Mapping

from ..modeling import ExactExpr, symbol
from .form_parsing import _parse_lattice_scalar, _parse_scalar
from .runtime import (
    _BRAVAIS_OPTION_CACHE,
    _MSG_FAMILY_CACHE,
    _WYCKOFF_OPTION_CACHE,
    _catalog,
)


def _data_exact_to_user(value: Mapping[str, Any]) -> Any:
    kind = value.get("kind")
    if kind == "integer":
        return int(value["value"])
    if kind == "rational":
        return Fraction(int(value["numerator"]), int(value["denominator"]))
    if kind == "symbol":
        return symbol(str(value["name"]))
    if kind in {"call", "root_of_unity"}:
        return ExactExpr(value)
    raise ValueError(f"unsupported exact Bravais expression kind {kind!r}")


def _bravais_matrix_to_user(value: Mapping[str, Any]) -> tuple[tuple[Any, ...], ...]:
    if value.get("kind") != "matrix" or int(value.get("rows", 0)) != 3 or int(
        value.get("columns", 0)
    ) != 3:
        raise ValueError("stable Bravais primitive vectors must be a 3 by 3 matrix")
    entries = value.get("entries")
    if not isinstance(entries, list) or len(entries) != 9:
        raise ValueError("stable Bravais primitive-vector entries are malformed")
    decoded = [_data_exact_to_user(entry) for entry in entries]
    return tuple(tuple(decoded[index : index + 3]) for index in range(0, 9, 3))


def _named_parameters(
    values: Mapping[str, Any], allowed_names: List[str], field_name: str
) -> Dict[str, Any]:
    provided = set(values)
    allowed = set(allowed_names)
    missing = [name for name in allowed_names if name not in provided]
    extra = sorted(provided - allowed)
    if missing or extra:
        parts = []
        if missing:
            parts.append(f"missing {', '.join(missing)}")
        if extra:
            parts.append(f"unexpected {', '.join(extra)}")
        raise ValueError(f"{field_name}: {'; '.join(parts)}")
    return {
        name: _parse_lattice_scalar(values[name], f"{field_name}.{name}")
        for name in allowed_names
    }


def _bravais_option(stable_id: str) -> Dict[str, Any]:
    cached = _BRAVAIS_OPTION_CACHE.get(stable_id)
    if cached is not None:
        return cached
    record = _catalog().bravais(stable_id)
    primitive = record["primitive_vectors"]
    conventional = record["conventional_vectors"]
    symbols: List[str] = []
    _collect_symbols(primitive, symbols)
    preferred_order = ("a", "b", "c", "α", "β", "γ")
    symbols = [name for name in preferred_order if name in symbols] + [
        name for name in symbols if name not in preferred_order
    ]

    def rows(matrix: Mapping[str, Any]) -> List[List[str]]:
        entries = matrix.get("entries")
        row_count = int(matrix.get("rows", 0))
        column_count = int(matrix.get("columns", 0))
        if not isinstance(entries, list) or len(entries) != row_count * column_count:
            raise ValueError("stable Bravais display matrix is malformed")
        text_entries = [_data_expression_text(entry) for entry in entries]
        return [
            text_entries[index : index + column_count]
            for index in range(0, len(text_entries), column_count)
        ]

    result = {
        "stable_id": stable_id,
        "parameter_symbols": symbols,
        "parameter_defaults": {
            name: "Pi/2" if name in {"α", "β", "γ"} else "1"
            for name in symbols
        },
        "primitive_vectors": rows(primitive),
        "conventional_vectors": rows(conventional),
    }
    _BRAVAIS_OPTION_CACHE[stable_id] = result
    return result


def _msg_family_options(space_group_number: int) -> List[Dict[str, Any]]:
    cached = _MSG_FAMILY_CACHE.get(space_group_number)
    if cached is not None:
        return cached
    catalog = _catalog()
    selections = [
        item
        for item in catalog.classification_maps["bns"]
        if int(item["source_key"][0]) == space_group_number
    ]
    result = []
    for selection in selections:
        record = catalog.msg(str(selection["msg_id"]))
        result.append(
            {
                "stable_id": record["stable_id"],
                "source_index": record["source_index"],
                "bns": record["bns_key"],
                "og": record["og_key"],
                "symbol": record["symbol"],
                "classification": record["classification"],
                "bravais_lattice_id": record["bravais_lattice_id"],
                "operation_count": len(record["operations"]),
            }
        )
    _MSG_FAMILY_CACHE[space_group_number] = result
    return result


def _collect_symbols(value: Any, result: List[str]) -> None:
    if isinstance(value, Mapping):
        if value.get("kind") == "symbol":
            name = str(value["name"]).split("`")[-1]
            if name not in result:
                result.append(name)
        for nested in value.values():
            _collect_symbols(nested, result)
    elif isinstance(value, list):
        for nested in value:
            _collect_symbols(nested, result)


def _data_expression_text(value: Any) -> str:
    if isinstance(value, int):
        return str(value)
    if not isinstance(value, Mapping):
        return str(value)
    kind = value.get("kind")
    if kind == "integer":
        return str(value["value"])
    if kind == "rational":
        return f"{value['numerator']}/{value['denominator']}"
    if kind == "symbol":
        return str(value["name"]).split("`")[-1]
    if kind == "root_of_unity":
        return _root_of_unity_text(int(value["order"]), int(value["power"]))
    if kind == "call":
        head = str(value.get("head", "Call")).split("`")[-1]
        arguments = [_data_expression_text(item) for item in value.get("arguments", [])]
        if head == "Plus":
            return " + ".join(arguments)
        if head == "Times":
            return " * ".join(arguments)
        if head == "Power" and len(arguments) == 2:
            return f"({arguments[0]})^({arguments[1]})"
        return f"{head}({', '.join(arguments)})"
    return str(value)


def _selected_wyckoff_option(
    msg_id: str, source_ordinal: int, letter: str
) -> Dict[str, Any]:
    for option in _wyckoff_options(msg_id):
        if option["source_ordinal"] == source_ordinal and option["letter"] == letter:
            return option
    raise ValueError(
        f"MSG {msg_id} has no Wyckoff entry {letter} at source ordinal {source_ordinal}"
    )


def _parameter_seed(
    values: Mapping[str, Any],
    allowed_names: List[str],
    seed_names: tuple[str, str, str],
    field_name: str,
) -> tuple[Any, Any, Any]:
    provided = set(values)
    allowed = set(allowed_names)
    missing = [name for name in allowed_names if name not in provided]
    extra = sorted(provided - allowed)
    if missing or extra:
        parts = []
        if missing:
            parts.append(f"missing {', '.join(missing)}")
        if extra:
            parts.append(f"unexpected {', '.join(extra)}")
        raise ValueError(f"{field_name}: {'; '.join(parts)}")
    unsupported = [name for name in allowed_names if name not in seed_names]
    if unsupported:
        raise ValueError(
            f"{field_name} contains unsupported stable symbols: {', '.join(unsupported)}"
        )
    return tuple(
        _parse_scalar(values[name], f"{field_name}.{name}") if name in allowed else 0
        for name in seed_names
    )  # type: ignore[return-value]


def _wyckoff_options(msg_id: str) -> List[Dict[str, Any]]:
    cached = _WYCKOFF_OPTION_CACHE.get(msg_id)
    if cached is not None:
        return cached
    entries = _catalog().wyckoff(msg_id)["entries"]
    result = []
    for entry in entries:
        coordinate_symbols: List[str] = []
        moment_symbols: List[str] = []
        for position in entry["positions"]:
            _collect_symbols(position["coordinates"], coordinate_symbols)
            _collect_symbols(position["moment"], moment_symbols)
        result.append(
            {
                "source_ordinal": entry["source_ordinal"],
                "multiplicity": entry["multiplicity"],
                "letter": entry["letter"],
                "coordinate_symbols": coordinate_symbols,
                "moment_symbols": moment_symbols,
                # Visible form defaults for a generic point. Denominator 29 is
                # outside crystallographic rotation/translation orders, and
                # distinct numerators avoid common special loci such as x=y or
                # x=-x modulo the lattice. Explicit user values are never
                # changed or repaired.
                "coordinate_parameter_defaults": {
                    name: f"{index + 1}/29"
                    for index, name in enumerate(coordinate_symbols)
                },
                "moment_parameter_defaults": {
                    name: "0" for name in moment_symbols
                },
                "positions": [
                    {
                        "coordinates": [
                            _data_expression_text(value)
                            for value in position["coordinates"]
                        ],
                        "moment": [
                            _data_expression_text(value) for value in position["moment"]
                        ],
                    }
                    for position in entry["positions"]
                ],
            }
        )
    _WYCKOFF_OPTION_CACHE[msg_id] = result
    return result
