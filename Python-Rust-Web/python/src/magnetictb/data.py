"""Data boundary for the Rust-backed MagneticTB Python API.

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

from .data_api import DataCatalog
from .errors import DataError

from .linear_algebra import (
    CyclotomicElement,
    ExactAtom,
    ExactExpr,
    _cyclotomic_to_user,
    _decode_element,
    _decode_matrix,
    _field_from_context,
)


@lru_cache(maxsize=1)
def _data_catalog() -> DataCatalog:
    """Return the immutable Rust-owned stable Data catalog."""

    return DataCatalog()


@dataclass(frozen=True)
class MSG:
    """A stable Magnetic Space Group selection owned by Rust DataCatalog."""

    stable_id: Optional[str] = None
    bns: Optional[Tuple[int, int]] = None

    def __post_init__(self) -> None:
        if (self.stable_id is None) == (self.bns is None):
            raise ValueError("MSG requires exactly one of stable_id or bns")
        if self.stable_id is not None and (
            not isinstance(self.stable_id, str) or not self.stable_id
        ):
            raise TypeError("stable_id must be a nonempty string")
        if self.bns is not None and (
            not isinstance(self.bns, tuple)
            or len(self.bns) != 2
            or any(isinstance(value, bool) or not isinstance(value, int) for value in self.bns)
        ):
            raise TypeError("bns must be a pair of integers")


@dataclass(frozen=True)
class MagneticGroupSelector:
    """An exact key in one stable Mathematica MSG Association."""

    classification: Literal[
        "gray", "typeI", "typeIII", "typeIV", "bns", "og", "og_number"
    ]
    source_key: Tuple[int, ...]

    def __post_init__(self) -> None:
        expected_lengths = {
            "gray": 1,
            "typeI": 2,
            "typeIII": 2,
            "typeIV": 2,
            "bns": 2,
            "og": 3,
            "og_number": 1,
        }
        if self.classification not in expected_lengths:
            raise ValueError("unsupported magnetic-group classification")
        object.__setattr__(self, "source_key", tuple(self.source_key))
        if len(self.source_key) != expected_lengths[self.classification] or any(
            isinstance(value, bool) or not isinstance(value, int) or value < 1
            for value in self.source_key
        ):
            raise TypeError(
                f"{self.classification} key must contain "
                f"{expected_lengths[self.classification]} positive integer(s)"
            )


@dataclass(frozen=True)
class _MagneticGroupDictionary(Mapping):
    """Read-only Python mapping of one stable Mathematica Association."""

    classification: str
    key_length: int

    def _selector(self, key: Any) -> MagneticGroupSelector:
        if self.key_length == 1 and isinstance(key, int) and not isinstance(key, bool):
            values = (key,)
        elif isinstance(key, Sequence) and not isinstance(key, (str, bytes)):
            values = tuple(key)
        else:
            raise TypeError(
                f"{self.classification} key must contain {self.key_length} positive integer(s)"
            )
        return MagneticGroupSelector(self.classification, values)  # type: ignore[arg-type]

    def __call__(self, key: Any) -> int:
        return self[key]

    def __getitem__(self, key: Any) -> int:
        selector = self._selector(key)
        return _data_catalog().resolve_msg_source_key(
            selector.classification, selector.source_key
        )

    def selector(self, key: Any) -> MagneticGroupSelector:
        """Return the advanced typed selector without resolving its MSGOP index."""

        return self._selector(key)

    def keys(self) -> Tuple[Any, ...]:
        catalog_name = {
            "typeI": "type_i",
            "typeIII": "type_iii",
            "typeIV": "type_iv",
        }.get(self.classification, self.classification)
        keys = tuple(
            tuple(int(value) for value in entry["source_key"])
            for entry in _data_catalog().classification_maps[catalog_name]
        )
        if self.key_length == 1:
            return tuple(key[0] for key in keys)
        return keys

    def __iter__(self):
        return iter(self.keys())

    def __len__(self) -> int:
        return len(self.keys())

    def items(self) -> Tuple[Tuple[Any, int], ...]:
        return tuple((key, self[key]) for key in self.keys())


gray = _MagneticGroupDictionary("gray", 1)


typeI = _MagneticGroupDictionary("typeI", 2)


typeIII = _MagneticGroupDictionary("typeIII", 2)


typeIV = _MagneticGroupDictionary("typeIV", 2)


bnsdict = _MagneticGroupDictionary("bns", 2)


ogdict = _MagneticGroupDictionary("og", 3)


ognumdict = _MagneticGroupDictionary("og_number", 1)


@dataclass(frozen=True)
class _SubperiodicGrayDictionary(Mapping):
    """Read-only stable ``grayrod``/``graylayer`` Association view."""

    kind: Literal["rod", "layer"]

    def __getitem__(self, key: Any) -> Tuple[int, int, int]:
        return _data_catalog().resolve_subperiodic_gray_key(self.kind, key)

    def __iter__(self):
        source_key = 1
        while True:
            try:
                self[source_key]
            except DataError as error:
                if error.tag == "UnknownSubperiodicGraySelector":
                    return
                raise
            yield source_key
            source_key += 1

    def __len__(self) -> int:
        return sum(1 for _ in self)


grayrod = _SubperiodicGrayDictionary("rod")


graylayer = _SubperiodicGrayDictionary("layer")


@dataclass(frozen=True)
class Wyckoff:
    letter: str
    position: Sequence[ExactAtom] = (0, 0, 0)
    moment: Sequence[ExactAtom] = (0, 0, 0)
    source_ordinal: Optional[int] = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "position", tuple(self.position))
        object.__setattr__(self, "moment", tuple(self.moment))
        if not isinstance(self.letter, str) or not self.letter:
            raise TypeError("Wyckoff letter must be a nonempty string")
        if self.source_ordinal is not None and (
            isinstance(self.source_ordinal, bool)
            or not isinstance(self.source_ordinal, int)
            or self.source_ordinal < 1
        ):
            raise TypeError("source_ordinal must be a positive integer or None")


def msgop(
    group: Union[int, MSG, MagneticGroupSelector],
) -> list[list[Any]]:
    """Return stable ``{label, rotation, translation, F/T}`` operation lists."""

    catalog = _data_catalog()
    if isinstance(group, bool):
        raise TypeError(
            "msgop group must be a one-based integer, MSG(...), or gray/type*/bnsdict/ogdict selector"
        )
    if isinstance(group, int):
        stable_id = catalog.resolve_msg_source_index(group)
    elif isinstance(group, MSG):
        stable_id = group.stable_id or catalog.resolve_msg_id(bns=group.bns)  # type: ignore[arg-type]
    elif isinstance(group, MagneticGroupSelector):
        stable_id = catalog.resolve_msg_key(group.classification, group.source_key)
    else:
        raise TypeError(
            "msgop group must be a one-based integer, MSG(...), or gray/type*/bnsdict/ogdict selector"
        )
    record = catalog.msg(stable_id)
    bravais = catalog.bravais(str(record["bravais_lattice_id"]))
    print("Magnetic space group (BNS):", list(record["bns_key"]))
    print("Lattice:", record["bravais_lattice_id"])
    print("Primitive Lattice Vactor:", _catalog_matrix_to_user(bravais["primitive_vectors"]))
    print(
        "Conventional Lattice Vactor:",
        _catalog_matrix_to_user(bravais["conventional_vectors"]),
    )
    compiled = catalog.compile_msg_group(stable_id)
    field_spec = _field_from_context(compiled["context"])
    return [
        [
            str(operation["label"]),
            [
                [_stable_record_scalar(value) for value in row]
                for row in _decode_matrix(operation["rotation"], field_spec).rows
            ],
            [
                _stable_record_scalar(_decode_element(value, field_spec))
                for value in operation["translation"]
            ],
            "T" if bool(operation["antiunitary"]) else "F",
        ]
        for operation in compiled["operations"]
    ]


def _subperiodic_operations(
    kind: Literal["rod", "layer"], group: Union[str, Sequence[int]]
) -> list[list[Any]]:
    catalog = _data_catalog()
    if isinstance(group, str):
        if not group:
            raise TypeError("group stable_id must be a nonempty string")
        stable_id = group
    else:
        stable_id = catalog.resolve_subperiodic_group_id(kind, group)
    record = (
        catalog.rod_group(stable_id)
        if kind == "rod"
        else catalog.layer_group(stable_id)
    )
    title = "Magnetic rod group:" if kind == "rod" else "Magnetic layer group:"
    print(title, [".".join(str(value) for value in record["og_key"]), record["symbol"]])
    print("Lattice:", record["bravais_lattice_id"])
    print(
        "Primitive Lattice Vactor:",
        _catalog_matrix_to_user(catalog.subperiodic_basic_vectors(kind, stable_id)),
    )
    return [
        [
            str(operation["label"]),
            _catalog_matrix_to_user(operation["rotation"]),
            [_catalog_exact_to_user(value) for value in operation["translation"]],
            "T" if bool(operation["antiunitary"]) else "F",
        ]
        for operation in record["operations"]
    ]


def mlgop(group: Union[str, Sequence[int]]) -> list[list[Any]]:
    """Return stable ordered magnetic layer-group operations."""

    return _subperiodic_operations("layer", group)


def mrgop(group: Union[str, Sequence[int]]) -> list[list[Any]]:
    """Return stable ordered magnetic rod-group operations."""

    return _subperiodic_operations("rod", group)


def _stable_record_scalar(value: CyclotomicElement) -> Any:
    """Decode one exact Rust scalar using ordinary stable-style Python atoms."""

    decoded = _cyclotomic_to_user(value)
    if isinstance(decoded, Fraction) and decoded.denominator == 1:
        return decoded.numerator
    return decoded


def _catalog_exact_to_user(value: Mapping[str, Any]) -> Any:
    kind = value.get("kind")
    if kind == "integer":
        return int(value["value"])
    if kind == "rational":
        return Fraction(int(value["numerator"]), int(value["denominator"]))
    if kind == "symbol":
        name = str(value["name"])
        return name.rsplit("`", 1)[-1]
    return ExactExpr(value)


def _catalog_matrix_to_user(value: Mapping[str, Any]) -> list[list[Any]]:
    if value.get("kind") != "matrix":
        raise RuntimeError("Rust Data returned a non-matrix Bravais record")
    rows = int(value["rows"])
    columns = int(value["columns"])
    entries = tuple(_catalog_exact_to_user(entry) for entry in value["entries"])
    if len(entries) != rows * columns:
        raise RuntimeError("Rust Data returned an invalid Bravais matrix")
    return [
        [entries[row * columns + column] for column in range(columns)]
        for row in range(rows)
    ]


def _normalise_wyckoff(value: Any) -> list[Tuple[Sequence[ExactAtom], Sequence[ExactAtom]]]:
    if isinstance(value, Wyckoff):
        value = [value]
    result = []
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise TypeError("wyckoffposition must be a Wyckoff or a sequence of {position,moment} pairs")
    for site in value:
        if isinstance(site, Wyckoff):
            position, moment = site.position, site.moment
        else:
            if not isinstance(site, Sequence) or len(site) != 2:
                raise TypeError("each Wyckoff seed must contain position and moment")
            position, moment = site
        if len(position) != 3 or len(moment) != 3:
            raise ValueError("Wyckoff position and moment must each have length 3")
        result.append((position, moment))
    if not result:
        raise ValueError("at least one Wyckoff seed is required")
    return result


def _normalise_data_wyckoff(value: Any) -> Tuple[Wyckoff, ...]:
    if isinstance(value, Wyckoff):
        return (value,)
    if (
        not isinstance(value, Sequence)
        or isinstance(value, (str, bytes))
        or not value
        or not all(isinstance(item, Wyckoff) for item in value)
    ):
        raise TypeError(
            "MSG input requires wyckoffposition=Wyckoff(...) or an ordered Wyckoff sequence"
        )
    return tuple(value)
