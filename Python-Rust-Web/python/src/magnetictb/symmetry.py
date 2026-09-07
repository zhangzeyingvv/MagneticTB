"""Symmetry boundary for the Rust-backed MagneticTB Python API.

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

from .linear_algebra import (
    ExactAtom,
    ExactExpr,
    ExactField,
    MatrixInput,
    _decode_element,
    _decode_input_exact,
    _decode_matrix,
    _encode_exact,
    _freeze_matrix_input,
    _matrix,
    _ordered,
    _tagged_list,
)


@dataclass(frozen=True)
class SymmetryOperation:
    label: str
    rotation: MatrixInput
    translation: Sequence[ExactAtom] = (0, 0, 0)
    antiunitary: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(self, "rotation", _freeze_matrix_input(self.rotation))
        object.__setattr__(self, "translation", tuple(self.translation))
        if not isinstance(self.antiunitary, bool):
            raise TypeError("antiunitary must be bool")

    @property
    def parity(self) -> Literal["F", "T"]:
        return "T" if self.antiunitary else "F"

    def to_tuple(self) -> Tuple[str, MatrixInput, Tuple[ExactAtom, ...], str]:
        return (self.label, self.rotation, tuple(self.translation), self.parity)

    def __len__(self) -> int:
        return 4

    def __iter__(self):
        return iter(self.to_tuple())

    def __getitem__(self, index: int) -> Any:
        return self.to_tuple()[index]


@dataclass(frozen=True)
class SpinSpaceOperation:
    label: str
    space_rotation: MatrixInput
    space_translation: Sequence[ExactAtom]
    spin_rotation: MatrixInput
    antiunitary: bool = False
    continuous_parameter: Optional[ExactExpr] = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "space_rotation", _freeze_matrix_input(self.space_rotation)
        )
        object.__setattr__(self, "space_translation", tuple(self.space_translation))
        object.__setattr__(
            self, "spin_rotation", _freeze_matrix_input(self.spin_rotation)
        )
        if not isinstance(self.antiunitary, bool):
            raise TypeError("antiunitary must be bool")
        if self.continuous_parameter is not None and not isinstance(
            self.continuous_parameter, ExactExpr
        ):
            raise TypeError("continuous_parameter must be an exact symbol AST or None")


def _encode_symmetry(value: Any) -> Dict[str, Any]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or not value:
        raise TypeError("symminformation must be a nonempty ordered sequence")
    spin_space = any(isinstance(operation, SpinSpaceOperation) for operation in value)
    if spin_space:
        entries = []
        for operation in value:
            if not isinstance(operation, SpinSpaceOperation):
                raise TypeError("spin-space symmetry entries cannot mix traditional records")
            fields = [
                (
                    "space",
                    _tagged_list(
                        [_matrix(operation.space_rotation), _tagged_list([_encode_exact(v) for v in operation.space_translation])]
                    ),
                ),
                (
                    "spin",
                    _tagged_list([_matrix(operation.spin_rotation), int(operation.antiunitary)]),
                ),
            ]
            if operation.continuous_parameter is not None:
                fields.append(
                    ("continuous", _encode_exact(operation.continuous_parameter))
                )
            entries.append(
                (
                    operation.label,
                    _ordered(fields),
                )
            )
        return _ordered(entries)
    encoded = []
    for operation in value:
        if isinstance(operation, SymmetryOperation):
            label = operation.label
            rotation = operation.rotation
            translation = operation.translation
            antiunitary = operation.antiunitary
        else:
            if not isinstance(operation, Sequence) or len(operation) != 4:
                raise TypeError("traditional symmetry records require label, rotation, translation, F/T")
            label, rotation, translation, parity = operation
            if parity not in ("F", "T"):
                raise ValueError("traditional symmetry parity must be 'F' or 'T'")
            antiunitary = parity == "T"
        if not isinstance(label, str) or not label:
            raise TypeError("symmetry label must be a nonempty string")
        encoded.append(
            _tagged_list(
                [label, _matrix(rotation), _tagged_list([_encode_exact(v) for v in translation]), "T" if antiunitary else "F"]
            )
        )
    return _tagged_list(encoded)


def _decode_broken_symmetry_operations(value: Mapping[str, Any]) -> Tuple[Any, ...]:
    kind = value.get("kind")
    if kind == "list":
        operations = []
        for operation in value["items"]:
            record = _decode_input_exact(operation)
            if not isinstance(record, tuple) or len(record) != 4:
                raise RuntimeError("Rust returned a malformed symmetry operation")
            operations.append(record)
        return tuple(operations)
    if kind != "ordered_association":
        raise RuntimeError("Rust returned malformed subgroup symmetry information")
    operations = []
    for entry in value["entries"]:
        fields = dict(_decode_input_exact(entry["value"]))
        space = fields.get("space")
        spin = fields.get("spin")
        if (
            not isinstance(space, tuple)
            or len(space) != 2
            or not isinstance(spin, tuple)
            or len(spin) != 2
        ):
            raise RuntimeError("Rust returned malformed spin-space subgroup data")
        operations.append(
            SpinSpaceOperation(
                label=str(entry["key"]),
                space_rotation=space[0],
                space_translation=space[1],
                spin_rotation=spin[0],
                antiunitary=bool(spin[1]),
                continuous_parameter=fields.get("continuous"),
            )
        )
    return tuple(operations)


def _decode_prepared_symmetry_operations(
    raw: Mapping[str, Any], field_spec: ExactField
) -> Tuple[Union[SymmetryOperation, SpinSpaceOperation], ...]:
    labels = tuple(str(label) for label in raw["operation_labels"])
    spatial = tuple(raw["spatial_actions"])
    antiunitary = tuple(bool(value) for value in raw["antiunitary_flags"])
    spin_rotations = raw.get("spin_rotations")
    if len(labels) != len(spatial) or len(spatial) != len(antiunitary):
        raise RuntimeError("Rust returned misaligned ordered symmetry operations")
    result = []
    for index, (label, operation, flag) in enumerate(
        zip(labels, spatial, antiunitary)
    ):
        rotation = _decode_matrix(operation["rotation"], field_spec).rows
        translation = tuple(
            _decode_element(value, field_spec) for value in operation["translation"]
        )
        if spin_rotations is None:
            result.append(
                SymmetryOperation(
                    label=label,
                    rotation=rotation,
                    translation=translation,
                    antiunitary=flag,
                )
            )
        else:
            result.append(
                SpinSpaceOperation(
                    label=label,
                    space_rotation=rotation,
                    space_translation=translation,
                    spin_rotation=_decode_matrix(
                        spin_rotations[index], field_spec
                    ).rows,
                    antiunitary=flag,
                )
            )
    return tuple(result)
