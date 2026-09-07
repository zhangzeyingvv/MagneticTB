"""Physical Representation boundary for the Rust-backed MagneticTB Python API.

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

from .core_bindings import geometry
from .errors import ExactError, ModelError

from .linear_algebra import (
    ExactField,
    MatrixInput,
    _decode_matrix,
    _encode_exact,
    _matrix,
    _matrix_to_user,
)

from .symmetry import (
    SymmetryOperation,
)

from .representation_theory import (
    ExplicitPolynomialBasis,
)


def pointMatrix(
    symmetry_information: Sequence[Any],
    orbit_basis: Union[Sequence[str], ExplicitPolynomialBasis],
    evaluated_lattice: MatrixInput,
    *,
    exact_field: Optional[ExactField] = None,
) -> list[list[list[Any]]]:
    """Compile stable exact point-operation matrices in the Rust core."""

    if (
        not isinstance(symmetry_information, Sequence)
        or isinstance(symmetry_information, (str, bytes))
        or not symmetry_information
    ):
        raise TypeError("symmetry_information must be a nonempty ordered sequence")
    spatial_actions = []
    antiunitary_flags = []
    for operation in symmetry_information:
        if isinstance(operation, SymmetryOperation):
            rotation = operation.rotation
            translation = operation.translation
            antiunitary = operation.antiunitary
        else:
            if (
                not isinstance(operation, Sequence)
                or isinstance(operation, (str, bytes))
                or len(operation) != 4
            ):
                raise TypeError(
                    "pointMatrix operations require label, rotation, translation, F/T"
                )
            _label, rotation, translation, parity = operation
            if parity not in ("F", "T"):
                raise ValueError("pointMatrix operation parity must be 'F' or 'T'")
            antiunitary = parity == "T"
        spatial_actions.append(
            {
                "rotation": _matrix(rotation),
                "translation": [_encode_exact(value) for value in translation],
            }
        )
        antiunitary_flags.append(antiunitary)
    if isinstance(orbit_basis, ExplicitPolynomialBasis):
        basis = _matrix(orbit_basis._matrix_rows())
    else:
        if (
            not isinstance(orbit_basis, Sequence)
            or isinstance(orbit_basis, (str, bytes))
            or not orbit_basis
            or not all(isinstance(label, str) and label for label in orbit_basis)
        ):
            raise TypeError(
                "orbit_basis must be catalog labels or ExplicitPolynomialBasis"
            )
        basis = list(orbit_basis)
    field_spec = exact_field or ExactField.crystallographic()
    try:
        raw = geometry(
            "point_matrix",
            context=field_spec._context(),
            spatial_actions=spatial_actions,
            antiunitary_flags=antiunitary_flags,
            basis=basis,
            lattice=_matrix(evaluated_lattice),
        )
    except ExactError as error:
        raise ModelError(error.tag, error.detail) from None
    return [
        _matrix_to_user(_decode_matrix(matrix, field_spec))
        for matrix in raw["matrices"]
    ]
