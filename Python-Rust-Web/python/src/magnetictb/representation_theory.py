"""Representation Theory boundary for the Rust-backed MagneticTB Python API.

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
    MatrixInput,
    _encode_exact,
    _freeze_matrix_input,
    _matrix,
    _ordered,
    _tagged_list,
)


def _reject_inexact_polynomial_value(value: Any, *, _inside_ast: bool = False) -> None:
    if isinstance(value, (float, complex)):
        if _inside_ast:
            raise TypeError("raw JSON numbers are invalid inside an ExactExpr payload")
        _encode_exact(value)
        return
    if isinstance(value, ExactExpr):
        _reject_inexact_polynomial_value(value._value, _inside_ast=True)
    elif isinstance(value, Mapping):
        for item in value.values():
            _reject_inexact_polynomial_value(item, _inside_ast=_inside_ast)
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        for item in value:
            _reject_inexact_polynomial_value(item, _inside_ast=_inside_ast)


@dataclass(frozen=True)
class ExplicitPolynomialBasis:
    """An ordered scalar or two-component exact polynomial basis for one orbit."""

    functions: Sequence[Any]
    labels: Optional[Sequence[str]] = None

    def __post_init__(self) -> None:
        if (
            not isinstance(self.functions, Sequence)
            or isinstance(self.functions, (str, bytes))
            or not self.functions
        ):
            raise TypeError("functions must be a nonempty ordered sequence")
        frozen_functions = []
        component_count: Optional[int] = None
        for function in self.functions:
            if isinstance(function, Sequence) and not isinstance(function, (str, bytes)):
                components = tuple(function)
                if len(components) != 2:
                    raise ValueError(
                        "spinor polynomial functions must contain exactly two components"
                    )
                columns = 2
                for component in components:
                    _reject_inexact_polynomial_value(component)
                    _encode_exact(component)
                frozen_functions.append(components)
            else:
                columns = 1
                _reject_inexact_polynomial_value(function)
                _encode_exact(function)
                frozen_functions.append(function)
            if component_count is None:
                component_count = columns
            elif component_count != columns:
                raise ValueError(
                    "scalar and two-component polynomial functions cannot mix within one orbit"
                )
        object.__setattr__(self, "functions", tuple(frozen_functions))
        if self.labels is not None:
            if isinstance(self.labels, (str, bytes)):
                raise TypeError("labels must be an ordered sequence of strings or None")
            labels = tuple(self.labels)
            if len(labels) != len(frozen_functions):
                raise ValueError("labels length must equal functions length")
            if not all(isinstance(label, str) and label for label in labels):
                raise TypeError("labels must contain nonempty strings")
            object.__setattr__(self, "labels", labels)

    @property
    def display_labels(self) -> Tuple[str, ...]:
        if self.labels is not None:
            return tuple(self.labels)
        return tuple(f"basis_{index}" for index in range(1, len(self.functions) + 1))

    def _matrix_rows(self) -> Tuple[Tuple[ExactAtom, ...], ...]:
        if isinstance(self.functions[0], tuple):
            return tuple(self.functions)  # type: ignore[return-value]
        return tuple((function,) for function in self.functions)


@dataclass(frozen=True)
class InducedOrbit:
    reference_site_index: int
    site_symmetry_operation_indices: Sequence[int]
    site_symmetry_matrices: Sequence[MatrixInput]

    def __post_init__(self) -> None:
        indices = tuple(self.site_symmetry_operation_indices)
        matrices = tuple(
            _freeze_matrix_input(matrix) for matrix in self.site_symmetry_matrices
        )
        object.__setattr__(
            self,
            "site_symmetry_operation_indices",
            indices,
        )
        object.__setattr__(
            self,
            "site_symmetry_matrices",
            matrices,
        )


@dataclass(frozen=True)
class InducedRepresentation:
    """Ordered induced orbits plus optional exact continuous site matrices."""

    orbits: Sequence[InducedOrbit]
    continuous_site_matrices: Optional[Sequence[Sequence[MatrixInput]]] = None

    def __post_init__(self) -> None:
        if (
            not isinstance(self.orbits, Sequence)
            or isinstance(self.orbits, (str, bytes))
            or not self.orbits
            or not all(isinstance(orbit, InducedOrbit) for orbit in self.orbits)
        ):
            raise TypeError("orbits must be a nonempty ordered InducedOrbit sequence")
        orbits = tuple(self.orbits)
        object.__setattr__(self, "orbits", orbits)
        if self.continuous_site_matrices is None:
            return
        if (
            not isinstance(self.continuous_site_matrices, Sequence)
            or isinstance(self.continuous_site_matrices, (str, bytes))
            or not self.continuous_site_matrices
        ):
            raise TypeError(
                "continuous_site_matrices must be a nonempty ordered orbit sequence"
            )
        continuous = []
        for orbit in self.continuous_site_matrices:
            if (
                not isinstance(orbit, Sequence)
                or isinstance(orbit, (str, bytes))
                or not orbit
            ):
                raise TypeError(
                    "each continuous induced orbit must contain a nonempty matrix sequence"
                )
            continuous.append(tuple(_freeze_matrix_input(matrix) for matrix in orbit))
        object.__setattr__(self, "continuous_site_matrices", tuple(continuous))


@dataclass(frozen=True)
class DirectProductRepresentation:
    """Ordered finite matrices plus optional exact continuous site matrices."""

    matrices_by_orbit: Sequence[Sequence[MatrixInput]]
    continuous_site_matrices: Optional[Sequence[Sequence[MatrixInput]]] = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "matrices_by_orbit",
            tuple(
                tuple(_freeze_matrix_input(matrix) for matrix in orbit)
                for orbit in self.matrices_by_orbit
            ),
        )
        if self.continuous_site_matrices is not None:
            object.__setattr__(
                self,
                "continuous_site_matrices",
                tuple(
                    tuple(_freeze_matrix_input(matrix) for matrix in orbit)
                    for orbit in self.continuous_site_matrices
                ),
            )


def _normalise_basis(
    values: Any, orbit_count: int
) -> Tuple[Union[Tuple[str, ...], ExplicitPolynomialBasis], ...]:
    if isinstance(values, ExplicitPolynomialBasis):
        if orbit_count != 1:
            raise ValueError(
                "a direct ExplicitPolynomialBasis value is valid only for one Wyckoff orbit"
            )
        return (values,)
    if not isinstance(values, Sequence) or isinstance(values, (str, bytes)):
        raise TypeError("basis_functions must be an ordered basis specification")
    if orbit_count == 1 and values and all(isinstance(value, str) for value in values):
        values = (values,)
    if len(values) != orbit_count:
        raise ValueError("basis_functions must contain one ordered basis list per Wyckoff orbit")
    result: list[Union[Tuple[str, ...], ExplicitPolynomialBasis]] = []
    for orbit in values:
        if isinstance(orbit, ExplicitPolynomialBasis):
            result.append(orbit)
            continue
        if not isinstance(orbit, Sequence) or isinstance(orbit, (str, bytes)) or not orbit:
            raise TypeError("each basis-function orbit must be a nonempty sequence")
        if not all(isinstance(label, str) and label for label in orbit):
            raise TypeError(
                "each orbit must use either catalog strings or one ExplicitPolynomialBasis"
            )
        result.append(tuple(orbit))
    return tuple(result)


def _encode_basis_functions(
    values: Sequence[Union[Tuple[str, ...], ExplicitPolynomialBasis]],
) -> Dict[str, Any]:
    encoded = []
    for orbit in values:
        if isinstance(orbit, ExplicitPolynomialBasis):
            encoded.append(_matrix(orbit._matrix_rows()))
        else:
            encoded.append(_tagged_list(list(orbit)))
    return _tagged_list(encoded)


def _basis_display_labels(
    values: Sequence[Union[Tuple[str, ...], ExplicitPolynomialBasis]],
) -> Tuple[Tuple[str, ...], ...]:
    return tuple(
        orbit.display_labels if isinstance(orbit, ExplicitPolynomialBasis) else orbit
        for orbit in values
    )


def _normalise_orbital_labels(
    values: Optional[Sequence[Sequence[str]]], orbit_count: int
) -> Tuple[Tuple[str, ...], ...]:
    if values is None:
        return ()
    if len(values) != orbit_count:
        raise ValueError(
            "orbital_labels must contain one ordered label list per Wyckoff orbit"
        )
    result = tuple(tuple(orbit) for orbit in values)
    if any(not orbit or not all(isinstance(label, str) and label for label in orbit) for orbit in result):
        raise TypeError("orbital_labels must contain nonempty string lists")
    return result


def _encode_representation_information(value: Any, mode: str) -> Any:
    if mode == "Induced":
        if isinstance(value, DirectProductRepresentation):
            raise TypeError(
                "DirectProductRepresentation cannot be used with representation_mode='Induced'"
            )
        continuous = None
        if isinstance(value, InducedRepresentation):
            orbits = value.orbits
            continuous = value.continuous_site_matrices
        else:
            orbits = value
        if not isinstance(orbits, Sequence) or not orbits:
            raise TypeError("Induced repinformation must be a nonempty sequence")
        entries = []
        for orbit in orbits:
            if not isinstance(orbit, InducedOrbit):
                raise TypeError("Induced repinformation entries must be InducedOrbit")
            entries.append(
                _ordered(
                    [
                        ("ReferenceSiteIndex", orbit.reference_site_index),
                        ("SiteSymmetryOperationIndices", _tagged_list(list(orbit.site_symmetry_operation_indices))),
                        ("SiteSymmetryMatrices", _tagged_list([_matrix(matrix) for matrix in orbit.site_symmetry_matrices])),
                    ]
                )
            )
        discrete = _tagged_list(entries)
        if continuous is None:
            return discrete
        return _ordered(
            [
                ("Discrete", discrete),
                (
                    "Continuous",
                    _ordered(
                        [
                            (
                                "SiteMatrices",
                                _tagged_list(
                                    [
                                        _tagged_list(
                                            [_matrix(matrix) for matrix in orbit]
                                        )
                                        for orbit in continuous
                                    ]
                                ),
                            )
                        ]
                    ),
                ),
            ]
        )
    if mode != "DirectProduct":
        raise ValueError("unsupported representation mode")
    if isinstance(value, InducedRepresentation):
        raise TypeError(
            "InducedRepresentation requires representation_mode='Induced'"
        )
    if isinstance(value, DirectProductRepresentation):
        discrete = _tagged_list(
            [
                _tagged_list([_matrix(matrix) for matrix in orbit])
                for orbit in value.matrices_by_orbit
            ]
        )
        if value.continuous_site_matrices is None:
            return discrete
        continuous = _tagged_list(
            [
                _tagged_list([_matrix(matrix) for matrix in orbit])
                for orbit in value.continuous_site_matrices
            ]
        )
        return _ordered(
            [
                ("Discrete", discrete),
                ("Continuous", _ordered([("SiteMatrices", continuous)])),
            ]
        )
    if not isinstance(value, Sequence) or not value:
        raise TypeError("DirectProduct repinformation must contain one matrix list per orbit")
    return _tagged_list(
        [_tagged_list([_matrix(matrix) for matrix in orbit]) for orbit in value]
    )
