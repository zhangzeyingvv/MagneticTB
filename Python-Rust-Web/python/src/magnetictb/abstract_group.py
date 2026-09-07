"""Abstract Group bindings backed by the Rust extension."""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from . import _core

from .errors import (
    _translate_group_error,
)


class GroupAlgebra:
    """Thin Python handle to a validated zero-based Rust Cayley table."""

    def __init__(self, multiplication_table: Sequence[Sequence[int]]) -> None:
        self._inner = _translate_group_error(_core.GroupAlgebra, multiplication_table)

    @classmethod
    def _from_core(cls, inner: Any) -> "GroupAlgebra":
        result = cls.__new__(cls)
        result._inner = inner
        return result

    @property
    def multiplication_table(self) -> List[List[int]]:
        return self._inner.multiplication_table

    @property
    def order(self) -> int:
        return self._inner.order

    @property
    def identity(self) -> int:
        return self._inner.identity

    @property
    def inverse_indices(self) -> List[int]:
        return self._inner.inverse_indices

    def product(self, left: int, right: int) -> int:
        return _translate_group_error(self._inner.product, left, right)

    def generated_subgroup(self, generators: Sequence[int]) -> List[int]:
        return _translate_group_error(self._inner.generated_subgroup, generators)

    def find_generator_indices(
        self, subgroup: Optional[Sequence[int]] = None
    ) -> List[int]:
        return _translate_group_error(self._inner.find_generator_indices, subgroup)

    def right_coset(
        self, subgroup: Sequence[int], representative: int
    ) -> List[int]:
        return _translate_group_error(
            self._inner.right_coset, subgroup, representative
        )

    def left_coset(
        self, subgroup: Sequence[int], representative: int
    ) -> List[int]:
        return _translate_group_error(self._inner.left_coset, subgroup, representative)

    def schreier_decomposition(
        self,
        subgroup: Sequence[int],
        representatives: Sequence[int],
        operation: int,
        source: int,
    ) -> Tuple[int, int]:
        return _translate_group_error(
            self._inner.schreier_decomposition,
            subgroup,
            representatives,
            operation,
            source,
        )

    def compile_action(self, action_table: Sequence[Sequence[int]]) -> "GroupAction":
        inner = _translate_group_error(self._inner.compile_action, action_table)
        return GroupAction._from_core(inner)


class GroupAction:
    """Thin Python handle to a Rust-validated operation-by-source action."""

    @classmethod
    def _from_core(cls, inner: Any) -> "GroupAction":
        result = cls.__new__(cls)
        result._inner = inner
        return result

    @property
    def action_table(self) -> List[List[int]]:
        return self._inner.action_table

    @property
    def object_count(self) -> int:
        return self._inner.object_count

    @property
    def faithful(self) -> bool:
        return self._inner.faithful

    def image(self, operation: int, source: int) -> int:
        return _translate_group_error(self._inner.image, operation, source)

    def orbit(self, source: int) -> List[int]:
        return _translate_group_error(self._inner.orbit, source)

    def orbits(self) -> List[List[int]]:
        return _translate_group_error(self._inner.orbits)

    def stabilizer(self, source: int) -> List[int]:
        return _translate_group_error(self._inner.stabilizer, source)

    def transporters(self, source: int, target: int) -> List[int]:
        return _translate_group_error(self._inner.transporters, source, target)


@dataclass(frozen=True)
class OrderedElement:
    """Stable metadata aligned with one Cayley-table index."""

    stable_id: str
    label: str
    antiunitary: bool = False


class OrderedFiniteGroup:
    """A Rust-validated table together with stable ordered element metadata."""

    def __init__(
        self,
        multiplication_table: Sequence[Sequence[int]],
        ordered_elements: Sequence[OrderedElement],
    ) -> None:
        packed = [
            (element.stable_id, element.label, element.antiunitary)
            for element in ordered_elements
        ]
        self._inner = _translate_group_error(
            _core.OrderedFiniteGroup, multiplication_table, packed
        )

    @property
    def algebra(self) -> GroupAlgebra:
        return GroupAlgebra._from_core(self._inner.algebra)

    @property
    def ordered_elements(self) -> List[OrderedElement]:
        return [OrderedElement(*values) for values in self._inner.ordered_elements]


def group_algebra_q(multiplication_table: Sequence[Sequence[int]]) -> bool:
    """Return whether Rust accepts a zero-based Cayley table."""

    return _core.GroupAlgebra.is_valid_table(multiplication_table)


def group_action_table_q(
    group: GroupAlgebra, action_table: Sequence[Sequence[int]]
) -> bool:
    """Return whether Rust accepts an operation-by-source action table."""

    return group._inner.action_table_is_valid(action_table)


def GenerateGroup(
    generators: Sequence[Any],
    identity_element: Any,
    multiply: Callable[[Any, Any], Any],
    *,
    SameTest: Optional[Callable[[Any, Any], bool]] = None,
) -> List[Any]:
    """Return the stable ordered finite closure using the Rust algorithm."""

    if isinstance(generators, (str, bytes)) or not isinstance(generators, Sequence):
        raise TypeError("generators must be an ordered sequence")
    if not callable(multiply):
        raise TypeError("multiply must be callable")
    if SameTest is not None and not callable(SameTest):
        raise TypeError("SameTest must be callable or None")
    return list(
        _translate_group_error(
            _core.generate_group_objects,
            list(generators),
            identity_element,
            multiply,
            SameTest,
        )
    )


def getGenerator(
    elements: Sequence[Any],
    identity_element: Any,
    multiply: Callable[[Any, Any], Any],
    *,
    SameTest: Optional[Callable[[Any, Any], bool]] = None,
) -> List[Any]:
    """Return stable concrete generators selected by the Rust algorithm."""

    if isinstance(elements, (str, bytes)) or not isinstance(elements, Sequence):
        raise TypeError("elements must be an ordered sequence")
    if not callable(multiply):
        raise TypeError("multiply must be callable")
    if SameTest is not None and not callable(SameTest):
        raise TypeError("SameTest must be callable or None")
    return list(
        _translate_group_error(
            _core.find_generator_objects,
            list(elements),
            identity_element,
            multiply,
            SameTest,
        )
    )
