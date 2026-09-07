"""Data Api bindings backed by the Rust extension."""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from . import _core

from .errors import (
    _translate_core_error,
    _translate_data_error,
    _translate_symmetry_error,
)

from .core_bindings import (
    _encode,
)


class DataCatalog:
    """Thin Python handle to the complete, Rust-owned MagneticTB data snapshot."""

    def __init__(self) -> None:
        self._inner = _translate_data_error(_core.DataCatalog)

    @staticmethod
    def _decode(payload: str) -> Dict[str, Any]:
        value = json.loads(payload)
        if not isinstance(value, dict):
            raise RuntimeError("Rust data core returned a non-object JSON result")
        return value

    @property
    def counts(self) -> Dict[str, int]:
        return self._decode(_translate_data_error(self._inner.counts_json))

    @property
    def manifest(self) -> Dict[str, Any]:
        return self._decode(_translate_data_error(self._inner.manifest_json))

    @property
    def bravais_ids(self) -> List[str]:
        return self._inner.bravais_ids()

    @property
    def msg_ids(self) -> List[str]:
        return self._inner.msg_ids()

    @property
    def rod_group_ids(self) -> List[str]:
        return self._inner.rod_group_ids()

    @property
    def layer_group_ids(self) -> List[str]:
        return self._inner.layer_group_ids()

    @property
    def classification_maps(self) -> Dict[str, Any]:
        return self._decode(
            _translate_data_error(self._inner.classification_maps_json)
        )

    def bravais(self, stable_id: str) -> Dict[str, Any]:
        return self._decode(
            _translate_data_error(self._inner.bravais_json, stable_id)
        )

    def msg(self, stable_id: str) -> Dict[str, Any]:
        return self._decode(_translate_data_error(self._inner.msg_json, stable_id))

    def resolve_msg_id(self, *, bns: Tuple[int, int]) -> str:
        """Resolve one exact BNS key in the Rust-owned catalog."""

        if (
            not isinstance(bns, tuple)
            or len(bns) != 2
            or any(isinstance(value, bool) or not isinstance(value, int) for value in bns)
        ):
            raise TypeError("bns must be a pair of integers")
        return str(
            _translate_data_error(
                self._inner.resolve_msg_id_by_bns,
                bns[0],
                bns[1],
            )
        )

    def resolve_msg_selector(self, classification: str, source_index: int) -> str:
        """Resolve a one-based position in a classification list in Rust.

        This is retained as an advanced compatibility endpoint.  Stable
        Mathematica ``gray/type*`` dictionaries are keyed associations, so
        ordinary public calls use :meth:`resolve_msg_key` instead.
        """

        if classification not in ("gray", "typeI", "typeIII", "typeIV"):
            raise ValueError("classification must be gray, typeI, typeIII, or typeIV")
        if (
            isinstance(source_index, bool)
            or not isinstance(source_index, int)
            or source_index < 1
        ):
            raise TypeError("source_index must be a positive integer")
        return str(
            _translate_data_error(
                self._inner.resolve_msg_id_by_classification,
                classification,
                source_index,
            )
        )

    def resolve_msg_source_index(self, source_index: int) -> str:
        """Resolve the one-based integer accepted directly by stable ``msgop``."""

        if (
            isinstance(source_index, bool)
            or not isinstance(source_index, int)
            or source_index < 1
        ):
            raise TypeError("source_index must be a positive integer")
        return str(
            _translate_data_error(
                self._inner.resolve_msg_id_by_source_index,
                source_index,
            )
        )

    def resolve_msg_key(self, classification: str, source_key: Sequence[int]) -> str:
        """Resolve an exact stable ``gray/type*/bns/og`` Association key."""

        expected_lengths = {
            "gray": 1,
            "typeI": 2,
            "typeIII": 2,
            "typeIV": 2,
            "bns": 2,
            "og": 3,
            "og_number": 1,
        }
        if classification not in expected_lengths:
            raise ValueError(
                "classification must be gray, typeI, typeIII, typeIV, bns, og, or og_number"
            )
        if isinstance(source_key, (str, bytes)) or not isinstance(source_key, Sequence):
            raise TypeError("source_key must be an ordered integer sequence")
        key = tuple(source_key)
        if len(key) != expected_lengths[classification] or any(
            isinstance(value, bool) or not isinstance(value, int) or value < 1
            for value in key
        ):
            raise TypeError(
                f"{classification} key must contain {expected_lengths[classification]} positive integer(s)"
            )
        return str(
            _translate_data_error(
                self._inner.resolve_msg_id_by_classification_key,
                classification,
                list(key),
            )
        )

    def resolve_msg_source_key(
        self, classification: str, source_key: Sequence[int]
    ) -> int:
        """Return the one-based ``MSGOP`` index stored by a stable Association."""

        expected_lengths = {
            "gray": 1,
            "typeI": 2,
            "typeIII": 2,
            "typeIV": 2,
            "bns": 2,
            "og": 3,
            "og_number": 1,
        }
        if classification not in expected_lengths:
            raise ValueError(
                "classification must be gray, typeI, typeIII, typeIV, bns, og, or og_number"
            )
        if isinstance(source_key, (str, bytes)) or not isinstance(source_key, Sequence):
            raise TypeError("source_key must be an ordered integer sequence")
        key = tuple(source_key)
        if len(key) != expected_lengths[classification] or any(
            isinstance(value, bool) or not isinstance(value, int) or value < 1
            for value in key
        ):
            raise TypeError(
                f"{classification} key must contain "
                f"{expected_lengths[classification]} positive integer(s)"
            )
        return int(
            _translate_data_error(
                self._inner.resolve_msg_source_index_by_classification_key,
                classification,
                list(key),
            )
        )

    def wyckoff(self, msg_id: str) -> Dict[str, Any]:
        return self._decode(_translate_data_error(self._inner.wyckoff_json, msg_id))

    def resolve_wyckoff(
        self,
        msg_id: str,
        letter: str,
        source_ordinal: Optional[int] = None,
    ) -> Dict[str, Any]:
        """Resolve one letter/ordinal selection without Python-side data logic."""

        if not isinstance(letter, str) or not letter:
            raise TypeError("letter must be a nonempty string")
        if source_ordinal is not None and (
            isinstance(source_ordinal, bool)
            or not isinstance(source_ordinal, int)
            or source_ordinal < 1
        ):
            raise TypeError("source_ordinal must be a positive integer or None")
        return self._decode(
            _translate_data_error(
                self._inner.resolve_wyckoff_json,
                msg_id,
                letter,
                source_ordinal,
            )
        )

    def rod_group(self, stable_id: str) -> Dict[str, Any]:
        return self._decode(
            _translate_data_error(self._inner.rod_group_json, stable_id)
        )

    def layer_group(self, stable_id: str) -> Dict[str, Any]:
        return self._decode(
            _translate_data_error(self._inner.layer_group_json, stable_id)
        )

    def resolve_subperiodic_group_id(
        self, kind: str, og_key: Sequence[int]
    ) -> str:
        """Resolve one stable magnetic rod/layer OG key in Rust."""

        if kind not in ("rod", "layer"):
            raise ValueError("kind must be rod or layer")
        if isinstance(og_key, (str, bytes)) or not isinstance(og_key, Sequence):
            raise TypeError("og_key must be a three-integer sequence")
        key = tuple(og_key)
        if len(key) != 3 or any(
            isinstance(value, bool) or not isinstance(value, int) or value < 1
            for value in key
        ):
            raise TypeError("og_key must be a three-positive-integer sequence")
        return str(
            _translate_data_error(
                self._inner.resolve_subperiodic_group_id_by_og,
                kind,
                list(key),
            )
        )

    def resolve_subperiodic_gray_key(self, kind: str, source_key: int) -> Tuple[int, int, int]:
        """Resolve one stable ``grayrod``/``graylayer`` Association key."""

        if kind not in ("rod", "layer"):
            raise ValueError("kind must be rod or layer")
        if isinstance(source_key, bool) or not isinstance(source_key, int) or source_key < 1:
            raise TypeError("source_key must be a positive integer")
        result = tuple(
            int(value)
            for value in _translate_data_error(
                self._inner.resolve_subperiodic_gray_key,
                kind,
                source_key,
            )
        )
        if len(result) != 3:
            raise RuntimeError("Rust Data returned a malformed subperiodic OG key")
        return result  # type: ignore[return-value]

    def subperiodic_basic_vectors(self, kind: str, stable_id: str) -> Any:
        """Return the stable basic-vector matrix for one rod/layer group."""

        if kind not in ("rod", "layer"):
            raise ValueError("kind must be rod or layer")
        return json.loads(
            _translate_data_error(
                self._inner.subperiodic_basic_vectors_json,
                kind,
                stable_id,
            )
        )

    def compile_msg_group(self, stable_id: str) -> Dict[str, Any]:
        """Compile one frozen MSG record through the Rust exact Seitz core."""

        return self._decode(
            _translate_symmetry_error(self._inner.compile_msg_group_json, stable_id)
        )

    def compile_msg_wyckoff_sites(
        self,
        msg_id: str,
        source_ordinal: int,
        context: Dict[str, Any],
        bindings: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Evaluate a frozen Wyckoff orbit and its exact Rust site action."""

        payload = json.dumps(
            {
                "msg_id": msg_id,
                "source_ordinal": source_ordinal,
                "context": context,
                "bindings": bindings,
            },
            separators=(",", ":"),
        )
        return self._decode(
            _translate_core_error(self._inner.compile_msg_wyckoff_sites_json, payload)
        )

    def compile_directed_bond_orbits(
        self,
        msg_id: str,
        source_ordinal: int,
        context: Dict[str, Any],
        bindings: Dict[str, Any],
        sites: Sequence[Sequence[Dict[str, Any]]],
        lattice: Sequence[Sequence[Dict[str, Any]]],
        translation_bound: int,
        requested_shells: int,
        shell_index: int,
        scalar_constraints: bool = False,
        solve_and_verify: bool = False,
    ) -> Dict[str, Any]:
        """Compile exact directed-bond orbits from Data-owned symmetry."""

        payload = json.dumps(
            {
                "operation": "compile_data_directed_bond_orbits",
                "msg_id": msg_id,
                "source_ordinal": source_ordinal,
                "context": context,
                "bindings": bindings,
                "sites": sites,
                "lattice": lattice,
                "translation_bound": translation_bound,
                "requested_shells": requested_shells,
                "shell_index": shell_index,
                "scalar_constraints": scalar_constraints,
                "solve_and_verify": solve_and_verify,
            },
            separators=(",", ":"),
        )
        return self._decode(
            _translate_core_error(self._inner.compile_directed_bond_orbits_json, payload)
        )

    def compile_scalar_physical_representation(
        self,
        msg_id: str,
        context: Dict[str, Any],
        site_orbits: Sequence[Sequence[Sequence[Dict[str, Any]]]],
    ) -> Dict[str, Any]:
        """Compile a Data-owned scalar physical representation in Rust."""

        payload = json.dumps(
            {
                "msg_id": msg_id,
                "context": context,
                "site_orbits": site_orbits,
            },
            separators=(",", ":"),
        )
        return self._decode(
            _translate_core_error(
                self._inner.compile_scalar_physical_representation_json,
                payload,
            )
        )

    def compile_scalar_model_closure(
        self,
        msg_id: str,
        context: Dict[str, Any],
        site_orbits: Sequence[Sequence[Sequence[Dict[str, Any]]]],
        sites: Sequence[Sequence[Dict[str, Any]]],
        lattice: Sequence[Sequence[Dict[str, Any]]],
        translation_bound: int,
        requested_shells: int,
        shell_index: int,
    ) -> Dict[str, Any]:
        """Compile and exactly verify a scalar Data model in the Rust core."""

        payload = json.dumps(
            {
                "operation": "compile_data_scalar_model_closure",
                "msg_id": msg_id,
                "context": context,
                "site_orbits": site_orbits,
                "sites": sites,
                "lattice": lattice,
                "translation_bound": translation_bound,
                "requested_shells": requested_shells,
                "shell_index": shell_index,
                "scalar_constraints": True,
                "solve_and_verify": True,
            },
            separators=(",", ":"),
        )
        return self._decode(
            _translate_core_error(self._inner.compile_directed_bond_orbits_json, payload)
        )

    def compile_model_input_closure(
        self,
        context: Dict[str, Any],
        compiler_input: Dict[str, Any],
        data_trace: Dict[str, Any],
        target_shell: int,
        *,
        hermitian: bool = True,
        kernel_method: str = "Iterative",
        validation_level: str = "Basic",
    ) -> Dict[str, Any]:
        """Compile Data, raw init input, and the exact Hamiltonian entirely in Rust."""

        payload = json.dumps(
            {
                "operation": "compile_data_model_input_closure",
                "context": context,
                "compiler_input": compiler_input,
                "data_trace": data_trace,
                "target_shell": target_shell,
                "hermitian": hermitian,
                "kernel_method": kernel_method,
                "validation_level": validation_level,
            },
            separators=(",", ":"),
        )
        return self._decode(
            _translate_core_error(self._inner.compile_directed_bond_orbits_json, payload)
        )

    def compile_model_input_closures(
        self,
        context: Dict[str, Any],
        compiler_input: Dict[str, Any],
        data_trace: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Compile every prepared shell while sharing the Rust model state."""

        payload = json.dumps(
            {
                "operation": "compile_data_model_input_closures",
                "context": context,
                "compiler_input": compiler_input,
                "data_trace": data_trace,
            },
            separators=(",", ":"),
        )
        return self._decode(
            _translate_core_error(self._inner.compile_directed_bond_orbits_json, payload)
        )

    def geometry(self, operation: str, **arguments: Any) -> Dict[str, Any]:
        """Call a Rust geometry operation with this frozen Data catalog."""

        if not isinstance(operation, str):
            raise TypeError("operation must be a string")
        payload = _encode({"operation": operation, **arguments})
        return self._decode(
            _translate_core_error(self._inner.compile_directed_bond_orbits_json, payload)
        )
