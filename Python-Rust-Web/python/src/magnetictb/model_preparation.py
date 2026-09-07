"""Model Preparation for the Rust-backed model lifecycle."""

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

from .model_dependencies import (
    DataCatalog,
    ExactAtom,
    ExactExpr,
    ExactField,
    MSG,
    MatrixInput,
    Wyckoff,
    _basis_display_labels,
    _encode_basis_functions,
    _encode_lattpar,
    _encode_representation_information,
    _encode_symmetry,
    _matrix,
    _normalise_basis,
    _normalise_data_wyckoff,
    _normalise_orbital_labels,
    _normalise_wyckoff,
    _ordered,
    _strict_bool,
    _tagged_list,
    geometry,
)

from .model_schema import (
    PreparedModel,
    _CompilationRecipe,
)


def prepare_model(
    *,
    lattice: MatrixInput,
    lattpar: Mapping[str | ExactExpr, ExactAtom] | Sequence[Tuple[str | ExactExpr, ExactAtom]] = (),
    wyckoffposition: Any,
    symminformation: Any,
    basis_functions: Any = ("s",),
    initial_bond_shells: int = 10,
    generate_symmetry_group: bool = False,
    representation_mode: str = "DirectProduct",
    exact_field: Optional[ExactField] = None,
    debug_q: bool = False,
) -> PreparedModel:
    """Prepare an exact model without installing the stable-style session.

    ``wyckoffposition`` is either one :class:`Wyckoff`, an ordered Wyckoff
    sequence for an :class:`MSG`, or an ordered sequence of explicit
    ``(position, moment)`` seeds.  For each seed/orbit, ``basis_functions``
    contains either an ordered catalog-label sequence or one
    :class:`ExplicitPolynomialBasis`.  ``symminformation`` is exclusively
    either an MSG selection or an ordered explicit-operation sequence.
    """

    debug_q = _strict_bool(debug_q, "debug_q")
    if isinstance(symminformation, MSG):
        wyckoff = _normalise_data_wyckoff(wyckoffposition)
        return _prepare_data_model(
            lattice=lattice,
            lattpar=lattpar,
            wyckoff=wyckoff,
            msg=symminformation,
            basis_functions=basis_functions,
            initial_bond_shells=initial_bond_shells,
            representation_mode=representation_mode,
            exact_field=exact_field or ExactField.crystallographic(),
            debug_q=debug_q,
        )
    return _prepare_explicit_model(
        lattice=lattice,
        lattpar=lattpar,
        wyckoffposition=wyckoffposition,
        symminformation=symminformation,
        basis_functions=basis_functions,
        initial_bond_shells=initial_bond_shells,
        generate_symmetry_group=generate_symmetry_group,
        representation_mode=representation_mode,
        exact_field=exact_field or ExactField.crystallographic(),
        debug_q=debug_q,
    )


def prepare_model_from_rep(
    *,
    lattice: MatrixInput,
    lattpar: Mapping[str | ExactExpr, ExactAtom] | Sequence[Tuple[str | ExactExpr, ExactAtom]] = (),
    wyckoffposition: Any,
    symminformation: Any,
    repinformation: Any,
    orbital_labels: Optional[Sequence[Sequence[str]]] = None,
    initial_bond_shells: int = 10,
    representation_mode: str = "DirectProduct",
    exact_field: Optional[ExactField] = None,
    debug_q: bool = False,
) -> PreparedModel:
    """Prepare an explicit-representation model without installing a session.

    ``repinformation`` and optional ``orbital_labels`` contain exactly one
    ordered item per Wyckoff seed/orbit.  DirectProduct consumes complete-group
    matrices in operation order; Induced consumes :class:`InducedOrbit`
    records whose explicit ``*_index`` values remain zero-based.
    """

    debug_q = _strict_bool(debug_q, "debug_q")
    field_spec = exact_field or ExactField.crystallographic()
    if isinstance(symminformation, MSG):
        wyckoff = _normalise_data_wyckoff(wyckoffposition)
        compiler_input = _base_compiler_input(
            lattice,
            lattpar,
            wyckoff,
            initial_bond_shells,
            representation_mode,
            debug_q,
        )
        compiler_input["entries"].extend(
            [
                {"key": "RepresentationSource", "value": "Matrices"},
                {
                    "key": "RepresentationInformation",
                    "value": _encode_representation_information(
                        repinformation, representation_mode
                    ),
                },
            ]
        )
        labels = _normalise_orbital_labels(orbital_labels, len(wyckoff))
        return _compile_data_prepared_model(
            field_spec=field_spec,
            compiler_input=compiler_input,
            msg=symminformation,
            wyckoff=wyckoff,
            initial_bond_shells=initial_bond_shells,
            basis_labels=labels,
        )
    compiler_input = _base_compiler_input(
        lattice,
        lattpar,
        wyckoffposition,
        initial_bond_shells,
        representation_mode,
        debug_q,
    )
    compiler_input["entries"].extend(
        [
            {"key": "SymmetryInformation", "value": _encode_symmetry(symminformation)},
            {"key": "RepresentationSource", "value": "Matrices"},
            {
                "key": "RepresentationInformation",
                "value": _encode_representation_information(
                    repinformation, representation_mode
                ),
            },
        ]
    )
    labels = _normalise_orbital_labels(
        orbital_labels, len(_normalise_wyckoff(wyckoffposition))
    )
    results = _compile_all_shells(field_spec, compiler_input, initial_bond_shells)
    return PreparedModel(
        field_spec=field_spec,
        shell_results=results,
        basis_labels=labels,
        compilation_recipe=_CompilationRecipe(
            "explicit", field_spec._context(), compiler_input
        ),
    )


def _prepare_explicit_model(
    *,
    lattice: MatrixInput,
    lattpar: Any,
    wyckoffposition: Any,
    symminformation: Any,
    basis_functions: Any,
    initial_bond_shells: int,
    generate_symmetry_group: bool,
    representation_mode: str,
    exact_field: ExactField,
    debug_q: bool,
) -> PreparedModel:
    compiler_input = _base_compiler_input(
        lattice,
        lattpar,
        wyckoffposition,
        initial_bond_shells,
        representation_mode,
        debug_q,
    )
    basis = _normalise_basis(basis_functions, len(_normalise_wyckoff(wyckoffposition)))
    compiler_input["entries"].extend(
        [
            {"key": "SymmetryInformation", "value": _encode_symmetry(symminformation)},
            {"key": "BasisFunctions", "value": _encode_basis_functions(basis)},
            {"key": "GenerateSymmetryGroup", "value": _strict_bool(generate_symmetry_group, "generate_symmetry_group")},
        ]
    )
    results = _compile_all_shells(exact_field, compiler_input, initial_bond_shells)
    return PreparedModel(
        field_spec=exact_field,
        shell_results=results,
        basis_labels=_basis_display_labels(basis),
        compilation_recipe=_CompilationRecipe(
            "explicit", exact_field._context(), compiler_input
        ),
    )


def _prepare_data_model(
    *,
    lattice: MatrixInput,
    lattpar: Any,
    wyckoff: Sequence[Wyckoff],
    msg: MSG,
    basis_functions: Any,
    initial_bond_shells: int,
    representation_mode: str,
    exact_field: ExactField,
    debug_q: bool,
) -> PreparedModel:
    compiler_input = _base_compiler_input(
        lattice,
        lattpar,
        wyckoff,
        initial_bond_shells,
        representation_mode,
        debug_q,
    )
    basis = _normalise_basis(basis_functions, len(wyckoff))
    compiler_input["entries"].append(
        {
            "key": "BasisFunctions",
            "value": _encode_basis_functions(basis),
        }
    )
    return _compile_data_prepared_model(
        field_spec=exact_field,
        compiler_input=compiler_input,
        msg=msg,
        wyckoff=wyckoff,
        initial_bond_shells=initial_bond_shells,
        basis_labels=_basis_display_labels(basis),
    )


def _compile_data_prepared_model(
    *,
    field_spec: ExactField,
    compiler_input: Dict[str, Any],
    msg: MSG,
    wyckoff: Sequence[Wyckoff],
    initial_bond_shells: int,
    basis_labels: Sequence[Sequence[str]],
) -> PreparedModel:
    catalog = DataCatalog()
    msg_id = msg.stable_id or catalog.resolve_msg_id(bns=msg.bns)  # type: ignore[arg-type]
    resolved = [
        catalog.resolve_wyckoff(
            msg_id,
            selection.letter,
            selection.source_ordinal,
        )
        for selection in wyckoff
    ]
    selections = [
        _ordered(
            [
                ("source_ordinal_1based", int(entry["source_ordinal"])),
                ("letter", selection.letter),
            ]
        )
        for selection, entry in zip(wyckoff, resolved)
    ]
    trace_entries: list[Tuple[str, Any]] = [
        ("msg_stable_id", msg_id),
        ("wyckoff_selections", _tagged_list(selections)),
    ]
    if len(selections) == 1:
        trace_entries.extend(
            [
                (
                    "wyckoff_source_ordinal_1based",
                    int(resolved[0]["source_ordinal"]),
                ),
                ("wyckoff_letter", wyckoff[0].letter),
            ]
        )
    trace = _ordered(trace_entries)
    payload = catalog.compile_model_input_closures(
        field_spec._context(), compiler_input, trace
    )
    results = payload.get("shell_results")
    if not isinstance(results, list) or len(results) != initial_bond_shells:
        raise RuntimeError("Rust returned an invalid batch of Data shell results")
    return PreparedModel(
        field_spec=field_spec,
        shell_results=results,
        basis_labels=basis_labels,
        compilation_recipe=_CompilationRecipe(
            "data", field_spec._context(), compiler_input, trace
        ),
    )


def _compile_all_shells(
    field_spec: ExactField,
    compiler_input: Dict[str, Any],
    shell_count: int,
) -> list[Dict[str, Any]]:
    payload = geometry(
        "compile_model_input_closures",
        context=field_spec._context(),
        compiler_input=compiler_input,
    )
    results = payload.get("shell_results")
    if not isinstance(results, list) or len(results) != shell_count:
        raise RuntimeError("Rust returned an invalid batch of prepared shell results")
    return results


def _base_compiler_input(
    lattice: MatrixInput,
    lattpar: Any,
    wyckoffposition: Any,
    initial_bond_shells: int,
    representation_mode: str,
    debug_q: bool,
) -> Dict[str, Any]:
    if (
        isinstance(initial_bond_shells, bool)
        or not isinstance(initial_bond_shells, int)
        or initial_bond_shells < 1
    ):
        raise TypeError("initial_bond_shells must be a positive integer")
    if representation_mode not in ("DirectProduct", "Induced"):
        raise ValueError("representation_mode must be 'DirectProduct' or 'Induced'")
    entries = [
        ("Lattice", _matrix(lattice)),
        ("LatticeParameters", _encode_lattpar(lattpar)),
        (
            "WyckoffPosition",
            _tagged_list([_matrix([site[0], site[1]]) for site in _normalise_wyckoff(wyckoffposition)]),
        ),
        ("Debug", debug_q),
        ("InitialBondShells", initial_bond_shells),
        ("RepresentationMode", representation_mode),
    ]
    if representation_mode == "Induced":
        # Stable 2.0.10 compiles the ordinary induced-basis route from the
        # internal Automatic marker.  It is intentionally not a public Python
        # or Web option; Rust remains the sole owner of site induction.
        entries.append(
            ("SiteLocalData", {"kind": "symbol", "name": "System`Automatic"})
        )
    return _ordered(entries)
