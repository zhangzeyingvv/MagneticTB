"""Model Session for the Rust-backed model lifecycle."""

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
    ExactAtom,
    ExactExpr,
    ExactField,
    Hamiltonian,
    HamiltonianExpression,
    HamiltonianSpace,
    MatrixInput,
    ModelError,
    sqrt,
    symbol,
)

from .model_schema import (
    PreparedModel,
    _broken_symmetry_rules_for_model,
    _symham_ii_for_model,
)

from .model_preparation import (
    prepare_model,
    prepare_model_from_rep,
)


_CURRENT_MODEL: ContextVar[Optional[PreparedModel]] = ContextVar(
    "magnetictb_current_model", default=None
)


def _stable_default_lattice() -> MatrixInput:
    a = symbol("a")
    c = symbol("c")
    return (
        (a, 0, 0),
        (-a / 2, sqrt(3) * a / 2, 0),
        (0, 0, c),
    )


def _stable_default_lattpar() -> Tuple[Tuple[str, ExactAtom], ...]:
    return (("a", 1), ("c", 3))


def _stable_default_wyckoff() -> Tuple[Tuple[Tuple[ExactAtom, ...], Tuple[ExactAtom, ...]], ...]:
    return (((Fraction(2, 3), Fraction(1, 3), 0), (0, 0, Fraction(1, 2))),)


def _stable_default_symmetry() -> Tuple[Tuple[Any, ...], ...]:
    return (("1", ((1, 0, 0), (0, 1, 0), (0, 0, 1)), (0, 0, 0), "F"),)


def _print_generators(model: PreparedModel) -> None:
    generators = [
        (
            label,
            "T" if antiunitary else "F",
        )
        for label, antiunitary in model._generator_display_records
    ]
    print("Generators:", generators)


def current_model() -> PreparedModel:
    """Return the current stable-style model session."""

    model = _CURRENT_MODEL.get()
    if model is None:
        raise ModelError("ModelNotInitialized", "call init or initfromrep first")
    return model


def CurrentModelSession() -> PreparedModel:
    """Return the read-only session prepared by ``init`` or ``initfromrep``.

    The returned :class:`PreparedModel` supports the stable 2.0.10 Association
    keys through ``model["ModelSpecification", "SiteOrbits"]`` while retaining
    the typed Python properties used by the ordinary API.
    """

    return current_model()


def CompileMagneticTBInput(input: Mapping[str, Any]) -> PreparedModel:
    """Compile one stable 2.0.10 init-input mapping without installing it.

    Mathematica returns an Association.  Python returns the corresponding
    immutable :class:`PreparedModel`, which exposes the same major compiled
    sections through stable string keys and keeps the Rust-generated canonical
    payload available through :meth:`PreparedModel.to_canonical_dict`.
    """

    if not isinstance(input, Mapping):
        raise ModelError(
            "InvalidCompilerInput",
            "CompileMagneticTBInput requires a mapping of stable init fields",
        )
    values = dict(input)
    allowed = {
        "Lattice",
        "LatticeParameters",
        "WyckoffPosition",
        "SymmetryInformation",
        "BasisFunctions",
        "RepresentationSource",
        "RepresentationInformation",
        "OrbitalLabels",
        "Debug",
        "InitialBondShells",
        "GenerateSymmetryGroup",
        "RepresentationMode",
        "SiteLocalData",
    }
    unknown = sorted(set(values) - allowed)
    if unknown:
        raise ModelError(
            "InvalidCompilerInput",
            "unknown stable compiler field(s): " + ", ".join(unknown),
        )
    missing = [
        name
        for name in (
            "Lattice",
            "LatticeParameters",
            "WyckoffPosition",
            "SymmetryInformation",
        )
        if name not in values
    ]
    if missing:
        raise ModelError(
            "InvalidCompilerInput",
            "missing stable compiler field(s): " + ", ".join(missing),
        )
    site_local_data = values.get("SiteLocalData", "Automatic")
    if site_local_data not in (None, "Automatic"):
        raise ModelError(
            "UnsupportedSiteLocalData",
            "explicit SiteLocalData is not part of the Python 2.0.10 input surface",
        )
    source = values.get("RepresentationSource", "BasisFunctions")
    common = {
        "lattice": values["Lattice"],
        "lattpar": values["LatticeParameters"],
        "wyckoffposition": values["WyckoffPosition"],
        "symminformation": values["SymmetryInformation"],
        "initial_bond_shells": values.get("InitialBondShells", 10),
        "representation_mode": values.get("RepresentationMode", "DirectProduct"),
        "debug_q": values.get("Debug", False),
    }
    if source == "BasisFunctions":
        if "BasisFunctions" not in values:
            raise ModelError(
                "InvalidCompilerInput",
                "BasisFunctions is required when RepresentationSource is BasisFunctions",
            )
        return prepare_model(
            **common,
            basis_functions=values["BasisFunctions"],
            generate_symmetry_group=values.get("GenerateSymmetryGroup", False),
        )
    if source == "Matrices":
        if "RepresentationInformation" not in values:
            raise ModelError(
                "InvalidCompilerInput",
                "RepresentationInformation is required when RepresentationSource is Matrices",
            )
        return prepare_model_from_rep(
            **common,
            repinformation=values["RepresentationInformation"],
            orbital_labels=values.get("OrbitalLabels"),
        )
    raise ModelError(
        "InvalidRepresentationSource",
        "RepresentationSource must be 'BasisFunctions' or 'Matrices'",
    )


def compile_magnetic_tb_input(input: Mapping[str, Any]) -> PreparedModel:
    """Snake-case alias of :func:`CompileMagneticTBInput`."""

    return CompileMagneticTBInput(input)


def _stable_keyword_alias(
    options: Dict[str, Any],
    stable_name: str,
    current: Any,
    default: Any,
    domain: str,
) -> Any:
    if stable_name not in options:
        return current
    value = options.pop(stable_name)
    if current != default and current != value:
        raise ModelError(
            f"Duplicate{domain}Option",
            f"both the Python spelling and stable spelling {stable_name} were supplied",
        )
    return value


def init(
    *,
    lattice: Optional[MatrixInput] = None,
    lattpar: Mapping[str | ExactExpr, ExactAtom]
    | Sequence[Tuple[str | ExactExpr, ExactAtom]]
    | None = None,
    wyckoffposition: Any = None,
    symminformation: Any = None,
    basis_functions: Any = ("s",),
    debug_q: bool = False,
    initial_bond_shells: int = 10,
    generate_symmetry_group: bool = False,
    representation_mode: str = "DirectProduct",
    exact_field: Optional[ExactField] = None,
    **stable_options: Any,
) -> None:
    """Prepare and install the current model, matching stable ``init`` lifecycle."""

    _CURRENT_MODEL.set(None)
    basis_functions = _stable_keyword_alias(
        stable_options, "basisFunctions", basis_functions, ("s",), "Init"
    )
    debug_q = _stable_keyword_alias(stable_options, "debugQ", debug_q, False, "Init")
    initial_bond_shells = _stable_keyword_alias(
        stable_options, "InitialBondShells", initial_bond_shells, 10, "Init"
    )
    generate_symmetry_group = _stable_keyword_alias(
        stable_options,
        "GenerateSymmetryGroup",
        generate_symmetry_group,
        False,
        "Init",
    )
    representation_mode = _stable_keyword_alias(
        stable_options,
        "RepresentationMode",
        representation_mode,
        "DirectProduct",
        "Init",
    )
    if stable_options:
        names = ", ".join(sorted(stable_options))
        raise ModelError("UnknownInitOption", f"unknown init option(s): {names}")
    model = prepare_model(
        lattice=_stable_default_lattice() if lattice is None else lattice,
        lattpar=_stable_default_lattpar() if lattpar is None else lattpar,
        wyckoffposition=_stable_default_wyckoff() if wyckoffposition is None else wyckoffposition,
        symminformation=_stable_default_symmetry() if symminformation is None else symminformation,
        basis_functions=basis_functions,
        debug_q=debug_q,
        initial_bond_shells=initial_bond_shells,
        generate_symmetry_group=generate_symmetry_group,
        representation_mode=representation_mode,
        exact_field=exact_field,
    )
    _CURRENT_MODEL.set(model)
    _print_generators(model)


def initfromrep(
    *,
    lattice: Optional[MatrixInput] = None,
    lattpar: Mapping[str | ExactExpr, ExactAtom]
    | Sequence[Tuple[str | ExactExpr, ExactAtom]]
    | None = None,
    wyckoffposition: Any = None,
    symminformation: Any = None,
    repinformation: Any = None,
    orbital_labels: Optional[Sequence[Sequence[str]]] = None,
    debug_q: bool = False,
    initial_bond_shells: int = 10,
    representation_mode: str = "DirectProduct",
    exact_field: Optional[ExactField] = None,
    **stable_options: Any,
) -> None:
    """Prepare and install an explicit-representation current model."""

    _CURRENT_MODEL.set(None)
    orbital_labels = _stable_keyword_alias(
        stable_options, "orbitalLabels", orbital_labels, None, "InitFromRep"
    )
    debug_q = _stable_keyword_alias(
        stable_options, "debugQ", debug_q, False, "InitFromRep"
    )
    initial_bond_shells = _stable_keyword_alias(
        stable_options,
        "InitialBondShells",
        initial_bond_shells,
        10,
        "InitFromRep",
    )
    representation_mode = _stable_keyword_alias(
        stable_options,
        "RepresentationMode",
        representation_mode,
        "DirectProduct",
        "InitFromRep",
    )
    if stable_options:
        names = ", ".join(sorted(stable_options))
        raise ModelError(
            "UnknownInitFromRepOption",
            "unknown initfromrep option(s): "
            f"{names}; basisFunctions, SiteLocalData, and GenerateSymmetryGroup "
            "are intentionally unsupported",
        )
    model = prepare_model_from_rep(
        lattice=_stable_default_lattice() if lattice is None else lattice,
        lattpar=_stable_default_lattpar() if lattpar is None else lattpar,
        wyckoffposition=_stable_default_wyckoff() if wyckoffposition is None else wyckoffposition,
        symminformation=_stable_default_symmetry() if symminformation is None else symminformation,
        repinformation=((((1,),),),) if repinformation is None else repinformation,
        orbital_labels=orbital_labels,
        debug_q=debug_q,
        initial_bond_shells=initial_bond_shells,
        representation_mode=representation_mode,
        exact_field=exact_field,
    )
    _CURRENT_MODEL.set(model)
    _print_generators(model)


def symham(
    shell: int,
    *,
    hermitian: bool = True,
    cartesian_coordinates: bool = False,
    validation_level: Literal["None", "Basic", "Full"] = "Basic",
    kernel_method: Literal["Iterative", "Stacked", "Cyclotomic"] = "Iterative",
    **stable_options: Any,
) -> list[list[HamiltonianExpression]]:
    """Return the current session's stable-style symbolic Hamiltonian matrix."""

    hamiltonian = hamiltonian_object(
        shell,
        hermitian=hermitian,
        cartesian_coordinates=cartesian_coordinates,
        validation_level=validation_level,
        kernel_method=kernel_method,
        **stable_options,
    )
    print("params:", list(hamiltonian.parameter_names))
    return hamiltonian.to_list()


def unsymham(shell: int) -> list[list[HamiltonianExpression]]:
    """Build the unconstrained symbolic Hamiltonian for the current session."""

    hamiltonian = current_model().unconstrained_hamiltonian(shell)
    print("params:", list(hamiltonian.parameter_names))
    return hamiltonian.to_list()


def symhamII(
    ham: Union[Hamiltonian, Sequence[Sequence[HamiltonianExpression]]],
    *,
    wcc: Optional[Sequence[Sequence[ExactAtom]]] = None,
) -> list[list[HamiltonianExpression]]:
    """Rewrite ``ham`` using stable convention-II Wannier-center phases."""

    return _symham_ii_for_model(current_model(), ham, wcc=wcc)


def brokenSymmetryInitRules(selector: Any) -> Dict[str, Any]:
    """Split current Wyckoff orbits under a retained stable subgroup."""

    return _broken_symmetry_rules_for_model(current_model(), selector)


def hamiltonian_object(
    shell: int,
    *,
    hermitian: bool = True,
    cartesian_coordinates: bool = False,
    validation_level: Literal["None", "Basic", "Full"] = "Basic",
    kernel_method: Literal["Iterative", "Stacked", "Cyclotomic"] = "Iterative",
    **stable_options: Any,
) -> Hamiltonian:
    """Return the advanced immutable Rust Hamiltonian object for one shell."""

    hermitian = _stable_keyword_alias(
        stable_options, "Hermitian", hermitian, True, "Symham"
    )
    cartesian_coordinates = _stable_keyword_alias(
        stable_options,
        "CartesianCoordinates",
        cartesian_coordinates,
        False,
        "Symham",
    )
    validation_level = _stable_keyword_alias(
        stable_options, "ValidationLevel", validation_level, "Basic", "Symham"
    )
    kernel_method = _stable_keyword_alias(
        stable_options, "KernelMethod", kernel_method, "Iterative", "Symham"
    )
    if stable_options:
        names = ", ".join(sorted(stable_options))
        raise ModelError("UnknownSymhamOption", f"unknown symham option(s): {names}")
    return current_model().hamiltonian(
        shell,
        hermitian=hermitian,
        cartesian_coordinates=cartesian_coordinates,
        validation_level=validation_level,
        kernel_method=kernel_method,
    )


def hamiltonian_space(
    shell: int,
    *,
    hermitian: bool = True,
    cartesian_coordinates: bool = False,
    validation_level: Literal["None", "Basic", "Full"] = "Basic",
    kernel_method: Literal["Iterative", "Stacked", "Cyclotomic"] = "Iterative",
) -> HamiltonianSpace:
    return current_model().hamiltonian_space(
        shell,
        hermitian=hermitian,
        cartesian_coordinates=cartesian_coordinates,
        validation_level=validation_level,
        kernel_method=kernel_method,
    )


def _clear_current_model_for_tests() -> None:
    _CURRENT_MODEL.set(None)
