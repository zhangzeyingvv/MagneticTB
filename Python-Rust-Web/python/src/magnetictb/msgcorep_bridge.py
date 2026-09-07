"""Stable-shaped bridge to the separately installed Python MSGCorep package.

MSGCorep remains the authority for its database and corepresentation labels.
MagneticTB's Rust core remains the authority for the tight-binding model,
Hamiltonian evaluation, eigenspaces, little group, and character traces.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from fractions import Fraction
from importlib import import_module
from math import isfinite
from numbers import Number
from typing import Any

from .api import ExactError, ModelError, geometry
from .linear_algebra import (
    ExactExpr,
    _encode_evaluation_scalar,
    _encode_exact,
    _matrix,
    _thaw_json,
)
from .model import PreparedModel, current_model
from .tight_binding import Hamiltonian


def _msgcorep() -> Any:
    try:
        return import_module("msgcorep")
    except ImportError as error:
        raise ModelError(
            "MSGCorepNotInstalled",
            "install the Python MSGCorep package before calling this bridge",
        ) from error


def _bns_identifier(value: Any) -> tuple[int, int]:
    if (
        not isinstance(value, Sequence)
        or isinstance(value, (str, bytes))
        or len(value) != 2
        or any(isinstance(item, bool) or not isinstance(item, int) for item in value)
    ):
        raise TypeError("MSG must be a two-integer BNS identifier")
    return int(value[0]), int(value[1])


def _msgcorep_failure(module: Any, operation: str, call: Any) -> Any:
    try:
        return call()
    except getattr(module, "MSGCorepError", ()) as error:
        raise ModelError("MSGCorepFailure", f"{operation}: {error}") from error


def _raw_operations(
    module: Any, identifier: tuple[int, int]
) -> tuple[str, tuple[tuple[str, tuple[tuple[int, ...], ...], tuple[Any, ...], bool], ...]]:
    bravais = _msgcorep_failure(
        module, "getSGLatt", lambda: module.getSGLatt(identifier[0])
    )
    elements = _msgcorep_failure(
        module, "getMSGElem", lambda: module.getMSGElem(identifier)
    )
    result = []
    for label, translation, antiunitary in elements:
        rotation = _msgcorep_failure(
            module,
            "getRotMat",
            lambda label=label: module.getRotMat(
                bravais, str(label).removeprefix("bar")
            ),
        )
        if (
            not isinstance(rotation, Sequence)
            or len(rotation) != 3
            or any(
                not isinstance(row, Sequence)
                or len(row) != 3
                or any(isinstance(value, bool) or not isinstance(value, int) for value in row)
                for row in rotation
            )
        ):
            raise ModelError(
                "InvalidMSGCorepRotation",
                f"MSGCorep rotation {label} is not an exact integer 3 by 3 matrix",
            )
        result.append(
            (
                str(label),
                tuple(tuple(value for value in row) for row in rotation),
                tuple(translation),
                bool(antiunitary),
            )
        )
    return bravais, tuple(result)


def getMSGElemFromMSGCorep(MSG: Any) -> list[list[Any]]:
    """Return MSGCorep operations in the original stable MagneticTB record shape."""

    identifier = _bns_identifier(MSG)
    module = _msgcorep()
    bravais, operations = _raw_operations(module, identifier)
    primitive = _msgcorep_failure(
        module, "BasicVectors", lambda: module.BasicVectors[bravais]
    )
    print("Primitive Lattice Vactor:", primitive)
    return [
        [
            label,
            [list(row) for row in rotation],
            list(translation),
            "T" if antiunitary else "F",
        ]
        for label, rotation, translation, antiunitary in operations
    ]


def _parameter_mapping(param: Any) -> dict[str, Any]:
    if isinstance(param, Mapping):
        items = tuple(param.items())
    elif isinstance(param, Sequence) and not isinstance(param, (str, bytes)):
        items = tuple(param)
        if any(
            not isinstance(item, Sequence)
            or isinstance(item, (str, bytes))
            or len(item) != 2
            for item in items
        ):
            raise TypeError("param must be a mapping or an ordered sequence of pairs")
    else:
        raise TypeError("param must be a mapping or an ordered sequence of pairs")
    result: dict[str, Any] = {}
    for key, value in items:
        name = str(key)
        if name in result:
            raise ModelError(
                "DuplicateHamiltonianParameter", f"parameter {name} occurs more than once"
            )
        result[name] = _encode_evaluation_scalar(value, f"parameter {name}")
    return result


def _kset(value: Any) -> tuple[tuple[Any, Any, Any], ...]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or not value:
        raise TypeError("kset must be a nonempty ordered sequence of three-vectors")
    points = []
    for index, point in enumerate(value, start=1):
        if (
            not isinstance(point, Sequence)
            or isinstance(point, (str, bytes))
            or len(point) != 3
        ):
            raise TypeError(f"kset point {index} must contain exactly three coordinates")
        if any(
            isinstance(coordinate, (bool, complex))
            or not isinstance(coordinate, Number)
            for coordinate in point
        ):
            raise TypeError(
                f"kset point {index} must use real Python numeric coordinates"
            )
        encoded = tuple(
            _encode_evaluation_scalar(coordinate, f"kset point {index}")
            for coordinate in point
        )
        points.append((encoded[0], encoded[1], encoded[2]))
    return tuple(points)


def _tracked_hamiltonian(model: PreparedModel, ham: Any) -> Hamiltonian:
    if isinstance(ham, Hamiltonian):
        if ham.model_identity_sha256 != model.model_identity_sha256:
            raise ModelError(
                "BandCorepModelMismatch",
                "Hamiltonian was not generated by the current prepared model",
            )
        return ham
    if isinstance(ham, Sequence) and not isinstance(ham, (str, bytes)):
        for shell in range(1, model.initial_bond_shells + 1):
            candidate = model.hamiltonian(shell)
            if candidate.to_list() == ham:
                return candidate
        raise ModelError(
            "UntrackedHamiltonianMatrix",
            "the stable-style matrix must be an unchanged symham result from the current model",
        )
    raise TypeError("ham must be a Hamiltonian or the matrix returned by symham")


def _numeric_complex(value: Any, evaluate_numeric: Any) -> complex:
    if isinstance(value, Number) and not isinstance(value, bool):
        result = complex(value)
    else:
        result = complex(evaluate_numeric(value))
    if not (isfinite(result.real) and isfinite(result.imag)):
        raise ModelError("NonFiniteMSGCorepSpinRotation", "spin rotation is not finite")
    return result


def _canonical_spin_rotations(module: Any, operations: Sequence[tuple[Any, ...]]) -> list[Any]:
    try:
        evaluate_numeric = import_module("msgcorep.core").evaluate_numeric
    except (ImportError, AttributeError) as error:
        raise ModelError(
            "InvalidMSGCorepInterface",
            "MSGCorep does not expose its exact numeric evaluator",
        ) from error
    result = []
    for label, _, _, _ in operations:
        try:
            matrix = module.getSpinRotOp[str(label).removeprefix("bar")][0]
        except (KeyError, IndexError, TypeError) as error:
            raise ModelError(
                "MSGCorepSpinRotationUnavailable",
                f"MSGCorep has no spin rotation for {label}",
            ) from error
        numeric = [
            [_numeric_complex(entry, evaluate_numeric) for entry in row] for row in matrix
        ]
        if len(numeric) != 2 or any(len(row) != 2 for row in numeric):
            raise ModelError(
                "InvalidMSGCorepSpinRotation", f"spin rotation {label} is not 2 by 2"
            )
        # Stable Corepresentation.wl uses PauliMatrix[1].spin.PauliMatrix[1].
        swapped = ((numeric[1][1], numeric[1][0]), (numeric[0][1], numeric[0][0]))
        result.append(
            [
                [
                    {"real": value.real, "imaginary": value.imag}
                    for value in row
                ]
                for row in swapped
            ]
        )
    return result


def _plain_exact(value: Any) -> Any:
    if isinstance(value, (int, Fraction, ExactExpr)) and not isinstance(value, bool):
        return _encode_exact(value)
    raise TypeError(f"unsupported MSGCorep exact scalar {type(value).__name__}")


def _complex_trace(value: Mapping[str, Any]) -> complex:
    return complex(float(value["real"]), float(value["imaginary"]))


def getTBBandCorep(MSG: Any, ham: Any, param: Any, kset: Any) -> Any:
    """Identify MSGCorep band corepresentations for a MagneticTB Hamiltonian."""

    identifier = _bns_identifier(MSG)
    module = _msgcorep()
    _, operations = _raw_operations(module, identifier)
    model = current_model()
    hamiltonian = _tracked_hamiltonian(model, ham)
    encoded_kset = _kset(kset)
    flags = model.model_identity_payload.get("spinor_basis_flags")
    if flags is None:
        raise ModelError(
            "BandCorepSpinMetadataUnavailable",
            "getTBBandCorep requires a model prepared from catalog or explicit polynomial basis functions",
        )
    flags = tuple(bool(value) for value in flags)
    if not flags or any(value != flags[0] for value in flags):
        raise ModelError(
            "MixedBandCorepSpinBasis",
            "all Wyckoff orbits must consistently use scalar or spinor bases",
        )
    soc = flags[0]

    operation_payload = [
        {
            "rotation": _matrix(rotation),
            "translation": [_plain_exact(value) for value in translation],
            "antiunitary": antiunitary,
        }
        for _, rotation, translation, antiunitary in operations
    ]
    arguments = {
        "hamiltonian": hamiltonian.to_canonical_dict(),
        "parameters": _parameter_mapping(param),
        "kset": [list(point) for point in encoded_kset],
        "model_identity": _thaw_json(model.model_identity_payload),
        "msgcorep_operations": operation_payload,
        "soc": soc,
    }
    if soc:
        arguments["msgcorep_spin_rotations"] = _canonical_spin_rotations(
            module, operations
        )
    try:
        trace = geometry("band_corep_trace", **arguments)
    except ExactError as error:
        raise ModelError(error.tag, error.detail) from None

    rotations = [[list(row) for row in operation[1]] for operation in operations]
    translations = [list(operation[2]) for operation in operations]
    trace_data = {
        "nelec": hamiltonian.shape[0],
        "soc": int(soc),
        "nsym": len(operations),
        "rot": rotations,
        "trans": translations,
        "srot": [
            [
                [_complex_trace(value) for value in row]
                for row in matrix
            ]
            for matrix in trace["spin_rotations"]
        ],
        "unitary": [1 if not operation[3] else -1 for operation in operations],
        "nk": len(encoded_kset),
        "kpt": [tuple(point) for point in kset],
        "nband": hamiltonian.shape[0],
        "ene": trace["energies"],
        "deg": trace["degeneracies"],
        "knsym": trace["little_symmetry_counts"],
        "kisym": trace["little_symmetry_indices"],
        "trace": [
            [
                [_complex_trace(value) for value in band]
                for band in point
            ]
            for point in trace["traces"]
        ],
    }
    return _msgcorep_failure(
        module,
        "getBandCorep",
        lambda: module.getBandCorep(identifier, trace_data),
    )


__all__ = ["getMSGElemFromMSGCorep", "getTBBandCorep"]
