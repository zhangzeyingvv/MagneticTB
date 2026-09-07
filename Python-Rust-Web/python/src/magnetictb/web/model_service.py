"""Prepared-model and Hamiltonian services shared by the HTTP routes."""

from __future__ import annotations

import io
from contextlib import redirect_stdout
from typing import Any, Dict, List, Mapping

from ..errors import ModelError
from ..model import (
    PreparedModel,
    _CURRENT_MODEL,
    current_model,
    init,
    initfromrep,
    prepare_model,
    prepare_model_from_rep,
)
from ..properties import hoppingData
from ..tight_binding import Hamiltonian
from ..utilities import texOutput
from .codec import decode_tagged, guard_json_tree
from .model_schemas import HamiltonianOptions, PrepareRequest
from .properties_schemas import SolvedShellRequest
from .registry import ModelEntry


def _model_summary(entry: ModelEntry, *, include_records: bool = False) -> Dict[str, Any]:
    model = entry.model
    result: Dict[str, Any] = {
        "model_id": entry.model_id,
        "label": entry.label,
        "model_identity_sha256": model.model_identity_sha256,
        "initial_bond_shells": model.initial_bond_shells,
        "group_order": int(model.group["order"]),
        "operation_count": len(model.symmetry_operations),
        "generator_indices": list(model.generator_indices),
        "orbital_count": len(model.orbital_table),
        "bond_shell_counts": [len(shell) for shell in model.bond_shells],
    }
    if include_records:
        result.update(
            {
                "group": model.group,
                "site_orbits": model.site_orbits,
                "site_actions": model.site_actions,
                "cell_translations": model.cell_translations,
                "orbitals": model.orbital_table,
            }
        )
    return result


def _prepare(request: PrepareRequest) -> tuple[PreparedModel, str]:
    guard_json_tree(request.arguments)
    arguments = decode_tagged(request.arguments)
    if not isinstance(arguments, dict):
        raise TypeError("arguments must decode to an object")
    output = io.StringIO()
    token = _CURRENT_MODEL.set(None)
    try:
        with redirect_stdout(output):
            if request.entrypoint == "init":
                init(**arguments)
                model = current_model()
            elif request.entrypoint == "initfromrep":
                initfromrep(**arguments)
                model = current_model()
            elif request.entrypoint == "prepare_model":
                model = prepare_model(**arguments)
            else:
                model = prepare_model_from_rep(**arguments)
    finally:
        _CURRENT_MODEL.reset(token)
    return model, output.getvalue()


def _hamiltonian(model: PreparedModel, options: HamiltonianOptions) -> Hamiltonian:
    return model.hamiltonian(
        options.shell,
        hermitian=options.hermitian,
        cartesian_coordinates=options.cartesian_coordinates,
        validation_level=options.validation_level,
        kernel_method=options.kernel_method,
    )


def _cumulative_hamiltonian(
    model: PreparedModel, options: HamiltonianOptions
) -> Hamiltonian:
    values = [
        model.hamiltonian(
            shell,
            hermitian=options.hermitian,
            cartesian_coordinates=options.cartesian_coordinates,
            validation_level=options.validation_level,
            kernel_method=options.kernel_method,
        )
        for shell in range(1, options.shell + 1)
    ]
    result = values[0]
    for value in values[1:]:
        result = result + value
    return result


def _property_shells(request: SolvedShellRequest, model: PreparedModel) -> List[int]:
    shells: List[int] = []
    for shell in request.shells:
        if isinstance(shell, bool) or shell < 1 or shell > model.initial_bond_shells:
            raise ModelError(
                "ShellNotPrepared",
                f"shell {shell} was not prepared; available shells are 1..{model.initial_bond_shells}",
            )
        if shell in shells:
            raise ModelError(
                "DuplicateHamiltonianShell", "one shell may occur only once"
            )
        shells.append(shell)
    return shells


def _property_parameters(request: SolvedShellRequest) -> Dict[str, Any]:
    guard_json_tree(request.parameters)
    decoded = decode_tagged(request.parameters)
    if not isinstance(decoded, dict):
        raise TypeError("parameters must decode to an object")
    return decoded


def _property_hamiltonian(
    model: PreparedModel, request: SolvedShellRequest
) -> Hamiltonian:
    values = [
        model.hamiltonian(
            shell,
            hermitian=request.hermitian,
            validation_level=request.validation_level,
            kernel_method=request.kernel_method,
        )
        for shell in _property_shells(request, model)
    ]
    result = values[0]
    for value in values[1:]:
        result = result + value
    return result


def _property_hopping_data(
    model: PreparedModel, request: SolvedShellRequest
) -> Dict[str, Any]:
    shells = _property_shells(request, model)
    # The browser action explicitly evaluates the requested stable symham shells.
    # The public Python hoppingData function itself remains cache-only, matching 2.0.10.
    _property_hamiltonian(model, request)
    token = _CURRENT_MODEL.set(model)
    try:
        return hoppingData(
            shells,
            _property_parameters(request),
            Hermitian=request.hermitian,
            ValidationLevel=request.validation_level,
            KernelMethod=request.kernel_method,
        )
    finally:
        _CURRENT_MODEL.reset(token)


def _property_centers(data: Mapping[str, Any]) -> Any:
    centers = data.get("WannierCenters")
    if centers is None:
        raise ModelError(
            "MissingPositionData",
            "the selected model does not provide one Wannier center per orbital",
        )
    return centers


def _hamiltonian_payload(value: Hamiltonian) -> Dict[str, Any]:
    result = value.to_dict()
    result["matrix_text"] = [[str(item) for item in row] for row in value.matrix]
    for result_term, term in zip(result["terms"], value.terms):
        for result_record, (displacement, matrix) in zip(
            result_term["fourier_coefficients"], term.fourier_coefficients
        ):
            result_record["displacement_text"] = [
                str(component) for component in displacement
            ]
            result_record["matrix_text"] = [
                [str(coefficient) for coefficient in row] for row in matrix.rows
            ]
    result["canonical_sha256"] = value.model_identity_sha256
    result["tex"] = texOutput(value)
    return result
