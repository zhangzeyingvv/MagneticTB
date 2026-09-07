"""Model Routes for the local MagneticTB workbench."""

from __future__ import annotations

from fastapi import APIRouter

from .route_dependencies import (
    Any,
    CombineRequest,
    Dict,
    EvaluateRequest,
    FriendlyInitRequest,
    Hamiltonian,
    HamiltonianBasisReportRequest,
    HamiltonianOptions,
    HoppingReportRequest,
    JSONResponse,
    Optional,
    PrepareRequest,
    Query,
    ShellReportRequest,
    SolvedShellRequest,
    SymmetryReportRequest,
    _REGISTRY,
    _bounded,
    _compute,
    _cumulative_hamiltonian,
    _friendly_init,
    _hamiltonian,
    _hamiltonian_payload,
    _model_summary,
    _prepare,
    _property_hopping_data,
    bondTable,
    decode_tagged,
    guard_json_tree,
    orbitalTable,
    partial,
    showHamiltonianBasis,
    showHoppingParameters,
    showSymmetryRepresentations,
)


router = APIRouter()


@router.post("/api/models/init", status_code=201)
async def create_init_model(payload: FriendlyInitRequest) -> JSONResponse:
    """Ordinary form-oriented ``init`` endpoint without raw canonical JSON."""

    model, messages = await _compute(partial(_friendly_init, payload))
    entry = _REGISTRY.put(model, payload.label)
    result = _model_summary(entry, include_records=True)
    result["messages"] = [line for line in messages.splitlines() if line]
    result["entrypoint"] = "init"
    return _bounded(result, status_code=201)


@router.get("/api/models")
async def list_models() -> Dict[str, Any]:
    return {"models": [_model_summary(entry) for entry in _REGISTRY.list()]}


@router.post("/api/models", status_code=201)
async def create_model(payload: PrepareRequest) -> JSONResponse:
    model, messages = await _compute(partial(_prepare, payload))
    entry = _REGISTRY.put(model, payload.label)
    result = _model_summary(entry, include_records=True)
    result["messages"] = [line for line in messages.splitlines() if line]
    return _bounded(result, status_code=201)


@router.get("/api/models/{model_id}")
async def get_model(model_id: str, details: bool = False) -> JSONResponse:
    return _bounded(_model_summary(_REGISTRY.get(model_id), include_records=details))


@router.delete("/api/models/{model_id}")
async def delete_model(model_id: str) -> Dict[str, Any]:
    if not _REGISTRY.delete(model_id):
        raise ModelError("WebModelNotFound", "the model id is unknown or expired")
    return {"deleted": True, "model_id": model_id}


@router.get("/api/models/{model_id}/symmetry-operations")
async def model_symmetry_operations(model_id: str) -> JSONResponse:
    model = _REGISTRY.get(model_id).model
    return _bounded({"operations": model.symmetry_operations})


@router.get("/api/models/{model_id}/orbitals")
async def model_orbitals(model_id: str) -> JSONResponse:
    return _bounded({"orbitals": _REGISTRY.get(model_id).model.orbital_table})


@router.get("/api/models/{model_id}/reports/orbital-table")
async def model_orbital_table_report(model_id: str) -> JSONResponse:
    model = _REGISTRY.get(model_id).model
    return _bounded((await _compute(partial(orbitalTable, model=model))).to_dict())


@router.post("/api/models/{model_id}/reports/hamiltonian-basis")
async def model_hamiltonian_basis_report(
    model_id: str, payload: HamiltonianBasisReportRequest
) -> JSONResponse:
    if (payload.row is None) != (payload.column is None):
        raise ModelError(
            "InvalidHamiltonianBasisIndex", "row and column must be supplied together"
        )
    model = _REGISTRY.get(model_id).model
    report = await _compute(
        partial(showHamiltonianBasis, payload.row, payload.column, model=model)
    )
    return _bounded(report.to_dict())


@router.post("/api/models/{model_id}/reports/bonds")
async def model_bond_report(model_id: str, payload: ShellReportRequest) -> JSONResponse:
    model = _REGISTRY.get(model_id).model
    report = await _compute(partial(bondTable, payload.shell, model=model))
    return _bounded(report.to_dict())


@router.post("/api/models/{model_id}/reports/symmetry-representations")
async def model_symmetry_report(
    model_id: str, payload: SymmetryReportRequest
) -> JSONResponse:
    guard_json_tree(payload.selection)
    model = _REGISTRY.get(model_id).model
    report = await _compute(
        partial(showSymmetryRepresentations, payload.selection, model=model)
    )
    return _bounded(report.to_dict())


@router.post("/api/models/{model_id}/reports/hopping-parameters")
async def model_hopping_report(
    model_id: str, payload: HoppingReportRequest
) -> JSONResponse:
    model = _REGISTRY.get(model_id).model
    report = await _compute(
        partial(
            showHoppingParameters,
            payload.shell,
            payload.parameter,
            model=model,
            hermitian=payload.hermitian,
            kernel_method=payload.kernel_method,
            validation_level=payload.validation_level,
        )
    )
    return _bounded(report.to_dict())


@router.get("/api/models/{model_id}/representations")
async def model_representations(model_id: str) -> JSONResponse:
    model = _REGISTRY.get(model_id).model
    return _bounded({"representation_matrices": model.representation_matrices})


@router.get("/api/models/{model_id}/identity")
async def model_identity(model_id: str) -> JSONResponse:
    model = _REGISTRY.get(model_id).model
    return _bounded(
        {
            "model_identity_sha256": model.model_identity_sha256,
            "payload": model.model_identity_payload,
        }
    )


@router.get("/api/models/{model_id}/canonical")
async def model_canonical(
    model_id: str, shell: int = Query(default=1, ge=1, le=128)
) -> JSONResponse:
    model = _REGISTRY.get(model_id).model
    return _bounded(model.to_canonical_dict(shell))


@router.get("/api/models/{model_id}/bonds")
async def model_bonds(model_id: str, shell: Optional[int] = Query(default=None, ge=1, le=128)) -> JSONResponse:
    model = _REGISTRY.get(model_id).model
    if shell is None:
        result = model.bond_shells
    else:
        if shell > model.initial_bond_shells:
            raise ModelError("ShellNotPrepared", f"shell {shell} is not prepared")
        result = model.bond_shells[shell - 1]
    return _bounded({"shell": shell, "bonds": result})


@router.post("/api/models/{model_id}/symham")
async def model_symham(model_id: str, options: HamiltonianOptions) -> JSONResponse:
    model = _REGISTRY.get(model_id).model
    result = await _compute(partial(_hamiltonian, model, options))
    readable = result.to_dict()
    return _bounded(
        {
            "result_kind": "stable_matrix",
            "shape": list(result.shape),
            "shell": options.shell,
            "parameter_names": list(result.parameter_names),
            "matrix": readable["matrix"],
            "matrix_text": [[str(item) for item in row] for row in result.matrix],
        }
    )


@router.post("/api/models/{model_id}/hamiltonian")
async def model_hamiltonian(model_id: str, options: HamiltonianOptions) -> JSONResponse:
    """Return the advanced immutable Hamiltonian object structure."""

    model = _REGISTRY.get(model_id).model
    result = await _compute(partial(_hamiltonian, model, options))
    payload = _hamiltonian_payload(result)
    payload["result_kind"] = "hamiltonian_object"
    return _bounded(payload)


@router.post("/api/models/{model_id}/unsymham")
async def model_unsymham(model_id: str, payload: ShellReportRequest) -> JSONResponse:
    """Return the stable unconstrained shell matrix generated by Rust."""

    model = _REGISTRY.get(model_id).model
    result = await _compute(partial(model.unconstrained_hamiltonian, payload.shell))
    readable = result.to_dict()
    return _bounded(
        {
            "result_kind": "stable_unconstrained_matrix",
            "shape": list(result.shape),
            "shell": payload.shell,
            "parameter_names": list(result.parameter_names),
            "matrix": readable["matrix"],
            "matrix_text": [[str(item) for item in row] for row in result.matrix],
        }
    )


@router.post("/api/models/{model_id}/hamiltonian-space")
async def model_hamiltonian_space(model_id: str, options: HamiltonianOptions) -> JSONResponse:
    model = _REGISTRY.get(model_id).model
    result = await _compute(
        partial(
            model.hamiltonian_space,
            options.shell,
            hermitian=options.hermitian,
            cartesian_coordinates=options.cartesian_coordinates,
            validation_level=options.validation_level,
            kernel_method=options.kernel_method,
        )
    )
    return _bounded(result.to_dict())


@router.post("/api/models/{model_id}/evaluate")
async def evaluate_hamiltonian(model_id: str, payload: EvaluateRequest) -> JSONResponse:
    guard_json_tree(payload.parameters)
    guard_json_tree(payload.momentum)
    model = _REGISTRY.get(model_id).model
    options = HamiltonianOptions(
        shell=payload.shell,
        hermitian=payload.hermitian,
        cartesian_coordinates=payload.cartesian_coordinates,
        validation_level=payload.validation_level,
        kernel_method=payload.kernel_method,
    )
    hamiltonian = await _compute(
        partial(
            _cumulative_hamiltonian if payload.cumulative else _hamiltonian,
            model,
            options,
        )
    )
    parameters = decode_tagged(payload.parameters)
    momentum = decode_tagged(payload.momentum)
    result = await _compute(partial(hamiltonian.evaluate, parameters, momentum))
    return _bounded(
        {
            "shell": payload.shell,
            "shells": list(hamiltonian.shells),
            "matrix": result,
        }
    )


@router.post("/api/models/{model_id}/combine")
async def combine_hamiltonians(model_id: str, payload: CombineRequest) -> JSONResponse:
    model = _REGISTRY.get(model_id).model
    if len(set(payload.shells)) != len(payload.shells):
        raise ModelError(
            "DuplicateHamiltonianShell", "one shell may occur only once in a combination"
        )

    def combine() -> Hamiltonian:
        values = [
            model.hamiltonian(
                shell,
                hermitian=payload.hermitian,
                validation_level=payload.validation_level,
                kernel_method=payload.kernel_method,
            )
            for shell in payload.shells
        ]
        result = values[0]
        for value in values[1:]:
            result = result + value
        return result

    return _bounded(_hamiltonian_payload(await _compute(combine)))


@router.post("/api/models/{model_id}/properties/hopping-data")
async def model_hopping_data(
    model_id: str, payload: SolvedShellRequest
) -> JSONResponse:
    model = _REGISTRY.get(model_id).model
    return _bounded(await _compute(partial(_property_hopping_data, model, payload)))
