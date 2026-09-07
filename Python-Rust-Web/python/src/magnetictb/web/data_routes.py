"""Data Routes for the local MagneticTB workbench."""

from __future__ import annotations

from fastapi import APIRouter

from .route_dependencies import (
    Any,
    DataCatalog,
    Dict,
    JSONResponse,
    List,
    MagneticGroupSelector,
    Mapping,
    MsgopRequest,
    Query,
    SubperiodicOperationRequest,
    _bounded,
    _catalog,
    _compute,
    guard_json_tree,
    io,
    mlgop,
    mrgop,
    msgop,
    partial,
    redirect_stdout,
    showMSGWyckoff,
)


router = APIRouter()


def _decode_group_selector(value: Mapping[str, Any]) -> Any:
    kind = value.get("kind")
    if kind == "source_index":
        return int(value["value"])
    if kind == "msg":
        stable_id = value.get("stable_id")
        bns = value.get("bns")
        return MSG(
            stable_id=str(stable_id) if stable_id is not None else None,
            bns=tuple(int(item) for item in bns) if bns is not None else None,
        )
    if kind == "selector":
        return MagneticGroupSelector(
            str(value["classification"]), tuple(int(item) for item in value["source_key"])
        )
    raise ValueError("group kind must be source_index, msg, or selector")


@router.post("/api/msgop")
async def data_msgop(payload: MsgopRequest) -> JSONResponse:
    guard_json_tree(payload.group)

    def run() -> Dict[str, Any]:
        output = io.StringIO()
        with redirect_stdout(output):
            operations = msgop(_decode_group_selector(payload.group))
        return {
            "operations": operations,
            "count": len(operations),
            "messages": [line for line in output.getvalue().splitlines() if line],
        }

    return _bounded(await _compute(run))


@router.post("/api/data/subperiodic-operations")
async def data_subperiodic_operations(
    payload: SubperiodicOperationRequest,
) -> JSONResponse:
    def run() -> Dict[str, Any]:
        if payload.selector == "gray":
            if len(payload.key) != 1:
                raise ValueError("grayrod/graylayer selector requires one integer")
            og_key = list(
                _catalog().resolve_subperiodic_gray_key(payload.kind, payload.key[0])
            )
        else:
            og_key = payload.key
        output = io.StringIO()
        with redirect_stdout(output):
            operations = (
                mrgop(og_key)
                if payload.kind == "rod"
                else mlgop(og_key)
            )
        return {
            "kind": payload.kind,
            "selector": payload.selector,
            "key": payload.key,
            "og_key": og_key,
            "operations": operations,
            "count": len(operations),
            "messages": [line for line in output.getvalue().splitlines() if line],
        }

    return _bounded(await _compute(run))


@router.post("/api/data/msg-wyckoff-report")
async def data_msg_wyckoff_report(payload: MsgopRequest) -> JSONResponse:
    guard_json_tree(payload.group)
    report = await _compute(partial(showMSGWyckoff, _decode_group_selector(payload.group)))
    return _bounded(report.to_dict())


@router.get("/api/data/summary")
async def data_summary() -> JSONResponse:
    def load() -> Dict[str, Any]:
        catalog = _catalog()
        return {"counts": catalog.counts, "manifest": catalog.manifest}

    return _bounded(await _compute(load))


@router.get("/api/data-classifications")
async def data_classifications() -> JSONResponse:
    return _bounded(await _compute(lambda: _catalog().classification_maps))


def _family_ids(catalog: DataCatalog, family: str) -> List[str]:
    if family == "bravais":
        return catalog.bravais_ids
    if family == "msg":
        return catalog.msg_ids
    if family == "rod":
        return catalog.rod_group_ids
    if family == "layer":
        return catalog.layer_group_ids
    raise ValueError("family must be bravais, msg, rod, or layer")


@router.get("/api/data/{family}")
async def list_data(
    family: str,
    offset: int = Query(default=0, ge=0),
    limit: int = Query(default=50, ge=1, le=250),
    search: str = Query(default="", max_length=120),
) -> JSONResponse:
    def load() -> Dict[str, Any]:
        values = _family_ids(_catalog(), family)
        query = search.casefold().strip()
        if query:
            values = [value for value in values if query in value.casefold()]
        return {
            "family": family,
            "total": len(values),
            "offset": offset,
            "limit": limit,
            "ids": values[offset : offset + limit],
        }

    return _bounded(await _compute(load))


@router.get("/api/data/{family}/{stable_id}")
async def get_data_record(family: str, stable_id: str) -> JSONResponse:
    def load() -> Any:
        catalog = _catalog()
        if family == "bravais":
            return catalog.bravais(stable_id)
        if family == "msg":
            return catalog.msg(stable_id)
        if family == "rod":
            return catalog.rod_group(stable_id)
        if family == "layer":
            return catalog.layer_group(stable_id)
        if family == "wyckoff":
            return catalog.wyckoff(stable_id)
        raise ValueError("family must be bravais, msg, wyckoff, rod, or layer")

    return _bounded(await _compute(load))
