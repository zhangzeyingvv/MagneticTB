"""Shared bounded runtime state for the local MagneticTB web service."""

from __future__ import annotations

import asyncio
import json
import os
from typing import Any, Dict, List, Optional

from fastapi import HTTPException
from fastapi.responses import JSONResponse
from starlette.concurrency import run_in_threadpool

from ..data_api import DataCatalog
from .codec import jsonable
from .registry import ModelRegistry


MAX_REQUEST_BYTES = int(os.environ.get("MAGNETICTB_WEB_MAX_REQUEST_BYTES", "2097152"))
MAX_RESPONSE_BYTES = int(os.environ.get("MAGNETICTB_WEB_MAX_RESPONSE_BYTES", "8388608"))
MAX_CONCURRENT = int(os.environ.get("MAGNETICTB_WEB_MAX_CONCURRENT", "2"))
COMPUTE_TIMEOUT = float(os.environ.get("MAGNETICTB_WEB_TIMEOUT", "120"))
MODEL_LIMIT = int(os.environ.get("MAGNETICTB_WEB_MODEL_LIMIT", "8"))
MODEL_TTL = float(os.environ.get("MAGNETICTB_WEB_MODEL_TTL", "3600"))

_REGISTRY = ModelRegistry(maximum=MODEL_LIMIT, ttl_seconds=MODEL_TTL)
_COMPUTE_SEMAPHORE = asyncio.Semaphore(MAX_CONCURRENT)
_CATALOG: Optional[DataCatalog] = None
_MSG_FAMILY_CACHE: Dict[int, List[Dict[str, Any]]] = {}
_WYCKOFF_OPTION_CACHE: Dict[str, List[Dict[str, Any]]] = {}
_BRAVAIS_OPTION_CACHE: Dict[str, Dict[str, Any]] = {}


def _catalog() -> DataCatalog:
    global _CATALOG
    if _CATALOG is None:
        _CATALOG = DataCatalog()
    return _CATALOG


async def _compute(call: Any) -> Any:
    try:
        await asyncio.wait_for(_COMPUTE_SEMAPHORE.acquire(), timeout=5.0)
    except asyncio.TimeoutError as error:
        raise HTTPException(
            status_code=429,
            detail={"tag": "WebBusy", "detail": "all local computation slots are busy"},
        ) from error
    try:
        return await asyncio.wait_for(run_in_threadpool(call), timeout=COMPUTE_TIMEOUT)
    except asyncio.TimeoutError as error:
        raise HTTPException(
            status_code=408,
            detail={
                "tag": "WebComputationTimeout",
                "detail": "the configured response deadline was exceeded",
            },
        ) from error
    except (TypeError, ValueError, KeyError) as error:
        raise HTTPException(
            status_code=400,
            detail={"tag": "WebInvalidRequest", "detail": str(error)},
        ) from error
    finally:
        _COMPUTE_SEMAPHORE.release()


def _bounded(value: Any, *, status_code: int = 200) -> JSONResponse:
    payload = jsonable(value)
    encoded = json.dumps(
        payload, ensure_ascii=False, allow_nan=False, separators=(",", ":")
    ).encode("utf-8")
    if len(encoded) > MAX_RESPONSE_BYTES:
        raise HTTPException(
            status_code=413,
            detail={
                "tag": "WebResultTooLarge",
                "detail": "result exceeds the configured response limit; request a smaller page or section",
            },
        )
    return JSONResponse(content=payload, status_code=status_code)
