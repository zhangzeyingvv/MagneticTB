"""Web request schemas for Mathematica ``Data`` selectors."""

from __future__ import annotations

from typing import Any, Dict, List, Literal

from pydantic import BaseModel


class MsgopRequest(BaseModel):
    group: Dict[str, Any]


class SubperiodicOperationRequest(BaseModel):
    kind: Literal["rod", "layer"]
    selector: Literal["og", "gray"] = "og"
    key: List[int]
