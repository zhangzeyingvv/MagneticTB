"""Bounded in-process registry for immutable prepared models."""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
from threading import RLock
from time import monotonic
from typing import List
from uuid import uuid4

from ..api import ModelError
from ..modeling import PreparedModel


@dataclass
class ModelEntry:
    model_id: str
    label: str
    model: PreparedModel
    last_access: float


class ModelRegistry:
    """Small TTL/LRU registry; one web process owns one independent registry."""

    def __init__(self, *, maximum: int, ttl_seconds: float) -> None:
        if maximum < 1 or ttl_seconds <= 0:
            raise ValueError("registry limits must be positive")
        self._maximum = maximum
        self._ttl_seconds = ttl_seconds
        self._entries: "OrderedDict[str, ModelEntry]" = OrderedDict()
        self._lock = RLock()

    @property
    def capacity(self) -> int:
        return self._maximum

    def _expire(self, now: float) -> None:
        expired = [
            key
            for key, entry in self._entries.items()
            if now - entry.last_access > self._ttl_seconds
        ]
        for key in expired:
            self._entries.pop(key, None)

    def put(self, model: PreparedModel, label: str) -> ModelEntry:
        now = monotonic()
        with self._lock:
            self._expire(now)
            while len(self._entries) >= self._maximum:
                self._entries.popitem(last=False)
            identity = model.model_identity_sha256[:12]
            model_id = f"{identity}-{uuid4().hex[:12]}"
            entry = ModelEntry(model_id, label, model, now)
            self._entries[model_id] = entry
            return entry

    def get(self, model_id: str) -> ModelEntry:
        now = monotonic()
        with self._lock:
            self._expire(now)
            entry = self._entries.get(model_id)
            if entry is None:
                raise ModelError(
                    "WebModelNotFound",
                    "the model id is unknown or its local web session expired",
                )
            entry.last_access = now
            self._entries.move_to_end(model_id)
            return entry

    def delete(self, model_id: str) -> bool:
        with self._lock:
            return self._entries.pop(model_id, None) is not None

    def list(self) -> List[ModelEntry]:
        now = monotonic()
        with self._lock:
            self._expire(now)
            return list(reversed(self._entries.values()))

    def clear(self) -> None:
        with self._lock:
            self._entries.clear()
