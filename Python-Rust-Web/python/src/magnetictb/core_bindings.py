"""Core Bindings bindings backed by the Rust extension."""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from . import _core

from .errors import (
    _translate_core_error,
    _translate_fitting_error,
    _translate_properties_error,
)


def null_space_json(payload: str) -> str:
    """Return an exact single-matrix null-space result as canonical JSON."""

    if not isinstance(payload, str):
        raise TypeError("payload must be a JSON string")
    return _translate_core_error(_core.null_space_json, payload)


def common_kernel_json(payload: str) -> str:
    """Return an exact multi-matrix common-kernel result as canonical JSON."""

    if not isinstance(payload, str):
        raise TypeError("payload must be a JSON string")
    return _translate_core_error(_core.common_kernel_json, payload)


def cyclotomic_context(conductor: int, max_degree: int = 128) -> Dict[str, Any]:
    """Construct one exact Rust cyclotomic context without Python algebra."""

    if (
        isinstance(conductor, bool)
        or not isinstance(conductor, int)
        or conductor < 1
        or isinstance(max_degree, bool)
        or not isinstance(max_degree, int)
        or max_degree < 1
    ):
        raise TypeError("conductor and max_degree must be positive integers")
    return _decode_result(
        _translate_core_error(
            _core.cyclotomic_context_json,
            conductor,
            max_degree,
        )
    )


def linear_algebra_json(payload: str) -> str:
    """Dispatch an exact linear-algebra request to the Rust core."""

    if not isinstance(payload, str):
        raise TypeError("payload must be a JSON string")
    return _translate_core_error(_core.linear_algebra_json, payload)


def representation_json(payload: str) -> str:
    """Dispatch an exact representation request to the Rust core."""

    if not isinstance(payload, str):
        raise TypeError("payload must be a JSON string")
    return _translate_core_error(_core.representation_json, payload)


def tight_binding_json(payload: str) -> str:
    """Dispatch an exact tight-binding request to the Rust core."""

    if not isinstance(payload, str):
        raise TypeError("payload must be a JSON string")
    return _translate_core_error(_core.tight_binding_json, payload)


def geometry_json(payload: str) -> str:
    """Dispatch an exact crystal-geometry request to the Rust core."""

    if not isinstance(payload, str):
        raise TypeError("payload must be a JSON string")
    return _translate_core_error(_core.geometry_json, payload)


def properties_json(payload: str) -> str:
    """Dispatch a numerical post-processing request to the Rust core."""

    if not isinstance(payload, str):
        raise TypeError("payload must be a JSON string")
    return _translate_properties_error(_core.properties_json, payload)


def fitting_json(payload: str) -> str:
    """Dispatch a numerical band-fitting request to the Rust core."""

    if not isinstance(payload, str):
        raise TypeError("payload must be a JSON string")
    return _translate_fitting_error(_core.fitting_json, payload)


def gapless_points_json(
    hamiltonian: Callable[[Sequence[float]], Any], payload: str
) -> str:
    """Run the Rust-owned periodic gapless-point search over a Python callback."""

    if not callable(hamiltonian):
        raise TypeError("hamiltonian must be callable")
    if not isinstance(payload, str):
        raise TypeError("payload must be a JSON string")
    return _translate_properties_error(
        _core.gapless_points_json,
        hamiltonian,
        payload,
    )


def _encode(value: Dict[str, Any]) -> str:
    return json.dumps(value, ensure_ascii=True, allow_nan=False, separators=(",", ":"))


def _decode_result(payload: str) -> Dict[str, Any]:
    value = json.loads(payload)
    if not isinstance(value, dict):
        raise RuntimeError("Rust core returned a non-object JSON result")
    return value


def null_space(context: Dict[str, Any], matrix: Dict[str, Any]) -> Dict[str, Any]:
    """Compute an exact null space in Rust from canonical context/matrix data."""

    request = {"context": context, "matrix": matrix}
    return _decode_result(null_space_json(_encode(request)))


def common_kernel(problem: Dict[str, Any]) -> Dict[str, Any]:
    """Compute an exact common kernel in Rust from a canonical problem."""

    return _decode_result(common_kernel_json(_encode(problem)))


def linear_algebra(
    operation: str, context: Dict[str, Any], **arguments: Any
) -> Dict[str, Any]:
    """Call one exact Rust linear-algebra operation through typed JSON data."""

    if not isinstance(operation, str):
        raise TypeError("operation must be a string")
    request = {"operation": operation, "context": context, **arguments}
    return _decode_result(linear_algebra_json(_encode(request)))


def representation(operation: str, **arguments: Any) -> Dict[str, Any]:
    """Call one ordered-group or exact representation operation in Rust."""

    if not isinstance(operation, str):
        raise TypeError("operation must be a string")
    request = {"operation": operation, **arguments}
    return _decode_result(representation_json(_encode(request)))


def tight_binding(
    operation: str, context: Dict[str, Any], **arguments: Any
) -> Dict[str, Any]:
    """Call one exact Rust tight-binding reconstruction or assembly operation."""

    if not isinstance(operation, str):
        raise TypeError("operation must be a string")
    request = {"operation": operation, "context": context, **arguments}
    return _decode_result(tight_binding_json(_encode(request)))


def geometry(operation: str, **arguments: Any) -> Dict[str, Any]:
    """Call one exact Rust crystal-geometry operation."""

    if not isinstance(operation, str):
        raise TypeError("operation must be a string")
    request = {"operation": operation, **arguments}
    return _decode_result(geometry_json(_encode(request)))


def properties(operation: str, **arguments: Any) -> Dict[str, Any]:
    """Call one Rust numerical-properties operation through typed JSON data."""

    if not isinstance(operation, str):
        raise TypeError("operation must be a string")
    request = {"operation": operation, **arguments}
    return _decode_result(properties_json(_encode(request)))


def fitting(operation: str, **arguments: Any) -> Dict[str, Any]:
    """Call one Rust numerical-fitting operation through typed JSON data."""

    if not isinstance(operation, str):
        raise TypeError("operation must be a string")
    request = {"operation": operation, **arguments}
    return _decode_result(fitting_json(_encode(request)))
