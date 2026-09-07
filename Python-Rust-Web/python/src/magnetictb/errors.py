"""Errors bindings backed by the Rust extension."""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from . import _core


class MagneticTBError(RuntimeError):
    """Base class for stable failures exposed by MagneticTB."""

    def __init__(self, tag: str, detail: str) -> None:
        self.tag = tag
        self.detail = detail
        super().__init__(detail)


class ExactError(MagneticTBError):
    """An explicit failure reported by the Rust exact core."""

    def __init__(self, tag: str, detail: str) -> None:
        super().__init__(tag, detail)


class GroupError(MagneticTBError):
    """An explicit failure reported by the Rust finite-group core."""

    def __init__(self, tag: str, detail: str) -> None:
        super().__init__(tag, detail)


class DataError(MagneticTBError):
    """An explicit failure reported by the Rust data catalog."""

    def __init__(self, tag: str, detail: str) -> None:
        super().__init__(tag, detail)


class SymmetryError(MagneticTBError):
    """An explicit failure reported by the Rust symmetry core."""

    def __init__(self, tag: str, detail: str) -> None:
        super().__init__(tag, detail)


class ModelError(MagneticTBError):
    """A stable high-level model lifecycle or option failure."""


class PropertiesError(MagneticTBError):
    """An explicit failure reported by the Rust numerical-properties core."""


class FittingError(MagneticTBError):
    """An explicit failure reported by the Rust numerical-fitting core."""


def _translate_core_error(call: Callable[..., str], *arguments: Any) -> str:
    try:
        return call(*arguments)
    except _core.ExactError as error:
        arguments = error.args
        tag = str(arguments[0]) if arguments else "RustBindingFailure"
        detail = str(arguments[1]) if len(arguments) > 1 else str(error)
        raise ExactError(tag, detail) from None


def _translate_group_error(call: Callable[..., Any], *arguments: Any) -> Any:
    try:
        return call(*arguments)
    except _core.GroupError as error:
        values = error.args
        tag = str(values[0]) if values else "RustGroupBindingFailure"
        detail = str(values[1]) if len(values) > 1 else str(error)
        raise GroupError(tag, detail) from None


def _translate_data_error(call: Callable[..., Any], *arguments: Any) -> Any:
    try:
        return call(*arguments)
    except _core.DataError as error:
        values = error.args
        tag = str(values[0]) if values else "RustDataBindingFailure"
        detail = str(values[1]) if len(values) > 1 else str(error)
        raise DataError(tag, detail) from None


def _translate_symmetry_error(call: Callable[..., Any], *arguments: Any) -> Any:
    try:
        return call(*arguments)
    except _core.SymmetryError as error:
        values = error.args
        tag = str(values[0]) if values else "RustSymmetryBindingFailure"
        detail = str(values[1]) if len(values) > 1 else str(error)
        raise SymmetryError(tag, detail) from None


def _translate_properties_error(call: Callable[..., Any], *arguments: Any) -> Any:
    try:
        return call(*arguments)
    except _core.PropertiesError as error:
        values = error.args
        tag = str(values[0]) if values else "RustPropertiesBindingFailure"
        detail = str(values[1]) if len(values) > 1 else str(error)
        raise PropertiesError(tag, detail) from None


def _translate_fitting_error(call: Callable[..., Any], *arguments: Any) -> Any:
    try:
        return call(*arguments)
    except _core.FittingError as error:
        values = error.args
        tag = str(values[0]) if values else "RustFittingBindingFailure"
        detail = str(values[1]) if len(values) > 1 else str(error)
        raise FittingError(tag, detail) from None
