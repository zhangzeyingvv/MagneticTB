"""Web request schemas for Mathematica ``CrystalGeometry`` views."""

from __future__ import annotations

from typing import Any, Optional

from pydantic import BaseModel, Field


class CrystalStructurePlotRequest(BaseModel):
    cell_range: Any = 0
    moment_scale: Optional[float] = Field(default=None, gt=0)
    atom_radius: Optional[float] = Field(default=None, gt=0)
    show_atom_labels: bool = False


class BrillouinZonePlotRequest(BaseModel):
    k_path: Any = "Automatic"
    show_k_path: bool = True
    translation_range: int = Field(default=2, ge=1, le=16)
    tolerance: float = Field(default=1.0e-8, gt=0)
