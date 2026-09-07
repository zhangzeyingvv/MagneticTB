"""Web request schemas for Mathematica ``Properties`` APIs."""

from __future__ import annotations

from typing import Any, Dict, List, Literal, Optional

from pydantic import BaseModel, Field


class SolvedShellRequest(BaseModel):
    shells: List[int] = Field(default_factory=lambda: [1], min_length=1, max_length=128)
    parameters: Dict[str, Any] = Field(default_factory=dict)
    hermitian: bool = True
    validation_level: Literal["None", "Basic", "Full"] = "Basic"
    kernel_method: Literal["Iterative", "Stacked", "Cyclotomic"] = "Iterative"


class BlochRequest(SolvedShellRequest):
    momentum: List[Any] = Field(default_factory=lambda: [0, 0, 0], min_length=3, max_length=3)


class BandPlotRequest(SolvedShellRequest):
    npoint: int = Field(default=40, ge=1, le=2000)
    bravais_type: str = Field(default="Automatic", min_length=1, max_length=80)
    tolerance: float = Field(default=1.0e-6, gt=0)


class RealSpaceRequest(SolvedShellRequest):
    geometry: List[Any] = Field(default_factory=lambda: [4, 1, 1], min_length=1)
    boundary_conditions: List[Literal["Open", "Periodic"]] = Field(
        default_factory=lambda: ["Open", "Open", "Open"], min_length=3, max_length=3
    )


class SlabRequest(SolvedShellRequest):
    size: List[int] = Field(default_factory=lambda: [1, 1, 8], min_length=3, max_length=3)
    momentum: List[Any] = Field(default_factory=lambda: [0, 0], min_length=1, max_length=2)
    periodic_directions: List[int] = Field(default_factory=lambda: [1, 2], min_length=1, max_length=2)
    cell_matrix: Optional[List[List[int]]] = None


class SurfaceRequest(SolvedShellRequest):
    momentum: List[Any] = Field(default_factory=lambda: [0, 0], min_length=2, max_length=2)
    energy: Any = 0
    broadening: float = Field(default=1.0e-3, gt=0)
    tolerance: float = Field(default=1.0e-10, gt=0)
    max_iterations: int = Field(default=200, ge=1)
    surface: Literal["Positive", "Negative"] = "Positive"
    cell_matrix: Optional[List[List[int]]] = None


class WilsonRequest(SolvedShellRequest):
    occupied: int = Field(ge=1)
    start: List[Any]
    end: List[Any]
    path_subdivisions: int = Field(default=50, ge=1)


class BerryPhaseRequest(SolvedShellRequest):
    occupied: int = Field(ge=1)
    path: List[List[Any]] = Field(min_length=2)


class BerryCurvatureRequest(SolvedShellRequest):
    occupied: int = Field(ge=1)
    point: List[Any]
    directions: List[int] = Field(default_factory=lambda: [1, 2], min_length=2, max_length=2)
    step_size: Any = 1.0e-3


class PointChernRequest(SolvedShellRequest):
    occupied: int = Field(ge=1)
    point: List[Any] = Field(min_length=3, max_length=3)
    radius: Any
    surface_subdivisions: int = Field(default=8, ge=1)
    require_gapless_center: bool = True


class GaplessRequest(SolvedShellRequest):
    occupied: int = Field(ge=1)
    brillouin_zone: List[List[Any]] = Field(
        default_factory=lambda: [[-3.141592653589793, 3.141592653589793]] * 3,
        min_length=1,
        max_length=3,
    )
    grid_size: Any = 15
    candidate_count: int = Field(default=32, ge=1)
    gap_tolerance: float = Field(default=1.0e-7, gt=0)
    merge_tolerance: float = Field(default=1.0e-4, gt=0)
    max_iterations: int = Field(default=500, ge=1)
    refinement_method: Literal["PrincipalAxis", "QuasiNewton"] = "PrincipalAxis"


class WilsonPlotRequest(WilsonRequest):
    parameter_path: List[Any]
    parameter_subdivisions: int = Field(default=60, ge=1, le=1000)


class SurfaceSpectrumRequest(SolvedShellRequest):
    momentum_path: List[Any]
    energy_range: List[Any] = Field(min_length=2, max_length=2)
    momentum_subdivisions: int = Field(default=60, ge=1, le=1000)
    energy_points: int = Field(default=201, ge=2, le=2001)
    broadening: float = Field(default=1.0e-3, gt=0)


class BerryCurvature2DPlotRequest(SolvedShellRequest):
    occupied: int = Field(ge=1)
    ranges: List[List[Any]] = Field(min_length=2, max_length=2)
    directions: List[int] = Field(default_factory=lambda: [1, 2], min_length=2, max_length=2)
    fixed_coordinates: Any = "Automatic"
    grid_size: Any = 31
    step_size: Any = "Automatic"


class BerryCurvature3DPlotRequest(SolvedShellRequest):
    occupied: int = Field(ge=1)
    ranges: List[List[Any]] = Field(min_length=3, max_length=3)
    grid_size: Any = 7
    step_size: Any = "Automatic"
    magnitude_threshold: float = Field(default=0.0, ge=0)


class WavefunctionPlotRequest(SolvedShellRequest):
    geometry: List[Any] = Field(default_factory=lambda: [4, 1, 1], min_length=1)
    state: List[Any] = Field(min_length=1)
    aggregation: Literal["Atom", "Orbital"] = "Atom"
    normalize: bool = True
    weight_threshold: float = Field(default=1.0e-8, ge=0)
    phase_coloring: bool = False
