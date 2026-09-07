"""Web request schemas for Mathematica ``Model`` and ``TightBinding`` APIs."""

from __future__ import annotations

from typing import Any, Dict, List, Literal, Optional

from pydantic import BaseModel, Field


class PrepareRequest(BaseModel):
    entrypoint: Literal["init", "initfromrep", "prepare_model", "prepare_model_from_rep"] = "init"
    arguments: Dict[str, Any] = Field(default_factory=dict)
    label: str = Field(default="MagneticTB model", min_length=1, max_length=120)


class InitLatticeRequest(BaseModel):
    preset: Literal["msg_bravais", "stable_default", "simple_cubic", "custom"] = "msg_bravais"
    matrix: Optional[List[List[Any]]] = None
    parameters: Dict[str, Any] = Field(default_factory=dict)


class InitOperationRequest(BaseModel):
    kind: Literal["spatial", "spin_space"] = "spatial"
    label: str = Field(min_length=1, max_length=120)
    rotation: List[List[Any]]
    translation: List[Any] = Field(default_factory=lambda: [0, 0, 0], min_length=3, max_length=3)
    antiunitary: bool = False
    spin_rotation: Optional[List[List[Any]]] = None
    continuous_parameter: Optional[str] = Field(default=None, max_length=120)


class InitSymmetryRequest(BaseModel):
    source: Literal["msg", "explicit"] = "msg"
    msg_id: Optional[str] = None
    operations: List[InitOperationRequest] = Field(default_factory=list, max_length=512)


class InitOrbitRequest(BaseModel):
    letter: str = Field(default="a", min_length=1, max_length=12)
    source_ordinal: Optional[int] = Field(default=None, ge=1)
    coordinate_parameters: Dict[str, Any] = Field(default_factory=dict)
    moment_parameters: Dict[str, Any] = Field(default_factory=dict)
    position: Optional[List[Any]] = Field(default=None, min_length=3, max_length=3)
    moment: Optional[List[Any]] = Field(default=None, min_length=3, max_length=3)
    basis_functions: List[str] = Field(default_factory=lambda: ["s"], min_length=1, max_length=128)


class InitFieldRequest(BaseModel):
    preset: Literal["crystallographic", "rational", "gaussian", "custom"] = "crystallographic"
    conductor: Optional[int] = Field(default=None, ge=1, le=65536)
    cyclotomic_polynomial: Optional[List[int]] = Field(default=None, min_length=2, max_length=1025)


class InitDirectRepresentationRequest(BaseModel):
    matrices_by_orbit: List[List[List[List[Any]]]] = Field(min_length=1)
    continuous_site_matrices: Optional[List[List[List[List[Any]]]]] = None
    orbital_labels: List[List[str]] = Field(min_length=1)


class FriendlyInitRequest(BaseModel):
    label: str = Field(default="MagneticTB model", min_length=1, max_length=120)
    lattice: InitLatticeRequest = Field(default_factory=InitLatticeRequest)
    symmetry: InitSymmetryRequest = Field(default_factory=InitSymmetryRequest)
    orbits: List[InitOrbitRequest] = Field(default_factory=lambda: [InitOrbitRequest()], min_length=1, max_length=64)
    initial_bond_shells: int = Field(default=10, ge=1, le=128)
    generate_symmetry_group: bool = False
    representation_mode: Literal["DirectProduct", "Induced"] = "DirectProduct"
    exact_field: InitFieldRequest = Field(default_factory=InitFieldRequest)
    direct_representation: Optional[InitDirectRepresentationRequest] = None


class HamiltonianOptions(BaseModel):
    shell: int = Field(ge=1, le=128)
    hermitian: bool = True
    cartesian_coordinates: bool = False
    validation_level: Literal["None", "Basic", "Full"] = "Basic"
    kernel_method: Literal["Iterative", "Stacked", "Cyclotomic"] = "Iterative"

class ShellReportRequest(BaseModel):
    shell: int = Field(ge=1, le=128)


class HamiltonianBasisReportRequest(BaseModel):
    row: Optional[int] = Field(default=None, ge=1)
    column: Optional[int] = Field(default=None, ge=1)


class SymmetryReportRequest(BaseModel):
    selection: Any = "Automatic"


class HoppingReportRequest(HamiltonianOptions):
    parameter: Optional[str] = Field(default=None, min_length=1, max_length=120)


class EvaluateRequest(HamiltonianOptions):
    parameters: Dict[str, Any] = Field(default_factory=dict)
    momentum: List[Any] = Field(default_factory=lambda: [0, 0, 0], min_length=3, max_length=3)
    cumulative: bool = False


class CombineRequest(BaseModel):
    shells: List[int] = Field(min_length=2, max_length=128)
    hermitian: bool = True
    validation_level: Literal["None", "Basic", "Full"] = "Basic"
    kernel_method: Literal["Iterative", "Stacked", "Cyclotomic"] = "Iterative"
