"""Tight Binding boundary for the Rust-backed MagneticTB Python API.

This module performs only immutable input/output adaptation; mathematical
algorithms remain in the Rust core.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from contextvars import ContextVar
from dataclasses import dataclass, field
from fractions import Fraction
from functools import lru_cache
from math import isfinite
from threading import RLock
from types import MappingProxyType
from typing import Any, Dict, Literal, Optional, Tuple, Union

from .core_bindings import geometry
from .errors import ExactError, ModelError

from .linear_algebra import (
    CyclotomicElement,
    ExactAtom,
    ExactField,
    ExactMatrix,
    PhaseCoordinate,
    QuadraticNumber,
    _cyclotomic_to_user,
    _decode_element,
    _decode_matrix,
    _decode_phase_coordinate,
    _decode_quadratic,
    _encode_cyclotomic_element,
    _encode_evaluation_scalar,
    _encode_phase_coordinate,
    _field_from_context,
    _freeze_json,
    _matrix_to_user,
    _phase_coordinate_to_user,
    _quadratic_to_user,
    _thaw_json,
)


@dataclass(frozen=True)
class OrbitalRecord:
    orbital_id: int
    wyckoff_orbit: int
    equivalent_atom: int
    site_id: int
    local_orbital: int
    fractional_position: Tuple[CyclotomicElement, ...]
    cartesian_position: Tuple[CyclotomicElement, ...]
    basis_function: str
    basis_state: Any
    spatial_orbital: Any
    spin_structure: str
    basis_convention: Literal["DirectProduct", "Induced"]
    spin_state: Any
    reference_equivalent_atom: Optional[int]
    transport_operation: Optional[int]
    transport_operation_label: Optional[str]
    core_orbital_index: int = field(repr=False)

    @property
    def basis_function_input(self) -> str:
        """Stable input label, distinct from the transported actual state."""

        return self.basis_function


@dataclass(frozen=True)
class BondRecord:
    shell: int
    bond_id: int
    source_site: int
    target_site: int
    translation: Tuple[int, ...]
    displacement: Tuple[QuadraticNumber, ...]
    source_endpoint: Tuple[QuadraticNumber, ...]
    target_endpoint: Tuple[QuadraticNumber, ...]
    squared_distance: QuadraticNumber
    core_bond_index: int = field(repr=False)

    @property
    def row_site(self) -> int:
        """Stable ``RowSite`` (bra/destination) index."""

        return self.source_site

    @property
    def column_site(self) -> int:
        """Stable ``ColumnSite`` (ket/source) index."""

        return self.target_site

    @property
    def row_endpoint(self) -> Tuple[QuadraticNumber, ...]:
        return self.source_endpoint

    @property
    def column_endpoint(self) -> Tuple[QuadraticNumber, ...]:
        return self.target_endpoint


@dataclass(frozen=True)
class HamiltonianParameter:
    name: str
    parameter_number: int
    constraint_orbit: int
    orbit_parameter: int
    representative_bond: int
    core_parameter_index: int = field(repr=False)


@dataclass(frozen=True)
class HamiltonianBasisTerm:
    parameter: HamiltonianParameter
    representative_hoppings: Tuple[ExactMatrix, ...]
    fourier_coefficients: Tuple[Tuple[Tuple[QuadraticNumber, ...], ExactMatrix], ...]
    gamma_hamiltonian: ExactMatrix
    verified: bool


@dataclass(frozen=True)
class HamiltonianTerm:
    """One Rust-generated exact parameter contribution to ``H(kappa)``."""

    parameter_name: str
    shell: int
    parameter_number: int
    core_parameter_index: int
    fourier_coefficients: Tuple[
        Tuple[Tuple[PhaseCoordinate, ...], ExactMatrix], ...
    ]
    gamma_hamiltonian: ExactMatrix
    verified: bool


@dataclass(frozen=True)
class HamiltonianMatrixContribution:
    """One exact ``parameter*coefficient*exp(+i*k.d)`` term."""

    parameter_name: str
    displacement: Tuple[PhaseCoordinate, ...]
    coefficient: CyclotomicElement


def _format_quadratic(value: QuadraticNumber) -> str:
    return str(value).replace(" + -", " - ")


def _format_phase_coordinate(value: PhaseCoordinate) -> str:
    return _format_quadratic(value) if isinstance(value, QuadraticNumber) else str(value)


def _phase_coordinate_is_zero(value: PhaseCoordinate) -> bool:
    if isinstance(value, QuadraticNumber):
        return all(coefficient == 0 for coefficient in value.coefficients)
    return value == 0


def _format_cyclotomic(value: CyclotomicElement) -> str:
    if all(coefficient == 0 for coefficient in value.coefficients[1:]):
        return str(value.coefficients[0])
    rendered = str(value)
    if sum(coefficient != 0 for coefficient in value.coefficients) > 1:
        return f"({rendered})"
    return rendered


@dataclass(frozen=True)
class HamiltonianExpression:
    """One immutable Rust-generated symbolic matrix element of ``H(k)``."""

    contributions: Tuple[HamiltonianMatrixContribution, ...]

    def __bool__(self) -> bool:
        return bool(self.contributions)

    def __str__(self) -> str:
        if not self.contributions:
            return "0"
        rendered = []
        for contribution in self.contributions:
            factors = []
            coefficient = _format_cyclotomic(contribution.coefficient)
            if coefficient != "1":
                factors.append(coefficient)
            factors.append(contribution.parameter_name)
            phase_terms = []
            for index, displacement in enumerate(contribution.displacement):
                if _phase_coordinate_is_zero(displacement):
                    continue
                variable = ("kx", "ky", "kz")[index] if index < 3 else f"k{index + 1}"
                rendered_displacement = _format_phase_coordinate(displacement)
                if rendered_displacement == "1":
                    phase_terms.append(variable)
                elif rendered_displacement == "-1":
                    phase_terms.append(f"-{variable}")
                else:
                    phase_terms.append(f"({rendered_displacement})*{variable}")
            if phase_terms:
                phase = " + ".join(phase_terms).replace("+ -", "- ")
                factors.append(f"exp(I*({phase}))")
            rendered.append("*".join(factors))
        return " + ".join(rendered).replace("+ -", "- ")

    def __repr__(self) -> str:
        return str(self)

    def to_dict(self) -> list[dict[str, Any]]:
        return [
            {
                "parameter": contribution.parameter_name,
                "displacement": [
                    _phase_coordinate_to_user(value) for value in contribution.displacement
                ],
                "coefficient": _cyclotomic_to_user(contribution.coefficient),
            }
            for contribution in self.contributions
        ]


@dataclass(frozen=True)
class Hamiltonian:
    """Immutable, matrix-like Rust-generated stable ``symham`` result."""

    shape: Tuple[int, int]
    shells: Tuple[int, ...]
    terms: Tuple[HamiltonianTerm, ...]
    matrix: Tuple[Tuple[HamiltonianExpression, ...], ...]
    covariance_verified: bool
    model_identity_sha256: str
    field: ExactField
    gauge: Mapping[str, Any]
    _raw: Mapping[str, Any] = field(repr=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "shells", tuple(self.shells))
        object.__setattr__(self, "terms", tuple(self.terms))
        object.__setattr__(self, "matrix", tuple(tuple(row) for row in self.matrix))
        object.__setattr__(self, "gauge", _freeze_json(self.gauge))
        object.__setattr__(self, "_raw", _freeze_json(self._raw))

    def __len__(self) -> int:
        return self.shape[0]

    def __iter__(self):
        return iter(self.matrix)

    def __getitem__(self, index: int) -> Tuple[HamiltonianExpression, ...]:
        return self.matrix[index]

    def to_list(self) -> list[list[HamiltonianExpression]]:
        """Return ordinary matrix rows while preserving immutable expressions."""

        return [list(row) for row in self.matrix]

    @property
    def parameter_names(self) -> Tuple[str, ...]:
        return tuple(term.parameter_name for term in self.terms)

    @property
    def gamma_hamiltonians(self) -> Mapping[str, ExactMatrix]:
        return MappingProxyType(
            {term.parameter_name: term.gamma_hamiltonian for term in self.terms}
        )

    def __str__(self) -> str:
        return "[" + ",\n ".join(
            "[" + ", ".join(str(element) for element in row) + "]"
            for row in self.matrix
        ) + "]"

    def __add__(self, other: object) -> "Hamiltonian":
        if not isinstance(other, Hamiltonian):
            return NotImplemented
        try:
            raw = geometry(
                "combine_symbolic_hamiltonians",
                left=self.to_canonical_dict(),
                right=other.to_canonical_dict(),
            )
        except ExactError as error:
            raise ModelError(error.tag, error.detail) from None
        return _decode_hamiltonian(raw)

    def __radd__(self, other: object) -> "Hamiltonian":
        if other == 0:
            return self
        return self.__add__(other)

    def evaluate(
        self,
        parameters: Mapping[str, ExactAtom],
        momentum: Sequence[ExactAtom] = (0, 0, 0),
    ) -> Tuple[Tuple[complex, ...], ...]:
        """Numerically evaluate the Rust Hamiltonian at one ``(kx, ky, kz)``."""

        if not isinstance(parameters, Mapping):
            raise TypeError("parameters must be a mapping keyed by parameter name")
        if not isinstance(momentum, Sequence) or isinstance(momentum, (str, bytes)):
            raise TypeError("momentum must be a length-3 sequence")
        if len(momentum) != 3:
            raise ValueError("momentum must contain exactly kx, ky, kz")
        encoded_parameters = {
            str(name): _encode_evaluation_scalar(value, f"parameter {name}")
            for name, value in parameters.items()
        }
        encoded_momentum = [
            _encode_evaluation_scalar(value, name)
            for value, name in zip(momentum, ("kx", "ky", "kz"))
        ]
        try:
            evaluated = geometry(
                "evaluate_symbolic_hamiltonian",
                hamiltonian=self.to_canonical_dict(),
                parameters=encoded_parameters,
                momentum=encoded_momentum,
            )
        except ExactError as error:
            raise ModelError(error.tag, error.detail) from None
        return tuple(
            tuple(
                complex(float(value["real"]), float(value["imaginary"]))
                for value in row
            )
            for row in evaluated["matrix"]
        )

    def gamma(
        self, parameters: Mapping[str, ExactAtom]
    ) -> Tuple[Tuple[complex, ...], ...]:
        """Evaluate at ``kx=ky=kz=0`` in the recorded gauge."""

        return self.evaluate(parameters, (0, 0, 0))

    def to_dict(self) -> Dict[str, Any]:
        """Return the readable exact Fourier structure without canonical tags."""

        return {
            "shape": list(self.shape),
            "shells": list(self.shells),
            "parameter_names": list(self.parameter_names),
            "gauge": _thaw_json(self.gauge),
            "terms": [
                {
                    "parameter": term.parameter_name,
                    "shell": term.shell,
                    "parameter_number": term.parameter_number,
                    "fourier_coefficients": [
                        {
                            "displacement": [
                                _phase_coordinate_to_user(value) for value in displacement
                            ],
                            "matrix": _matrix_to_user(matrix),
                        }
                        for displacement, matrix in term.fourier_coefficients
                    ],
                    "gamma_hamiltonian": _matrix_to_user(term.gamma_hamiltonian),
                    "verified": term.verified,
                }
                for term in self.terms
            ],
            "matrix": [
                [
                    cell.to_dict()
                    for cell in row
                ]
                for row in self.matrix
            ],
            "covariance_verified": self.covariance_verified,
            "model_identity_sha256": self.model_identity_sha256,
        }

    def to_canonical_dict(self) -> Dict[str, Any]:
        """Return the Rust-owned canonical Hamiltonian AST."""

        return _thaw_json(self._raw)


@dataclass(frozen=True)
class HamiltonianSpace:
    shell: int
    parameters: Tuple[HamiltonianParameter, ...]
    bonds: Tuple[BondRecord, ...]
    basis: Tuple[HamiltonianBasisTerm, ...]
    covariance_verified: bool
    _raw: Mapping[str, Any] = field(repr=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "_raw", _freeze_json(self._raw))

    @property
    def parameter_names(self) -> Tuple[str, ...]:
        return tuple(parameter.name for parameter in self.parameters)

    @property
    def fourier_coefficients(
        self,
    ) -> Mapping[str, Tuple[Tuple[Tuple[QuadraticNumber, ...], ExactMatrix], ...]]:
        return MappingProxyType(
            {term.parameter.name: term.fourier_coefficients for term in self.basis}
        )

    @property
    def gamma_hamiltonians(self) -> Mapping[str, ExactMatrix]:
        return MappingProxyType(
            {term.parameter.name: term.gamma_hamiltonian for term in self.basis}
        )

    def __str__(self) -> str:
        parameters = ", ".join(self.parameter_names) or "no free parameters"
        return (
            f"HamiltonianSpace(shell={self.shell}, parameters=[{parameters}], "
            f"bonds={len(self.bonds)}, exact_covariance={self.covariance_verified})"
        )

    def to_dict(self) -> Dict[str, Any]:
        """Return a user-facing structure without binding serialization tags."""

        return {
            "shell": self.shell,
            "parameter_names": list(self.parameter_names),
            "parameters": [
                {
                    "name": parameter.name,
                    "number": parameter.parameter_number,
                    "constraint_orbit": parameter.constraint_orbit,
                    "orbit_parameter": parameter.orbit_parameter,
                    "representative_bond": parameter.representative_bond,
                }
                for parameter in self.parameters
            ],
            "bonds": [_bond_to_user(bond) for bond in self.bonds],
            "basis": [
                {
                    "parameter": term.parameter.name,
                    "representative_hoppings": [
                        _matrix_to_user(matrix)
                        for matrix in term.representative_hoppings
                    ],
                    "fourier_coefficients": [
                        {
                            "displacement": [
                                _quadratic_to_user(value) for value in displacement
                            ],
                            "matrix": _matrix_to_user(matrix),
                        }
                        for displacement, matrix in term.fourier_coefficients
                    ],
                    "gamma_hamiltonian": _matrix_to_user(term.gamma_hamiltonian),
                    "verified": term.verified,
                }
                for term in self.basis
            ],
            "covariance_verified": self.covariance_verified,
        }

    def to_canonical_dict(self) -> Dict[str, Any]:
        """Return the advanced Rust binding payload for reproducibility."""

        return _thaw_json(self._raw)


def _encode_symbolic_matrix_entries(
    matrix: Sequence[Sequence[HamiltonianExpression]],
) -> list[list[list[Dict[str, Any]]]]:
    return [
        [
            [
                {
                    "parameter": {"name": contribution.parameter_name},
                    "displacement": [
                        _encode_phase_coordinate(value)
                        for value in contribution.displacement
                    ],
                    "coefficient": _encode_cyclotomic_element(
                        contribution.coefficient
                    ),
                }
                for contribution in cell.contributions
            ]
            for cell in row
        ]
        for row in matrix
    ]


def _decode_bond(shell: int, index: int, value: Mapping[str, Any]) -> BondRecord:
    return BondRecord(
        shell=shell,
        bond_id=index,
        source_site=int(value["source_site"]) + 1,
        target_site=int(value["target_site"]) + 1,
        translation=tuple(int(item) for item in value["translation"]),
        displacement=tuple(_decode_quadratic(item) for item in value["displacement"]),
        source_endpoint=tuple(
            _decode_quadratic(item) for item in value["source_endpoint"]
        ),
        target_endpoint=tuple(
            _decode_quadratic(item) for item in value["target_endpoint"]
        ),
        squared_distance=_decode_quadratic(value["squared_distance"]),
        core_bond_index=index - 1,
    )


def _decode_hamiltonian(value: Mapping[str, Any]) -> Hamiltonian:
    if value.get("schema") != "magnetictb.symbolic_hamiltonian.v1":
        raise RuntimeError("Rust returned an unsupported Hamiltonian schema")
    field_spec = _field_from_context(value["field_context"])
    terms = []
    for raw_term in value["terms"]:
        parameter = raw_term["parameter"]
        terms.append(
            HamiltonianTerm(
                parameter_name=str(parameter["name"]),
                shell=int(parameter["shell"]),
                parameter_number=int(parameter["parameter_number"]),
                core_parameter_index=int(parameter["core_parameter_index"]),
                fourier_coefficients=tuple(
                    (
                        tuple(
                            _decode_phase_coordinate(item)
                            for item in coefficient["displacement"]
                        ),
                        _decode_matrix(coefficient["matrix"], field_spec),
                    )
                    for coefficient in raw_term["fourier_coefficients"]
                ),
                gamma_hamiltonian=_decode_matrix(
                    raw_term["gamma_hamiltonian"], field_spec
                ),
                verified=bool(raw_term["verified"]),
            )
        )
    matrix = tuple(
        tuple(
            HamiltonianExpression(
                tuple(
                    HamiltonianMatrixContribution(
                        parameter_name=str(contribution["parameter"]["name"]),
                        displacement=tuple(
                            _decode_phase_coordinate(item)
                            for item in contribution["displacement"]
                        ),
                        coefficient=_decode_element(
                            contribution["coefficient"], field_spec
                        ),
                    )
                    for contribution in cell
                )
            )
            for cell in row
        )
        for row in value["matrix_entries"]
    )
    shape = tuple(int(item) for item in value["shape"])
    if len(shape) != 2:
        raise RuntimeError("Rust returned an invalid Hamiltonian shape")
    return Hamiltonian(
        shape=shape,  # type: ignore[arg-type]
        shells=tuple(int(item) for item in value["shells"]),
        terms=tuple(terms),
        matrix=matrix,
        covariance_verified=bool(value["covariance_verified"]),
        model_identity_sha256=str(value["model_identity_sha256"]),
        field=field_spec,
        gauge=value["gauge"],
        _raw=value,
    )


def _bond_to_user(value: BondRecord) -> Dict[str, Any]:
    return {
        "shell": value.shell,
        "bond_id": value.bond_id,
        "source_site": value.source_site,
        "target_site": value.target_site,
        "translation": list(value.translation),
        "displacement": [_quadratic_to_user(item) for item in value.displacement],
        "squared_distance": _quadratic_to_user(value.squared_distance),
    }


def _parameter_name(shell: int, number: int) -> str:
    prefix = {1: "e", 2: "t", 3: "r", 4: "s"}.get(shell)
    return f"{prefix}{number}" if prefix else f"p{shell}n{number}"
