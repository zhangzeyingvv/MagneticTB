"""Model Schema for the Rust-backed model lifecycle."""

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

from .model_dependencies import (
    BondRecord,
    DataCatalog,
    ExactAtom,
    ExactError,
    ExactField,
    ExactMatrix,
    ExplicitPolynomialBasis,
    Hamiltonian,
    HamiltonianBasisTerm,
    HamiltonianExpression,
    HamiltonianMatrixContribution,
    HamiltonianParameter,
    HamiltonianSpace,
    ModelError,
    OrbitalRecord,
    SpinSpaceOperation,
    SymmetryOperation,
    _decode_bond,
    _decode_broken_symmetry_operations,
    _decode_element,
    _decode_hamiltonian,
    _decode_input_exact,
    _decode_lattice_parameters,
    _decode_matrix,
    _decode_nested_elements,
    _decode_phase_coordinate,
    _decode_prepared_symmetry_operations,
    _decode_quadratic,
    _encode_cyclotomic_element,
    _encode_exact,
    _encode_symbolic_matrix_entries,
    _field_from_context,
    _freeze_json,
    _parameter_name,
    _thaw_json,
    geometry,
)


def _format_basis_state_payload(value: Any) -> Any:
    """Render Rust-owned exact polynomial states without exposing JSON tags."""

    if not isinstance(value, Mapping):
        return value
    kind = value.get("kind")
    if kind == "list":
        return tuple(_format_basis_state_payload(item) for item in value["items"])
    if kind != "call":
        return str(_decode_input_exact(value))
    head = str(value["head"]).rsplit("`", 1)[-1]
    arguments = tuple(_format_basis_state_payload(item) for item in value["arguments"])
    if head == "Times":
        if len(arguments) >= 2 and arguments[0] == "-1":
            product = "*".join(str(item) for item in arguments[1:])
            return f"-{product}"
        return "*".join(str(item) for item in arguments)
    if head == "Plus":
        return " + ".join(str(item) for item in arguments).replace("+ -", "- ")
    if head == "Power" and len(arguments) == 2:
        return f"{arguments[0]}^{arguments[1]}"
    return str(_decode_input_exact(value))


@dataclass(frozen=True)
class _CompilationRecipe:
    kind: Literal["explicit", "data"]
    context: Mapping[str, Any]
    compiler_input: Mapping[str, Any]
    data_trace: Optional[Mapping[str, Any]] = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "context", _freeze_json(self.context))
        object.__setattr__(self, "compiler_input", _freeze_json(self.compiler_input))
        if self.data_trace is not None:
            object.__setattr__(self, "data_trace", _freeze_json(self.data_trace))


@dataclass(frozen=True, init=False)
class PreparedModel:
    """One fully prepared, exact Rust model with immutable per-shell results."""

    _field: ExactField = field(repr=False)
    _shell_results: Tuple[Mapping[str, Any], ...] = field(repr=False)
    _basis_labels: Tuple[Tuple[str, ...], ...] = field(repr=False)
    _model_identity_payload: Mapping[str, Any] = field(repr=False)
    _generator_display_records: Tuple[Tuple[str, bool], ...] = field(repr=False)
    _compilation_recipe: Optional[_CompilationRecipe] = field(repr=False)
    _mode_cache: Dict[Tuple[int, bool, str, str], Mapping[str, Any]] = field(repr=False)
    _evaluated_shell_modes: set[Tuple[int, bool, str, str]] = field(repr=False)
    _mode_cache_lock: RLock = field(repr=False)
    group: Mapping[str, Any]
    generator_indices: Tuple[int, ...]
    symmetry_operations: Tuple[Union[SymmetryOperation, SpinSpaceOperation], ...]
    site_orbits: Any
    site_actions: Any
    cell_translations: Any
    representation_matrices: Tuple[ExactMatrix, ...]
    spin_actions: Tuple[ExactMatrix, ...]
    representation_mode: Literal["DirectProduct", "Induced"]
    representation_method: str
    model_lattice: ExactMatrix
    orbital_table: Tuple[OrbitalRecord, ...]
    bond_shells: Tuple[Tuple[BondRecord, ...], ...]

    def __init__(
        self,
        *,
        field_spec: ExactField,
        shell_results: Sequence[Mapping[str, Any]],
        basis_labels: Sequence[Sequence[str]],
        compilation_recipe: Optional[_CompilationRecipe] = None,
    ) -> None:
        if not shell_results:
            raise ValueError("a prepared model requires at least one shell")
        frozen_results = tuple(_freeze_json(result) for result in shell_results)
        frozen_labels = tuple(tuple(labels) for labels in basis_labels)
        object.__setattr__(self, "_field", field_spec)
        object.__setattr__(self, "_shell_results", frozen_results)
        object.__setattr__(self, "_basis_labels", frozen_labels)
        object.__setattr__(self, "_compilation_recipe", compilation_recipe)
        object.__setattr__(self, "_mode_cache", {})
        object.__setattr__(self, "_evaluated_shell_modes", set())
        object.__setattr__(self, "_mode_cache_lock", RLock())
        first = frozen_results[0]
        object.__setattr__(
            self,
            "group",
            _freeze_json(
                {
                    "order": first["group_order"],
                    "multiplication_table": first["multiplication_table"],
                    "antiunitary_flags": first["antiunitary_flags"],
                }
            ),
        )
        object.__setattr__(
            self,
            "generator_indices",
            tuple(int(index) for index in first["generator_indices"]),
        )
        display_records = first.get("generator_display_records")
        if display_records is None:
            display_records = tuple(
                {
                    "label": first["operation_labels"][index],
                    "antiunitary": first["antiunitary_flags"][index],
                }
                for index in self.generator_indices
            )
        object.__setattr__(
            self,
            "_generator_display_records",
            tuple(
                (str(record["label"]), bool(record["antiunitary"]))
                for record in display_records
            ),
        )
        operations = _decode_prepared_symmetry_operations(first, field_spec)
        object.__setattr__(self, "symmetry_operations", operations)
        object.__setattr__(
            self,
            "_model_identity_payload",
            _freeze_json(
                {
                    "schema": "magnetictb.model_identity.v1",
                    "field_context": first["field_context"],
                    "gauge": first["hamiltonian"]["gauge"],
                    "lattice": first["model_lattice"],
                    "multiplication_table": first["multiplication_table"],
                    "antiunitary_flags": first["antiunitary_flags"],
                    "operation_labels": first["operation_labels"],
                    "spatial_actions": first["spatial_actions"],
                    "spin_rotations": first["spin_rotations"],
                    "image_site_indices": first["image_site_indices"],
                    "cell_translations": first["cell_translations"],
                    "local_dimensions": first["local_dimensions"],
                    "spinor_basis_flags": first.get("spinor_basis_flags"),
                    "site_dimensions": first["site_dimensions"],
                    "orbital_layout": first["orbital_layout"],
                    "representation_matrices": first["representation_matrices"],
                }
            ),
        )
        object.__setattr__(
            self,
            "site_orbits",
            _decode_nested_elements(first["site_orbits"], field_spec),
        )
        object.__setattr__(self, "site_actions", _freeze_json(first["image_site_indices"]))
        object.__setattr__(
            self,
            "cell_translations",
            _decode_nested_elements(first["cell_translations"], field_spec),
        )
        object.__setattr__(
            self,
            "representation_matrices",
            tuple(
                _decode_matrix(matrix, field_spec)
                for matrix in first["representation_matrices"]
            ),
        )
        object.__setattr__(
            self,
            "spin_actions",
            tuple(_decode_matrix(matrix, field_spec) for matrix in first["spin_actions"]),
        )
        object.__setattr__(
            self,
            "representation_mode",
            str(first["representation_mode"]),
        )
        object.__setattr__(
            self,
            "representation_method",
            str(first["representation_method"]),
        )
        object.__setattr__(
            self,
            "model_lattice",
            _decode_matrix(first["evaluated_model_lattice"], field_spec),
        )
        object.__setattr__(
            self,
            "orbital_table",
            self._make_orbital_table(
                first["orbital_layout"],
                representation_source=str(first.get("representation_source", "Matrices")),
            ),
        )
        object.__setattr__(
            self,
            "bond_shells",
            tuple(
                tuple(
                    _decode_bond(shell, index + 1, bond)
                    for index, bond in enumerate(result["bonds"])
                )
                for shell, result in enumerate(frozen_results, start=1)
            ),
        )

    @property
    def initial_bond_shells(self) -> int:
        return len(self._shell_results)

    @property
    def model_identity_sha256(self) -> str:
        return str(self._shell_results[0]["model_identity_sha256"])

    @property
    def model_identity_payload(self) -> Mapping[str, Any]:
        """Return the immutable payload audited by the Rust identity hash."""

        return self._model_identity_payload

    def symmetry_operation(
        self, number: int
    ) -> Union[SymmetryOperation, SpinSpaceOperation]:
        if (
            isinstance(number, bool)
            or not isinstance(number, int)
            or number < 1
            or number > len(self.symmetry_operations)
        ):
            raise ModelError(
                "InvalidSymmetryOperationIndex",
                f"operation {number} is outside 1..{len(self.symmetry_operations)}",
            )
        return self.symmetry_operations[number - 1]

    def to_canonical_dict(self, shell: int = 1) -> Dict[str, Any]:
        """Return one advanced, zero-based Rust result for exact comparison."""

        if (
            isinstance(shell, bool)
            or not isinstance(shell, int)
            or shell < 1
            or shell > len(self._shell_results)
        ):
            raise ModelError(
                "ShellNotPrepared",
                f"shell {shell} was not prepared; available shells are 1..{len(self._shell_results)}",
            )
        return _thaw_json(self._shell_results[shell - 1])

    def _stable_session_view(self) -> Mapping[str, Any]:
        """Return the immutable Python view of the stable 2.0.10 Association."""

        first = self._shell_results[0]
        source = "BasisFunctions"
        recipe = self._compilation_recipe
        if recipe is not None and recipe.compiler_input.get("kind") == "ordered_association":
            for entry in recipe.compiler_input["entries"]:
                if entry["key"] == "RepresentationSource":
                    source = str(entry["value"])
                    break
        model_specification = MappingProxyType(
            {
                "Lattice": self.model_lattice,
                "SiteOrbits": self.site_orbits,
                "LocalDimensions": tuple(int(value) for value in first["local_dimensions"]),
                "LocalBases": self._basis_labels,
                "BondClasses": self.bond_shells,
            }
        )
        full_representation = MappingProxyType(
            {
                "Method": self.representation_method,
                "Convention": "row-vector-site-action",
                "RepresentationMatrices": self.representation_matrices,
                "Dimension": len(self.orbital_table),
                "UnitaryVerified": bool(
                    first.get("parameter_space", {}).get(
                        "full_parameter_space_verified", True
                    )
                ),
            }
        )
        compiled_shells = tuple(
            MappingProxyType(
                {
                    "Shell": shell,
                    "Bonds": bonds,
                    "BondCount": len(bonds),
                }
            )
            for shell, bonds in enumerate(self.bond_shells, start=1)
        )
        return MappingProxyType(
            {
                "Schema": "MagneticTBModelSession",
                "SchemaVersion": 7,
                "RepresentationMode": self.representation_mode,
                "RepresentationSource": source,
                "BasisOrderingData": self.orbital_table,
                "ModelSpecification": model_specification,
                "FullRepresentation": full_representation,
                "GeneratorIndices": tuple(index + 1 for index in self.generator_indices),
                "BondSearch": MappingProxyType(
                    {"PreparedShellCount": self.initial_bond_shells}
                ),
                "CompiledBondShells": compiled_shells,
                "BondConstraintCache": MappingProxyType({}),
                "SolvedShellCache": MappingProxyType({}),
                "RealSpaceShellCache": MappingProxyType({}),
            }
        )

    def __getitem__(self, key: Union[str, Tuple[str, ...]]) -> Any:
        """Read stable Association keys without exposing mutable session state."""

        keys = key if isinstance(key, tuple) else (key,)
        value: Any = self._stable_session_view()
        for name in keys:
            if not isinstance(name, str):
                raise TypeError("stable session keys must be strings")
            value = value[name]
        return value

    def to_dict(self) -> Dict[str, Any]:
        """Return the ordinary Python form of the stable session/compiler view."""

        return _thaw_json(self._stable_session_view())

    def _make_orbital_table(
        self,
        layout: Sequence[Mapping[str, Any]],
        *,
        representation_source: str,
    ) -> Tuple[OrbitalRecord, ...]:
        records = []
        for item in layout:
            orbit = int(item["orbit_index"])
            local = int(item["local_orbital_index"])
            labels = self._basis_labels[orbit] if orbit < len(self._basis_labels) else ()
            label = labels[local] if local < len(labels) else f"orb{local + 1}"
            reference = item.get("reference_site_index")
            transport = item.get("transport_operation_index")
            reference_equivalent_atom = (
                None if reference is None else int(reference) + 1
            )
            transport_operation = None if transport is None else int(transport) + 1
            transport_label = item.get("transport_operation_label")
            if "basis_state" in item:
                basis_state = _format_basis_state_payload(item["basis_state"])
            elif self.representation_mode == "Induced" and representation_source == "Matrices":
                if transport_operation is None:
                    raise RuntimeError(
                        "Rust omitted induced transport data for an explicit representation"
                    )
                basis_state = MappingProxyType(
                    {
                        "AbstractOrbitalLabel": label,
                        "TransportOperationIndex": transport_operation,
                    }
                )
            else:
                basis_state = label
            spatial_orbital = (
                None
                if item.get("spatial_orbital") is None
                else _format_basis_state_payload(item["spatial_orbital"])
            )
            spin_state = (
                None
                if item.get("spin_state") is None
                else _format_basis_state_payload(item["spin_state"])
            )
            spin_structure = str(
                item.get(
                    "spin_structure",
                    "Spinless" if spin_state is None else "GeneralInternalState",
                )
            )
            records.append(
                OrbitalRecord(
                    orbital_id=int(item["orbital_index"]) + 1,
                    wyckoff_orbit=orbit + 1,
                    equivalent_atom=int(item["site_in_orbit"]) + 1,
                    site_id=int(item["site_index"]) + 1,
                    local_orbital=local + 1,
                    fractional_position=tuple(
                        _decode_element(value, self._field)
                        for value in item["fractional_position"]
                    ),
                    cartesian_position=tuple(
                        _decode_element(value, self._field)
                        for value in item["cartesian_position"]
                    ),
                    basis_function=label,
                    basis_state=basis_state,
                    spatial_orbital=spatial_orbital,
                    spin_structure=spin_structure,
                    basis_convention=self.representation_mode,
                    spin_state=spin_state,
                    reference_equivalent_atom=reference_equivalent_atom,
                    transport_operation=transport_operation,
                    transport_operation_label=(
                        None if transport_label is None else str(transport_label)
                    ),
                    core_orbital_index=int(item["orbital_index"]),
                )
            )
        return tuple(records)

    def hamiltonian_space(
        self,
        shell: int,
        *,
        hermitian: bool = True,
        cartesian_coordinates: bool = False,
        validation_level: Literal["None", "Basic", "Full"] = "Basic",
        kernel_method: Literal["Iterative", "Stacked", "Cyclotomic"] = "Iterative",
    ) -> HamiltonianSpace:
        self._validate_symham_options(
            shell,
            hermitian=hermitian,
            cartesian_coordinates=cartesian_coordinates,
            validation_level=validation_level,
            kernel_method=kernel_method,
        )
        if cartesian_coordinates:
            raise ModelError(
                "UnsupportedCoordinateMode",
                "CartesianCoordinates changes the symbolic symham matrix only; "
                "HamiltonianSpace retains exact fractional bond diagnostics",
            )
        raw = self._shell_result_for_options(
            shell,
            hermitian=hermitian,
            validation_level=validation_level,
            kernel_method=kernel_method,
        )
        parameter_records = []
        for item in raw["parameter_order"]:
            index = int(item["parameter_index"])
            parameter_records.append(
                HamiltonianParameter(
                    name=_parameter_name(shell, index + 1),
                    parameter_number=index + 1,
                    constraint_orbit=int(item["constraint_orbit_index"]) + 1,
                    orbit_parameter=int(item["orbit_parameter_index"]) + 1,
                    representative_bond=int(item["representative_bond_index"]) + 1,
                    core_parameter_index=index,
                )
            )
        parameters = tuple(parameter_records)
        solutions = raw["parameter_space"]["parameter_basis_solutions"]
        basis = tuple(
            HamiltonianBasisTerm(
                parameter=parameters[int(solution["parameter_index"])],
                representative_hoppings=tuple(
                    _decode_matrix(matrix, self._field)
                    for matrix in solution["representative_hoppings"]
                ),
                fourier_coefficients=tuple(
                    (
                        tuple(
                            _decode_quadratic(value)
                            for value in coefficient["displacement"]
                        ),
                        _decode_matrix(coefficient["matrix"], self._field),
                    )
                    for coefficient in solution["fourier_coefficients"]
                ),
                gamma_hamiltonian=_decode_matrix(
                    solution["gamma_hamiltonian"], self._field
                ),
                verified=bool(solution["verified"]),
            )
            for solution in solutions
        )
        verified = bool(raw["parameter_space"]["full_parameter_space_verified"])
        return HamiltonianSpace(
            shell=shell,
            parameters=parameters,
            bonds=self.bond_shells[shell - 1],
            basis=basis,
            covariance_verified=verified and all(item.verified for item in basis),
            _raw=raw,
        )

    def hamiltonian(
        self,
        shell: int,
        *,
        hermitian: bool = True,
        cartesian_coordinates: bool = False,
        validation_level: Literal["None", "Basic", "Full"] = "Basic",
        kernel_method: Literal["Iterative", "Stacked", "Cyclotomic"] = "Iterative",
    ) -> Hamiltonian:
        """Return the advanced Rust-owned single-shell Hamiltonian object."""

        self._validate_symham_options(
            shell,
            hermitian=hermitian,
            cartesian_coordinates=cartesian_coordinates,
            validation_level=validation_level,
            kernel_method=kernel_method,
        )
        raw = self._shell_result_for_options(
            shell,
            hermitian=hermitian,
            validation_level=validation_level,
            kernel_method=kernel_method,
        )
        with self._mode_cache_lock:
            self._evaluated_shell_modes.add(
                (shell, hermitian, kernel_method, validation_level)
            )
        raw_hamiltonian = raw["hamiltonian"]
        if cartesian_coordinates:
            raw_hamiltonian = geometry(
                "cartesian_symbolic_hamiltonian",
                hamiltonian=_thaw_json(raw_hamiltonian),
                lattice=_thaw_json(raw["model_lattice"]),
                symbolic_lattice=self._symbolic_lattice(raw),
            )
        return _decode_hamiltonian(raw_hamiltonian)

    def _cached_hopping_shell_results(
        self,
        shells: Sequence[int],
        *,
        hermitian: bool,
        kernel_method: str,
        validation_level: str,
    ) -> list[dict[str, Any]]:
        keys = [
            (shell, hermitian, kernel_method, validation_level)
            for shell in shells
        ]
        with self._mode_cache_lock:
            missing = [
                shell
                for shell, key in zip(shells, keys)
                if key not in self._evaluated_shell_modes
            ]
        if missing:
            raise ModelError(
                "UnsolvedShells",
                "evaluate symham for shell(s) "
                + ", ".join(str(shell) for shell in missing)
                + " with the same Hermitian, KernelMethod, and ValidationLevel first",
            )
        return [
            _thaw_json(
                self._shell_result_for_options(
                    shell,
                    hermitian=hermitian,
                    kernel_method=kernel_method,
                    validation_level=validation_level,
                )
            )
            for shell in shells
        ]

    def symham(
        self,
        shell: int,
        *,
        hermitian: bool = True,
        cartesian_coordinates: bool = False,
        validation_level: Literal["None", "Basic", "Full"] = "Basic",
        kernel_method: Literal["Iterative", "Stacked", "Cyclotomic"] = "Iterative",
    ) -> list[list[HamiltonianExpression]]:
        """Return the stable-style symbolic Hamiltonian matrix for one shell."""

        return self.hamiltonian(
            shell,
            hermitian=hermitian,
            cartesian_coordinates=cartesian_coordinates,
            validation_level=validation_level,
            kernel_method=kernel_method,
        ).to_list()

    def unconstrained_hamiltonian(self, shell: int) -> Hamiltonian:
        """Return the Rust-built unconstrained stable ``unsymham`` object."""

        self._validate_symham_options(
            shell,
            hermitian=True,
            cartesian_coordinates=False,
            validation_level="Basic",
            kernel_method="Iterative",
        )
        try:
            raw = geometry(
                "unconstrained_symbolic_hamiltonian",
                model_result=self.to_canonical_dict(shell),
                shell=shell,
            )
        except ExactError as error:
            raise ModelError(error.tag, error.detail) from None
        return _decode_hamiltonian(raw)

    def unsymham(self, shell: int) -> list[list[HamiltonianExpression]]:
        """Return the stable unconstrained symbolic matrix for one shell."""

        return self.unconstrained_hamiltonian(shell).to_list()

    def symhamII(
        self,
        ham: Union[Hamiltonian, Sequence[Sequence[HamiltonianExpression]]],
        *,
        wcc: Optional[Sequence[Sequence[ExactAtom]]] = None,
    ) -> list[list[HamiltonianExpression]]:
        """Rewrite one symbolic Hamiltonian from stable convention I to II."""

        return _symham_ii_for_model(self, ham, wcc=wcc)

    def brokenSymmetryInitRules(self, selector: Any) -> Dict[str, Any]:
        """Return init-ready Python options for the retained stable subgroup."""

        return _broken_symmetry_rules_for_model(self, selector)

    def _symbolic_lattice(self, raw: Mapping[str, Any]) -> Any:
        recipe = self._compilation_recipe
        if recipe is not None:
            compiler_input = recipe.compiler_input
            if compiler_input.get("kind") == "ordered_association":
                for entry in compiler_input["entries"]:
                    if entry["key"] == "Lattice":
                        return _thaw_json(entry["value"])
        return _thaw_json(raw["model_lattice"])

    def _validate_symham_options(
        self,
        shell: int,
        *,
        hermitian: bool,
        cartesian_coordinates: bool,
        validation_level: str,
        kernel_method: str,
    ) -> None:
        if (
            isinstance(shell, bool)
            or not isinstance(shell, int)
            or shell < 1
            or shell > len(self._shell_results)
        ):
            raise ModelError(
                "ShellNotPrepared",
                f"shell {shell} was not prepared; available shells are 1..{len(self._shell_results)}",
            )
        if not isinstance(hermitian, bool):
            raise ModelError("InvalidHermitianOption", "hermitian must be bool")
        if not isinstance(cartesian_coordinates, bool):
            raise ModelError(
                "InvalidCartesianCoordinates",
                "cartesian_coordinates must be bool",
            )
        if validation_level not in ("None", "Basic", "Full"):
            raise ModelError(
                "InvalidValidationLevel",
                "validation_level must be 'None', 'Basic', or 'Full'",
            )
        if kernel_method not in ("Iterative", "Stacked", "Cyclotomic"):
            raise ModelError(
                "InvalidKernelMethod",
                "kernel_method must be 'Iterative', 'Stacked', or 'Cyclotomic'",
            )

    def _shell_result_for_options(
        self,
        shell: int,
        *,
        hermitian: bool,
        validation_level: str,
        kernel_method: str,
    ) -> Mapping[str, Any]:
        key = (shell, hermitian, kernel_method, validation_level)
        if key == (shell, True, "Iterative", "Basic"):
            return self._shell_results[shell - 1]
        with self._mode_cache_lock:
            cached = self._mode_cache.get(key)
            if cached is not None:
                return cached
            recipe = self._compilation_recipe
            if recipe is None:
                raise ModelError(
                    "ModelRecompilationUnavailable",
                    "this advanced PreparedModel has no stable input recipe for option recompilation",
                )
            try:
                if recipe.kind == "explicit":
                    result = geometry(
                        "compile_model_input_closure",
                        context=_thaw_json(recipe.context),
                        compiler_input=_thaw_json(recipe.compiler_input),
                        target_shell=shell,
                        hermitian=hermitian,
                        kernel_method=kernel_method,
                        validation_level=validation_level,
                    )
                else:
                    if recipe.data_trace is None:
                        raise RuntimeError("Data compilation recipe is missing data_trace")
                    result = DataCatalog().compile_model_input_closure(
                        _thaw_json(recipe.context),
                        _thaw_json(recipe.compiler_input),
                        _thaw_json(recipe.data_trace),
                        shell,
                        hermitian=hermitian,
                        kernel_method=kernel_method,
                        validation_level=validation_level,
                    )
            except ExactError as error:
                raise ModelError(error.tag, error.detail) from None
            frozen = _freeze_json(result)
            self._mode_cache[key] = frozen
            return frozen


def _symham_ii_for_model(
    model: PreparedModel,
    ham: Union[Hamiltonian, Sequence[Sequence[HamiltonianExpression]]],
    *,
    wcc: Optional[Sequence[Sequence[ExactAtom]]],
) -> list[list[HamiltonianExpression]]:
    if isinstance(ham, Hamiltonian):
        raw = ham.to_canonical_dict()
        context = raw["field_context"]
        field_spec = ham.field
        matrix_entries = raw["matrix_entries"]
        dimension = ham.shape[0]
    else:
        if (
            isinstance(ham, (str, bytes))
            or not isinstance(ham, Sequence)
            or not ham
        ):
            raise TypeError("ham must be a nonempty symbolic Hamiltonian matrix")
        rows = tuple(tuple(row) for row in ham)
        dimension = len(rows)
        if any(len(row) != dimension for row in rows):
            raise ValueError("symhamII requires a square Hamiltonian matrix")
        if any(
            not isinstance(cell, HamiltonianExpression)
            for row in rows
            for cell in row
        ):
            raise TypeError(
                "symhamII matrix entries must be HamiltonianExpression values from Rust"
            )
        field_spec = model._field
        context = field_spec._context()
        matrix_entries = _encode_symbolic_matrix_entries(rows)
    if dimension != len(model.orbital_table) and wcc is None:
        raise ModelError(
            "InvalidWannierCenters",
            "the current model does not provide one Wannier center per Hamiltonian row",
        )
    if wcc is None:
        centers = [
            [_encode_cyclotomic_element(value) for value in orbital.fractional_position]
            for orbital in model.orbital_table
        ]
    else:
        if isinstance(wcc, (str, bytes)) or not isinstance(wcc, Sequence):
            raise TypeError("wcc must be one ordered fractional three-vector per row")
        center_rows = tuple(tuple(center) for center in wcc)
        if len(center_rows) != dimension or any(len(center) != 3 for center in center_rows):
            raise ValueError(
                "wcc must contain one fractional three-vector per Hamiltonian row"
            )
        centers = [[_encode_exact(value) for value in center] for center in center_rows]
    try:
        transformed = geometry(
            "symham_ii",
            context=context,
            matrix_entries=matrix_entries,
            centers=centers,
        )
    except ExactError as error:
        raise ModelError(error.tag, error.detail) from None
    output_field = _field_from_context(transformed["field_context"])
    return [
        [
            HamiltonianExpression(
                tuple(
                    HamiltonianMatrixContribution(
                        parameter_name=str(contribution["parameter"]["name"]),
                        displacement=tuple(
                            _decode_phase_coordinate(value)
                            for value in contribution["displacement"]
                        ),
                        coefficient=_decode_element(
                            contribution["coefficient"], output_field
                        ),
                    )
                    for contribution in cell
                )
            )
            for cell in row
        ]
        for row in transformed["matrix_entries"]
    ]


def _decode_broken_symmetry_basis(
    value: Mapping[str, Any],
    source_orbit_indices: Sequence[int],
    model: PreparedModel,
) -> Tuple[Any, ...]:
    if value.get("kind") != "list" or len(value["items"]) != len(source_orbit_indices):
        raise RuntimeError("Rust returned misaligned broken-symmetry basis data")
    result = []
    for raw, source_orbit in zip(value["items"], source_orbit_indices):
        if raw.get("kind") == "list":
            result.append(tuple(str(label) for label in raw["items"]))
            continue
        if raw.get("kind") != "matrix":
            raise RuntimeError("Rust returned an unsupported basis specification")
        rows = _decode_input_exact(raw)
        columns = len(rows[0])
        functions = tuple(row[0] if columns == 1 else tuple(row) for row in rows)
        result.append(
            ExplicitPolynomialBasis(
                functions,
                labels=model._basis_labels[int(source_orbit)],
            )
        )
    return tuple(result)


def _broken_symmetry_rules_for_model(model: PreparedModel, selector: Any) -> Dict[str, Any]:
    recipe = model._compilation_recipe
    if recipe is None:
        raise ModelError(
            "ModelRecompilationUnavailable",
            "brokenSymmetryInitRules requires the original stable init recipe",
        )
    if selector == "All":
        selected: list[Any] = list(range(len(model.symmetry_operations)))
    else:
        items = selector if isinstance(selector, Sequence) and not isinstance(selector, (str, bytes)) else (selector,)
        selected = []
        for item in items:
            if isinstance(item, bool):
                raise TypeError("retained generators must be 1-based indices or labels")
            if isinstance(item, int):
                if item < 1:
                    raise ValueError("retained generator indices are 1-based positive integers")
                selected.append(item - 1)
            elif isinstance(item, str) and item:
                selected.append(item)
            else:
                raise TypeError("retained generators must be 1-based indices or labels")
    try:
        raw = geometry(
            "broken_symmetry_init_rules",
            model_result=model.to_canonical_dict(1),
            compiler_input=_thaw_json(recipe.compiler_input),
            selector=selected,
        )
    except ExactError as error:
        raise ModelError(error.tag, error.detail) from None
    compiler_input = raw["compiler_input"]
    entries = {entry["key"]: entry["value"] for entry in compiler_input["entries"]}
    wyckoff = _decode_input_exact(entries["WyckoffPosition"])
    source_orbits = tuple(int(index) for index in raw["source_orbit_indices"])
    return {
        "lattice": _decode_input_exact(entries["Lattice"]),
        "lattpar": _decode_lattice_parameters(entries["LatticeParameters"]),
        "wyckoffposition": wyckoff,
        "symminformation": _decode_broken_symmetry_operations(
            entries["SymmetryInformation"]
        ),
        "basis_functions": _decode_broken_symmetry_basis(
            entries["BasisFunctions"], source_orbits, model
        ),
        "debug_q": bool(entries["Debug"]),
        "initial_bond_shells": int(entries["InitialBondShells"]),
        "generate_symmetry_group": False,
        "representation_mode": str(entries["RepresentationMode"]),
        "exact_field": model._field,
    }
