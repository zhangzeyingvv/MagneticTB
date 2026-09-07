"""Internal dependency surface shared by the split Model modules."""

from .core_bindings import geometry
from .data_api import DataCatalog
from .errors import ExactError, ModelError
from .linear_algebra import (
    ExactAtom,
    ExactExpr,
    ExactField,
    ExactMatrix,
    MatrixInput,
    _decode_element,
    _decode_input_exact,
    _decode_lattice_parameters,
    _decode_matrix,
    _decode_nested_elements,
    _decode_phase_coordinate,
    _decode_quadratic,
    _encode_cyclotomic_element,
    _encode_exact,
    _encode_lattpar,
    _field_from_context,
    _freeze_json,
    _matrix,
    _ordered,
    _strict_bool,
    _tagged_list,
    _thaw_json,
    sqrt,
    symbol,
)
from .symmetry import (
    SpinSpaceOperation,
    SymmetryOperation,
    _decode_broken_symmetry_operations,
    _decode_prepared_symmetry_operations,
    _encode_symmetry,
)
from .data import MSG, Wyckoff, _normalise_data_wyckoff, _normalise_wyckoff
from .representation_theory import (
    ExplicitPolynomialBasis,
    _basis_display_labels,
    _encode_basis_functions,
    _encode_representation_information,
    _normalise_basis,
    _normalise_orbital_labels,
)
from .tight_binding import (
    BondRecord,
    Hamiltonian,
    HamiltonianBasisTerm,
    HamiltonianExpression,
    HamiltonianMatrixContribution,
    HamiltonianParameter,
    HamiltonianSpace,
    OrbitalRecord,
    _decode_bond,
    _decode_hamiltonian,
    _encode_symbolic_matrix_entries,
    _parameter_name,
)
