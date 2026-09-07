"""Compatibility facade for split browser form, Data-option, and init services."""

from .form_parsing import (
    _field_from_form,
    _operation_from_form,
    _parse_lattice_matrix,
    _parse_lattice_scalar,
    _parse_matrix,
    _parse_representation_matrix,
    _parse_scalar,
    _parse_vector,
)

from .data_options_service import (
    _bravais_matrix_to_user,
    _bravais_option,
    _collect_symbols,
    _data_exact_to_user,
    _data_expression_text,
    _msg_family_options,
    _named_parameters,
    _parameter_seed,
    _selected_wyckoff_option,
    _wyckoff_options,
)

from .model_init_service import (
    _friendly_init,
    _missing_explicit_seed,
)
