"""Compatibility facade over domain-structured Rust extension bindings."""

from .errors import (
    DataError,
    ExactError,
    FittingError,
    GroupError,
    MagneticTBError,
    ModelError,
    PropertiesError,
    SymmetryError,
    _translate_core_error,
    _translate_fitting_error,
    _translate_data_error,
    _translate_group_error,
    _translate_properties_error,
    _translate_symmetry_error,
)

from .core_bindings import (
    _decode_result,
    _encode,
    common_kernel,
    common_kernel_json,
    cyclotomic_context,
    fitting,
    fitting_json,
    gapless_points_json,
    geometry,
    geometry_json,
    linear_algebra,
    linear_algebra_json,
    null_space,
    null_space_json,
    properties,
    properties_json,
    representation,
    representation_json,
    tight_binding,
    tight_binding_json,
)

from .data_api import (
    DataCatalog,
)

from .abstract_group import (
    GenerateGroup,
    GroupAction,
    GroupAlgebra,
    OrderedElement,
    OrderedFiniteGroup,
    getGenerator,
    group_action_table_q,
    group_algebra_q,
)
