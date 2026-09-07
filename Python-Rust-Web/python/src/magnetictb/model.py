"""Compatibility facade for the split Model schema, preparation, and session APIs."""

from .model_schema import (
    PreparedModel,
    _CompilationRecipe,
    _broken_symmetry_rules_for_model,
    _decode_broken_symmetry_basis,
    _symham_ii_for_model,
)

from .model_preparation import (
    _base_compiler_input,
    _compile_all_shells,
    _compile_data_prepared_model,
    _prepare_data_model,
    _prepare_explicit_model,
    prepare_model,
    prepare_model_from_rep,
)

from .model_session import (
    CompileMagneticTBInput,
    CurrentModelSession,
    _CURRENT_MODEL,
    _clear_current_model_for_tests,
    _print_generators,
    _stable_default_lattice,
    _stable_default_lattpar,
    _stable_default_symmetry,
    _stable_default_wyckoff,
    _stable_keyword_alias,
    brokenSymmetryInitRules,
    compile_magnetic_tb_input,
    current_model,
    hamiltonian_object,
    hamiltonian_space,
    init,
    initfromrep,
    symham,
    symhamII,
    unsymham,
)
