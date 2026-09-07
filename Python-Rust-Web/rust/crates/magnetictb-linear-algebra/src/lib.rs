#![allow(clippy::missing_errors_doc, clippy::missing_panics_doc)]

mod basis;
mod constraints;
mod matrix;

pub use basis::{
    action_matrix_from_images, action_matrix_in_basis, basis_gram_matrix, coordinates_in_basis,
    dagger_coordinate_matrix, hermitian_matrix_basis, matrix_to_real_coordinates,
    real_coordinates_to_matrix, real_hs_inner_product, rectangular_matrix_units,
    rectangular_real_basis,
};
pub use constraints::{
    ConstraintKernelMethod, ConstraintKernelResult, ConstraintTarget, ConstraintValidationLevel,
    RectangularConstraintOperation, SolveBasisStabilizerResult, assemble_constraint_matrix,
    compile_rectangular_constraint_blocks, invariant_basis,
    rectangular_generator_constraint_matrix, rectangular_transport_matrix, solve_basis_stabilizer,
    solve_basis_stabilizer_actions, solve_constraint_kernel, solve_constraint_kernel_with_method,
    solve_constraint_kernel_with_options, stabilizer_constraint_matrix,
};
pub use matrix::{
    add, block_matrix_2x2, conjugate, conjugate_transpose, exact_hermitian_matrix,
    exact_matrix_equal, exact_square_matrix, exact_unitary_matrix, integer_vector, kronecker,
    real_and_imaginary_parts, scale, subtract, transpose,
};

pub use cyclotomic_nullspace::{ExactError, ExactResult};
