use crate::basis::{action_matrix_in_basis, dagger_coordinate_matrix};
use crate::matrix::{
    block_matrix_2x2, conjugate, exact_hermitian_matrix, exact_square_matrix, kronecker,
    real_and_imaginary_parts, scale, subtract, transpose,
};
use cyclotomic_nullspace::{
    CommonKernelResult, CompiledProblem, Element, ExactError, ExactMatrix, ExactResult,
    common_kernel, null_space, vertical_stack,
};
use std::sync::Arc;

fn dimension_error(detail: impl Into<String>) -> ExactError {
    ExactError::new("DimensionMismatch", detail)
}

fn constraint_error(detail: impl Into<String>) -> ExactError {
    ExactError::new("InvalidConstraintData", detail)
}

fn negate(matrix: &ExactMatrix) -> ExactResult<ExactMatrix> {
    scale(matrix, &Element::one(matrix.context())?.negate()?)
}

fn real_imaginary_matrices(matrix: &ExactMatrix) -> ExactResult<(ExactMatrix, ExactMatrix)> {
    let parts = matrix
        .entries()
        .iter()
        .map(real_and_imaginary_parts)
        .collect::<ExactResult<Vec<_>>>()?;
    let real = ExactMatrix::new(
        Arc::clone(matrix.context()),
        matrix.rows(),
        matrix.columns(),
        parts.iter().map(|(real, _)| real.clone()).collect(),
    )?;
    let imaginary = ExactMatrix::new(
        Arc::clone(matrix.context()),
        matrix.rows(),
        matrix.columns(),
        parts.into_iter().map(|(_, imaginary)| imaginary).collect(),
    )?;
    Ok((real, imaginary))
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ConstraintTarget {
    Same,
    Reverse,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RectangularConstraintOperation {
    pub left: ExactMatrix,
    pub right: ExactMatrix,
    pub antiunitary: bool,
    pub target: ConstraintTarget,
}

pub fn rectangular_transport_matrix(
    left: &ExactMatrix,
    right: &ExactMatrix,
    antiunitary: bool,
) -> ExactResult<ExactMatrix> {
    if !exact_square_matrix(left) || !exact_square_matrix(right) {
        return Err(dimension_error(
            "rectangular transport representation matrices must be nonempty and square",
        ));
    }
    if left.context() != right.context() {
        return Err(ExactError::new(
            "ConductorMismatch",
            "rectangular transport contexts differ",
        ));
    }
    let complex_action = kronecker(left, &conjugate(right)?)?;
    let (real, imaginary) = real_imaginary_matrices(&complex_action)?;
    if antiunitary {
        block_matrix_2x2(&real, &imaginary, &imaginary, &negate(&real)?)
    } else {
        block_matrix_2x2(&real, &negate(&imaginary)?, &imaginary, &real)
    }
}

pub fn rectangular_generator_constraint_matrix(
    left: &ExactMatrix,
    right: &ExactMatrix,
) -> ExactResult<ExactMatrix> {
    if !exact_square_matrix(left) || !exact_square_matrix(right) {
        return Err(dimension_error(
            "rectangular generator matrices must be nonempty and square",
        ));
    }
    if left.context() != right.context() {
        return Err(ExactError::new(
            "ConductorMismatch",
            "rectangular generator contexts differ",
        ));
    }
    let left_term = kronecker(left, &ExactMatrix::identity(left.context(), right.rows())?)?;
    let right_term = kronecker(
        &ExactMatrix::identity(left.context(), left.rows())?,
        &transpose(right)?,
    )?;
    let complex_constraint = subtract(&left_term, &right_term)?;
    let (real, imaginary) = real_imaginary_matrices(&complex_constraint)?;
    block_matrix_2x2(&real, &negate(&imaginary)?, &imaginary, &real)
}

pub fn compile_rectangular_constraint_blocks(
    rows: usize,
    columns: usize,
    operations: &[RectangularConstraintOperation],
    continuous_pair: Option<(&ExactMatrix, &ExactMatrix)>,
) -> ExactResult<Vec<ExactMatrix>> {
    if rows == 0 || columns == 0 {
        return Err(dimension_error(
            "rectangular constraint dimensions must be positive",
        ));
    }
    let coordinate_dimension = 2usize
        .checked_mul(rows)
        .and_then(|value| value.checked_mul(columns))
        .ok_or_else(|| dimension_error("coordinate dimension overflows"))?;
    let context = operations
        .first()
        .map(|operation| Arc::clone(operation.left.context()))
        .or_else(|| continuous_pair.map(|pair| Arc::clone(pair.0.context())))
        .ok_or_else(|| {
            constraint_error("at least one finite operation or a continuous pair is required")
        })?;
    let identity = ExactMatrix::identity(&context, coordinate_dimension)?;
    let dagger = dagger_coordinate_matrix(&context, rows, columns)?;
    let mut blocks = Vec::with_capacity(operations.len() + usize::from(continuous_pair.is_some()));
    for operation in operations {
        if operation.left.context() != &context
            || operation.right.context() != &context
            || operation.left.rows() != rows
            || operation.left.columns() != rows
            || operation.right.rows() != columns
            || operation.right.columns() != columns
            || (operation.target == ConstraintTarget::Reverse && rows != columns)
        {
            return Err(constraint_error(
                "finite rectangular constraint operation has incompatible data",
            ));
        }
        let action =
            rectangular_transport_matrix(&operation.left, &operation.right, operation.antiunitary)?;
        blocks.push(match operation.target {
            ConstraintTarget::Same => subtract(&action, &identity)?,
            ConstraintTarget::Reverse => subtract(&action, &dagger)?,
        });
    }
    if let Some((left, right)) = continuous_pair {
        if left.context() != &context
            || right.context() != &context
            || left.rows() != rows
            || left.columns() != rows
            || right.rows() != columns
            || right.columns() != columns
            || !exact_hermitian_matrix(left)?
            || !exact_hermitian_matrix(right)?
        {
            return Err(constraint_error(
                "continuous rectangular generator pair has incompatible data",
            ));
        }
        blocks.insert(0, rectangular_generator_constraint_matrix(left, right)?);
    }
    Ok(blocks)
}

pub fn assemble_constraint_matrix(
    context: &Arc<cyclotomic_nullspace::CyclotomicContext>,
    blocks: &[ExactMatrix],
    coordinate_dimension: usize,
) -> ExactResult<ExactMatrix> {
    if coordinate_dimension == 0 {
        return Err(dimension_error("coordinate dimension must be positive"));
    }
    vertical_stack(blocks, context, coordinate_dimension)
}

// These independent flags mirror the stable Mathematica diagnostic fields;
// collapsing them into one state would lose ValidationLevel semantics.
#[allow(clippy::struct_excessive_bools)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ConstraintKernelResult {
    pub coordinate_dimension: usize,
    pub iteration_nullities: Vec<usize>,
    pub constraint_blocks: Vec<ExactMatrix>,
    pub constraint_matrix: ExactMatrix,
    pub nullspace_rows: ExactMatrix,
    pub basis_matrix: ExactMatrix,
    pub rank: usize,
    pub nullity: usize,
    pub residual: ExactMatrix,
    pub exact_residual_verified: bool,
    pub independent_verified: bool,
    pub rank_nullity_verified: bool,
    pub method: ConstraintKernelMethod,
    pub validation_level: ConstraintValidationLevel,
    pub residual_computed: bool,
    pub independent_checked: bool,
    pub rank_nullity_checked: bool,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ConstraintKernelMethod {
    Iterative,
    Stacked,
    Cyclotomic,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ConstraintValidationLevel {
    None,
    Basic,
    Full,
}

fn constraint_result(
    blocks: &[ExactMatrix],
    coordinate_dimension: usize,
    result: CommonKernelResult,
    method: ConstraintKernelMethod,
    validation_level: ConstraintValidationLevel,
) -> ExactResult<ConstraintKernelResult> {
    let context = Arc::clone(result.basis_matrix.context());
    let constraint_matrix = assemble_constraint_matrix(&context, blocks, coordinate_dimension)?;
    let residual_computed = validation_level != ConstraintValidationLevel::None;
    let residual = if residual_computed {
        cyclotomic_nullspace::multiply(&constraint_matrix, &result.basis_matrix)?
    } else {
        ExactMatrix::zero(&context, 0, result.nullity)?
    };
    let exact_residual_verified =
        residual_computed && residual.is_zero() && result.exact_residual_verified;
    let independent_checked = validation_level == ConstraintValidationLevel::Full;
    let independent_verified = independent_checked
        && (result.nullity == 0
            || cyclotomic_nullspace::rref(&result.nullspace_rows)?.rank == result.nullity);
    let rank_nullity_checked = validation_level == ConstraintValidationLevel::Full;
    let independent_rank = if rank_nullity_checked {
        cyclotomic_nullspace::rref(&constraint_matrix)?.rank
    } else {
        result.rank
    };
    let rank_nullity_verified =
        rank_nullity_checked && independent_rank + result.nullity == coordinate_dimension;
    let verified = match validation_level {
        ConstraintValidationLevel::None => true,
        ConstraintValidationLevel::Basic => exact_residual_verified,
        ConstraintValidationLevel::Full => {
            exact_residual_verified && independent_verified && rank_nullity_verified
        }
    };
    if !verified {
        return Err(ExactError::new(
            "VerificationFailure",
            "exact constraint-kernel certification failed",
        ));
    }
    Ok(ConstraintKernelResult {
        coordinate_dimension,
        iteration_nullities: result.iteration_nullities,
        constraint_blocks: blocks.to_vec(),
        constraint_matrix,
        nullspace_rows: result.nullspace_rows,
        basis_matrix: result.basis_matrix,
        rank: independent_rank,
        nullity: result.nullity,
        residual,
        exact_residual_verified,
        independent_verified,
        rank_nullity_verified,
        method,
        validation_level,
        residual_computed,
        independent_checked,
        rank_nullity_checked,
    })
}

pub fn solve_constraint_kernel(
    context: &Arc<cyclotomic_nullspace::CyclotomicContext>,
    blocks: &[ExactMatrix],
    coordinate_dimension: usize,
) -> ExactResult<ConstraintKernelResult> {
    solve_constraint_kernel_with_method(
        context,
        blocks,
        coordinate_dimension,
        ConstraintKernelMethod::Iterative,
    )
}

pub fn solve_constraint_kernel_with_method(
    context: &Arc<cyclotomic_nullspace::CyclotomicContext>,
    blocks: &[ExactMatrix],
    coordinate_dimension: usize,
    method: ConstraintKernelMethod,
) -> ExactResult<ConstraintKernelResult> {
    solve_constraint_kernel_with_options(
        context,
        blocks,
        coordinate_dimension,
        method,
        ConstraintValidationLevel::Full,
    )
}

pub fn solve_constraint_kernel_with_options(
    context: &Arc<cyclotomic_nullspace::CyclotomicContext>,
    blocks: &[ExactMatrix],
    coordinate_dimension: usize,
    method: ConstraintKernelMethod,
    validation_level: ConstraintValidationLevel,
) -> ExactResult<ConstraintKernelResult> {
    if coordinate_dimension == 0 {
        return Err(dimension_error("coordinate dimension must be positive"));
    }
    for block in blocks {
        if block.context() != context || block.columns() != coordinate_dimension {
            return Err(constraint_error(
                "constraint block has incompatible context or column count",
            ));
        }
    }
    let problem = CompiledProblem {
        context: Arc::clone(context),
        coordinate_dimension,
        constraints: blocks.to_vec(),
    };
    let result = match method {
        ConstraintKernelMethod::Iterative => common_kernel(&problem)?,
        ConstraintKernelMethod::Cyclotomic | ConstraintKernelMethod::Stacked
            if blocks.is_empty() =>
        {
            common_kernel(&problem)?
        }
        ConstraintKernelMethod::Cyclotomic => {
            let mut result = common_kernel(&problem)?;
            result.iteration_nullities = vec![coordinate_dimension, result.nullity];
            result
        }
        ConstraintKernelMethod::Stacked => {
            let stacked = assemble_constraint_matrix(context, blocks, coordinate_dimension)?;
            let kernel = null_space(&stacked)?;
            CommonKernelResult {
                basis_matrix: kernel.basis_matrix,
                nullspace_rows: kernel.nullspace_rows,
                iteration_nullities: vec![coordinate_dimension, kernel.nullity],
                rank: kernel.rank,
                nullity: kernel.nullity,
                exact_residual_verified: kernel.exact_residual_verified,
            }
        }
    };
    constraint_result(
        blocks,
        coordinate_dimension,
        result,
        method,
        validation_level,
    )
}

fn validate_actions(actions: &[ExactMatrix]) -> ExactResult<usize> {
    let Some(first) = actions.first() else {
        return Err(constraint_error("action list must be nonempty"));
    };
    if !exact_square_matrix(first)
        || actions.iter().any(|action| {
            !exact_square_matrix(action)
                || action.rows() != first.rows()
                || action.context() != first.context()
        })
    {
        return Err(constraint_error(
            "actions must be same-sized nonempty square matrices in one context",
        ));
    }
    Ok(first.rows())
}

fn stabilizer_blocks(actions: &[ExactMatrix]) -> ExactResult<Vec<ExactMatrix>> {
    let dimension = validate_actions(actions)?;
    let identity = ExactMatrix::identity(actions[0].context(), dimension)?;
    actions
        .iter()
        .map(|action| subtract(action, &identity))
        .collect()
}

pub fn stabilizer_constraint_matrix(actions: &[ExactMatrix]) -> ExactResult<ExactMatrix> {
    let dimension = validate_actions(actions)?;
    let blocks = stabilizer_blocks(actions)?;
    assemble_constraint_matrix(actions[0].context(), &blocks, dimension)
}

pub fn invariant_basis(actions: &[ExactMatrix]) -> ExactResult<ExactMatrix> {
    let dimension = validate_actions(actions)?;
    let blocks = stabilizer_blocks(actions)?;
    Ok(solve_constraint_kernel(actions[0].context(), &blocks, dimension)?.basis_matrix)
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SolveBasisStabilizerResult {
    pub operator_basis: Vec<ExactMatrix>,
    pub action_matrices: Vec<ExactMatrix>,
    pub kernel: ConstraintKernelResult,
}

pub fn solve_basis_stabilizer<F>(
    basis: &[ExactMatrix],
    transforms: &mut [F],
) -> ExactResult<SolveBasisStabilizerResult>
where
    F: FnMut(&ExactMatrix) -> ExactResult<ExactMatrix>,
{
    let Some(first) = basis.first() else {
        return Err(constraint_error("operator basis must be nonempty"));
    };
    let action_matrices = transforms
        .iter_mut()
        .map(|transform| action_matrix_in_basis(basis, transform))
        .collect::<ExactResult<Vec<_>>>()?;
    let blocks = stabilizer_blocks(&action_matrices)?;
    let kernel = solve_constraint_kernel(first.context(), &blocks, basis.len())?;
    Ok(SolveBasisStabilizerResult {
        operator_basis: basis.to_vec(),
        action_matrices,
        kernel,
    })
}

pub fn solve_basis_stabilizer_actions(
    basis: &[ExactMatrix],
    action_matrices: &[ExactMatrix],
) -> ExactResult<SolveBasisStabilizerResult> {
    let Some(first) = basis.first() else {
        return Err(constraint_error("operator basis must be nonempty"));
    };
    if action_matrices.is_empty() {
        return Err(constraint_error("action list must be nonempty"));
    }
    let blocks = stabilizer_blocks(action_matrices)?;
    if action_matrices[0].rows() != basis.len() || action_matrices[0].context() != first.context() {
        return Err(constraint_error(
            "action matrices do not act on the supplied operator basis",
        ));
    }
    let kernel = solve_constraint_kernel(first.context(), &blocks, basis.len())?;
    Ok(SolveBasisStabilizerResult {
        operator_basis: basis.to_vec(),
        action_matrices: action_matrices.to_vec(),
        kernel,
    })
}

#[cfg(test)]
mod tests {
    use super::{
        ConstraintKernelMethod, ConstraintTarget, ConstraintValidationLevel,
        RectangularConstraintOperation, compile_rectangular_constraint_blocks, invariant_basis,
        rectangular_transport_matrix, solve_constraint_kernel, solve_constraint_kernel_with_method,
        solve_constraint_kernel_with_options,
    };
    use cyclotomic_nullspace::{Element, ExactMatrix, Rational, make_context};
    use std::sync::Arc;

    fn rational_matrix(
        context: &Arc<cyclotomic_nullspace::CyclotomicContext>,
        rows: usize,
        columns: usize,
        values: &[i64],
    ) -> ExactMatrix {
        ExactMatrix::new(
            Arc::clone(context),
            rows,
            columns,
            values
                .iter()
                .map(|value| {
                    Element::from_polynomial(Arc::clone(context), &[Rational::from_i64(*value)])
                        .expect("rational element")
                })
                .collect(),
        )
        .expect("matrix")
    }

    #[test]
    fn fixed_space_and_empty_constraints_use_exact_common_kernel() {
        let context = make_context(1, 8).expect("context");
        let action = rational_matrix(&context, 2, 2, &[1, 0, 0, -1]);
        let basis = invariant_basis(&[action]).expect("invariant basis");
        assert_eq!((basis.rows(), basis.columns()), (2, 1));
        assert!(!basis.entry(0, 0).expect("entry").is_zero());
        assert!(basis.entry(1, 0).expect("entry").is_zero());

        let unconstrained = solve_constraint_kernel(&context, &[], 3).expect("kernel");
        assert_eq!(unconstrained.nullity, 3);
        assert_eq!(unconstrained.iteration_nullities, vec![3]);
    }

    #[test]
    fn iterative_and_stacked_methods_return_the_same_exact_subspace() {
        let context = make_context(1, 8).expect("context");
        let blocks = vec![
            rational_matrix(&context, 1, 3, &[1, 1, 0]),
            rational_matrix(&context, 1, 3, &[0, 1, 1]),
        ];
        let iterative = solve_constraint_kernel_with_method(
            &context,
            &blocks,
            3,
            ConstraintKernelMethod::Iterative,
        )
        .expect("iterative kernel");
        let stacked = solve_constraint_kernel_with_method(
            &context,
            &blocks,
            3,
            ConstraintKernelMethod::Stacked,
        )
        .expect("stacked kernel");
        let cyclotomic = solve_constraint_kernel_with_method(
            &context,
            &blocks,
            3,
            ConstraintKernelMethod::Cyclotomic,
        )
        .expect("cyclotomic kernel");
        assert_eq!(iterative.basis_matrix, stacked.basis_matrix);
        assert_eq!(iterative.basis_matrix, cyclotomic.basis_matrix);
        assert_eq!(iterative.nullity, 1);
        assert_eq!(iterative.iteration_nullities, vec![3, 2, 1]);
        assert_eq!(stacked.iteration_nullities, vec![3, 1]);
        assert_eq!(cyclotomic.iteration_nullities, vec![3, 1]);
    }

    #[test]
    fn validation_levels_match_stable_certification_boundaries() {
        let context = make_context(1, 8).expect("context");
        let blocks = vec![rational_matrix(&context, 1, 2, &[1, 0])];
        let none = solve_constraint_kernel_with_options(
            &context,
            &blocks,
            2,
            ConstraintKernelMethod::Iterative,
            ConstraintValidationLevel::None,
        )
        .expect("unvalidated kernel");
        assert!(!none.residual_computed);
        assert!(!none.independent_checked);
        assert!(!none.rank_nullity_checked);

        let basic = solve_constraint_kernel_with_options(
            &context,
            &blocks,
            2,
            ConstraintKernelMethod::Iterative,
            ConstraintValidationLevel::Basic,
        )
        .expect("residual-validated kernel");
        assert!(basic.residual_computed);
        assert!(basic.exact_residual_verified);
        assert!(!basic.independent_checked);
        assert!(!basic.rank_nullity_checked);

        let full = solve_constraint_kernel_with_options(
            &context,
            &blocks,
            2,
            ConstraintKernelMethod::Iterative,
            ConstraintValidationLevel::Full,
        )
        .expect("fully validated kernel");
        assert!(full.exact_residual_verified);
        assert!(full.independent_verified);
        assert!(full.rank_nullity_verified);
        assert!(full.independent_checked);
        assert!(full.rank_nullity_checked);
        assert_eq!(none.basis_matrix, basic.basis_matrix);
        assert_eq!(basic.basis_matrix, full.basis_matrix);
    }

    #[test]
    fn rectangular_operations_match_real_coordinate_contract() {
        let context = make_context(4, 8).expect("context");
        let identity = ExactMatrix::identity(&context, 1).expect("identity");
        let unitary =
            rectangular_transport_matrix(&identity, &identity, false).expect("unitary transport");
        assert_eq!(
            unitary,
            ExactMatrix::identity(&context, 2).expect("real identity")
        );
        let antiunitary = rectangular_transport_matrix(&identity, &identity, true)
            .expect("antiunitary transport");
        assert!(!antiunitary.entry(0, 0).expect("entry").is_zero());
        assert_eq!(
            antiunitary.entry(1, 1).expect("entry"),
            &Element::one(&context).expect("1").negate().expect("-1")
        );

        let blocks = compile_rectangular_constraint_blocks(
            1,
            1,
            &[RectangularConstraintOperation {
                left: identity.clone(),
                right: identity,
                antiunitary: true,
                target: ConstraintTarget::Same,
            }],
            None,
        )
        .expect("blocks");
        let kernel = solve_constraint_kernel(&context, &blocks, 2).expect("kernel");
        assert_eq!(kernel.nullity, 1);
    }

    #[test]
    fn incompatible_reverse_shape_fails_explicitly() {
        let context = make_context(4, 8).expect("context");
        let error = compile_rectangular_constraint_blocks(
            1,
            2,
            &[RectangularConstraintOperation {
                left: ExactMatrix::identity(&context, 1).expect("left"),
                right: ExactMatrix::identity(&context, 2).expect("right"),
                antiunitary: false,
                target: ConstraintTarget::Reverse,
            }],
            None,
        )
        .expect_err("rectangular reverse must fail");
        assert_eq!(error.tag(), "InvalidConstraintData");
    }
}
