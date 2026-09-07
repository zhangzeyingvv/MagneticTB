//! Direct Rust port of `ExactMatrix.wl`, `ExactNullSpace.wl`, and `CommonKernel.wl`.

use crate::error::{ExactResult, fail};
use crate::exact::{CyclotomicContext, Element};
use crate::native_backend::{
    coefficient_field_from_cyclotomic, coefficient_raw_common_kernel, coefficient_raw_from_exact,
    coefficient_raw_multiply, coefficient_raw_null_space, coefficient_raw_rref,
    coefficient_raw_to_exact, coefficient_work_plan, coefficient_work_plan_lift_basis,
    rational_raw_common_kernel, rational_raw_from_exact, rational_raw_multiply,
    rational_raw_null_space, rational_raw_rref, rational_raw_to_exact,
};
use std::sync::Arc;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ExactMatrix {
    context: Arc<CyclotomicContext>,
    rows: usize,
    columns: usize,
    entries: Vec<Element>,
}

impl ExactMatrix {
    /// Port of `cqMatrixCreate`.
    pub fn new(
        context: Arc<CyclotomicContext>,
        rows: usize,
        columns: usize,
        entries: Vec<Element>,
    ) -> ExactResult<Self> {
        let expected = rows.checked_mul(columns).ok_or_else(|| {
            crate::ExactError::new("MalformedSerialization", "matrix shape overflows")
        })?;
        if entries.len() != expected {
            return fail(
                "MalformedSerialization",
                "matrix entry count does not match shape",
            );
        }
        if entries.iter().any(|entry| entry.context() != &context) {
            return fail("ConductorMismatch", "matrix entry has a different context");
        }
        Ok(Self {
            context,
            rows,
            columns,
            entries,
        })
    }

    pub fn zero(
        context: &Arc<CyclotomicContext>,
        rows: usize,
        columns: usize,
    ) -> ExactResult<Self> {
        let count = rows.checked_mul(columns).ok_or_else(|| {
            crate::ExactError::new("MalformedSerialization", "matrix shape overflows")
        })?;
        Self::new(
            Arc::clone(context),
            rows,
            columns,
            vec![Element::zero(context)?; count],
        )
    }

    pub fn identity(context: &Arc<CyclotomicContext>, dimension: usize) -> ExactResult<Self> {
        let count = dimension.checked_mul(dimension).ok_or_else(|| {
            crate::ExactError::new("MalformedSerialization", "matrix shape overflows")
        })?;
        let mut entries = vec![Element::zero(context)?; count];
        for index in 0..dimension {
            entries[index * dimension + index] = Element::one(context)?;
        }
        Self::new(Arc::clone(context), dimension, dimension, entries)
    }

    #[must_use]
    pub fn context(&self) -> &Arc<CyclotomicContext> {
        &self.context
    }

    #[must_use]
    pub fn rows(&self) -> usize {
        self.rows
    }

    #[must_use]
    pub fn columns(&self) -> usize {
        self.columns
    }

    #[must_use]
    pub fn entries(&self) -> &[Element] {
        &self.entries
    }

    fn offset(&self, row: usize, column: usize) -> ExactResult<usize> {
        if row >= self.rows || column >= self.columns {
            return fail("DimensionMismatch", "matrix index is out of range");
        }
        Ok(row * self.columns + column)
    }

    pub fn entry(&self, row: usize, column: usize) -> ExactResult<&Element> {
        Ok(&self.entries[self.offset(row, column)?])
    }

    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.entries.iter().all(Element::is_zero)
    }

    /// Port of `cqMatrixTranspose`, preserving explicit zero-dimensional shapes.
    pub fn transpose(&self) -> ExactResult<Self> {
        let mut entries = Vec::with_capacity(self.entries.len());
        for column in 0..self.columns {
            for row in 0..self.rows {
                entries.push(self.entry(row, column)?.clone());
            }
        }
        Self::new(Arc::clone(&self.context), self.columns, self.rows, entries)
    }
}

/// Port of `cqMatrixMultiply`.
pub fn multiply(left: &ExactMatrix, right: &ExactMatrix) -> ExactResult<ExactMatrix> {
    if left.context != right.context {
        return fail("ConductorMismatch", "matrix multiplication contexts differ");
    }
    if left.columns != right.rows {
        return fail(
            "DimensionMismatch",
            "matrix multiplication dimensions do not match",
        );
    }
    if left.context.conductor == 1 {
        return rational_raw_to_exact(
            &left.context,
            &rational_raw_multiply(
                &rational_raw_from_exact(left)?,
                &rational_raw_from_exact(right)?,
            )?,
        );
    }
    let field = coefficient_field_from_cyclotomic(&left.context);
    coefficient_raw_to_exact(
        &left.context,
        &coefficient_raw_multiply(
            field,
            &coefficient_raw_from_exact(left)?,
            &coefficient_raw_from_exact(right)?,
        )?,
    )
}

/// Port of `cqMatrixVerticalStack`.
pub fn vertical_stack(
    matrices: &[ExactMatrix],
    context: &Arc<CyclotomicContext>,
    columns: usize,
) -> ExactResult<ExactMatrix> {
    let mut rows = 0usize;
    let mut entries = Vec::new();
    for matrix in matrices {
        if matrix.context != *context {
            return fail("ConductorMismatch", "stacked matrix context differs");
        }
        if matrix.columns != columns {
            return fail("DimensionMismatch", "stacked matrix column count differs");
        }
        rows = rows.checked_add(matrix.rows).ok_or_else(|| {
            crate::ExactError::new("MalformedSerialization", "stacked row count overflows")
        })?;
        entries.extend(matrix.entries.iter().cloned());
    }
    ExactMatrix::new(Arc::clone(context), rows, columns, entries)
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RrefResult {
    pub reduced_matrix: ExactMatrix,
    pub pivot_columns: Vec<usize>,
    pub free_columns: Vec<usize>,
    pub rank: usize,
    pub nullity: usize,
}

/// Direct port of `CyclotomicRREF`.
pub fn rref(matrix: &ExactMatrix) -> ExactResult<RrefResult> {
    if matrix.context.conductor == 1 {
        let raw = rational_raw_rref(&rational_raw_from_exact(matrix)?)?;
        return Ok(RrefResult {
            reduced_matrix: rational_raw_to_exact(&matrix.context, &raw.reduced)?,
            pivot_columns: raw.pivot_columns,
            free_columns: raw.free_columns,
            rank: raw.rank,
            nullity: raw.nullity,
        });
    }
    let raw_matrix = coefficient_raw_from_exact(matrix)?;
    let raw = coefficient_raw_rref(
        coefficient_field_from_cyclotomic(&matrix.context),
        &raw_matrix,
    )?;
    Ok(RrefResult {
        reduced_matrix: coefficient_raw_to_exact(&matrix.context, &raw.reduced)?,
        pivot_columns: raw.pivot_columns,
        free_columns: raw.free_columns,
        rank: raw.rank,
        nullity: raw.nullity,
    })
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct KernelResult {
    pub reduced_matrix: ExactMatrix,
    pub nullspace_rows: ExactMatrix,
    pub basis_matrix: ExactMatrix,
    pub pivot_columns: Vec<usize>,
    pub free_columns: Vec<usize>,
    pub rank: usize,
    pub nullity: usize,
    pub exact_residual_verified: bool,
}

/// Direct port of `CyclotomicNullSpace`.
pub fn null_space(matrix: &ExactMatrix) -> ExactResult<KernelResult> {
    if matrix.context.conductor == 1 {
        let raw = rational_raw_null_space(&rational_raw_from_exact(matrix)?)?;
        return Ok(KernelResult {
            reduced_matrix: rational_raw_to_exact(&matrix.context, &raw.rref.reduced)?,
            nullspace_rows: rational_raw_to_exact(&matrix.context, &raw.nullspace_rows)?,
            basis_matrix: rational_raw_to_exact(&matrix.context, &raw.basis_matrix)?,
            pivot_columns: raw.rref.pivot_columns,
            free_columns: raw.rref.free_columns,
            rank: raw.rref.rank,
            nullity: raw.rref.nullity,
            exact_residual_verified: true,
        });
    }
    let raw_matrix = coefficient_raw_from_exact(matrix)?;
    let raw = coefficient_raw_null_space(
        coefficient_field_from_cyclotomic(&matrix.context),
        &raw_matrix,
    )?;
    Ok(KernelResult {
        reduced_matrix: coefficient_raw_to_exact(&matrix.context, &raw.rref.reduced)?,
        nullspace_rows: coefficient_raw_to_exact(&matrix.context, &raw.nullspace_rows)?,
        basis_matrix: coefficient_raw_to_exact(&matrix.context, &raw.basis_matrix)?,
        pivot_columns: raw.rref.pivot_columns,
        free_columns: raw.rref.free_columns,
        rank: raw.rref.rank,
        nullity: raw.rref.nullity,
        exact_residual_verified: true,
    })
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CompiledProblem {
    pub context: Arc<CyclotomicContext>,
    pub coordinate_dimension: usize,
    pub constraints: Vec<ExactMatrix>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CommonKernelResult {
    pub basis_matrix: ExactMatrix,
    pub nullspace_rows: ExactMatrix,
    pub iteration_nullities: Vec<usize>,
    pub rank: usize,
    pub nullity: usize,
    pub exact_residual_verified: bool,
}

#[must_use]
pub fn verify_common_kernel(constraints: &[ExactMatrix], basis_matrix: &ExactMatrix) -> bool {
    constraints.iter().all(|constraint| {
        multiply(constraint, basis_matrix).is_ok_and(|residual| residual.is_zero())
    })
}

fn common_kernel_finalize(
    problem: &CompiledProblem,
    basis_matrix: ExactMatrix,
    iteration_nullities: Vec<usize>,
    rank: usize,
    nullity: usize,
    exact_residual_verified: bool,
) -> ExactResult<CommonKernelResult> {
    let monotone = iteration_nullities
        .windows(2)
        .all(|pair| pair[0] >= pair[1]);
    if !exact_residual_verified
        || basis_matrix.rows != problem.coordinate_dimension
        || basis_matrix.columns != nullity
        || rank + nullity != problem.coordinate_dimension
        || iteration_nullities.len() != problem.constraints.len() + 1
        || iteration_nullities.first() != Some(&problem.coordinate_dimension)
        || iteration_nullities.last() != Some(&nullity)
        || !monotone
    {
        return fail(
            "MalformedSerialization",
            "raw common-kernel result metadata differs",
        );
    }
    let nullspace_rows = basis_matrix.transpose()?;
    Ok(CommonKernelResult {
        basis_matrix,
        nullspace_rows,
        iteration_nullities,
        rank,
        nullity,
        exact_residual_verified: true,
    })
}

/// Direct port of `CyclotomicCommonKernel`.
pub fn common_kernel(problem: &CompiledProblem) -> ExactResult<CommonKernelResult> {
    for constraint in &problem.constraints {
        if constraint.context != problem.context {
            return fail("ConductorMismatch", "constraint context differs");
        }
        if constraint.columns != problem.coordinate_dimension {
            return fail(
                "DimensionMismatch",
                "constraint column count differs from coordinate dimension",
            );
        }
    }
    if problem.context.conductor == 1 {
        let constraints = problem
            .constraints
            .iter()
            .map(rational_raw_from_exact)
            .collect::<ExactResult<Vec<_>>>()?;
        let raw = rational_raw_common_kernel(problem.coordinate_dimension, &constraints)?;
        let basis = rational_raw_to_exact(&problem.context, &raw.basis_matrix)?;
        return common_kernel_finalize(
            problem,
            basis,
            raw.iteration_nullities,
            raw.rank,
            raw.nullity,
            raw.exact_residual_verified,
        );
    }
    let full_constraints = problem
        .constraints
        .iter()
        .map(coefficient_raw_from_exact)
        .collect::<ExactResult<Vec<_>>>()?;
    let plan = coefficient_work_plan(
        &problem.context,
        problem.coordinate_dimension,
        &full_constraints,
    )?;
    let raw = coefficient_raw_common_kernel(
        plan.coordinate_dimension,
        &plan.raw_constraints,
        Arc::clone(&plan.field_context),
    )?;
    let lifted = coefficient_work_plan_lift_basis(&plan, &raw.basis_matrix)?;
    let basis = coefficient_raw_to_exact(&problem.context, &lifted)?;
    common_kernel_finalize(
        problem,
        basis,
        raw.iteration_nullities,
        raw.rank,
        raw.nullity,
        raw.exact_residual_verified,
    )
}
