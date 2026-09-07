#![forbid(unsafe_code)]
#![allow(clippy::missing_errors_doc)]

mod band_trace;
mod bond_constraint;
mod bond_orbit;
mod solved_model;

pub use band_trace::{
    BandPointTrace, BandTraceInput, NumericSymmetryOperation, compute_band_trace,
};
pub use bond_constraint::{
    BondConstraintData, BondConstraintOrbit, compile_bond_constraints,
    compile_bond_constraints_with_continuous, compile_hermitian_bond_constraints,
    compile_hermitian_bond_constraints_with_continuous,
};
pub use bond_orbit::{
    DirectedBondOrbit, DirectedBondOrbitData, DirectedBondOrbitMember, compile_directed_bond_orbits,
};
pub use solved_model::{
    FourierCoefficient, HamiltonianSymmetryVerification, SolvedBondModel, SolvedBondTerm,
    fourier_coefficients, reconstruct_solved_bond_model, verify_hamiltonian_symmetry,
};

use std::sync::Arc;

use cyclotomic_nullspace::{Element, ExactError, ExactMatrix, ExactResult, multiply};
use magnetictb_linear_algebra::{
    add, conjugate, conjugate_transpose, real_coordinates_to_matrix, scale,
};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BlochTerm {
    pub bond_index: usize,
    pub orbit_index: usize,
    pub row_block: usize,
    pub column_block: usize,
    pub displacement: Vec<i64>,
    pub phase_sign: i64,
    pub matrix: ExactMatrix,
}

pub fn reconstruct_rectangular_hopping(
    basis_matrix: &ExactMatrix,
    parameters: &[Element],
    rows: usize,
    columns: usize,
) -> ExactResult<ExactMatrix> {
    let coordinate_dimension = rows
        .checked_mul(columns)
        .and_then(|value| value.checked_mul(2))
        .ok_or_else(|| ExactError::new("DimensionMismatch", "hopping dimension overflows"))?;
    if rows == 0
        || columns == 0
        || basis_matrix.rows() != coordinate_dimension
        || basis_matrix.columns() != parameters.len()
        || parameters
            .iter()
            .any(|parameter| parameter.context() != basis_matrix.context())
    {
        return Err(ExactError::new(
            "DimensionMismatch",
            "basis matrix must have 2*m*n rows and one column per parameter",
        ));
    }
    let parameter_matrix = ExactMatrix::new(
        Arc::clone(basis_matrix.context()),
        parameters.len(),
        1,
        parameters.to_vec(),
    )?;
    let coordinates = multiply(basis_matrix, &parameter_matrix)?;
    real_coordinates_to_matrix(basis_matrix.context(), coordinates.entries(), rows, columns)
}

pub fn reconstruct_from_basis(
    basis_matrix: &ExactMatrix,
    parameters: &[Element],
    operator_basis: &[ExactMatrix],
) -> ExactResult<ExactMatrix> {
    let Some(first) = operator_basis.first() else {
        return Err(ExactError::new(
            "InvalidOperatorBasis",
            "operator basis must be nonempty",
        ));
    };
    if basis_matrix.rows() != operator_basis.len()
        || basis_matrix.columns() != parameters.len()
        || basis_matrix.context() != first.context()
        || parameters
            .iter()
            .any(|parameter| parameter.context() != first.context())
        || operator_basis.iter().any(|matrix| {
            matrix.context() != first.context()
                || matrix.rows() != first.rows()
                || matrix.columns() != first.columns()
        })
    {
        return Err(ExactError::new(
            "DimensionMismatch",
            "basis coordinates, parameters, and operator basis do not align",
        ));
    }
    let parameter_matrix = ExactMatrix::new(
        Arc::clone(first.context()),
        parameters.len(),
        1,
        parameters.to_vec(),
    )?;
    let coordinates = multiply(basis_matrix, &parameter_matrix)?;
    let mut result = ExactMatrix::zero(first.context(), first.rows(), first.columns())?;
    for (coordinate, basis) in coordinates.entries().iter().zip(operator_basis) {
        result = add(&result, &scale(basis, coordinate)?)?;
    }
    Ok(result)
}

pub fn propagate_hopping(
    hopping: &ExactMatrix,
    left: &ExactMatrix,
    right: &ExactMatrix,
    antiunitary: bool,
) -> ExactResult<ExactMatrix> {
    if left.rows() == 0
        || left.rows() != left.columns()
        || right.rows() == 0
        || right.rows() != right.columns()
        || hopping.rows() != left.rows()
        || hopping.columns() != right.rows()
        || hopping.context() != left.context()
        || hopping.context() != right.context()
    {
        return Err(ExactError::new(
            "DimensionMismatch",
            "hopping and representation matrices have incompatible dimensions",
        ));
    }
    let transformed = if antiunitary {
        conjugate(hopping)?
    } else {
        hopping.clone()
    };
    multiply(&multiply(left, &transformed)?, &conjugate_transpose(right)?)
}

pub fn assemble_bloch_hamiltonian(
    block_dimensions: &[usize],
    terms: &[BlochTerm],
    phase_generators: &[Element],
) -> ExactResult<ExactMatrix> {
    if block_dimensions.is_empty() || block_dimensions.contains(&0) {
        return Err(ExactError::new(
            "DimensionMismatch",
            "block dimensions must be positive",
        ));
    }
    let context = terms
        .first()
        .map(|term| Arc::clone(term.matrix.context()))
        .or_else(|| {
            phase_generators
                .first()
                .map(|phase| Arc::clone(phase.context()))
        })
        .ok_or_else(|| {
            ExactError::new(
                "MissingExactContext",
                "at least one term or phase generator is required",
            )
        })?;
    if phase_generators
        .iter()
        .any(|phase| phase.context() != &context)
    {
        return Err(ExactError::new(
            "ConductorMismatch",
            "phase generators use different exact fields",
        ));
    }
    let mut blocks = block_dimensions
        .iter()
        .map(|&rows| {
            block_dimensions
                .iter()
                .map(|&columns| ExactMatrix::zero(&context, rows, columns))
                .collect::<ExactResult<Vec<_>>>()
        })
        .collect::<ExactResult<Vec<_>>>()?;
    for term in terms {
        if term.row_block >= block_dimensions.len()
            || term.column_block >= block_dimensions.len()
            || term.displacement.len() != phase_generators.len()
            || !matches!(term.phase_sign, -1 | 1)
            || term.matrix.context() != &context
            || term.matrix.rows() != block_dimensions[term.row_block]
            || term.matrix.columns() != block_dimensions[term.column_block]
        {
            return Err(ExactError::new(
                "MalformedBlochTerm",
                "a Bloch term has incompatible block, displacement, or matrix data",
            ));
        }
        let mut phase = Element::one(&context)?;
        for (generator, &displacement) in phase_generators.iter().zip(&term.displacement) {
            let exponent = term.phase_sign.checked_mul(displacement).ok_or_else(|| {
                ExactError::new("DimensionMismatch", "Bloch phase exponent overflows")
            })?;
            phase = phase.multiply(&generator.power(exponent)?)?;
        }
        let contribution = scale(&term.matrix, &phase)?;
        blocks[term.row_block][term.column_block] =
            add(&blocks[term.row_block][term.column_block], &contribution)?;
    }
    flatten_blocks(&context, &blocks, block_dimensions)
}

fn flatten_blocks(
    context: &Arc<cyclotomic_nullspace::CyclotomicContext>,
    blocks: &[Vec<ExactMatrix>],
    dimensions: &[usize],
) -> ExactResult<ExactMatrix> {
    let total: usize = dimensions.iter().sum();
    let mut entries = vec![Element::zero(context)?; total * total];
    let mut row_offset = 0;
    for (row, &row_dimension) in dimensions.iter().enumerate() {
        let mut column_offset = 0;
        for (column, &column_dimension) in dimensions.iter().enumerate() {
            let block = &blocks[row][column];
            for local_row in 0..row_dimension {
                for local_column in 0..column_dimension {
                    entries[(row_offset + local_row) * total + column_offset + local_column] =
                        block.entry(local_row, local_column)?.clone();
                }
            }
            column_offset += column_dimension;
        }
        row_offset += row_dimension;
    }
    ExactMatrix::new(Arc::clone(context), total, total, entries)
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use cyclotomic_nullspace::{Element, ExactMatrix, Rational, make_context};

    use super::{BlochTerm, assemble_bloch_hamiltonian, reconstruct_rectangular_hopping};

    fn integer(context: &Arc<cyclotomic_nullspace::CyclotomicContext>, value: i64) -> Element {
        Element::from_polynomial(Arc::clone(context), &[Rational::from_i64(value)])
            .expect("integer")
    }

    #[test]
    fn trivial_s_onsite_reconstructs_and_assembles_exactly() {
        let context = make_context(1, 8).expect("rational context");
        let basis = ExactMatrix::new(
            Arc::clone(&context),
            2,
            1,
            vec![integer(&context, 1), integer(&context, 0)],
        )
        .expect("onsite basis");
        let hopping = reconstruct_rectangular_hopping(&basis, &[integer(&context, 3)], 1, 1)
            .expect("onsite hopping");
        assert_eq!(hopping.entry(0, 0).expect("entry"), &integer(&context, 3));
        let term = BlochTerm {
            bond_index: 0,
            orbit_index: 0,
            row_block: 0,
            column_block: 0,
            displacement: vec![0, 0, 0],
            phase_sign: 1,
            matrix: hopping,
        };
        let phases = vec![Element::one(&context).expect("one"); 3];
        let hamiltonian = assemble_bloch_hamiltonian(&[1], &[term], &phases).expect("Hamiltonian");
        assert_eq!(
            hamiltonian.entry(0, 0).expect("entry"),
            &integer(&context, 3)
        );
    }

    #[test]
    fn malformed_bloch_term_fails_explicitly() {
        let context = make_context(1, 8).expect("context");
        let term = BlochTerm {
            bond_index: 0,
            orbit_index: 0,
            row_block: 1,
            column_block: 0,
            displacement: vec![],
            phase_sign: 1,
            matrix: ExactMatrix::identity(&context, 1).expect("matrix"),
        };
        let error =
            assemble_bloch_hamiltonian(&[1], &[term], &[]).expect_err("bad block index must fail");
        assert_eq!(error.tag(), "MalformedBlochTerm");
    }
}
