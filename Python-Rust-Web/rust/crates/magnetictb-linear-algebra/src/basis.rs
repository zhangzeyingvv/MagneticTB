use crate::matrix::{conjugate_transpose, real_and_imaginary_parts, scale, subtract};
use cyclotomic_nullspace::{Element, ExactError, ExactMatrix, ExactResult, multiply, rref};
use std::sync::Arc;

fn dimension_error(detail: impl Into<String>) -> ExactError {
    ExactError::new("DimensionMismatch", detail)
}

fn basis_error(detail: impl Into<String>) -> ExactError {
    ExactError::new("InvalidOperatorBasis", detail)
}

fn checked_entry_count(rows: usize, columns: usize) -> ExactResult<usize> {
    if rows == 0 || columns == 0 {
        return Err(dimension_error("matrix dimensions must be positive"));
    }
    rows.checked_mul(columns)
        .ok_or_else(|| dimension_error("matrix dimensions overflow"))
}

fn imaginary_unit(context: &Arc<cyclotomic_nullspace::CyclotomicContext>) -> ExactResult<Element> {
    if context.conductor % 4 != 0 {
        return Err(ExactError::new(
            "FieldRepresentationInsufficient",
            "the supplied cyclotomic context does not contain the imaginary unit",
        ));
    }
    Element::root_of_unity(context, 4, 1)
}

fn matrix_unit(
    context: &Arc<cyclotomic_nullspace::CyclotomicContext>,
    rows: usize,
    columns: usize,
    row: usize,
    column: usize,
) -> ExactResult<ExactMatrix> {
    let count = checked_entry_count(rows, columns)?;
    if row >= rows || column >= columns {
        return Err(dimension_error("matrix-unit index is out of range"));
    }
    let mut entries = vec![Element::zero(context)?; count];
    entries[row * columns + column] = Element::one(context)?;
    ExactMatrix::new(Arc::clone(context), rows, columns, entries)
}

pub fn rectangular_matrix_units(
    context: &Arc<cyclotomic_nullspace::CyclotomicContext>,
    rows: usize,
    columns: usize,
) -> ExactResult<Vec<ExactMatrix>> {
    let count = checked_entry_count(rows, columns)?;
    let mut basis = Vec::with_capacity(count);
    for row in 0..rows {
        for column in 0..columns {
            basis.push(matrix_unit(context, rows, columns, row, column)?);
        }
    }
    Ok(basis)
}

pub fn rectangular_real_basis(
    context: &Arc<cyclotomic_nullspace::CyclotomicContext>,
    rows: usize,
    columns: usize,
) -> ExactResult<Vec<ExactMatrix>> {
    let units = rectangular_matrix_units(context, rows, columns)?;
    let imaginary_unit = imaginary_unit(context)?;
    let imaginary_units = units
        .iter()
        .map(|unit| scale(unit, &imaginary_unit))
        .collect::<ExactResult<Vec<_>>>()?;
    Ok(units.into_iter().chain(imaginary_units).collect())
}

pub fn hermitian_matrix_basis(
    context: &Arc<cyclotomic_nullspace::CyclotomicContext>,
    dimension: usize,
) -> ExactResult<Vec<ExactMatrix>> {
    if dimension == 0 {
        return Err(dimension_error(
            "Hermitian matrix dimension must be positive",
        ));
    }
    let mut diagonal = Vec::with_capacity(dimension);
    let mut symmetric = Vec::with_capacity(dimension.saturating_mul(dimension - 1) / 2);
    let mut antisymmetric = Vec::with_capacity(symmetric.capacity());
    for row in 0..dimension {
        diagonal.push(matrix_unit(context, dimension, dimension, row, row)?);
        for column in row + 1..dimension {
            let upper = matrix_unit(context, dimension, dimension, row, column)?;
            let lower = matrix_unit(context, dimension, dimension, column, row)?;
            symmetric.push(crate::matrix::add(&upper, &lower)?);
            antisymmetric.push(scale(
                &subtract(&upper, &lower)?,
                &imaginary_unit(context)?,
            )?);
        }
    }
    diagonal.extend(symmetric);
    diagonal.extend(antisymmetric);
    Ok(diagonal)
}

pub fn matrix_to_real_coordinates(matrix: &ExactMatrix) -> ExactResult<Vec<Element>> {
    if matrix.rows() == 0 || matrix.columns() == 0 {
        return Err(dimension_error("matrix must be nonempty"));
    }
    let parts = matrix
        .entries()
        .iter()
        .map(real_and_imaginary_parts)
        .collect::<ExactResult<Vec<_>>>()?;
    let mut coordinates = Vec::with_capacity(parts.len() * 2);
    coordinates.extend(parts.iter().map(|(real, _)| real.clone()));
    coordinates.extend(parts.into_iter().map(|(_, imaginary)| imaginary));
    Ok(coordinates)
}

pub fn real_coordinates_to_matrix(
    context: &Arc<cyclotomic_nullspace::CyclotomicContext>,
    coordinates: &[Element],
    rows: usize,
    columns: usize,
) -> ExactResult<ExactMatrix> {
    let entry_count = checked_entry_count(rows, columns)?;
    if coordinates.len() != 2 * entry_count {
        return Err(dimension_error(
            "real coordinate count does not match matrix dimensions",
        ));
    }
    if coordinates.iter().any(|entry| entry.context() != context) {
        return Err(ExactError::new(
            "ConductorMismatch",
            "coordinate context differs",
        ));
    }
    for coordinate in coordinates {
        if coordinate.conjugate()? != *coordinate {
            return Err(ExactError::new(
                "InvalidRealCoordinate",
                "a real operator-space coordinate is not exactly real",
            ));
        }
    }
    let imaginary_coordinates = &coordinates[entry_count..];
    let imaginary_unit = if imaginary_coordinates.iter().any(|value| !value.is_zero()) {
        Some(imaginary_unit(context)?)
    } else {
        None
    };
    let entries = coordinates[..entry_count]
        .iter()
        .zip(imaginary_coordinates)
        .map(|(real, imaginary)| match &imaginary_unit {
            Some(unit) => real.add(&imaginary.multiply(unit)?),
            None => Ok(real.clone()),
        })
        .collect::<ExactResult<Vec<_>>>()?;
    ExactMatrix::new(Arc::clone(context), rows, columns, entries)
}

pub fn real_hs_inner_product(left: &ExactMatrix, right: &ExactMatrix) -> ExactResult<Element> {
    if left.rows() == 0
        || left.columns() == 0
        || left.rows() != right.rows()
        || left.columns() != right.columns()
        || left.context() != right.context()
    {
        return Err(dimension_error(
            "Hilbert--Schmidt arguments must be nonempty with matching shape and context",
        ));
    }
    let product = multiply(&conjugate_transpose(left)?, right)?;
    let mut trace = Element::zero(left.context())?;
    for index in 0..product.rows() {
        trace = trace.add(product.entry(index, index)?)?;
    }
    Ok(real_and_imaginary_parts(&trace)?.0)
}

fn validate_basis(basis: &[ExactMatrix]) -> ExactResult<(usize, usize)> {
    let Some(first) = basis.first() else {
        return Err(basis_error("operator basis must be nonempty"));
    };
    if first.rows() == 0 || first.columns() == 0 {
        return Err(basis_error("operator basis matrices must be nonempty"));
    }
    if basis.iter().any(|matrix| {
        matrix.context() != first.context()
            || matrix.rows() != first.rows()
            || matrix.columns() != first.columns()
    }) {
        return Err(basis_error(
            "operator basis matrices must share shape and context",
        ));
    }
    Ok((first.rows(), first.columns()))
}

pub fn basis_gram_matrix(basis: &[ExactMatrix]) -> ExactResult<ExactMatrix> {
    validate_basis(basis)?;
    let mut entries = Vec::with_capacity(basis.len() * basis.len());
    for left in basis {
        for right in basis {
            entries.push(real_hs_inner_product(left, right)?);
        }
    }
    ExactMatrix::new(
        Arc::clone(basis[0].context()),
        basis.len(),
        basis.len(),
        entries,
    )
}

fn solve_unique_square(matrix: &ExactMatrix, right: &[Element]) -> ExactResult<Vec<Element>> {
    if !crate::matrix::exact_square_matrix(matrix) || right.len() != matrix.rows() {
        return Err(dimension_error("exact linear system dimensions differ"));
    }
    let dimension = matrix.rows();
    let mut entries = Vec::with_capacity(dimension * (dimension + 1));
    for (row, right_entry) in right.iter().enumerate() {
        entries.extend_from_slice(&matrix.entries()[row * dimension..(row + 1) * dimension]);
        entries.push(right_entry.clone());
    }
    let augmented = ExactMatrix::new(
        Arc::clone(matrix.context()),
        dimension,
        dimension + 1,
        entries,
    )?;
    let reduced = rref(&augmented)?;
    if reduced.rank != dimension || reduced.pivot_columns != (0..dimension).collect::<Vec<_>>() {
        return Err(basis_error("operator basis Gram matrix is singular"));
    }
    (0..dimension)
        .map(|row| reduced.reduced_matrix.entry(row, dimension).cloned())
        .collect()
}

pub fn coordinates_in_basis(
    matrix: &ExactMatrix,
    basis: &[ExactMatrix],
) -> ExactResult<Vec<Element>> {
    let (rows, columns) = validate_basis(basis)?;
    if matrix.rows() != rows
        || matrix.columns() != columns
        || matrix.context() != basis[0].context()
    {
        return Err(basis_error(
            "matrix and operator basis must share shape and context",
        ));
    }
    let gram = basis_gram_matrix(basis)?;
    let right = basis
        .iter()
        .map(|element| real_hs_inner_product(element, matrix))
        .collect::<ExactResult<Vec<_>>>()?;
    let coordinates = solve_unique_square(&gram, &right)?;
    let reconstruction = linear_combination(basis, &coordinates)?;
    if !crate::matrix::exact_matrix_equal(matrix, &reconstruction) {
        return Err(basis_error("matrix is outside the supplied operator span"));
    }
    Ok(coordinates)
}

fn linear_combination(basis: &[ExactMatrix], coefficients: &[Element]) -> ExactResult<ExactMatrix> {
    validate_basis(basis)?;
    if basis.len() != coefficients.len() {
        return Err(dimension_error("basis and coefficient counts differ"));
    }
    let mut result = ExactMatrix::zero(basis[0].context(), basis[0].rows(), basis[0].columns())?;
    for (element, coefficient) in basis.iter().zip(coefficients) {
        result = crate::matrix::add(&result, &scale(element, coefficient)?)?;
    }
    Ok(result)
}

pub fn action_matrix_in_basis<F>(
    basis: &[ExactMatrix],
    mut transform: F,
) -> ExactResult<ExactMatrix>
where
    F: FnMut(&ExactMatrix) -> ExactResult<ExactMatrix>,
{
    validate_basis(basis)?;
    let images = basis
        .iter()
        .map(&mut transform)
        .collect::<ExactResult<Vec<_>>>()?;
    action_matrix_from_images(basis, &images)
}

pub fn action_matrix_from_images(
    basis: &[ExactMatrix],
    images: &[ExactMatrix],
) -> ExactResult<ExactMatrix> {
    validate_basis(basis)?;
    if images.len() != basis.len() {
        return Err(basis_error(
            "operator-basis image count differs from basis length",
        ));
    }
    let columns = images
        .iter()
        .map(|image| coordinates_in_basis(image, basis))
        .collect::<ExactResult<Vec<_>>>()?;
    let mut entries = Vec::with_capacity(basis.len() * basis.len());
    for row in 0..basis.len() {
        for column in &columns {
            entries.push(column[row].clone());
        }
    }
    ExactMatrix::new(
        Arc::clone(basis[0].context()),
        basis.len(),
        basis.len(),
        entries,
    )
}

pub fn dagger_coordinate_matrix(
    context: &Arc<cyclotomic_nullspace::CyclotomicContext>,
    rows: usize,
    columns: usize,
) -> ExactResult<ExactMatrix> {
    let entry_count = checked_entry_count(rows, columns)?;
    let mut entries = vec![Element::zero(context)?; 4 * entry_count * entry_count];
    let one = Element::one(context)?;
    let minus_one = one.negate()?;
    let dimension = 2 * entry_count;
    for row in 0..rows {
        for column in 0..columns {
            let source = row * columns + column;
            let target = column * rows + row;
            entries[target * dimension + source] = one.clone();
            entries[(entry_count + target) * dimension + entry_count + source] = minus_one.clone();
        }
    }
    ExactMatrix::new(Arc::clone(context), dimension, dimension, entries)
}

#[cfg(test)]
mod tests {
    use super::{
        dagger_coordinate_matrix, hermitian_matrix_basis, matrix_to_real_coordinates,
        real_coordinates_to_matrix, rectangular_matrix_units,
    };
    use crate::matrix::{exact_hermitian_matrix, exact_matrix_equal};
    use cyclotomic_nullspace::{Element, ExactMatrix, make_context, multiply};
    use std::sync::Arc;

    #[test]
    fn bases_have_mathematica_order_and_exact_shapes() {
        let context = make_context(4, 8).expect("context");
        let units = rectangular_matrix_units(&context, 2, 3).expect("units");
        assert_eq!(units.len(), 6);
        assert!(!units[0].entry(0, 0).expect("entry").is_zero());
        assert!(!units[1].entry(0, 1).expect("entry").is_zero());
        assert!(!units[3].entry(1, 0).expect("entry").is_zero());

        let hermitian = hermitian_matrix_basis(&context, 2).expect("Hermitian basis");
        assert_eq!(hermitian.len(), 4);
        assert!(
            hermitian
                .iter()
                .all(|matrix| exact_hermitian_matrix(matrix).expect("predicate"))
        );
    }

    #[test]
    fn real_coordinates_and_dagger_are_exact_inverses() {
        let context = make_context(4, 8).expect("context");
        let i = Element::root_of_unity(&context, 4, 1).expect("i");
        let matrix = ExactMatrix::new(
            Arc::clone(&context),
            2,
            1,
            vec![Element::one(&context).expect("1"), i],
        )
        .expect("matrix");
        let coordinates = matrix_to_real_coordinates(&matrix).expect("coordinates");
        let round_trip = real_coordinates_to_matrix(&context, &coordinates, 2, 1)
            .expect("matrix from coordinates");
        assert!(exact_matrix_equal(&matrix, &round_trip));

        let coordinate_column =
            ExactMatrix::new(Arc::clone(&context), coordinates.len(), 1, coordinates)
                .expect("coordinate column");
        let dagger_coordinates = multiply(
            &dagger_coordinate_matrix(&context, 2, 1).expect("dagger action"),
            &coordinate_column,
        )
        .expect("dagger coordinates");
        let dagger = real_coordinates_to_matrix(&context, dagger_coordinates.entries(), 1, 2)
            .expect("dagger matrix");
        assert!(exact_matrix_equal(
            &dagger,
            &crate::conjugate_transpose(&matrix).expect("conjugate transpose")
        ));
    }
}
