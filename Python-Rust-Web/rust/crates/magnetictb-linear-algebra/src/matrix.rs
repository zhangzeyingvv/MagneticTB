use cyclotomic_nullspace::{Element, ExactError, ExactMatrix, ExactResult, Rational, multiply};
use std::sync::Arc;

fn dimension_error(detail: impl Into<String>) -> ExactError {
    ExactError::new("DimensionMismatch", detail)
}

fn representation_error(detail: impl Into<String>) -> ExactError {
    ExactError::new("FieldRepresentationInsufficient", detail)
}

fn require_same_context(left: &ExactMatrix, right: &ExactMatrix) -> ExactResult<()> {
    if left.context() != right.context() {
        return Err(ExactError::new(
            "ConductorMismatch",
            "matrix contexts differ",
        ));
    }
    Ok(())
}

#[must_use]
pub fn exact_square_matrix(matrix: &ExactMatrix) -> bool {
    matrix.rows() > 0 && matrix.rows() == matrix.columns()
}

#[must_use]
pub fn exact_matrix_equal(left: &ExactMatrix, right: &ExactMatrix) -> bool {
    left.context() == right.context()
        && left.rows() == right.rows()
        && left.columns() == right.columns()
        && left.entries() == right.entries()
}

pub fn integer_vector(vector: &[Element]) -> ExactResult<bool> {
    for entry in vector {
        let conjugate = entry.conjugate()?;
        if &conjugate != entry {
            return Ok(false);
        }
        let coefficients = entry.coefficients();
        if coefficients.first().is_none_or(|value| !value.is_integer())
            || coefficients.iter().skip(1).any(|value| !value.is_zero())
        {
            return Ok(false);
        }
    }
    Ok(true)
}

pub fn transpose(matrix: &ExactMatrix) -> ExactResult<ExactMatrix> {
    matrix.transpose()
}

pub fn conjugate(matrix: &ExactMatrix) -> ExactResult<ExactMatrix> {
    let entries = matrix
        .entries()
        .iter()
        .map(Element::conjugate)
        .collect::<ExactResult<Vec<_>>>()?;
    ExactMatrix::new(
        Arc::clone(matrix.context()),
        matrix.rows(),
        matrix.columns(),
        entries,
    )
}

pub fn conjugate_transpose(matrix: &ExactMatrix) -> ExactResult<ExactMatrix> {
    conjugate(matrix)?.transpose()
}

pub fn add(left: &ExactMatrix, right: &ExactMatrix) -> ExactResult<ExactMatrix> {
    require_same_context(left, right)?;
    if left.rows() != right.rows() || left.columns() != right.columns() {
        return Err(dimension_error("matrix addition dimensions differ"));
    }
    let entries = left
        .entries()
        .iter()
        .zip(right.entries())
        .map(|(left, right)| left.add(right))
        .collect::<ExactResult<Vec<_>>>()?;
    ExactMatrix::new(
        Arc::clone(left.context()),
        left.rows(),
        left.columns(),
        entries,
    )
}

pub fn subtract(left: &ExactMatrix, right: &ExactMatrix) -> ExactResult<ExactMatrix> {
    require_same_context(left, right)?;
    if left.rows() != right.rows() || left.columns() != right.columns() {
        return Err(dimension_error("matrix subtraction dimensions differ"));
    }
    let entries = left
        .entries()
        .iter()
        .zip(right.entries())
        .map(|(left, right)| left.subtract(right))
        .collect::<ExactResult<Vec<_>>>()?;
    ExactMatrix::new(
        Arc::clone(left.context()),
        left.rows(),
        left.columns(),
        entries,
    )
}

pub fn scale(matrix: &ExactMatrix, scalar: &Element) -> ExactResult<ExactMatrix> {
    if matrix.context() != scalar.context() {
        return Err(ExactError::new(
            "ConductorMismatch",
            "matrix and scalar contexts differ",
        ));
    }
    let entries = matrix
        .entries()
        .iter()
        .map(|entry| entry.multiply(scalar))
        .collect::<ExactResult<Vec<_>>>()?;
    ExactMatrix::new(
        Arc::clone(matrix.context()),
        matrix.rows(),
        matrix.columns(),
        entries,
    )
}

pub fn kronecker(left: &ExactMatrix, right: &ExactMatrix) -> ExactResult<ExactMatrix> {
    require_same_context(left, right)?;
    let rows = left
        .rows()
        .checked_mul(right.rows())
        .ok_or_else(|| dimension_error("Kronecker row count overflows"))?;
    let columns = left
        .columns()
        .checked_mul(right.columns())
        .ok_or_else(|| dimension_error("Kronecker column count overflows"))?;
    let mut entries = Vec::with_capacity(
        rows.checked_mul(columns)
            .ok_or_else(|| dimension_error("Kronecker matrix size overflows"))?,
    );
    for left_row in 0..left.rows() {
        for right_row in 0..right.rows() {
            for left_column in 0..left.columns() {
                for right_column in 0..right.columns() {
                    entries.push(
                        left.entry(left_row, left_column)?
                            .multiply(right.entry(right_row, right_column)?)?,
                    );
                }
            }
        }
    }
    ExactMatrix::new(Arc::clone(left.context()), rows, columns, entries)
}

pub fn block_matrix_2x2(
    top_left: &ExactMatrix,
    top_right: &ExactMatrix,
    bottom_left: &ExactMatrix,
    bottom_right: &ExactMatrix,
) -> ExactResult<ExactMatrix> {
    require_same_context(top_left, top_right)?;
    require_same_context(top_left, bottom_left)?;
    require_same_context(top_left, bottom_right)?;
    if top_left.rows() != top_right.rows()
        || bottom_left.rows() != bottom_right.rows()
        || top_left.columns() != bottom_left.columns()
        || top_right.columns() != bottom_right.columns()
    {
        return Err(dimension_error("2 by 2 block dimensions do not align"));
    }
    let rows = top_left.rows() + bottom_left.rows();
    let columns = top_left.columns() + top_right.columns();
    let mut entries = Vec::with_capacity(rows * columns);
    for row in 0..top_left.rows() {
        entries.extend_from_slice(
            &top_left.entries()[row * top_left.columns()..(row + 1) * top_left.columns()],
        );
        entries.extend_from_slice(
            &top_right.entries()[row * top_right.columns()..(row + 1) * top_right.columns()],
        );
    }
    for row in 0..bottom_left.rows() {
        entries.extend_from_slice(
            &bottom_left.entries()[row * bottom_left.columns()..(row + 1) * bottom_left.columns()],
        );
        entries.extend_from_slice(
            &bottom_right.entries()
                [row * bottom_right.columns()..(row + 1) * bottom_right.columns()],
        );
    }
    ExactMatrix::new(Arc::clone(top_left.context()), rows, columns, entries)
}

pub fn real_and_imaginary_parts(value: &Element) -> ExactResult<(Element, Element)> {
    let conjugate = value.conjugate()?;
    let half = Rational::parse("1", "2")?;
    let real = value.add(&conjugate)?.scale(&half)?;
    let difference = value.subtract(&conjugate)?;
    if difference.is_zero() {
        return Ok((real, Element::zero(value.context())?));
    }
    if value.context().conductor % 4 != 0 {
        return Err(representation_error(
            "the supplied cyclotomic context cannot represent exact real and imaginary parts",
        ));
    }
    let imaginary_unit = Element::root_of_unity(value.context(), 4, 1)?;
    let twice_imaginary_unit = imaginary_unit.scale(&Rational::from_i64(2))?;
    let imaginary = difference.divide(&twice_imaginary_unit)?;
    Ok((real, imaginary))
}

pub fn exact_hermitian_matrix(matrix: &ExactMatrix) -> ExactResult<bool> {
    if !exact_square_matrix(matrix) {
        return Ok(false);
    }
    Ok(exact_matrix_equal(matrix, &conjugate_transpose(matrix)?))
}

pub fn exact_unitary_matrix(matrix: &ExactMatrix) -> ExactResult<bool> {
    if !exact_square_matrix(matrix) {
        return Ok(false);
    }
    let product = multiply(&conjugate_transpose(matrix)?, matrix)?;
    Ok(exact_matrix_equal(
        &product,
        &ExactMatrix::identity(matrix.context(), matrix.rows())?,
    ))
}

#[cfg(test)]
mod tests {
    use super::{exact_hermitian_matrix, exact_unitary_matrix, integer_vector};
    use cyclotomic_nullspace::{Element, ExactMatrix, Rational, make_context};
    use std::sync::Arc;

    #[test]
    fn exact_predicates_do_not_use_tolerances() {
        let context = make_context(4, 8).expect("context");
        let i = Element::root_of_unity(&context, 4, 1).expect("i");
        let minus_i = i.negate().expect("-i");
        let matrix = ExactMatrix::new(
            Arc::clone(&context),
            2,
            2,
            vec![
                Element::zero(&context).expect("0"),
                i,
                minus_i,
                Element::zero(&context).expect("0"),
            ],
        )
        .expect("matrix");
        assert!(exact_hermitian_matrix(&matrix).expect("Hermitian predicate"));
        assert!(exact_unitary_matrix(&matrix).expect("unitary predicate"));

        let integer = Element::from_polynomial(Arc::clone(&context), &[Rational::from_i64(-3)])
            .expect("integer");
        assert!(integer_vector(&[integer]).expect("integer predicate"));
        assert!(
            !integer_vector(&[Element::root_of_unity(&context, 4, 1).expect("i")])
                .expect("integer predicate")
        );
    }
}
