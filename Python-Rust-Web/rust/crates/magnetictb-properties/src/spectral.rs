use nalgebra::linalg::{SVD, SymmetricEigen};
use nalgebra::{Complex, DMatrix};

use crate::{PropertiesError, PropertiesResult};

#[derive(Clone, Debug)]
pub(crate) struct OccupiedFrame {
    pub eigenvalues: Vec<f64>,
    pub vectors: DMatrix<Complex<f64>>,
}

pub(crate) fn finite_complex(value: Complex<f64>) -> bool {
    value.re.is_finite() && value.im.is_finite()
}

pub(crate) fn validate_nonnegative(value: f64, name: &str) -> PropertiesResult<()> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(PropertiesError::new(
            "InvalidPropertiesOption",
            format!("{name} must be a finite nonnegative real number"),
        ))
    }
}

pub(crate) fn max_norm(matrix: &DMatrix<Complex<f64>>) -> f64 {
    matrix.iter().map(|value| value.norm()).fold(0.0, f64::max)
}

pub(crate) fn occupied_frame(
    matrix: &DMatrix<Complex<f64>>,
    occupied: usize,
    hermitian_tolerance: f64,
) -> PropertiesResult<OccupiedFrame> {
    if matrix.nrows() == 0
        || matrix.nrows() != matrix.ncols()
        || matrix.iter().copied().any(|value| !finite_complex(value))
    {
        return Err(PropertiesError::new(
            "InvalidHamiltonian",
            "Hamiltonian must be a finite nonempty square numeric matrix",
        ));
    }
    if occupied == 0 || occupied > matrix.nrows() {
        return Err(PropertiesError::new(
            "InvalidOccupiedBands",
            format!(
                "occupied-band count must be between 1 and {}; received {occupied}",
                matrix.nrows()
            ),
        ));
    }
    let residual = max_norm(&(matrix - matrix.adjoint()));
    if residual > hermitian_tolerance {
        return Err(PropertiesError::new(
            "NonHermitianHamiltonian",
            format!(
                "Hamiltonian Hermitian residual {residual} exceeds tolerance {hermitian_tolerance}"
            ),
        ));
    }
    let decomposition = SymmetricEigen::new(matrix.clone());
    let mut order = (0..matrix.nrows()).collect::<Vec<_>>();
    order.sort_by(|&left, &right| {
        decomposition.eigenvalues[left].total_cmp(&decomposition.eigenvalues[right])
    });
    let eigenvalues = order
        .iter()
        .map(|&index| decomposition.eigenvalues[index])
        .collect::<Vec<_>>();
    if eigenvalues.iter().any(|value| !value.is_finite()) {
        return Err(PropertiesError::new(
            "EigensystemFailure",
            "Hamiltonian eigensystem contains a non-finite eigenvalue",
        ));
    }
    let vectors = DMatrix::from_fn(matrix.nrows(), occupied, |row, column| {
        decomposition.eigenvectors[(row, order[column])]
    });
    Ok(OccupiedFrame {
        eigenvalues,
        vectors,
    })
}

pub fn hermitian_eigenvalues(
    matrix: &DMatrix<Complex<f64>>,
    hermitian_tolerance: f64,
) -> PropertiesResult<Vec<f64>> {
    validate_nonnegative(hermitian_tolerance, "HermitianTolerance")?;
    occupied_frame(matrix, matrix.nrows(), hermitian_tolerance).map(|frame| frame.eigenvalues)
}

pub(crate) fn unitary_part(
    matrix: &DMatrix<Complex<f64>>,
    overlap_tolerance: f64,
) -> PropertiesResult<(DMatrix<Complex<f64>>, f64)> {
    if matrix.nrows() == 0 || matrix.nrows() != matrix.ncols() {
        return Err(PropertiesError::new(
            "SingularOccupiedOverlap",
            "occupied-subspace overlap must be a nonempty square matrix",
        ));
    }
    let decomposition = SVD::new(matrix.clone(), true, true);
    let minimum = decomposition
        .singular_values
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    if !minimum.is_finite() || minimum <= overlap_tolerance {
        return Err(PropertiesError::new(
            "SingularOccupiedOverlap",
            format!(
                "minimum overlap singular value {minimum} does not exceed tolerance {overlap_tolerance}"
            ),
        ));
    }
    let left = decomposition.u.ok_or_else(|| {
        PropertiesError::new("SvdFailure", "SVD did not return left singular vectors")
    })?;
    let right_adjoint = decomposition.v_t.ok_or_else(|| {
        PropertiesError::new("SvdFailure", "SVD did not return right singular vectors")
    })?;
    Ok((left * right_adjoint, minimum))
}
