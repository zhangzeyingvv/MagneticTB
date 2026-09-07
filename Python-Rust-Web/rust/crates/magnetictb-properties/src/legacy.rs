use std::f64::consts::{PI, TAU};

use nalgebra::linalg::{SVD, Schur};
use nalgebra::{Complex, DMatrix};

use crate::spectral::occupied_frame;
use crate::{PropertiesError, PropertiesResult};

fn reverse_unitary_part(matrix: &DMatrix<Complex<f64>>) -> PropertiesResult<DMatrix<Complex<f64>>> {
    let decomposition = SVD::new(matrix.clone(), true, true);
    let left = decomposition.u.ok_or_else(|| {
        PropertiesError::new("SvdFailure", "SVD did not return left singular vectors")
    })?;
    let right_adjoint = decomposition.v_t.ok_or_else(|| {
        PropertiesError::new("SvdFailure", "SVD did not return right singular vectors")
    })?;
    Ok(right_adjoint.adjoint() * left.adjoint())
}

fn phases(matrix: DMatrix<Complex<f64>>, scale: f64) -> Vec<f64> {
    let (_, triangular) = Schur::new(matrix).unpack();
    triangular
        .diagonal()
        .iter()
        .map(|value| value.arg() / scale)
        .collect()
}

pub fn legacy_z2_path(
    matrices: &[DMatrix<Complex<f64>>],
    occupied: usize,
    hermitian_tolerance: f64,
) -> PropertiesResult<Vec<f64>> {
    if matrices.is_empty() {
        return Err(PropertiesError::new(
            "InvalidWilsonPath",
            "z2path requires at least two input path points",
        ));
    }
    let frames = matrices
        .iter()
        .map(|matrix| occupied_frame(matrix, occupied, hermitian_tolerance))
        .collect::<PropertiesResult<Vec<_>>>()?;
    let mut product = DMatrix::<Complex<f64>>::identity(occupied, occupied);
    for index in 0..frames.len() {
        let next = (index + 1) % frames.len();
        let overlap = frames[index].vectors.adjoint() * &frames[next].vectors;
        product = reverse_unitary_part(&overlap)? * product;
    }
    Ok(phases(product, -TAU)
        .into_iter()
        .map(|value| {
            let reduced = value.rem_euclid(1.0);
            if reduced.abs() <= hermitian_tolerance || (1.0 - reduced).abs() <= hermitian_tolerance
            {
                0.0
            } else {
                reduced
            }
        })
        .collect())
}

pub fn legacy_berryph(
    matrices: &[DMatrix<Complex<f64>>],
    occupied: usize,
    hermitian_tolerance: f64,
) -> PropertiesResult<Vec<f64>> {
    if matrices.len() < 2 {
        return Err(PropertiesError::new(
            "InvalidWilsonPath",
            "berryph requires at least two path points",
        ));
    }
    let mut frames = matrices
        .iter()
        .map(|matrix| occupied_frame(matrix, occupied, hermitian_tolerance))
        .collect::<PropertiesResult<Vec<_>>>()?;
    for index in 0..frames.len() - 1 {
        let overlap = frames[index].vectors.adjoint() * &frames[index + 1].vectors;
        let transport = reverse_unitary_part(&overlap)?;
        frames[index + 1].vectors = frames[index + 1].vectors.clone() * transport.transpose();
    }
    let endpoint = frames[0].vectors.adjoint() * &frames[frames.len() - 1].vectors;
    Ok(phases(endpoint, PI)
        .into_iter()
        .map(|value| {
            if value.abs() <= hermitian_tolerance {
                0.0
            } else {
                value
            }
        })
        .collect())
}

pub fn legacy_wloop(
    matrix_grid: &[Vec<DMatrix<Complex<f64>>>],
    occupied: usize,
    hermitian_tolerance: f64,
) -> PropertiesResult<Vec<Vec<f64>>> {
    if matrix_grid.is_empty() || matrix_grid.iter().any(|path| path.len() < 2) {
        return Err(PropertiesError::new(
            "InvalidWilsonPath",
            "wLoop requires a nonempty offset grid and at least two points per loop path",
        ));
    }
    matrix_grid
        .iter()
        .map(|matrices| {
            let frames = matrices
                .iter()
                .map(|matrix| occupied_frame(matrix, occupied, hermitian_tolerance))
                .collect::<PropertiesResult<Vec<_>>>()?;
            let mut product = DMatrix::<Complex<f64>>::identity(occupied, occupied);
            for pair in frames.windows(2) {
                let overlap = pair[0].vectors.adjoint() * &pair[1].vectors;
                product *= overlap;
            }
            Ok(phases(product, 1.0))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use nalgebra::{Complex, DMatrix};

    use super::{legacy_berryph, legacy_wloop, legacy_z2_path};

    fn scalar(value: f64) -> DMatrix<Complex<f64>> {
        DMatrix::from_element(1, 1, Complex::new(value, 0.0))
    }

    #[test]
    fn atomic_legacy_paths_match_stable_zero_phases() {
        let samples = vec![scalar(-1.0), scalar(-1.0), scalar(-1.0)];
        assert_eq!(
            legacy_z2_path(&samples[..2], 1, 1.0e-10).unwrap(),
            vec![0.0]
        );
        assert_eq!(legacy_berryph(&samples, 1, 1.0e-10).unwrap(), vec![0.0]);
        assert_eq!(
            legacy_wloop(&[samples], 1, 1.0e-10).unwrap(),
            vec![vec![0.0]]
        );
    }
}
