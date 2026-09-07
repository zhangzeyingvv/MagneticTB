use std::f64::consts::{PI, TAU};

use nalgebra::linalg::Schur;
use nalgebra::{Complex, DMatrix};
use serde::Serialize;

use crate::spectral::{max_norm, occupied_frame, unitary_part, validate_nonnegative};
use crate::{PropertiesError, PropertiesResult};

#[derive(Clone, Copy, Debug)]
pub struct WilsonLoopOptions {
    pub hermitian_tolerance: f64,
    pub gap_tolerance: f64,
    pub covariance_tolerance: f64,
    pub overlap_tolerance: f64,
}

impl Default for WilsonLoopOptions {
    fn default() -> Self {
        Self {
            hermitian_tolerance: 1.0e-10,
            gap_tolerance: 1.0e-9,
            covariance_tolerance: 1.0e-8,
            overlap_tolerance: 1.0e-10,
        }
    }
}

impl WilsonLoopOptions {
    pub(crate) fn validate(self) -> PropertiesResult<Self> {
        validate_nonnegative(self.hermitian_tolerance, "HermitianTolerance")?;
        validate_nonnegative(self.gap_tolerance, "GapTolerance")?;
        validate_nonnegative(self.covariance_tolerance, "CovarianceTolerance")?;
        validate_nonnegative(self.overlap_tolerance, "OverlapTolerance")?;
        Ok(self)
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct WilsonLoopData {
    pub eigenvalues: Vec<(f64, f64)>,
    pub phases_over_pi: Vec<f64>,
    pub wannier_centers: Vec<f64>,
    #[serde(skip)]
    pub wilson_matrix: DMatrix<Complex<f64>>,
    pub unitarity_residual: f64,
    pub minimum_link_singular_value: f64,
    pub minimum_direct_gap: f64,
    pub path: Vec<Vec<f64>>,
    pub path_subdivisions: usize,
    pub occupied_bands: usize,
    pub closure_vector: Vec<f64>,
    pub reciprocal_coordinates: Vec<i64>,
    pub endpoint_covariance_residual: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct BerryPhaseData {
    pub phase: f64,
    pub phase_over_pi: f64,
    pub wilson_determinant: (f64, f64),
    #[serde(skip)]
    pub wilson_matrix: DMatrix<Complex<f64>>,
    pub unitarity_residual: f64,
    pub minimum_link_singular_value: f64,
    pub minimum_direct_gap: f64,
    pub path: Vec<Vec<f64>>,
    pub occupied_bands: usize,
    pub closure_vector: Vec<f64>,
    pub reciprocal_coordinates: Vec<i64>,
    pub endpoint_covariance_residual: f64,
}

fn finite_real_vector(vector: &[f64]) -> bool {
    !vector.is_empty() && vector.iter().all(|value| value.is_finite())
}

fn reciprocal_closure(path: &[Vec<f64>], tolerance: f64) -> PropertiesResult<(Vec<f64>, Vec<i64>)> {
    let Some(first) = path.first() else {
        return Err(PropertiesError::new(
            "InvalidWilsonPath",
            "path must contain at least two momentum points",
        ));
    };
    if path.len() < 2
        || !finite_real_vector(first)
        || path
            .iter()
            .any(|point| point.len() != first.len() || !finite_real_vector(point))
    {
        return Err(PropertiesError::new(
            "InvalidWilsonPath",
            "path must contain at least two aligned finite real momentum vectors",
        ));
    }
    let last = path.last().expect("nonempty path");
    let closure = last
        .iter()
        .zip(first)
        .map(|(end, start)| end - start)
        .collect::<Vec<_>>();
    let reciprocal = closure
        .iter()
        .map(|value| (value / TAU).round())
        .collect::<Vec<_>>();
    if closure.iter().zip(&reciprocal).any(|(value, rounded)| {
        (value / TAU - rounded).abs() > tolerance
            || *rounded < i64::MIN as f64
            || *rounded > i64::MAX as f64
    }) {
        return Err(PropertiesError::new(
            "InvalidWilsonPath",
            "path endpoint must differ from its start by 2 Pi times integers",
        ));
    }
    Ok((
        closure,
        reciprocal.into_iter().map(|value| value as i64).collect(),
    ))
}

fn endpoint_sewing(centers: &[Vec<f64>], closure: &[f64]) -> DMatrix<Complex<f64>> {
    DMatrix::from_diagonal(&nalgebra::DVector::from_iterator(
        centers.len(),
        centers.iter().map(|center| {
            Complex::from_polar(
                1.0,
                closure
                    .iter()
                    .zip(center)
                    .map(|(left, right)| left * right)
                    .sum(),
            )
        }),
    ))
}

pub(crate) struct WilsonComputation {
    pub wilson_matrix: DMatrix<Complex<f64>>,
    pub minimum_link_singular_value: f64,
    pub minimum_direct_gap: f64,
    pub unitarity_residual: f64,
    pub endpoint_covariance_residual: f64,
    pub closure_vector: Vec<f64>,
    pub reciprocal_coordinates: Vec<i64>,
}

pub(crate) fn compute_wilson(
    matrices: &[DMatrix<Complex<f64>>],
    path: &[Vec<f64>],
    centers: &[Vec<f64>],
    occupied: usize,
    options: WilsonLoopOptions,
) -> PropertiesResult<WilsonComputation> {
    let options = options.validate()?;
    let (closure, reciprocal) = reciprocal_closure(path, options.covariance_tolerance)?;
    if matrices.len() != path.len() || matrices.is_empty() {
        return Err(PropertiesError::new(
            "InvalidHamiltonianSamples",
            "one Hamiltonian matrix is required for every path point",
        ));
    }
    let dimension = matrices[0].nrows();
    if dimension == 0
        || matrices
            .iter()
            .any(|matrix| matrix.nrows() != dimension || matrix.ncols() != dimension)
    {
        return Err(PropertiesError::new(
            "InvalidHamiltonian",
            "Hamiltonian dimension must remain finite, square, nonzero, and constant",
        ));
    }
    if centers.len() != dimension
        || centers.iter().any(|center| {
            center.len() != path[0].len() || !center.iter().all(|value| value.is_finite())
        })
    {
        return Err(PropertiesError::new(
            "InvalidWannierCenters",
            "one finite coordinate vector per Hamiltonian row is required",
        ));
    }
    let frames = matrices
        .iter()
        .map(|matrix| occupied_frame(matrix, occupied, options.hermitian_tolerance))
        .collect::<PropertiesResult<Vec<_>>>()?;
    let minimum_gap = if occupied < dimension {
        frames
            .iter()
            .map(|frame| frame.eigenvalues[occupied] - frame.eigenvalues[occupied - 1])
            .fold(f64::INFINITY, f64::min)
    } else {
        f64::INFINITY
    };
    if minimum_gap <= options.gap_tolerance {
        return Err(PropertiesError::new(
            "OccupiedSubspaceGapClosed",
            format!(
                "minimum direct gap {minimum_gap} does not exceed tolerance {}",
                options.gap_tolerance
            ),
        ));
    }
    let sewing = endpoint_sewing(centers, &closure);
    let covariance_residual = max_norm(
        &(matrices.last().expect("nonempty matrices") - sewing.adjoint() * &matrices[0] * &sewing),
    );
    if covariance_residual > options.covariance_tolerance {
        return Err(PropertiesError::new(
            "EndpointCovarianceFailure",
            format!(
                "endpoint covariance residual {covariance_residual} exceeds tolerance {}",
                options.covariance_tolerance
            ),
        ));
    }
    let mut links = Vec::with_capacity(frames.len());
    let mut minimum_singular = f64::INFINITY;
    for pair in frames.windows(2) {
        let overlap = pair[1].vectors.adjoint() * &pair[0].vectors;
        let (unitary, singular) = unitary_part(&overlap, options.overlap_tolerance)?;
        minimum_singular = minimum_singular.min(singular);
        links.push(unitary);
    }
    let closure_overlap =
        frames[0].vectors.adjoint() * sewing * &frames.last().expect("nonempty frames").vectors;
    let (closure_unitary, closure_singular) =
        unitary_part(&closure_overlap, options.overlap_tolerance)?;
    minimum_singular = minimum_singular.min(closure_singular);
    let mut wilson = closure_unitary;
    for link in links.iter().rev() {
        wilson *= link;
    }
    for value in wilson.iter_mut() {
        if value.norm() <= options.hermitian_tolerance {
            *value = Complex::new(0.0, 0.0);
        }
    }
    let identity = DMatrix::<Complex<f64>>::identity(occupied, occupied);
    let unitarity_residual = max_norm(&(wilson.adjoint() * &wilson - identity));
    Ok(WilsonComputation {
        wilson_matrix: wilson,
        minimum_link_singular_value: minimum_singular,
        minimum_direct_gap: minimum_gap,
        unitarity_residual,
        endpoint_covariance_residual: covariance_residual,
        closure_vector: closure,
        reciprocal_coordinates: reciprocal,
    })
}

pub fn wilson_loop(
    matrices: &[DMatrix<Complex<f64>>],
    path: &[Vec<f64>],
    centers: &[Vec<f64>],
    occupied: usize,
    options: WilsonLoopOptions,
) -> PropertiesResult<WilsonLoopData> {
    let computation = compute_wilson(matrices, path, centers, occupied, options)?;
    let (_, triangular) = Schur::new(computation.wilson_matrix.clone()).unpack();
    let mut eigenvalues = triangular.diagonal().iter().copied().collect::<Vec<_>>();
    eigenvalues.sort_by(|left, right| (left.arg() / PI).total_cmp(&(right.arg() / PI)));
    let phases = eigenvalues
        .iter()
        .map(|value| {
            let phase = value.arg() / PI;
            if phase.abs() <= options.hermitian_tolerance {
                0.0
            } else {
                phase
            }
        })
        .collect::<Vec<_>>();
    let centers_modulo_one = phases
        .iter()
        .map(|phase| {
            let value = phase / 2.0;
            let chopped = if value.abs() <= options.hermitian_tolerance {
                0.0
            } else {
                value
            };
            chopped.rem_euclid(1.0)
        })
        .collect();
    Ok(WilsonLoopData {
        eigenvalues: eigenvalues
            .iter()
            .map(|value| (value.re, value.im))
            .collect(),
        phases_over_pi: phases,
        wannier_centers: centers_modulo_one,
        wilson_matrix: computation.wilson_matrix,
        unitarity_residual: computation.unitarity_residual,
        minimum_link_singular_value: computation.minimum_link_singular_value,
        minimum_direct_gap: computation.minimum_direct_gap,
        path: path.to_vec(),
        path_subdivisions: path.len() - 1,
        occupied_bands: occupied,
        closure_vector: computation.closure_vector,
        reciprocal_coordinates: computation.reciprocal_coordinates,
        endpoint_covariance_residual: computation.endpoint_covariance_residual,
    })
}

pub fn berry_phase_from_samples(
    matrices: &[DMatrix<Complex<f64>>],
    path: &[Vec<f64>],
    centers: &[Vec<f64>],
    occupied: usize,
    options: WilsonLoopOptions,
) -> PropertiesResult<BerryPhaseData> {
    let computation = compute_wilson(matrices, path, centers, occupied, options)?;
    let determinant = computation.wilson_matrix.determinant();
    let mut phase = determinant.arg();
    if phase.abs() <= options.hermitian_tolerance {
        phase = 0.0;
    }
    Ok(BerryPhaseData {
        phase,
        phase_over_pi: phase / PI,
        wilson_determinant: (determinant.re, determinant.im),
        wilson_matrix: computation.wilson_matrix,
        unitarity_residual: computation.unitarity_residual,
        minimum_link_singular_value: computation.minimum_link_singular_value,
        minimum_direct_gap: computation.minimum_direct_gap,
        path: path.to_vec(),
        occupied_bands: occupied,
        closure_vector: computation.closure_vector,
        reciprocal_coordinates: computation.reciprocal_coordinates,
        endpoint_covariance_residual: computation.endpoint_covariance_residual,
    })
}
