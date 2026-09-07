use std::f64::consts::TAU;

use cyclotomic_nullspace::{ExactError, ExactResult};
use nalgebra::linalg::SymmetricEigen;
use nalgebra::{Complex, DMatrix, Matrix3, Vector3};

const LITTLE_GROUP_TOLERANCE: f64 = 1.0e-10;
const HERMITIAN_TOLERANCE: f64 = 1.0e-9;
const STABLE_DEGENERACY_STEP: f64 = 1.0e-4;

#[derive(Clone, Debug)]
pub struct NumericSymmetryOperation {
    pub rotation: Matrix3<f64>,
    pub translation: Vector3<f64>,
    pub antiunitary: bool,
    pub representation: DMatrix<Complex<f64>>,
    pub trace_sign: f64,
}

#[derive(Clone, Debug)]
pub struct BandTraceInput<'a> {
    pub hamiltonians: &'a [DMatrix<Complex<f64>>],
    pub reduced_kpoints: &'a [Vector3<f64>],
    pub orbital_positions: &'a [Vector3<f64>],
    pub operations: &'a [NumericSymmetryOperation],
}

#[derive(Clone, Debug, PartialEq)]
pub struct BandPointTrace {
    pub energies: Vec<f64>,
    pub degeneracies: Vec<usize>,
    pub little_symmetry_indices: Vec<usize>,
    pub traces: Vec<Vec<Complex<f64>>>,
}

fn band_error(tag: &str, detail: impl Into<String>) -> ExactError {
    ExactError::new(tag, detail)
}

fn finite_complex(value: Complex<f64>) -> bool {
    value.re.is_finite() && value.im.is_finite()
}

fn validate_input(input: &BandTraceInput<'_>) -> ExactResult<usize> {
    if input.hamiltonians.is_empty()
        || input.hamiltonians.len() != input.reduced_kpoints.len()
        || input.orbital_positions.is_empty()
        || input.operations.is_empty()
    {
        return Err(band_error(
            "InvalidBandCorepInput",
            "Hamiltonians, k points, orbitals, and symmetry operations must be nonempty and aligned",
        ));
    }
    let dimension = input.orbital_positions.len();
    if input.hamiltonians.iter().any(|matrix| {
        matrix.nrows() != dimension
            || matrix.ncols() != dimension
            || matrix.iter().any(|value| !finite_complex(*value))
    }) {
        return Err(band_error(
            "InvalidBandCorepHamiltonian",
            "every numerical Hamiltonian must be finite, square, and aligned with the orbital table",
        ));
    }
    if input
        .reduced_kpoints
        .iter()
        .any(|point| !point.iter().all(|value| value.is_finite()))
        || input
            .orbital_positions
            .iter()
            .any(|point| !point.iter().all(|value| value.is_finite()))
    {
        return Err(band_error(
            "NonFiniteBandCorepInput",
            "k points and orbital positions must be finite",
        ));
    }
    for operation in input.operations {
        if operation.representation.nrows() != dimension
            || operation.representation.ncols() != dimension
            || operation
                .representation
                .iter()
                .any(|value| !finite_complex(*value))
            || !operation.rotation.iter().all(|value| value.is_finite())
            || !operation.translation.iter().all(|value| value.is_finite())
            || !matches!(operation.trace_sign, -1.0 | 1.0)
        {
            return Err(band_error(
                "InvalidBandCorepSymmetry",
                "symmetry matrices, translations, representations, and trace signs must be finite and aligned",
            ));
        }
    }
    Ok(dimension)
}

fn modulo_one(value: f64) -> f64 {
    value - value.floor()
}

fn same_reduced_k(left: &Vector3<f64>, right: &Vector3<f64>) -> bool {
    left.iter().zip(right.iter()).all(|(left, right)| {
        let difference = left - right;
        (difference - difference.round()).abs() <= LITTLE_GROUP_TOLERANCE
    })
}

fn little_group(
    operations: &[NumericSymmetryOperation],
    reduced_kpoint: &Vector3<f64>,
) -> ExactResult<Vec<(usize, Vector3<f64>)>> {
    let reduced = reduced_kpoint.map(modulo_one);
    let mut result = Vec::new();
    for (index, operation) in operations.iter().enumerate() {
        if operation.antiunitary {
            continue;
        }
        let reciprocal = operation
            .rotation
            .transpose()
            .try_inverse()
            .ok_or_else(|| {
                band_error(
                    "InvalidBandCorepSymmetry",
                    format!("operation {} has a singular spatial rotation", index + 1),
                )
            })?;
        let image = reciprocal * reduced;
        if same_reduced_k(&image, &reduced) {
            result.push((index, image));
        }
    }
    if result.is_empty() {
        return Err(band_error(
            "InvalidBandCorepLittleGroup",
            "the requested k point has no unitary little-group operation",
        ));
    }
    Ok(result)
}

fn stable_gauge_hamiltonian(
    hamiltonian: &DMatrix<Complex<f64>>,
    reduced_kpoint: &Vector3<f64>,
    orbital_positions: &[Vector3<f64>],
) -> DMatrix<Complex<f64>> {
    let phases = orbital_positions
        .iter()
        .map(|position| Complex::from_polar(1.0, -TAU * reduced_kpoint.dot(position)))
        .collect::<Vec<_>>();
    DMatrix::from_fn(hamiltonian.nrows(), hamiltonian.ncols(), |row, column| {
        phases[row].conj() * hamiltonian[(row, column)] * phases[column]
    })
}

fn validate_hermitian(matrix: &DMatrix<Complex<f64>>) -> ExactResult<()> {
    for row in 0..matrix.nrows() {
        for column in 0..matrix.ncols() {
            if (matrix[(row, column)] - matrix[(column, row)].conj()).norm() > HERMITIAN_TOLERANCE {
                return Err(band_error(
                    "NonHermitianBandCorepHamiltonian",
                    "getTBBandCorep requires a Hermitian numerical Hamiltonian",
                ));
            }
        }
    }
    Ok(())
}

fn bloch_operation(
    operation: &NumericSymmetryOperation,
    transformed_kpoint: &Vector3<f64>,
    orbital_positions: &[Vector3<f64>],
) -> DMatrix<Complex<f64>> {
    let translation_phase = -TAU * operation.translation.dot(transformed_kpoint);
    DMatrix::from_fn(
        operation.representation.nrows(),
        operation.representation.ncols(),
        |row, column| {
            let displacement =
                orbital_positions[row] - operation.rotation * orbital_positions[column];
            let phase = translation_phase + TAU * transformed_kpoint.dot(&displacement);
            operation.representation[(row, column)] * Complex::from_polar(1.0, phase)
        },
    )
}

fn stable_energy_groups(energies: &[f64]) -> Vec<Vec<usize>> {
    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut previous_key = None;
    for (index, energy) in energies.iter().enumerate() {
        let key = (energy / STABLE_DEGENERACY_STEP).round();
        if previous_key.is_some_and(|previous: f64| previous.total_cmp(&key).is_eq()) {
            groups
                .last_mut()
                .expect("a previous key has a group")
                .push(index);
        } else {
            groups.push(vec![index]);
            previous_key = Some(key);
        }
    }
    groups
}

fn subspace_trace(
    eigenvectors: &DMatrix<Complex<f64>>,
    eigenvector_indices: &[usize],
    operation: &DMatrix<Complex<f64>>,
) -> Complex<f64> {
    eigenvector_indices
        .iter()
        .fold(Complex::new(0.0, 0.0), |total, &band| {
            let vector = eigenvectors.column(band);
            let transformed = operation * vector;
            total + vector.dotc(&transformed)
        })
}

pub fn compute_band_trace(input: &BandTraceInput<'_>) -> ExactResult<Vec<BandPointTrace>> {
    let dimension = validate_input(input)?;
    input
        .hamiltonians
        .iter()
        .zip(input.reduced_kpoints)
        .map(|(hamiltonian, reduced_kpoint)| {
            let gauged =
                stable_gauge_hamiltonian(hamiltonian, reduced_kpoint, input.orbital_positions);
            validate_hermitian(&gauged)?;
            let decomposition = SymmetricEigen::new(gauged);
            let mut order = (0..dimension).collect::<Vec<_>>();
            order.sort_by(|&left, &right| {
                decomposition.eigenvalues[left].total_cmp(&decomposition.eigenvalues[right])
            });
            let energies = order
                .iter()
                .map(|&index| decomposition.eigenvalues[index])
                .collect::<Vec<_>>();
            let eigenvectors = DMatrix::from_fn(dimension, dimension, |row, column| {
                decomposition.eigenvectors[(row, order[column])]
            });
            let groups = stable_energy_groups(&energies);
            let little = little_group(input.operations, reduced_kpoint)?;
            let operators = little
                .iter()
                .map(|(index, transformed)| {
                    (
                        *index,
                        bloch_operation(
                            &input.operations[*index],
                            transformed,
                            input.orbital_positions,
                        ),
                    )
                })
                .collect::<Vec<_>>();
            let mut degeneracies = vec![0; dimension];
            let mut traces = vec![Vec::new(); dimension];
            for group in groups {
                let characters = operators
                    .iter()
                    .map(|(index, operation)| {
                        subspace_trace(&eigenvectors, &group, operation)
                            * input.operations[*index].trace_sign
                    })
                    .collect::<Vec<_>>();
                for &band in &group {
                    degeneracies[band] = group.len();
                    traces[band].clone_from(&characters);
                }
            }
            if traces.iter().flatten().any(|value| !finite_complex(*value)) {
                return Err(band_error(
                    "NonFiniteBandCorepTrace",
                    "a band-character trace is not finite",
                ));
            }
            Ok(BandPointTrace {
                energies,
                degeneracies,
                little_symmetry_indices: operators.iter().map(|(index, _)| index + 1).collect(),
                traces,
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::{
        BandTraceInput, Complex, DMatrix, Matrix3, NumericSymmetryOperation, Vector3,
        compute_band_trace,
    };

    fn identity_operation(dimension: usize) -> NumericSymmetryOperation {
        NumericSymmetryOperation {
            rotation: Matrix3::identity(),
            translation: Vector3::zeros(),
            antiunitary: false,
            representation: DMatrix::identity(dimension, dimension),
            trace_sign: 1.0,
        }
    }

    #[test]
    fn one_band_trace_matches_the_stable_identity_character() {
        let hamiltonians = [DMatrix::from_element(1, 1, Complex::new(2.0, 0.0))];
        let kpoints = [Vector3::new(0.25, 0.0, 0.0)];
        let positions = [Vector3::zeros()];
        let operations = [identity_operation(1)];
        let result = compute_band_trace(&BandTraceInput {
            hamiltonians: &hamiltonians,
            reduced_kpoints: &kpoints,
            orbital_positions: &positions,
            operations: &operations,
        })
        .expect("band trace");
        assert_eq!(result[0].energies, [2.0]);
        assert_eq!(result[0].degeneracies, [1]);
        assert_eq!(result[0].little_symmetry_indices, [1]);
        assert!((result[0].traces[0][0] - Complex::new(1.0, 0.0)).norm() < 1.0e-12);
    }

    #[test]
    fn stable_degeneracy_is_repeated_for_every_band() {
        let hamiltonians = [DMatrix::zeros(2, 2)];
        let kpoints = [Vector3::zeros()];
        let positions = [Vector3::zeros(), Vector3::zeros()];
        let operations = [identity_operation(2)];
        let result = compute_band_trace(&BandTraceInput {
            hamiltonians: &hamiltonians,
            reduced_kpoints: &kpoints,
            orbital_positions: &positions,
            operations: &operations,
        })
        .expect("degenerate trace");
        assert_eq!(result[0].degeneracies, [2, 2]);
        for row in &result[0].traces {
            assert!((row[0] - Complex::new(2.0, 0.0)).norm() < 1.0e-12);
        }
    }

    #[test]
    fn nonhermitian_input_fails_explicitly() {
        let hamiltonians = [DMatrix::from_row_slice(
            2,
            2,
            &[
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
        )];
        let kpoints = [Vector3::zeros()];
        let positions = [Vector3::zeros(), Vector3::zeros()];
        let operations = [identity_operation(2)];
        let error = compute_band_trace(&BandTraceInput {
            hamiltonians: &hamiltonians,
            reduced_kpoints: &kpoints,
            orbital_positions: &positions,
            operations: &operations,
        })
        .expect_err("non-Hermitian input must fail");
        assert_eq!(error.tag(), "NonHermitianBandCorepHamiltonian");
    }
}
