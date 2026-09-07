use std::time::{Duration, Instant};

use magnetictb_properties::hermitian_eigenvalues;
use nalgebra::{Complex, DMatrix, DVector};

use crate::{FittingError, FittingResult};

#[derive(Clone, Debug)]
pub struct AffineBandModel {
    pub constant_matrices: Vec<DMatrix<Complex<f64>>>,
    pub coefficient_matrices: Vec<Vec<DMatrix<Complex<f64>>>>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum KPointNeighborhood {
    All,
    Radius {
        center: [f64; 3],
        radius: f64,
        periodic: bool,
    },
    Range {
        center: [f64; 3],
        range: [f64; 3],
        periodic: bool,
    },
}

#[derive(Clone, Debug)]
pub struct FitOptions {
    pub max_iterations: usize,
    pub time_constraint: Option<Duration>,
    pub energy_window: Option<(f64, f64)>,
    pub k_point_neighborhood: KPointNeighborhood,
    pub band_selection: Vec<usize>,
    pub residual_weights: Option<Vec<Vec<f64>>>,
    pub finite_difference_step: f64,
    pub initial_damping: f64,
    pub fit_tolerance: f64,
    pub hermitian_tolerance: f64,
}

#[derive(Clone, Debug)]
pub struct SelectionSummary {
    pub candidate_k_point_count: usize,
    pub used_k_point_indices: Vec<usize>,
    pub unique_k_point_count: usize,
    pub residual_count: usize,
    pub selection_restricted: bool,
    pub selected_band_indices_by_k_point: Vec<Vec<usize>>,
    pub band_selection: Vec<usize>,
    pub energy_window: Option<(f64, f64)>,
    pub k_point_neighborhood: KPointNeighborhood,
    pub explicit_weights: bool,
}

#[derive(Clone, Debug)]
pub struct FitOutput {
    pub values: Vec<f64>,
    pub loss: f64,
    pub initial_loss: f64,
    pub iterations: usize,
    pub converged: bool,
    pub final_damping: f64,
    pub initial_bands: Vec<Vec<f64>>,
    pub fitted_bands: Vec<Vec<f64>>,
    pub reference_bands: Vec<Vec<f64>>,
    pub selection_mask: Vec<Vec<bool>>,
    pub selection: SelectionSummary,
}

#[derive(Clone, Debug)]
struct SelectedProblem {
    model: AffineBandModel,
    reference_bands: Vec<Vec<f64>>,
    weights: Vec<Vec<f64>>,
    mask: Vec<Vec<bool>>,
    summary: SelectionSummary,
}

#[derive(Clone, Debug)]
struct OptimizerOutput {
    values: Vec<f64>,
    loss: f64,
    iterations: usize,
    converged: bool,
    final_damping: f64,
}

fn invalid(tag: &str, detail: impl Into<String>) -> FittingError {
    FittingError::new(tag, detail)
}

fn finite_complex(value: Complex<f64>) -> bool {
    value.re.is_finite() && value.im.is_finite()
}

fn validate_model(
    model: &AffineBandModel,
    k_points: &[[f64; 3]],
    reference_bands: &[Vec<f64>],
    parameter_count: usize,
) -> FittingResult<usize> {
    if model.constant_matrices.is_empty()
        || model.constant_matrices.len() != k_points.len()
        || reference_bands.len() != k_points.len()
        || model.coefficient_matrices.len() != parameter_count
        || model
            .coefficient_matrices
            .iter()
            .any(|matrices| matrices.len() != k_points.len())
    {
        return Err(invalid(
            "InvalidAffineBandModel",
            "constant and coefficient matrix samples must match all reference k points",
        ));
    }
    let dimension = model.constant_matrices[0].nrows();
    if dimension == 0
        || model.constant_matrices.iter().any(|matrix| {
            matrix.nrows() != dimension
                || matrix.ncols() != dimension
                || matrix.iter().copied().any(|value| !finite_complex(value))
        })
        || model.coefficient_matrices.iter().flatten().any(|matrix| {
            matrix.nrows() != dimension
                || matrix.ncols() != dimension
                || matrix.iter().copied().any(|value| !finite_complex(value))
        })
    {
        return Err(invalid(
            "InvalidHamiltonian",
            "Hamiltonian samples must be finite nonempty square matrices of one dimension",
        ));
    }
    if k_points.iter().flatten().any(|value| !value.is_finite())
        || reference_bands
            .iter()
            .any(|bands| bands.len() != dimension || bands.iter().any(|value| !value.is_finite()))
    {
        return Err(invalid(
            "InvalidReferenceBands",
            "reference data must contain finite three-component k points and one energy per band",
        ));
    }
    Ok(dimension)
}

fn validate_options(options: &FitOptions, dimension: usize) -> FittingResult<()> {
    if options.max_iterations == 0
        || options
            .time_constraint
            .is_some_and(|duration| duration.is_zero())
        || !options.finite_difference_step.is_finite()
        || options.finite_difference_step <= 0.0
        || !options.initial_damping.is_finite()
        || options.initial_damping <= 0.0
        || !options.fit_tolerance.is_finite()
        || options.fit_tolerance <= 0.0
        || !options.hermitian_tolerance.is_finite()
        || options.hermitian_tolerance < 0.0
    {
        return Err(invalid(
            "InvalidFittingOption",
            "iteration, time, finite-difference, damping, fit, and Hermitian options are invalid",
        ));
    }
    if let Some((minimum, maximum)) = options.energy_window
        && (!minimum.is_finite() || !maximum.is_finite() || minimum > maximum)
    {
        return Err(invalid(
            "InvalidEnergyWindow",
            "EnergyWindow must be an inclusive finite interval",
        ));
    }
    if options.band_selection.is_empty()
        || options
            .band_selection
            .iter()
            .any(|&index| index >= dimension)
        || options
            .band_selection
            .windows(2)
            .any(|pair| pair[0] >= pair[1])
    {
        return Err(invalid(
            "InvalidBandSelection",
            "band indices must be unique, strictly increasing, and in range",
        ));
    }
    match &options.k_point_neighborhood {
        KPointNeighborhood::All => {}
        KPointNeighborhood::Radius { center, radius, .. } => {
            if center.iter().any(|value| !value.is_finite()) || !radius.is_finite() || *radius < 0.0
            {
                return Err(invalid(
                    "InvalidKPointNeighborhood",
                    "radius neighborhoods require a finite center and nonnegative radius",
                ));
            }
        }
        KPointNeighborhood::Range { center, range, .. } => {
            if center.iter().any(|value| !value.is_finite())
                || range.iter().any(|value| !value.is_finite() || *value < 0.0)
            {
                return Err(invalid(
                    "InvalidKPointNeighborhood",
                    "component neighborhoods require a finite center and nonnegative range",
                ));
            }
        }
    }
    Ok(())
}

fn point_selected(point: [f64; 3], neighborhood: &KPointNeighborhood) -> bool {
    match neighborhood {
        KPointNeighborhood::All => true,
        KPointNeighborhood::Radius {
            center,
            radius,
            periodic,
        } => {
            let mut delta = [0.0; 3];
            for axis in 0..3 {
                delta[axis] = point[axis] - center[axis];
                if *periodic {
                    delta[axis] -= delta[axis].round();
                }
            }
            delta.iter().map(|value| value * value).sum::<f64>().sqrt() <= *radius
        }
        KPointNeighborhood::Range {
            center,
            range,
            periodic,
        } => (0..3).all(|axis| {
            let mut delta = point[axis] - center[axis];
            if *periodic {
                delta -= delta.round();
            }
            delta.abs() <= range[axis]
        }),
    }
}

fn select_problem(
    model: &AffineBandModel,
    k_points: &[[f64; 3]],
    reference_bands: &[Vec<f64>],
    k_range: &[usize],
    options: &FitOptions,
    parameter_count: usize,
) -> FittingResult<SelectedProblem> {
    let dimension = validate_model(model, k_points, reference_bands, parameter_count)?;
    validate_options(options, dimension)?;
    if k_range.is_empty() || k_range.iter().any(|&index| index >= k_points.len()) {
        return Err(invalid(
            "InvalidKRange",
            "kRange must be a nonempty list of in-range reference-data indices",
        ));
    }
    let all_weights = match &options.residual_weights {
        None => vec![vec![1.0; dimension]; k_points.len()],
        Some(weights)
            if weights.len() == k_points.len()
                && weights.iter().all(|row| {
                    row.len() == dimension
                        && row.iter().all(|value| value.is_finite() && *value >= 0.0)
                }) =>
        {
            weights.clone()
        }
        Some(_) => {
            return Err(invalid(
                "InvalidResidualWeights",
                "ResidualWeights must be a finite nonnegative reference-by-band matrix",
            ));
        }
    };
    let candidate = k_range
        .iter()
        .copied()
        .filter(|&index| point_selected(k_points[index], &options.k_point_neighborhood))
        .collect::<Vec<_>>();
    if candidate.is_empty() {
        return Err(invalid(
            "EmptyFittingSelection",
            "the combined selection contains no residuals",
        ));
    }
    let band_allowed = (0..dimension)
        .map(|index| options.band_selection.contains(&index))
        .collect::<Vec<_>>();
    let mut used_indices = Vec::new();
    let mut selected_reference = Vec::new();
    let mut selected_weights = Vec::new();
    let mut masks = Vec::new();
    for &index in &candidate {
        let mask = (0..dimension)
            .map(|band| {
                band_allowed[band]
                    && all_weights[index][band] > 0.0
                    && options.energy_window.is_none_or(|(minimum, maximum)| {
                        reference_bands[index][band] >= minimum
                            && reference_bands[index][band] <= maximum
                    })
            })
            .collect::<Vec<_>>();
        if mask.contains(&true) {
            used_indices.push(index);
            selected_reference.push(reference_bands[index].clone());
            selected_weights.push(all_weights[index].clone());
            masks.push(mask);
        }
    }
    let residual_count = masks.iter().flatten().filter(|&&value| value).count();
    if residual_count == 0 {
        return Err(invalid(
            "EmptyFittingSelection",
            "the combined selection contains no residuals",
        ));
    }
    let selection_restricted = options.energy_window.is_some()
        || !matches!(options.k_point_neighborhood, KPointNeighborhood::All)
        || options.band_selection != (0..dimension).collect::<Vec<_>>()
        || options.residual_weights.is_some();
    if selection_restricted && residual_count < parameter_count {
        return Err(invalid(
            "InsufficientFittingSelection",
            format!(
                "the selection contains {residual_count} residuals for {parameter_count} parameters"
            ),
        ));
    }
    let selected_band_indices_by_k_point = masks
        .iter()
        .map(|mask| {
            mask.iter()
                .enumerate()
                .filter_map(|(index, &selected)| selected.then_some(index))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let selected_model = AffineBandModel {
        constant_matrices: used_indices
            .iter()
            .map(|&index| model.constant_matrices[index].clone())
            .collect(),
        coefficient_matrices: model
            .coefficient_matrices
            .iter()
            .map(|matrices| {
                used_indices
                    .iter()
                    .map(|&index| matrices[index].clone())
                    .collect()
            })
            .collect(),
    };
    let mut unique = used_indices.clone();
    unique.sort_unstable();
    unique.dedup();
    Ok(SelectedProblem {
        model: selected_model,
        reference_bands: selected_reference,
        weights: selected_weights,
        mask: masks,
        summary: SelectionSummary {
            candidate_k_point_count: candidate.len(),
            used_k_point_indices: used_indices,
            unique_k_point_count: unique.len(),
            residual_count,
            selection_restricted,
            selected_band_indices_by_k_point,
            band_selection: options.band_selection.clone(),
            energy_window: options.energy_window,
            k_point_neighborhood: options.k_point_neighborhood.clone(),
            explicit_weights: options.residual_weights.is_some(),
        },
    })
}

fn numeric_hamiltonians(model: &AffineBandModel, values: &[f64]) -> Vec<DMatrix<Complex<f64>>> {
    model
        .constant_matrices
        .iter()
        .enumerate()
        .map(|(point, constant)| {
            model
                .coefficient_matrices
                .iter()
                .zip(values)
                .fold(constant.clone(), |matrix, (samples, &value)| {
                    matrix + &samples[point] * Complex::new(value, 0.0)
                })
        })
        .collect()
}

fn bands(
    model: &AffineBandModel,
    values: &[f64],
    hermitian_tolerance: f64,
) -> FittingResult<Vec<Vec<f64>>> {
    numeric_hamiltonians(model, values)
        .iter()
        .map(|matrix| {
            hermitian_eigenvalues(matrix, hermitian_tolerance)
                .map_err(|error| invalid(error.tag(), error.detail()))
        })
        .collect()
}

fn residual(
    problem: &SelectedProblem,
    values: &[f64],
    hermitian_tolerance: f64,
) -> FittingResult<Vec<f64>> {
    let model_bands = bands(&problem.model, values, hermitian_tolerance)?;
    let mut result = Vec::with_capacity(problem.summary.residual_count);
    for (point, model) in model_bands.iter().enumerate() {
        for (band, &model_energy) in model.iter().enumerate() {
            if problem.mask[point][band] {
                result.push(
                    problem.weights[point][band].sqrt()
                        * (model_energy - problem.reference_bands[point][band]),
                );
            }
        }
    }
    if result.iter().any(|value| !value.is_finite()) {
        return Err(invalid(
            "InvalidResidual",
            "band residual contains a non-finite value",
        ));
    }
    Ok(result)
}

fn squared_norm(values: &[f64]) -> f64 {
    values.iter().map(|value| value * value).sum()
}

fn timed_out(start: Instant, constraint: Option<Duration>) -> bool {
    constraint.is_some_and(|duration| start.elapsed() >= duration)
}

fn levenberg_marquardt(
    problem: &SelectedProblem,
    initial_values: &[f64],
    options: &FitOptions,
) -> FittingResult<OptimizerOutput> {
    let mut values = initial_values.to_vec();
    let mut residual_values = residual(problem, &values, options.hermitian_tolerance)?;
    let mut loss = squared_norm(&residual_values);
    let mut damping = options.initial_damping;
    let start = Instant::now();
    let parameter_count = values.len();
    let residual_count = residual_values.len();
    let mut converged = false;
    let mut completed_iteration = 0;

    for iteration in 1..=options.max_iterations {
        completed_iteration = iteration;
        if timed_out(start, options.time_constraint) {
            return Err(invalid(
                "FittingTimeConstraint",
                "the fit exceeded TimeConstraint",
            ));
        }
        let step_sizes = values
            .iter()
            .map(|value| options.finite_difference_step * value.abs().max(1.0))
            .collect::<Vec<_>>();
        let mut jacobian = DMatrix::<f64>::zeros(residual_count, parameter_count);
        for parameter in 0..parameter_count {
            if timed_out(start, options.time_constraint) {
                return Err(invalid(
                    "FittingTimeConstraint",
                    "the fit exceeded TimeConstraint",
                ));
            }
            let mut plus = values.clone();
            let mut minus = values.clone();
            plus[parameter] += step_sizes[parameter];
            minus[parameter] -= step_sizes[parameter];
            let plus_residual = residual(problem, &plus, options.hermitian_tolerance)?;
            let minus_residual = residual(problem, &minus, options.hermitian_tolerance)?;
            if plus_residual.len() != residual_count || minus_residual.len() != residual_count {
                return Err(invalid(
                    "InvalidResidual",
                    "finite-difference residual length changed",
                ));
            }
            for row in 0..residual_count {
                jacobian[(row, parameter)] =
                    (plus_residual[row] - minus_residual[row]) / (2.0 * step_sizes[parameter]);
            }
        }
        let residual_vector = DVector::from_column_slice(&residual_values);
        let normal = jacobian.transpose() * &jacobian;
        let gradient = jacobian.transpose() * residual_vector;
        if gradient.amax() <= options.fit_tolerance {
            converged = true;
            break;
        }
        let mut damped = normal.clone();
        for index in 0..parameter_count {
            damped[(index, index)] += damping * normal[(index, index)].max(1.0);
        }
        let Some(delta) = damped.lu().solve(&(-gradient)) else {
            damping *= 10.0;
            continue;
        };
        if delta.iter().any(|value| !value.is_finite()) {
            damping *= 10.0;
            continue;
        }
        let candidate_values = values
            .iter()
            .zip(delta.iter())
            .map(|(value, change)| value + change)
            .collect::<Vec<_>>();
        let candidate_residual = residual(problem, &candidate_values, options.hermitian_tolerance)?;
        let candidate_loss = squared_norm(&candidate_residual);
        if candidate_loss < loss {
            let improvement = loss - candidate_loss;
            let delta_norm = delta.norm();
            values = candidate_values;
            residual_values = candidate_residual;
            loss = candidate_loss;
            damping = (damping / 3.0).max(1.0e-15);
            let values_norm = DVector::from_column_slice(&values).norm();
            if delta_norm <= options.fit_tolerance * (1.0 + values_norm)
                || improvement <= options.fit_tolerance * loss.max(1.0)
            {
                converged = true;
                break;
            }
        } else {
            damping = (damping * 10.0).min(1.0e15);
        }
    }
    Ok(OptimizerOutput {
        values,
        loss,
        iterations: completed_iteration,
        converged,
        final_damping: damping,
    })
}

pub fn fit_bands(
    model: &AffineBandModel,
    k_points: &[[f64; 3]],
    reference_bands: &[Vec<f64>],
    k_range: &[usize],
    initial_values: &[f64],
    options: &FitOptions,
) -> FittingResult<FitOutput> {
    if initial_values.iter().any(|value| !value.is_finite()) {
        return Err(invalid(
            "InvalidInitialParameters",
            "every fitting parameter requires a finite real initial value",
        ));
    }
    let selected = select_problem(
        model,
        k_points,
        reference_bands,
        k_range,
        options,
        initial_values.len(),
    )?;
    let initial_bands = bands(&selected.model, initial_values, options.hermitian_tolerance)?;
    let initial_residual = residual(&selected, initial_values, options.hermitian_tolerance)?;
    let initial_loss = squared_norm(&initial_residual);
    let optimizer = if initial_values.is_empty() {
        OptimizerOutput {
            values: Vec::new(),
            loss: initial_loss,
            iterations: 0,
            converged: true,
            final_damping: options.initial_damping,
        }
    } else {
        levenberg_marquardt(&selected, initial_values, options)?
    };
    let fitted_bands = bands(
        &selected.model,
        &optimizer.values,
        options.hermitian_tolerance,
    )?;
    Ok(FitOutput {
        values: optimizer.values,
        loss: optimizer.loss,
        initial_loss,
        iterations: optimizer.iterations,
        converged: optimizer.converged,
        final_damping: optimizer.final_damping,
        initial_bands,
        fitted_bands,
        reference_bands: selected.reference_bands,
        selection_mask: selected.mask,
        selection: selected.summary,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scalar(value: f64) -> DMatrix<Complex<f64>> {
        DMatrix::from_element(1, 1, Complex::new(value, 0.0))
    }

    fn default_options() -> FitOptions {
        FitOptions {
            max_iterations: 100,
            time_constraint: Some(Duration::from_secs(60)),
            energy_window: None,
            k_point_neighborhood: KPointNeighborhood::All,
            band_selection: vec![0],
            residual_weights: None,
            finite_difference_step: 1.0e-5,
            initial_damping: 1.0e-3,
            fit_tolerance: 1.0e-8,
            hermitian_tolerance: 1.0e-10,
        }
    }

    #[test]
    fn affine_fit_recovers_parameter_and_reports_selection() {
        let points = vec![[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]];
        let model = AffineBandModel {
            constant_matrices: vec![scalar(0.0), scalar(0.0)],
            coefficient_matrices: vec![vec![scalar(1.0), scalar(-1.0)]],
        };
        let output = fit_bands(
            &model,
            &points,
            &[vec![0.3], vec![-0.3]],
            &[0, 1],
            &[0.0],
            &default_options(),
        )
        .expect("fit must succeed");
        assert!((output.values[0] - 0.3).abs() < 1.0e-8);
        assert!(output.loss < 1.0e-12);
        assert_eq!(output.selection.used_k_point_indices, vec![0, 1]);
        assert_eq!(output.selection.residual_count, 2);
    }

    #[test]
    fn periodic_radius_uses_minimum_image() {
        let points = vec![
            [15.0 / 16.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0 / 16.0, 0.0, 0.0],
            [0.25, 0.0, 0.0],
        ];
        let model = AffineBandModel {
            constant_matrices: vec![scalar(0.0); 4],
            coefficient_matrices: vec![vec![scalar(1.0); 4]],
        };
        let mut options = default_options();
        options.k_point_neighborhood = KPointNeighborhood::Radius {
            center: [0.0; 3],
            radius: 0.1,
            periodic: true,
        };
        let output = fit_bands(
            &model,
            &points,
            &vec![vec![0.3]; 4],
            &[0, 1, 2, 3],
            &[0.0],
            &options,
        )
        .expect("fit must succeed");
        assert_eq!(output.selection.used_k_point_indices, vec![0, 1, 2]);
    }

    #[test]
    fn restricted_underdetermined_selection_fails() {
        let model = AffineBandModel {
            constant_matrices: vec![scalar(0.0)],
            coefficient_matrices: vec![vec![scalar(1.0)], vec![scalar(1.0)]],
        };
        let mut options = default_options();
        options.band_selection = vec![0];
        options.energy_window = Some((-1.0, 1.0));
        let error = fit_bands(
            &model,
            &[[0.0; 3]],
            &[vec![0.0]],
            &[0],
            &[0.0, 0.0],
            &options,
        )
        .expect_err("restricted underdetermined fit must fail");
        assert_eq!(error.tag(), "InsufficientFittingSelection");
    }
}
