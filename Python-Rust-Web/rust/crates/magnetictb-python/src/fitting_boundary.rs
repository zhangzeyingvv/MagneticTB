use std::time::Duration;

use magnetictb_fitting::{
    AffineBandModel, Complex, DMatrix, FitOptions, FittingError, FittingResult, KPointNeighborhood,
    fit_bands, parse_vasp_eigenval,
};
use serde_json::{Map, Value, json};

fn malformed(detail: impl Into<String>) -> FittingError {
    FittingError::new("MalformedFittingInput", detail)
}

fn object(value: &Value) -> FittingResult<&Map<String, Value>> {
    value
        .as_object()
        .ok_or_else(|| malformed("expected a JSON object"))
}

fn array<'a>(value: &'a Value, field: &str) -> FittingResult<&'a [Value]> {
    value
        .as_array()
        .map(Vec::as_slice)
        .ok_or_else(|| malformed(format!("{field} must be an array")))
}

fn required<'a>(record: &'a Map<String, Value>, field: &str) -> FittingResult<&'a Value> {
    record
        .get(field)
        .ok_or_else(|| malformed(format!("missing required field {field}")))
}

fn finite_real(value: &Value, field: &str) -> FittingResult<f64> {
    let result = value
        .as_f64()
        .ok_or_else(|| malformed(format!("{field} must be a real number")))?;
    if result.is_finite() {
        Ok(result)
    } else {
        Err(malformed(format!("{field} must be finite")))
    }
}

fn usize_value(value: &Value, field: &str) -> FittingResult<usize> {
    value
        .as_u64()
        .and_then(|number| usize::try_from(number).ok())
        .ok_or_else(|| malformed(format!("{field} must be a nonnegative integer")))
}

fn real_vector(value: &Value, field: &str) -> FittingResult<Vec<f64>> {
    array(value, field)?
        .iter()
        .enumerate()
        .map(|(index, value)| finite_real(value, &format!("{field}[{index}]")))
        .collect()
}

fn real_three(value: &Value, field: &str) -> FittingResult<[f64; 3]> {
    let values = real_vector(value, field)?;
    if values.len() != 3 {
        return Err(malformed(format!("{field} must have length three")));
    }
    Ok([values[0], values[1], values[2]])
}

fn complex_value(value: &Value, field: &str) -> FittingResult<Complex<f64>> {
    if value.is_number() {
        return finite_real(value, field).map(|real| Complex::new(real, 0.0));
    }
    let record = object(value)?;
    Ok(Complex::new(
        finite_real(required(record, "real")?, &format!("{field}.real"))?,
        finite_real(
            required(record, "imaginary")?,
            &format!("{field}.imaginary"),
        )?,
    ))
}

fn complex_matrix(value: &Value, field: &str) -> FittingResult<DMatrix<Complex<f64>>> {
    let rows = array(value, field)?;
    if rows.is_empty() {
        return Err(malformed(format!("{field} must be nonempty")));
    }
    let columns = array(&rows[0], &format!("{field}[0]"))?.len();
    if columns == 0 {
        return Err(malformed(format!("{field} rows must be nonempty")));
    }
    let mut entries = Vec::with_capacity(rows.len() * columns);
    for (row_index, row) in rows.iter().enumerate() {
        let row = array(row, &format!("{field}[{row_index}]"))?;
        if row.len() != columns {
            return Err(malformed(format!("{field} rows must have equal length")));
        }
        for (column_index, value) in row.iter().enumerate() {
            entries.push(complex_value(
                value,
                &format!("{field}[{row_index}][{column_index}]"),
            )?);
        }
    }
    Ok(DMatrix::from_row_slice(rows.len(), columns, &entries))
}

fn complex_matrices(value: &Value, field: &str) -> FittingResult<Vec<DMatrix<Complex<f64>>>> {
    array(value, field)?
        .iter()
        .enumerate()
        .map(|(index, value)| complex_matrix(value, &format!("{field}[{index}]")))
        .collect()
}

fn coefficient_matrices(
    value: &Value,
    constants: &[DMatrix<Complex<f64>>],
) -> FittingResult<Vec<Vec<DMatrix<Complex<f64>>>>> {
    array(value, "unit_parameter_matrices")?
        .iter()
        .enumerate()
        .map(|(parameter, value)| {
            let units = complex_matrices(value, &format!("unit_parameter_matrices[{parameter}]"))?;
            if units.len() != constants.len() {
                return Err(malformed(
                    "each unit-parameter sample set must match constant_matrices",
                ));
            }
            Ok(units
                .into_iter()
                .zip(constants)
                .map(|(unit, constant)| unit - constant)
                .collect())
        })
        .collect()
}

fn neighborhood(value: &Value) -> FittingResult<KPointNeighborhood> {
    let record = object(value)?;
    let mode = required(record, "mode")?
        .as_str()
        .ok_or_else(|| malformed("k_point_neighborhood.mode must be a string"))?;
    if mode == "All" {
        return Ok(KPointNeighborhood::All);
    }
    let center = real_three(required(record, "center")?, "k_point_neighborhood.center")?;
    let periodic = required(record, "periodic")?
        .as_bool()
        .ok_or_else(|| malformed("k_point_neighborhood.periodic must be Boolean"))?;
    match mode {
        "Radius" => Ok(KPointNeighborhood::Radius {
            center,
            radius: finite_real(required(record, "radius")?, "k_point_neighborhood.radius")?,
            periodic,
        }),
        "Range" => Ok(KPointNeighborhood::Range {
            center,
            range: real_three(required(record, "range")?, "k_point_neighborhood.range")?,
            periodic,
        }),
        _ => Err(FittingError::new(
            "InvalidKPointNeighborhood",
            format!("unsupported KPointNeighborhood mode {mode}"),
        )),
    }
}

fn neighborhood_json(value: &KPointNeighborhood) -> Value {
    match value {
        KPointNeighborhood::All => json!({"Mode": "All"}),
        KPointNeighborhood::Radius {
            center,
            radius,
            periodic,
        } => json!({
            "Mode": "Radius", "Center": center, "Radius": radius, "Periodic": periodic
        }),
        KPointNeighborhood::Range {
            center,
            range,
            periodic,
        } => json!({
            "Mode": "Range", "Center": center, "Range": range, "Periodic": periodic
        }),
    }
}

fn parse_options(record: &Map<String, Value>, dimension: usize) -> FittingResult<FitOptions> {
    let options = object(required(record, "options")?)?;
    let energy_window = options
        .get("energy_window")
        .filter(|value| !value.is_null())
        .map(|value| {
            let values = real_vector(value, "energy_window")?;
            if values.len() != 2 {
                return Err(malformed("energy_window must contain two values"));
            }
            Ok((values[0], values[1]))
        })
        .transpose()?;
    let weights = options
        .get("residual_weights")
        .filter(|value| !value.is_null())
        .map(|value| {
            array(value, "residual_weights")?
                .iter()
                .enumerate()
                .map(|(index, row)| real_vector(row, &format!("residual_weights[{index}]")))
                .collect::<FittingResult<Vec<_>>>()
        })
        .transpose()?;
    let time_constraint = options
        .get("time_constraint")
        .filter(|value| !value.is_null())
        .map(|value| {
            finite_real(value, "time_constraint").and_then(|seconds| {
                if seconds > 0.0 {
                    Ok(Duration::from_secs_f64(seconds))
                } else {
                    Err(FittingError::new(
                        "InvalidFittingOption",
                        "TimeConstraint must be positive or None",
                    ))
                }
            })
        })
        .transpose()?;
    let band_selection = options
        .get("band_selection")
        .map(|value| {
            array(value, "band_selection")?
                .iter()
                .enumerate()
                .map(|(index, value)| usize_value(value, &format!("band_selection[{index}]")))
                .collect::<FittingResult<Vec<_>>>()
        })
        .transpose()?
        .unwrap_or_else(|| (0..dimension).collect());
    Ok(FitOptions {
        max_iterations: usize_value(required(options, "max_iterations")?, "max_iterations")?,
        time_constraint,
        energy_window,
        k_point_neighborhood: neighborhood(required(options, "k_point_neighborhood")?)?,
        band_selection,
        residual_weights: weights,
        finite_difference_step: finite_real(
            required(options, "finite_difference_step")?,
            "finite_difference_step",
        )?,
        initial_damping: finite_real(required(options, "initial_damping")?, "initial_damping")?,
        fit_tolerance: finite_real(required(options, "fit_tolerance")?, "fit_tolerance")?,
        hermitian_tolerance: finite_real(
            required(options, "hermitian_tolerance")?,
            "hermitian_tolerance",
        )?,
    })
}

#[allow(clippy::too_many_lines)]
fn compute_value(payload: &str) -> FittingResult<Value> {
    let value: Value = serde_json::from_str(payload)
        .map_err(|error| malformed(format!("invalid JSON: {error}")))?;
    let request = object(&value)?;
    let operation = required(request, "operation")?
        .as_str()
        .ok_or_else(|| malformed("operation must be a string"))?;
    if operation == "parse_vasp_eigenval" {
        let path = required(request, "path")?
            .as_str()
            .ok_or_else(|| malformed("path must be a string"))?;
        let text = std::fs::read_to_string(path).map_err(|error| {
            FittingError::new(
                "VaspEigenvalFile",
                format!("could not read EIGENVAL file {path}: {error}"),
            )
        })?;
        let records = parse_vasp_eigenval(
            &text,
            finite_real(required(request, "fermi_energy")?, "fermi_energy")?,
            usize_value(required(request, "spin")?, "spin")?,
            usize_value(required(request, "start_band")?, "start_band")?,
            usize_value(required(request, "end_band")?, "end_band")?,
        )?;
        return Ok(json!({
            "Schema": "MagneticTBVaspEigenvalues",
            "SchemaVersion": 1,
            "Records": records.into_iter().map(|record| {
                json!([record.k_point, record.energies])
            }).collect::<Vec<_>>()
        }));
    }
    if operation != "fit_bands" {
        return Err(FittingError::new(
            "UnsupportedOperation",
            format!("unsupported fitting operation: {operation}"),
        ));
    }
    let parameter_names = array(required(request, "parameter_names")?, "parameter_names")?
        .iter()
        .map(|value| {
            value
                .as_str()
                .map(str::to_owned)
                .ok_or_else(|| malformed("parameter_names must contain strings"))
        })
        .collect::<FittingResult<Vec<_>>>()?;
    let constants = complex_matrices(required(request, "constant_matrices")?, "constant_matrices")?;
    let model = AffineBandModel {
        coefficient_matrices: coefficient_matrices(
            required(request, "unit_parameter_matrices")?,
            &constants,
        )?,
        constant_matrices: constants,
    };
    let k_points = array(required(request, "k_points")?, "k_points")?
        .iter()
        .enumerate()
        .map(|(index, value)| real_three(value, &format!("k_points[{index}]")))
        .collect::<FittingResult<Vec<_>>>()?;
    let reference_bands = array(required(request, "reference_bands")?, "reference_bands")?
        .iter()
        .enumerate()
        .map(|(index, value)| real_vector(value, &format!("reference_bands[{index}]")))
        .collect::<FittingResult<Vec<_>>>()?;
    let k_range = array(required(request, "k_range")?, "k_range")?
        .iter()
        .enumerate()
        .map(|(index, value)| usize_value(value, &format!("k_range[{index}]")))
        .collect::<FittingResult<Vec<_>>>()?;
    let initial_values = real_vector(required(request, "initial_values")?, "initial_values")?;
    if initial_values.len() != parameter_names.len() {
        return Err(FittingError::new(
            "InvalidInitialParameters",
            "initial_values must assign every fitting parameter",
        ));
    }
    let dimension = model.constant_matrices.first().map_or(0, DMatrix::nrows);
    let options = parse_options(request, dimension)?;
    let output = fit_bands(
        &model,
        &k_points,
        &reference_bands,
        &k_range,
        &initial_values,
        &options,
    )?;
    let selection = &output.selection;
    let energy_window = selection.energy_window.map_or_else(
        || json!("All"),
        |(minimum, maximum)| json!([minimum, maximum]),
    );
    Ok(json!({
        "Schema": "MagneticTBFittingResult",
        "SchemaVersion": 1,
        "FittedParams": parameter_names.iter().zip(&output.values).map(|(name, value)| {
            json!({"name": name, "value": value})
        }).collect::<Vec<_>>(),
        "Objective": if selection.explicit_weights {
            "WeightedSumSquaredBandResiduals"
        } else {
            "SumSquaredBandResiduals"
        },
        "Optimizer": if parameter_names.is_empty() { "None" } else { "LevenbergMarquardt" },
        "ParameterModel": if parameter_names.is_empty() {
            "Constant"
        } else {
            "CachedAffineCoefficientMatrices"
        },
        "Converged": output.converged,
        "Iterations": output.iterations,
        "FinalDamping": output.final_damping,
        "LSQForInitParams": output.initial_loss,
        "LSQForFittedParams": output.loss,
        "CandidateKPointCount": selection.candidate_k_point_count,
        "UsedKPointCount": selection.used_k_point_indices.len(),
        "UniqueKPointCount": selection.unique_k_point_count,
        "UsedKPointIndices": selection.used_k_point_indices.iter().map(|index| index + 1).collect::<Vec<_>>(),
        "ResidualCount": selection.residual_count,
        "SelectionRestricted": selection.selection_restricted,
        "SelectedBandIndicesByKPoint": selection.selected_band_indices_by_k_point.iter().map(|indices| {
            indices.iter().map(|index| index + 1).collect::<Vec<_>>()
        }).collect::<Vec<_>>(),
        "BandSelection": selection.band_selection.iter().map(|index| index + 1).collect::<Vec<_>>(),
        "EnergyWindow": energy_window,
        "KPointNeighborhood": neighborhood_json(&selection.k_point_neighborhood),
        "Weighting": if selection.explicit_weights { "ExplicitMatrix" } else { "Uniform" },
        "BandPairing": "SortedModelEigenvalueToReferencePosition",
        "InitialBands": output.initial_bands,
        "FittedBands": output.fitted_bands,
        "ReferenceBands": output.reference_bands,
        "SelectionMask": output.selection_mask,
    }))
}

pub fn compute(payload: &str) -> FittingResult<String> {
    serde_json::to_string(&compute_value(payload)?)
        .map_err(|error| malformed(format!("result serialization failed: {error}")))
}
