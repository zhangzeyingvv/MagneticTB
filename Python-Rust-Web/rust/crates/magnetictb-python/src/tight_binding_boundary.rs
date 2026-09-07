use std::sync::Arc;

use cyclotomic_nullspace::fixture::{
    array, bool_field, context_from_json, encoded_element_from_json, encoded_matrix_from_json,
    matrix_to_json, object, required, size_field, string_field,
};
use cyclotomic_nullspace::{Element, ExactError, ExactResult};
use magnetictb_tight_binding::{
    BlochTerm, assemble_bloch_hamiltonian, propagate_hopping, reconstruct_rectangular_hopping,
};
use serde::de::DeserializeOwned;
use serde_json::{Value, json};

fn malformed(detail: impl Into<String>) -> ExactError {
    ExactError::new("MalformedSerialization", detail)
}

fn from_value<T: DeserializeOwned>(value: &Value, field: &str) -> ExactResult<T> {
    serde_json::from_value(value.clone())
        .map_err(|error| malformed(format!("invalid {field}: {error}")))
}

fn elements_from_json(
    context: &Arc<cyclotomic_nullspace::CyclotomicContext>,
    value: &Value,
) -> ExactResult<Vec<Element>> {
    array(value)?
        .iter()
        .map(|element| encoded_element_from_json(context, element))
        .collect()
}

fn term_from_json(
    context: &Arc<cyclotomic_nullspace::CyclotomicContext>,
    value: &Value,
) -> ExactResult<BlochTerm> {
    let record = object(value)?;
    let phase_sign = record
        .get("phase_sign")
        .map_or(Ok(1), |value| from_value::<i64>(value, "phase_sign"))?;
    Ok(BlochTerm {
        bond_index: record
            .get("bond_index")
            .map_or(Ok(0), |value| from_value::<usize>(value, "bond_index"))?,
        orbit_index: record
            .get("orbit_index")
            .map_or(Ok(0), |value| from_value::<usize>(value, "orbit_index"))?,
        row_block: size_field(record, "row_block")?,
        column_block: size_field(record, "column_block")?,
        displacement: from_value(required(record, "displacement")?, "displacement")?,
        phase_sign,
        matrix: encoded_matrix_from_json(context, required(record, "matrix")?)?,
    })
}

fn evaluate(request: &serde_json::Map<String, Value>) -> ExactResult<Value> {
    let context = context_from_json(required(request, "context")?, 128)?;
    match string_field(request, "operation")? {
        "validate_static_bond_data" => {
            let value = object(required(request, "static")?)?;
            for field in ["bonds", "directed_orbits", "site_dimensions"] {
                required(value, field)?;
            }
            Ok(json!({"verified": true}))
        }
        "validate_constraint_problem" => {
            let value = object(required(request, "problem")?)?;
            for field in ["coordinate_dimension", "constraint_blocks"] {
                required(value, field)?;
            }
            Ok(json!({"verified": true}))
        }
        "validate_solved_bond_layers" => {
            for (layer, fields) in [
                ("static", &["bonds", "directed_orbits"][..]),
                (
                    "problem",
                    &["coordinate_dimension", "constraint_blocks"][..],
                ),
                ("solved", &["terms", "parameters"][..]),
            ] {
                let value = object(required(request, layer)?)?;
                for field in fields {
                    required(value, field)?;
                }
            }
            Ok(json!({"verified": true}))
        }
        "reconstruct_rectangular_hopping" => {
            let basis = encoded_matrix_from_json(&context, required(request, "basis_matrix")?)?;
            let parameters = elements_from_json(&context, required(request, "parameters")?)?;
            let rows = size_field(request, "rows")?;
            let columns = size_field(request, "columns")?;
            Ok(json!({
                "matrix": matrix_to_json(&reconstruct_rectangular_hopping(
                    &basis,
                    &parameters,
                    rows,
                    columns,
                )?)
            }))
        }
        "propagate_hopping" => {
            let hopping = encoded_matrix_from_json(&context, required(request, "hopping")?)?;
            let left = encoded_matrix_from_json(&context, required(request, "left")?)?;
            let right = encoded_matrix_from_json(&context, required(request, "right")?)?;
            let antiunitary = bool_field(request, "antiunitary")?;
            Ok(json!({
                "matrix": matrix_to_json(&propagate_hopping(
                    &hopping,
                    &left,
                    &right,
                    antiunitary,
                )?)
            }))
        }
        "assemble_bloch_hamiltonian" => {
            let block_dimensions: Vec<usize> =
                from_value(required(request, "block_dimensions")?, "block_dimensions")?;
            let terms = array(required(request, "terms")?)?
                .iter()
                .map(|term| term_from_json(&context, term))
                .collect::<ExactResult<Vec<_>>>()?;
            let phases = elements_from_json(&context, required(request, "phase_generators")?)?;
            Ok(json!({
                "matrix": matrix_to_json(&assemble_bloch_hamiltonian(
                    &block_dimensions,
                    &terms,
                    &phases,
                )?)
            }))
        }
        operation => Err(ExactError::new(
            "UnsupportedOperation",
            format!("unsupported tight-binding operation: {operation}"),
        )),
    }
}

pub fn compute(payload: &str) -> ExactResult<String> {
    let value: Value = serde_json::from_str(payload)
        .map_err(|error| malformed(format!("invalid JSON: {error}")))?;
    let result = evaluate(object(&value)?)?;
    serde_json::to_string(&result)
        .map_err(|error| malformed(format!("result serialization failed: {error}")))
}
