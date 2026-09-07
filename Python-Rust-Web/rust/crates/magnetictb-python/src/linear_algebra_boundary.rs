use cyclotomic_nullspace::fixture::{
    array, bool_field, common_result_to_json, context_from_json, element_from_json,
    element_to_json, encoded_element_from_json, encoded_matrix_from_json, matrix_from_json,
    matrix_to_json, object, required, size_field, string_field,
};
use cyclotomic_nullspace::{
    CompiledProblem, Element, ExactError, ExactMatrix, ExactResult, common_kernel, rref,
};
use magnetictb_linear_algebra::{
    ConstraintKernelMethod, ConstraintKernelResult, ConstraintTarget, ConstraintValidationLevel,
    RectangularConstraintOperation, action_matrix_from_images, assemble_constraint_matrix,
    basis_gram_matrix, compile_rectangular_constraint_blocks, coordinates_in_basis,
    dagger_coordinate_matrix, exact_hermitian_matrix, exact_square_matrix, exact_unitary_matrix,
    hermitian_matrix_basis, integer_vector, invariant_basis, matrix_to_real_coordinates,
    real_coordinates_to_matrix, real_hs_inner_product, rectangular_generator_constraint_matrix,
    rectangular_matrix_units, rectangular_real_basis, rectangular_transport_matrix,
    solve_basis_stabilizer_actions, solve_constraint_kernel_with_options,
    stabilizer_constraint_matrix,
};
use serde_json::{Map, Value, json};
use std::sync::Arc;

fn malformed(detail: impl Into<String>) -> ExactError {
    ExactError::new("MalformedSerialization", detail)
}

fn matrices_from_json(
    context: &Arc<cyclotomic_nullspace::CyclotomicContext>,
    value: &Value,
) -> ExactResult<Vec<ExactMatrix>> {
    array(value)?
        .iter()
        .map(|matrix| encoded_matrix_from_json(context, matrix))
        .collect()
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

fn matrices_to_json(matrices: &[ExactMatrix]) -> Value {
    Value::Array(matrices.iter().map(matrix_to_json).collect())
}

fn stable_compile_error(error: ExactError) -> ExactError {
    if error.tag() == "MalformedSerialization"
        && error
            .to_string()
            .contains("exact bare JSON number must use integer syntax")
    {
        return ExactError::new("NonExactInput", error.to_string());
    }
    if error.tag() == "UnsupportedExactExpression" {
        return ExactError::new("UnsupportedScalarInput", error.to_string());
    }
    error
}

fn elements_to_json(elements: &[Element]) -> Value {
    Value::Array(elements.iter().map(element_to_json).collect())
}

fn dimensions(request: &Map<String, Value>) -> ExactResult<(usize, usize)> {
    let rows = size_field(request, "rows")?;
    let columns = size_field(request, "columns")?;
    if rows == 0 || columns == 0 {
        return Err(malformed("rows and columns must be positive"));
    }
    Ok((rows, columns))
}

fn kernel_to_json(result: &ConstraintKernelResult) -> Value {
    let method = match result.method {
        ConstraintKernelMethod::Iterative => "Iterative",
        ConstraintKernelMethod::Stacked => "Stacked",
        ConstraintKernelMethod::Cyclotomic => "Cyclotomic",
    };
    let validation = match result.validation_level {
        ConstraintValidationLevel::None => "None",
        ConstraintValidationLevel::Basic => "Basic",
        ConstraintValidationLevel::Full => "Full",
    };
    let verified = match result.validation_level {
        ConstraintValidationLevel::None => None,
        ConstraintValidationLevel::Basic => Some(result.exact_residual_verified),
        ConstraintValidationLevel::Full => Some(
            result.exact_residual_verified
                && result.independent_verified
                && result.rank_nullity_verified,
        ),
    };
    json!({
        "coordinate_dimension": result.coordinate_dimension,
        "iteration_nullities": result.iteration_nullities,
        "constraint_blocks": matrices_to_json(&result.constraint_blocks),
        "constraint_matrix": matrix_to_json(&result.constraint_matrix),
        "nullspace_rows": matrix_to_json(&result.nullspace_rows),
        "basis_matrix": matrix_to_json(&result.basis_matrix),
        "rank": result.rank,
        "nullity": result.nullity,
        "method": method,
        "validation_level": validation,
        "residual": result.residual_computed.then(|| matrix_to_json(&result.residual)),
        "exact_residual_verified": result.residual_computed
            .then_some(result.exact_residual_verified),
        "independent_verified": result.independent_checked
            .then_some(result.independent_verified),
        "rank_nullity_verified": result.rank_nullity_checked
            .then_some(result.rank_nullity_verified),
        "verified": verified
    })
}

fn operation_records(
    context: &Arc<cyclotomic_nullspace::CyclotomicContext>,
    value: &Value,
) -> ExactResult<Vec<RectangularConstraintOperation>> {
    array(value)?
        .iter()
        .map(|value| {
            let record = object(value)?;
            let target = match string_field(record, "target")? {
                "Same" => ConstraintTarget::Same,
                "Reverse" => ConstraintTarget::Reverse,
                _ => return Err(malformed("constraint target must be Same or Reverse")),
            };
            Ok(RectangularConstraintOperation {
                left: encoded_matrix_from_json(context, required(record, "left")?)?,
                right: encoded_matrix_from_json(context, required(record, "right")?)?,
                antiunitary: bool_field(record, "antiunitary")?,
                target,
            })
        })
        .collect()
}

fn optional_continuous_pair(
    context: &Arc<cyclotomic_nullspace::CyclotomicContext>,
    request: &Map<String, Value>,
) -> ExactResult<Option<(ExactMatrix, ExactMatrix)>> {
    let Some(value) = request.get("continuous_pair") else {
        return Ok(None);
    };
    if value.is_null() {
        return Ok(None);
    }
    let pair = matrices_from_json(context, value)?;
    if pair.len() != 2 {
        return Err(malformed(
            "continuous_pair must contain exactly two matrices",
        ));
    }
    Ok(Some((pair[0].clone(), pair[1].clone())))
}

#[allow(clippy::too_many_lines)]
fn evaluate(request: &Map<String, Value>) -> ExactResult<Value> {
    let context = context_from_json(required(request, "context")?, 128)?;
    let operation = string_field(request, "operation")?;
    match operation {
        "compile_exact_matrices" => {
            let matrices = matrices_from_json(&context, required(request, "matrices")?)
                .map_err(stable_compile_error)?;
            if let Some(value) = request.get("coordinate_dimension") {
                let coordinate_dimension = value
                    .as_u64()
                    .and_then(|value| usize::try_from(value).ok())
                    .ok_or_else(|| {
                        malformed("coordinate_dimension must be a nonnegative integer")
                    })?;
                if matrices
                    .iter()
                    .any(|matrix| matrix.columns() != coordinate_dimension)
                {
                    return Err(ExactError::new(
                        "DimensionMismatch",
                        "constraint columns do not match coordinate dimension",
                    ));
                }
            }
            Ok(json!({"matrices": matrices_to_json(&matrices)}))
        }
        "common_null_space" => {
            let encoded = array(required(request, "matrices")?)?;
            let coordinate_dimension = if let Some(value) = request.get("coordinate_dimension") {
                value
                    .as_u64()
                    .and_then(|value| usize::try_from(value).ok())
                    .ok_or_else(|| {
                        malformed("coordinate_dimension must be a nonnegative integer")
                    })?
            } else if let Some(first) = encoded.first() {
                let record = object(first)?;
                size_field(record, "columns")?
            } else {
                return Err(ExactError::new(
                    "CoordinateDimensionRequired",
                    "empty constraints require an explicit coordinate dimension",
                ));
            };
            let matrices = matrices_from_json(&context, required(request, "matrices")?)
                .map_err(stable_compile_error)?;
            if matrices
                .iter()
                .any(|matrix| matrix.columns() != coordinate_dimension)
            {
                return Err(ExactError::new(
                    "DimensionMismatch",
                    "constraint columns do not match coordinate dimension",
                ));
            }
            let result = common_kernel(&CompiledProblem {
                context: Arc::clone(&context),
                coordinate_dimension,
                constraints: matrices,
            })?;
            Ok(common_result_to_json(&result))
        }
        "canonical_exact_element" => {
            let element = element_from_json(&context, required(request, "element")?)?;
            Ok(json!({"element": element_to_json(&element)}))
        }
        "rref" => {
            let matrix = matrix_from_json(&context, required(request, "matrix")?)?;
            let result = rref(&matrix)?;
            Ok(json!({
                "reduced_matrix": matrix_to_json(&result.reduced_matrix),
                "pivot_columns": result.pivot_columns,
                "free_columns": result.free_columns,
                "rank": result.rank,
                "nullity": result.nullity
            }))
        }
        "matrix_predicates" => {
            let matrix = encoded_matrix_from_json(&context, required(request, "matrix")?)?;
            Ok(json!({
                "square": exact_square_matrix(&matrix),
                "zero": matrix.is_zero(),
                "hermitian": exact_hermitian_matrix(&matrix)?,
                "unitary": exact_unitary_matrix(&matrix)?
            }))
        }
        "integer_vector" => {
            let vector = elements_from_json(&context, required(request, "vector")?)?;
            Ok(json!({"integer": integer_vector(&vector)?}))
        }
        "rectangular_matrix_units" => {
            let (rows, columns) = dimensions(request)?;
            Ok(json!({"basis": matrices_to_json(&rectangular_matrix_units(
                &context, rows, columns
            )?)}))
        }
        "rectangular_real_basis" => {
            let (rows, columns) = dimensions(request)?;
            Ok(json!({"basis": matrices_to_json(&rectangular_real_basis(
                &context, rows, columns
            )?)}))
        }
        "hermitian_matrix_basis" => {
            let dimension = size_field(request, "dimension")?;
            Ok(json!({"basis": matrices_to_json(&hermitian_matrix_basis(
                &context, dimension
            )?)}))
        }
        "matrix_to_real_coordinates" => {
            let matrix = encoded_matrix_from_json(&context, required(request, "matrix")?)?;
            Ok(json!({"coordinates": elements_to_json(
                &matrix_to_real_coordinates(&matrix)?
            )}))
        }
        "real_coordinates_to_matrix" => {
            let (rows, columns) = dimensions(request)?;
            let coordinates = elements_from_json(&context, required(request, "coordinates")?)?;
            Ok(json!({"matrix": matrix_to_json(&real_coordinates_to_matrix(
                &context, &coordinates, rows, columns
            )?)}))
        }
        "real_hs_inner_product" => {
            let left = encoded_matrix_from_json(&context, required(request, "left")?)?;
            let right = encoded_matrix_from_json(&context, required(request, "right")?)?;
            Ok(json!({"value": element_to_json(&real_hs_inner_product(
                &left, &right
            )?)}))
        }
        "basis_gram_matrix" => {
            let basis = matrices_from_json(&context, required(request, "basis")?)?;
            Ok(json!({"matrix": matrix_to_json(&basis_gram_matrix(&basis)?)}))
        }
        "coordinates_in_basis" => {
            let matrix = encoded_matrix_from_json(&context, required(request, "matrix")?)?;
            let basis = matrices_from_json(&context, required(request, "basis")?)?;
            Ok(
                json!({"coordinates": elements_to_json(&coordinates_in_basis(
                &matrix, &basis
            )?)}),
            )
        }
        "action_matrix_in_basis" => {
            let basis = matrices_from_json(&context, required(request, "basis")?)?;
            let images = matrices_from_json(&context, required(request, "images")?)?;
            Ok(json!({"matrix": matrix_to_json(&action_matrix_from_images(
                &basis, &images
            )?)}))
        }
        "dagger_coordinate_matrix" => {
            let (rows, columns) = dimensions(request)?;
            Ok(json!({"matrix": matrix_to_json(&dagger_coordinate_matrix(
                &context, rows, columns
            )?)}))
        }
        "rectangular_transport_matrix" => {
            let left = encoded_matrix_from_json(&context, required(request, "left")?)?;
            let right = encoded_matrix_from_json(&context, required(request, "right")?)?;
            let antiunitary = bool_field(request, "antiunitary")?;
            Ok(
                json!({"matrix": matrix_to_json(&rectangular_transport_matrix(
                &left, &right, antiunitary
            )?)}),
            )
        }
        "rectangular_generator_constraint_matrix" => {
            let left = encoded_matrix_from_json(&context, required(request, "left")?)?;
            let right = encoded_matrix_from_json(&context, required(request, "right")?)?;
            Ok(json!({"matrix": matrix_to_json(
                &rectangular_generator_constraint_matrix(&left, &right)?
            )}))
        }
        "compile_rectangular_constraint_blocks" => {
            let (rows, columns) = dimensions(request)?;
            let operations = operation_records(&context, required(request, "operations")?)?;
            let continuous = optional_continuous_pair(&context, request)?;
            let continuous_refs = continuous.as_ref().map(|(left, right)| (left, right));
            Ok(json!({"blocks": matrices_to_json(
                &compile_rectangular_constraint_blocks(
                    rows,
                    columns,
                    &operations,
                    continuous_refs,
                )?
            )}))
        }
        "assemble_constraint_matrix" => {
            let coordinate_dimension = size_field(request, "coordinate_dimension")?;
            let blocks = matrices_from_json(&context, required(request, "blocks")?)?;
            Ok(json!({"matrix": matrix_to_json(&assemble_constraint_matrix(
                &context,
                &blocks,
                coordinate_dimension,
            )?)}))
        }
        "solve_constraint_kernel" => {
            let method = match request
                .get("method")
                .and_then(Value::as_str)
                .unwrap_or("Iterative")
            {
                "Iterative" => ConstraintKernelMethod::Iterative,
                "Stacked" => ConstraintKernelMethod::Stacked,
                "Cyclotomic" | "CyclotomicExact" => ConstraintKernelMethod::Cyclotomic,
                method => {
                    return Err(ExactError::new(
                        "InvalidKernelMethod",
                        format!("unsupported constraint-kernel method {method}"),
                    ));
                }
            };
            let validation = match request
                .get("validation_level")
                .and_then(Value::as_str)
                .unwrap_or("Basic")
            {
                "None" => ConstraintValidationLevel::None,
                "Basic" => ConstraintValidationLevel::Basic,
                "Full" => ConstraintValidationLevel::Full,
                level => {
                    return Err(ExactError::new(
                        "InvalidValidationLevel",
                        format!("unsupported validation level {level}"),
                    ));
                }
            };
            let coordinate_dimension = size_field(request, "coordinate_dimension")?;
            let blocks = matrices_from_json(&context, required(request, "blocks")?)?;
            Ok(kernel_to_json(&solve_constraint_kernel_with_options(
                &context,
                &blocks,
                coordinate_dimension,
                method,
                validation,
            )?))
        }
        "stabilizer_constraint_matrix" => {
            let actions = matrices_from_json(&context, required(request, "actions")?)?;
            Ok(
                json!({"matrix": matrix_to_json(&stabilizer_constraint_matrix(
                &actions
            )?)}),
            )
        }
        "invariant_basis" => {
            let actions = matrices_from_json(&context, required(request, "actions")?)?;
            Ok(json!({"basis_matrix": matrix_to_json(&invariant_basis(
                &actions
            )?)}))
        }
        "solve_basis_stabilizer" => {
            let basis = matrices_from_json(&context, required(request, "basis")?)?;
            let actions = matrices_from_json(&context, required(request, "actions")?)?;
            let result = solve_basis_stabilizer_actions(&basis, &actions)?;
            Ok(json!({
                "operator_basis": matrices_to_json(&result.operator_basis),
                "action_matrices": matrices_to_json(&result.action_matrices),
                "kernel": kernel_to_json(&result.kernel)
            }))
        }
        _ => Err(ExactError::new(
            "UnsupportedOperation",
            format!("unsupported exact linear-algebra operation: {operation}"),
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

#[cfg(test)]
mod tests {
    use super::compute;
    use serde_json::json;

    fn rational_context() -> serde_json::Value {
        json!({
            "conductor": 1,
            "degree": 1,
            "cyclotomic_polynomial": [
                {"numerator": "-1", "denominator": "1"},
                {"numerator": "1", "denominator": "1"}
            ]
        })
    }

    #[test]
    fn dispatcher_calls_rust_linear_algebra_and_rejects_unknown_operations() {
        let request = json!({
            "operation": "hermitian_matrix_basis",
            "context": rational_context(),
            "dimension": 1
        });
        let result: serde_json::Value =
            serde_json::from_str(&compute(&request.to_string()).expect("Hermitian basis"))
                .expect("result JSON");
        assert_eq!(result["basis"].as_array().expect("basis").len(), 1);

        let error =
            compute(&json!({"operation": "unknown", "context": rational_context()}).to_string())
                .expect_err("unknown operation must fail");
        assert_eq!(error.tag(), "UnsupportedOperation");
    }
}
