use std::sync::Arc;

use cyclotomic_nullspace::fixture::{
    array, context_from_json, element_to_json, encoded_element_from_json, encoded_matrix_from_json,
    matrix_to_json, object, required, string_field,
};
use cyclotomic_nullspace::{Element, ExactError, ExactMatrix, ExactResult};
use magnetictb_abstract_group::{GroupAction, GroupAlgebra};
use magnetictb_crystal_geometry::{SpatialOperation, compile_site_permutations};
use magnetictb_linear_algebra::exact_matrix_equal;
use magnetictb_representation::{
    InducedOrbitSpec, RepresentationError, compile_direct_product, compile_induced,
    ordered_permutation_group,
};
use serde::de::DeserializeOwned;
use serde_json::{Value, json};

fn malformed(detail: impl Into<String>) -> ExactError {
    ExactError::new("MalformedSerialization", detail)
}

fn representation_error(error: &RepresentationError) -> ExactError {
    ExactError::new(error.tag(), error.to_string())
}

fn from_value<T: DeserializeOwned>(value: &Value, field: &str) -> ExactResult<T> {
    serde_json::from_value(value.clone())
        .map_err(|error| malformed(format!("invalid {field}: {error}")))
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

fn matrix_grid_from_json(
    context: &Arc<cyclotomic_nullspace::CyclotomicContext>,
    value: &Value,
) -> ExactResult<Vec<Vec<ExactMatrix>>> {
    array(value)?
        .iter()
        .map(|matrices| matrices_from_json(context, matrices))
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

fn site_orbits_from_json(
    context: &Arc<cyclotomic_nullspace::CyclotomicContext>,
    value: &Value,
) -> ExactResult<Vec<Vec<Vec<Element>>>> {
    array(value)?
        .iter()
        .map(|orbit| {
            array(orbit)?
                .iter()
                .map(|site| elements_from_json(context, site))
                .collect()
        })
        .collect()
}

fn translations_to_json(translations: &[Vec<Vec<Vec<Element>>>]) -> Value {
    Value::Array(
        translations
            .iter()
            .map(|orbit| {
                Value::Array(
                    orbit
                        .iter()
                        .map(|operation| {
                            Value::Array(
                                operation
                                    .iter()
                                    .map(|translation| {
                                        Value::Array(
                                            translation.iter().map(element_to_json).collect(),
                                        )
                                    })
                                    .collect(),
                            )
                        })
                        .collect(),
                )
            })
            .collect(),
    )
}

fn matrix_grid_to_json(grid: &[Vec<ExactMatrix>]) -> Value {
    Value::Array(
        grid.iter()
            .map(|matrices| Value::Array(matrices.iter().map(matrix_to_json).collect()))
            .collect(),
    )
}

fn matrix_cube_to_json(cube: &[Vec<Vec<ExactMatrix>>]) -> Value {
    Value::Array(cube.iter().map(|grid| matrix_grid_to_json(grid)).collect())
}

fn group_and_actions(
    request: &serde_json::Map<String, Value>,
) -> ExactResult<(GroupAlgebra, Vec<GroupAction>)> {
    let table: Vec<Vec<usize>> = from_value(required(request, "multiplication_table")?, "table")?;
    let group = GroupAlgebra::new(table)
        .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
    let action_tables: Vec<Vec<Vec<usize>>> =
        from_value(required(request, "action_tables")?, "action_tables")?;
    let actions = action_tables
        .into_iter()
        .map(|table| {
            GroupAction::compile(&group, table)
                .map_err(|error| ExactError::new(error.tag(), error.to_string()))
        })
        .collect::<ExactResult<Vec<_>>>()?;
    Ok((group, actions))
}

#[allow(clippy::too_many_lines)]
fn evaluate(request: &serde_json::Map<String, Value>) -> ExactResult<Value> {
    match string_field(request, "operation")? {
        "ordered_permutation_group" => {
            let permutations: Vec<Vec<usize>> =
                from_value(required(request, "permutations")?, "permutations")?;
            let group = ordered_permutation_group(&permutations)
                .map_err(|error| representation_error(&error))?;
            Ok(json!({
                "multiplication_table": group.multiplication_table(),
                "identity": group.identity(),
                "inverse_indices": group.inverse_indices()
            }))
        }
        "compile_site_permutations" => {
            let context = context_from_json(required(request, "context")?, 128)?;
            let table: Vec<Vec<usize>> =
                from_value(required(request, "multiplication_table")?, "table")?;
            let group = GroupAlgebra::new(table)
                .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
            let site_orbits = site_orbits_from_json(&context, required(request, "site_orbits")?)?;
            let operations = array(required(request, "spatial_actions")?)?
                .iter()
                .map(|value| {
                    let record = object(value)?;
                    SpatialOperation::new(
                        encoded_matrix_from_json(&context, required(record, "rotation")?)?,
                        elements_from_json(&context, required(record, "translation")?)?,
                    )
                    .map_err(|error| ExactError::new(error.tag(), error.to_string()))
                })
                .collect::<ExactResult<Vec<_>>>()?;
            let data = compile_site_permutations(&group, site_orbits, operations)
                .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
            let translations = data.cell_translations();
            Ok(json!({
                "image_site_indices": data.image_site_indices(),
                "cell_translations": translations_to_json(&translations),
                "site_symmetry_operations": data.actions().iter().enumerate().map(
                    |(orbit, action)| (0..action.object_count()).map(|site| {
                        data.site_symmetry_operations(orbit, site)
                    }).collect::<Result<Vec<_>, _>>()
                ).collect::<Result<Vec<_>, _>>()
                .map_err(|error| ExactError::new(error.tag(), error.to_string()))?
            }))
        }
        "compile_direct_product" => {
            let context = context_from_json(required(request, "context")?, 128)?;
            let (group, actions) = group_and_actions(request)?;
            let flags: Vec<bool> = from_value(required(request, "antiunitary_flags")?, "flags")?;
            let local =
                matrix_grid_from_json(&context, required(request, "local_representations")?)?;
            let result = compile_direct_product(&group, &actions, &local, &flags)
                .map_err(|error| representation_error(&error))?;
            Ok(json!({
                "method": result.method(),
                "local_dimensions": result.local_dimensions(),
                "block_dimensions": result.block_dimensions(),
                "dimension": result.dimension(),
                "orbit_representation_matrices": matrix_grid_to_json(
                    result.orbit_representation_matrices()
                ),
                "representation_matrices": Value::Array(
                    result.representation_matrices().iter().map(matrix_to_json).collect()
                )
            }))
        }
        "compile_induced" => {
            let context = context_from_json(required(request, "context")?, 128)?;
            let (group, actions) = group_and_actions(request)?;
            let flags: Vec<bool> = from_value(required(request, "antiunitary_flags")?, "flags")?;
            let specifications = array(required(request, "specifications")?)?
                .iter()
                .map(|value| {
                    let record = object(value)?;
                    Ok(InducedOrbitSpec {
                        subgroup_indices: from_value(
                            required(record, "subgroup_indices")?,
                            "subgroup_indices",
                        )?,
                        site_symmetry_matrices: matrices_from_json(
                            &context,
                            required(record, "site_symmetry_matrices")?,
                        )?,
                        coset_representatives: from_value(
                            required(record, "coset_representatives")?,
                            "coset_representatives",
                        )?,
                    })
                })
                .collect::<ExactResult<Vec<_>>>()?;
            let result = compile_induced(&group, &actions, &specifications, &flags)
                .map_err(|error| representation_error(&error))?;
            let mathematical = &result.representation;
            Ok(json!({
                "method": mathematical.method(),
                "local_dimensions": mathematical.local_dimensions(),
                "block_dimensions": mathematical.block_dimensions(),
                "dimension": mathematical.dimension(),
                "local_blocks": matrix_cube_to_json(mathematical.local_blocks()),
                "orbit_representation_matrices": matrix_grid_to_json(
                    mathematical.orbit_representation_matrices()
                ),
                "representation_matrices": Value::Array(
                    mathematical.representation_matrices().iter().map(matrix_to_json).collect()
                ),
                "site_symmetry_data": result.site_symmetry_data.iter().map(|orbit| json!({
                    "subgroup_indices": orbit.subgroup_indices,
                    "coset_representatives": orbit.coset_representatives,
                    "schreier_records": orbit.schreier_records,
                    "image_site_indices": orbit.image_site_indices,
                    "local_blocks": matrix_grid_to_json(&orbit.local_blocks),
                    "local_dimension": orbit.local_dimension
                })).collect::<Vec<_>>()
            }))
        }
        "encoded_matrices_equal" => {
            let context = context_from_json(required(request, "context")?, 128)?;
            let left = array(required(request, "left")?)?;
            let right = array(required(request, "right")?)?;
            if left.len() != right.len() {
                return Ok(json!({"equal": false}));
            }
            let equal = left
                .iter()
                .zip(right)
                .try_fold(true, |equal, (left, right)| {
                    Ok::<bool, ExactError>(
                        equal
                            && exact_matrix_equal(
                                &encoded_matrix_from_json(&context, left)?,
                                &encoded_matrix_from_json(&context, right)?,
                            ),
                    )
                })?;
            Ok(json!({"equal": equal}))
        }
        "representation_site_action" => {
            let (_, actions) = group_and_actions(request)?;
            let orbit: usize = from_value(required(request, "orbit")?, "orbit")?;
            let source: usize = from_value(required(request, "source")?, "source")?;
            let operation: usize =
                from_value(required(request, "group_operation")?, "group_operation")?;
            let action = actions.get(orbit).ok_or_else(|| {
                ExactError::new(
                    "SiteActionIndexOutOfRange",
                    "site-action orbit index is outside the compiled orbit list",
                )
            })?;
            Ok(json!({
                "target": action.image(operation, source).map_err(|error| {
                    ExactError::new(error.tag(), error.to_string())
                })?
            }))
        }
        operation => Err(ExactError::new(
            "UnsupportedOperation",
            format!("unsupported representation operation: {operation}"),
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
