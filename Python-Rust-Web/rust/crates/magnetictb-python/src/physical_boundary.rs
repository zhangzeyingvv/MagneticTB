use std::sync::Arc;

use cyclotomic_nullspace::fixture::{
    array, context_from_json, element_from_json, element_to_json, matrix_to_json, object, required,
};
use cyclotomic_nullspace::{Element, ExactError, ExactMatrix, ExactResult};
use magnetictb_crystal_geometry::compile_site_permutations;
use magnetictb_data::DataCatalog;
use magnetictb_representation::compile_direct_product;
use magnetictb_symmetry::compile_msg_group;
use serde_json::{Value, json};

fn malformed(detail: impl Into<String>) -> ExactError {
    ExactError::new("MalformedSerialization", detail)
}

fn elements_from_json(
    context: &Arc<cyclotomic_nullspace::CyclotomicContext>,
    value: &Value,
) -> ExactResult<Vec<Element>> {
    array(value)?
        .iter()
        .map(|element| element_from_json(context, element))
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

/// Compiles the scalar `DirectProduct` physical representation of one Data-owned MSG.
///
/// All group, geometry, and representation work remains in Rust.  The caller only
/// supplies the exact site orbit selected by the model fixture.
pub fn compute_with_catalog(payload: &str, catalog: &DataCatalog) -> ExactResult<String> {
    let value: Value = serde_json::from_str(payload)
        .map_err(|error| malformed(format!("invalid JSON: {error}")))?;
    let request = object(&value)?;
    let context = context_from_json(required(request, "context")?, 128)?;
    let msg_id = required(request, "msg_id")?
        .as_str()
        .ok_or_else(|| malformed("msg_id must be a string"))?;
    let source = catalog
        .msg(msg_id)
        .ok_or_else(|| ExactError::new("UnknownStableId", format!("unknown MSG ID {msg_id}")))?;
    let symmetry = compile_msg_group(source)
        .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
    let site_orbits = site_orbits_from_json(&context, required(request, "site_orbits")?)?;
    let spatial_actions = symmetry
        .operations()
        .iter()
        .map(|operation| operation.spatial().clone())
        .collect();
    let sites = compile_site_permutations(symmetry.group(), site_orbits, spatial_actions)
        .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
    let identity = ExactMatrix::identity(&context, 1)?;
    let local_representations = sites
        .actions()
        .iter()
        .map(|_| vec![identity.clone(); symmetry.group().order()])
        .collect::<Vec<_>>();
    let flags = symmetry.antiunitary_flags();
    let representation = compile_direct_product(
        symmetry.group(),
        sites.actions(),
        &local_representations,
        &flags,
    )
    .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;

    let result = json!({
        "msg_id": msg_id,
        "multiplication_table": symmetry.group().multiplication_table(),
        "antiunitary_flags": flags,
        "image_site_indices": sites.image_site_indices(),
        "cell_translations": translations_to_json(&sites.cell_translations()),
        "method": representation.method(),
        "local_dimensions": representation.local_dimensions(),
        "block_dimensions": representation.block_dimensions(),
        "dimension": representation.dimension(),
        "representation_matrices": representation
            .representation_matrices()
            .iter()
            .map(matrix_to_json)
            .collect::<Vec<_>>(),
        "orbit_representation_matrices": representation
            .orbit_representation_matrices()
            .iter()
            .map(|matrices| matrices.iter().map(matrix_to_json).collect::<Vec<_>>())
            .collect::<Vec<_>>(),
        "verified": true
    });
    serde_json::to_string(&result)
        .map_err(|error| malformed(format!("result serialization failed: {error}")))
}
