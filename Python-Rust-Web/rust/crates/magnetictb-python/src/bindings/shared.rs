fn python_error(error: &CoreExactError) -> PyErr {
    ExactError::new_err((error.tag().to_owned(), error.to_string()))
}
fn python_group_error(error: &CoreGroupError) -> PyErr {
    GroupError::new_err((error.tag().to_owned(), error.to_string()))
}

fn python_data_error(error: &CoreDataError) -> PyErr {
    DataError::new_err((error.tag().to_owned(), error.to_string()))
}

fn python_symmetry_error(error: &CoreSymmetryError) -> PyErr {
    SymmetryError::new_err((error.tag().to_owned(), error.to_string()))
}

fn python_properties_error(error: &CorePropertiesError) -> PyErr {
    PropertiesError::new_err((error.tag().to_owned(), error.to_string()))
}

fn python_fitting_error(error: &CoreFittingError) -> PyErr {
    FittingError::new_err((error.tag().to_owned(), error.to_string()))
}

fn serialize_data<T: Serialize>(value: &T) -> PyResult<String> {
    serde_json::to_string(value).map_err(|error| {
        python_data_error(&CoreDataError::new(
            "MalformedSerialization",
            error.to_string(),
        ))
    })
}

fn unknown_data_id(kind: &str, stable_id: &str) -> PyErr {
    python_data_error(&CoreDataError::new(
        "UnknownStableId",
        format!("unknown {kind} stable ID {stable_id}"),
    ))
}
