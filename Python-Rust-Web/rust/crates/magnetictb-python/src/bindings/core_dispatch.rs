fn parse_json(payload: &str) -> Result<Value, CoreExactError> {
    serde_json::from_str(payload)
        .map_err(|error| CoreExactError::new("MalformedSerialization", error.to_string()))
}
fn serialize_json(value: &Value) -> Result<String, CoreExactError> {
    serde_json::to_string(value)
        .map_err(|error| CoreExactError::new("MalformedSerialization", error.to_string()))
}

fn compute_null_space(payload: &str) -> Result<String, CoreExactError> {
    let request_value = parse_json(payload)?;
    let request = object(&request_value)?;
    let context = context_from_json(required(request, "context")?, 128)?;
    let matrix = matrix_from_json(&context, required(request, "matrix")?)?;
    let result = null_space(&matrix)?;
    serialize_json(&kernel_result_to_json(&result))
}

fn compute_common_kernel(payload: &str) -> Result<String, CoreExactError> {
    let request_value = parse_json(payload)?;
    let problem = problem_from_json(&request_value)?;
    let result = common_kernel(&problem)?;
    serialize_json(&common_result_to_json(&result))
}

#[pyfunction]
fn null_space_json(payload: &str) -> PyResult<String> {
    compute_null_space(payload).map_err(|error| python_error(&error))
}

#[pyfunction]
fn common_kernel_json(payload: &str) -> PyResult<String> {
    compute_common_kernel(payload).map_err(|error| python_error(&error))
}

#[pyfunction]
fn cyclotomic_context_json(conductor: usize, max_degree: usize) -> PyResult<String> {
    let context = make_context(conductor, max_degree).map_err(|error| python_error(&error))?;
    serialize_json(&context_to_json(&context)).map_err(|error| python_error(&error))
}

#[pyfunction]
fn linear_algebra_json(payload: &str) -> PyResult<String> {
    linear_algebra_boundary::compute(payload).map_err(|error| python_error(&error))
}

#[pyfunction]
fn representation_json(payload: &str) -> PyResult<String> {
    representation_boundary::compute(payload).map_err(|error| python_error(&error))
}

#[pyfunction]
fn tight_binding_json(payload: &str) -> PyResult<String> {
    tight_binding_boundary::compute(payload).map_err(|error| python_error(&error))
}

#[pyfunction]
fn geometry_json(payload: &str) -> PyResult<String> {
    geometry_boundary::compute(payload).map_err(|error| python_error(&error))
}

#[pyfunction]
fn properties_json(payload: &str) -> PyResult<String> {
    properties_boundary::compute(payload).map_err(|error| python_properties_error(&error))
}

#[pyfunction]
fn fitting_json(payload: &str) -> PyResult<String> {
    fitting_boundary::compute(payload).map_err(|error| python_fitting_error(&error))
}

#[pyfunction]
fn gapless_points_json(hamiltonian: &Bound<'_, PyAny>, payload: &str) -> PyResult<String> {
    let (occupied, options) = properties_boundary::parse_gapless_request(payload)
        .map_err(|error| python_properties_error(&error))?;
    let result = core_find_gapless_points(
        |point| {
            let value = hamiltonian.call1((point.to_vec(),)).map_err(|error| {
                CorePropertiesError::new(
                    "PythonHamiltonianFailure",
                    format!("Hamiltonian callback failed at {point:?}: {error}"),
                )
            })?;
            let rows = value
                .extract::<Vec<Vec<NumericComplex<f64>>>>()
                .map_err(|error| {
                    CorePropertiesError::new(
                        "InvalidHamiltonian",
                        format!("Hamiltonian callback must return a numeric matrix: {error}"),
                    )
                })?;
            if rows.is_empty()
                || rows[0].is_empty()
                || rows.iter().any(|row| row.len() != rows[0].len())
            {
                return Err(CorePropertiesError::new(
                    "InvalidHamiltonian",
                    "Hamiltonian callback returned an empty or ragged matrix",
                ));
            }
            let columns = rows[0].len();
            let entries = rows.into_iter().flatten().collect::<Vec<_>>();
            Ok(NumericMatrix::from_row_slice(
                entries.len() / columns,
                columns,
                &entries,
            ))
        },
        occupied,
        options,
    )
    .map_err(|error| python_properties_error(&error))?;
    properties_boundary::gapless_result_json(&result)
        .map_err(|error| python_properties_error(&error))
}
