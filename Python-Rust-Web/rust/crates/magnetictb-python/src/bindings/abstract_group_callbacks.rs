fn python_same_test<'py>(
    left: &Bound<'py, PyAny>,
    right: &Bound<'py, PyAny>,
    same_test: Option<&Bound<'py, PyAny>>,
    error: &RefCell<Option<String>>,
) -> bool {
    if error.borrow().is_some() {
        return false;
    }
    let comparison = if let Some(same_test) = same_test {
        same_test
            .call1((left, right))
            .and_then(|value| value.is_truthy())
    } else {
        left.eq(right)
    };
    match comparison {
        Ok(value) => value,
        Err(value) => {
            *error.borrow_mut() = Some(value.to_string());
            false
        }
    }
}
#[pyfunction(signature = (generators, identity, multiply, same_test=None))]
#[allow(clippy::needless_pass_by_value)]
fn generate_group_objects(
    generators: Vec<Bound<'_, PyAny>>,
    identity: Bound<'_, PyAny>,
    multiply: Bound<'_, PyAny>,
    same_test: Option<Bound<'_, PyAny>>,
) -> PyResult<Vec<Py<PyAny>>> {
    let same_error = RefCell::new(None);
    let result = generate_group_by(
        &generators,
        &identity,
        |left, right| {
            multiply
                .call1((left, right))
                .map_err(|error| error.to_string())
        },
        |left, right| python_same_test(left, right, same_test.as_ref(), &same_error),
    )
    .map_err(|error| python_group_error(&error))?;
    if let Some(error) = same_error.into_inner() {
        return Err(python_group_error(&CoreGroupError::new(
            "same_test_failed",
            "generate_group",
            error,
        )));
    }
    Ok(result.into_iter().map(Bound::unbind).collect())
}

#[pyfunction(signature = (elements, identity, multiply, same_test=None))]
#[allow(clippy::needless_pass_by_value)]
fn find_generator_objects(
    elements: Vec<Bound<'_, PyAny>>,
    identity: Bound<'_, PyAny>,
    multiply: Bound<'_, PyAny>,
    same_test: Option<Bound<'_, PyAny>>,
) -> PyResult<Vec<Py<PyAny>>> {
    let same_error = RefCell::new(None);
    let result = find_concrete_generators_by(
        &elements,
        &identity,
        |left, right| {
            multiply
                .call1((left, right))
                .map_err(|error| error.to_string())
        },
        |left, right| python_same_test(left, right, same_test.as_ref(), &same_error),
    )
    .map_err(|error| python_group_error(&error))?;
    if let Some(error) = same_error.into_inner() {
        return Err(python_group_error(&CoreGroupError::new(
            "same_test_failed",
            "find_concrete_generators",
            error,
        )));
    }
    Ok(result.into_iter().map(Bound::unbind).collect())
}
