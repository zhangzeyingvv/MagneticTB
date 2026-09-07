struct ExplicitSymmetry {
    group: CompiledSeitzGroup,
    spatial: Vec<SpatialOperation>,
    antiunitary_flags: Vec<bool>,
    spin_rotations: Option<Vec<ExactMatrix>>,
    operation_labels: Vec<String>,
}
fn spin_space_operations(
    context: &Arc<CyclotomicContext>,
    values: &Value,
) -> ExactResult<Vec<SpinSpaceOperation>> {
    array(values)?
        .iter()
        .map(|value| {
            let record = object(value)?;
            let spatial_record = object(required(record, "spatial")?)?;
            SpinSpaceOperation::new(
                SpatialOperation::new(
                    encoded_matrix_from_json(context, required(spatial_record, "rotation")?)?,
                    exact_elements(context, required(spatial_record, "translation")?)?,
                )
                .map_err(|error| ExactError::new(error.tag(), error.to_string()))?,
                encoded_matrix_from_json(context, required(record, "spin_rotation")?)?,
                from_value(required(record, "antiunitary")?, "antiunitary")?,
            )
            .map_err(|error| ExactError::new(error.tag(), error.to_string()))
        })
        .collect()
}

#[allow(clippy::too_many_lines)]
fn explicit_symmetry(
    context: &Arc<CyclotomicContext>,
    request: &serde_json::Map<String, Value>,
) -> ExactResult<ExplicitSymmetry> {
    if request.contains_key("spin_space_generators") || request.contains_key("spin_space_elements")
    {
        let compiled = if let Some(raw_generators) = request.get("spin_space_generators") {
            generate_spin_space_group(&spin_space_operations(context, raw_generators)?)
        } else {
            compile_spin_space_group(spin_space_operations(
                context,
                required(request, "spin_space_elements")?,
            )?)
        }
        .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
        let spatial = compiled
            .elements()
            .iter()
            .map(|element| element.spatial().clone())
            .collect::<Vec<_>>();
        let flags = compiled
            .elements()
            .iter()
            .map(SpinSpaceOperation::antiunitary)
            .collect::<Vec<_>>();
        let spin_rotations = compiled.spin_rotations();
        let supplied_labels = request
            .get("operation_labels")
            .map(|value| from_value::<Vec<String>>(value, "operation_labels"))
            .transpose()?
            .unwrap_or_default();
        let operation_labels = if supplied_labels.len() == compiled.seitz_group().group().order() {
            supplied_labels
        } else {
            compiled
                .seitz_group()
                .operations()
                .iter()
                .map(|operation| operation.label().to_owned())
                .collect()
        };
        return Ok(ExplicitSymmetry {
            group: compiled.seitz_group().clone(),
            spatial,
            antiunitary_flags: flags,
            spin_rotations: Some(spin_rotations),
            operation_labels,
        });
    }
    if let Some(raw_generators) = request.get("seitz_generators") {
        let generators = array(raw_generators)?
            .iter()
            .enumerate()
            .map(|(index, value)| {
                let record = object(value)?;
                let spatial_record = object(required(record, "spatial")?)?;
                Ok(SeitzOperation::new(
                    format!("generator-{index}"),
                    record
                        .get("label")
                        .and_then(Value::as_str)
                        .map_or_else(|| format!("g{index}"), str::to_owned),
                    SpatialOperation::new(
                        encoded_matrix_from_json(context, required(spatial_record, "rotation")?)?,
                        exact_elements(context, required(spatial_record, "translation")?)?,
                    )
                    .map_err(|error| ExactError::new(error.tag(), error.to_string()))?,
                    from_value(required(record, "antiunitary")?, "antiunitary")?,
                ))
            })
            .collect::<ExactResult<Vec<_>>>()?;
        let group = generate_seitz_group(&generators)
            .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
        let spatial = group
            .operations()
            .iter()
            .map(|operation| operation.spatial().clone())
            .collect::<Vec<_>>();
        let antiunitary_flags = group.antiunitary_flags();
        let operation_labels = group
            .operations()
            .iter()
            .map(|operation| operation.label().to_owned())
            .collect();
        return Ok(ExplicitSymmetry {
            group,
            spatial,
            antiunitary_flags,
            spin_rotations: None,
            operation_labels,
        });
    }
    let flags: Vec<bool> =
        from_value(required(request, "antiunitary_flags")?, "antiunitary_flags")?;
    let spatial = spatial_operations(context, required(request, "spatial_actions")?)?;
    if spatial.len() != flags.len() {
        return Err(ExactError::new(
            "InvalidHamiltonianSymmetryInput",
            "spatial actions and antiunitary flags must align",
        ));
    }
    let supplied_labels = request
        .get("operation_labels")
        .map(|value| from_value::<Vec<String>>(value, "operation_labels"))
        .transpose()?
        .unwrap_or_default();
    if !supplied_labels.is_empty() && supplied_labels.len() != spatial.len() {
        return Err(malformed(
            "operation_labels must align with the complete ordered operations",
        ));
    }
    let operations = spatial
        .iter()
        .cloned()
        .zip(&flags)
        .enumerate()
        .map(|(index, (operation, &antiunitary))| {
            SeitzOperation::new(
                format!("explicit-{index}"),
                supplied_labels
                    .get(index)
                    .cloned()
                    .unwrap_or_else(|| format!("g{index}")),
                operation,
                antiunitary,
            )
        })
        .collect();
    let symmetry = if let Some(table) = request.get("multiplication_table") {
        compile_seitz_group_with_table(operations, from_value(table, "multiplication_table")?)
    } else {
        compile_seitz_group(operations)
    }
    .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
    let spin_rotations = request
        .get("spin_rotations")
        .map(|value| exact_matrices(context, value))
        .transpose()?;
    let operation_labels = symmetry
        .operations()
        .iter()
        .map(|operation| operation.label().to_owned())
        .collect();
    Ok(ExplicitSymmetry {
        group: symmetry,
        spatial,
        antiunitary_flags: flags,
        spin_rotations,
        operation_labels,
    })
}
