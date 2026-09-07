fn quadratic_matrix_rows(matrix: &ExactMatrix) -> ExactResult<Vec<Vec<Value>>> {
    (0..matrix.rows())
        .map(|row| {
            (0..matrix.columns())
                .map(|column| quadratic_expression_from_element(matrix.entry(row, column)?))
                .collect()
        })
        .collect()
}
fn site_orbits_from_seeds(
    context: &Arc<CyclotomicContext>,
    wyckoff_seeds: &Value,
    spatial: &[SpatialOperation],
) -> ExactResult<Vec<Vec<Vec<Element>>>> {
    let mut result = Vec::new();
    for seed in tagged_items(wyckoff_seeds)? {
        let seed = encoded_matrix_from_json(context, seed)?;
        if seed.rows() != 2 || seed.columns() != 3 {
            return Err(ExactError::new(
                "InvalidWyckoffSeed",
                "each Wyckoff seed must be an exact 2 by 3 matrix",
            ));
        }
        let coordinate = (0..3)
            .map(|column| seed.entry(0, column).cloned())
            .collect::<ExactResult<Vec<_>>>()?;
        let mut orbit = Vec::new();
        for operation in spatial {
            let mut image = Vec::with_capacity(3);
            for row in 0..3 {
                let mut entry = operation.translation()[row].clone();
                for (column, coordinate) in coordinate.iter().enumerate() {
                    entry = entry.add(
                        &operation
                            .rotation()
                            .entry(row, column)?
                            .multiply(coordinate)?,
                    )?;
                }
                image.push(
                    fractional_part(&entry)
                        .map_err(|error| ExactError::new(error.tag(), error.to_string()))?,
                );
            }
            if !orbit.contains(&image) {
                orbit.push(image);
            }
        }
        if orbit.is_empty() {
            return Err(ExactError::new(
                "InvalidWyckoffSeed",
                "a Wyckoff seed generated an empty orbit",
            ));
        }
        result.push(orbit);
    }
    if result.is_empty() {
        return Err(ExactError::new(
            "InvalidWyckoffSeed",
            "at least one Wyckoff seed is required",
        ));
    }
    Ok(result)
}

#[allow(clippy::too_many_lines)]
fn compiler_symmetry_arguments(
    compiler_input: &Value,
    output: &mut serde_json::Map<String, Value>,
) -> ExactResult<Option<String>> {
    let raw = ordered_lookup(compiler_input, "SymmetryInformation")?;
    let explicit_representation = optional_ordered_lookup(compiler_input, "RepresentationSource")
        .and_then(Value::as_str)
        == Some("Matrices");
    let generate = if explicit_representation {
        false
    } else {
        optional_ordered_lookup(compiler_input, "GenerateSymmetryGroup")
            .map_or(Ok(false), |value| {
                from_value(value, "GenerateSymmetryGroup")
            })?
    };
    if raw.get("kind").and_then(Value::as_str) == Some("list") {
        let mut operations = Vec::new();
        let mut flags = Vec::new();
        let mut labels = Vec::new();
        for value in tagged_items(raw)? {
            let fields = tagged_items(value)?;
            if fields.len() != 4 {
                return Err(malformed(
                    "a magnetic symmetry operation must have label, rotation, translation, and F/T parity",
                ));
            }
            let flag = match fields[3].as_str() {
                Some("F") => false,
                Some("T") => true,
                _ => return Err(malformed("magnetic symmetry parity must be F or T")),
            };
            operations.push(json!({
                "label": fields[0].clone(),
                "spatial": {
                    "rotation": fields[1].clone(),
                    "translation": decode_tagged_lists(&fields[2])?
                },
                "antiunitary": flag
            }));
            flags.push(flag);
            labels.push(fields[0].clone());
        }
        output.insert("operation_labels".to_owned(), Value::Array(labels));
        if generate {
            output.insert("seitz_generators".to_owned(), Value::Array(operations));
        } else {
            output.insert(
                "spatial_actions".to_owned(),
                Value::Array(
                    operations
                        .iter()
                        .map(|operation| operation["spatial"].clone())
                        .collect(),
                ),
            );
            output.insert("antiunitary_flags".to_owned(), json!(flags));
        }
        return Ok(None);
    }

    let association = object(raw)?;
    if association.get("kind").and_then(Value::as_str) != Some("ordered_association") {
        return Err(malformed(
            "SymmetryInformation must be a tagged list or ordered association",
        ));
    }
    let mut finite = Vec::new();
    let mut labels = Vec::new();
    let mut continuous_parameter = None;
    for entry in array(required(association, "entries")?)? {
        let entry = object(entry)?;
        let label = required(entry, "key")?
            .as_str()
            .ok_or_else(|| malformed("spin-space operation label must be a string"))?;
        let element = required(entry, "value")?;
        if let Some(parameter) = optional_ordered_lookup(element, "continuous") {
            let parameter = object(parameter)?;
            if string_field(parameter, "kind")? != "symbol" {
                return Err(malformed(
                    "a continuous symmetry parameter must be a symbol",
                ));
            }
            if continuous_parameter
                .replace(string_field(parameter, "name")?.to_owned())
                .is_some()
            {
                return Err(ExactError::new(
                    "InvalidContinuousSymmetry",
                    "only one continuous symmetry parameter is supported",
                ));
            }
            continue;
        }
        let space = tagged_items(ordered_lookup(element, "space")?)?;
        let spin = tagged_items(ordered_lookup(element, "spin")?)?;
        if space.len() != 2 || spin.len() != 2 {
            return Err(malformed(
                "spin-space records must contain two space and spin fields",
            ));
        }
        let antiunitary = match spin[1].as_i64() {
            Some(0) => false,
            Some(1) => true,
            _ => return Err(malformed("spin-space antiunitary parity must be 0 or 1")),
        };
        finite.push(json!({
            "label": label,
            "spatial": {
                "rotation": space[0].clone(),
                "translation": decode_tagged_lists(&space[1])?
            },
            "spin_rotation": spin[0].clone(),
            "antiunitary": antiunitary
        }));
        labels.push(Value::String(label.to_owned()));
    }
    output.insert("operation_labels".to_owned(), Value::Array(labels));
    output.insert(
        if generate {
            "spin_space_generators"
        } else {
            "spin_space_elements"
        }
        .to_owned(),
        Value::Array(finite),
    );
    Ok(continuous_parameter)
}

#[allow(clippy::too_many_lines)]
fn compiler_representation_arguments(
    compiler_input: &Value,
    site_orbits: &[Vec<Vec<Element>>],
    continuous_parameter: Option<&str>,
    output: &mut serde_json::Map<String, Value>,
) -> ExactResult<()> {
    let source = optional_ordered_lookup(compiler_input, "RepresentationSource")
        .and_then(Value::as_str)
        .unwrap_or("BasisFunctions");
    if source == "BasisFunctions" {
        output.insert(
            "basis_functions".to_owned(),
            decode_tagged_lists(ordered_lookup(compiler_input, "BasisFunctions")?)?,
        );
        let mode = optional_ordered_lookup(compiler_input, "RepresentationMode")
            .and_then(Value::as_str)
            .unwrap_or("DirectProduct");
        if mode == "Induced" {
            output.insert(
                "representation_mode".to_owned(),
                Value::String("Induced".to_owned()),
            );
            let local_data = ordered_lookup(compiler_input, "SiteLocalData")?;
            let automatic = local_data.as_object().is_some_and(|record| {
                record.get("kind").and_then(Value::as_str) == Some("symbol")
                    && record.get("name").and_then(Value::as_str) == Some("System`Automatic")
            });
            if !automatic {
                let specifications = tagged_items(local_data)?
                    .iter()
                    .map(|specification| {
                        let mut encoded = serde_json::Map::new();
                        encoded.insert(
                            "reference_site_index".to_owned(),
                            ordered_lookup(specification, "ReferenceSiteIndex")?.clone(),
                        );
                        encoded.insert(
                            "site_symmetry_operation_indices".to_owned(),
                            decode_tagged_lists(ordered_lookup(
                                specification,
                                "SiteSymmetryOperationIndices",
                            )?)?,
                        );
                        encoded.insert(
                            "site_symmetry_matrices".to_owned(),
                            decode_tagged_lists(ordered_lookup(
                                specification,
                                "SiteSymmetryMatrices",
                            )?)?,
                        );
                        if let Some(representatives) =
                            optional_ordered_lookup(specification, "CosetRepresentativeIndices")
                        {
                            encoded.insert(
                                "coset_representative_indices".to_owned(),
                                decode_tagged_lists(representatives)?,
                            );
                        }
                        Ok(Value::Object(encoded))
                    })
                    .collect::<ExactResult<Vec<_>>>()?;
                output.insert(
                    "induced_specifications".to_owned(),
                    Value::Array(specifications),
                );
            }
        } else if mode != "DirectProduct" {
            return Err(ExactError::new(
                "InvalidRepresentationMode",
                format!("unsupported representation mode {mode}"),
            ));
        }
        return Ok(());
    }
    if source != "Matrices" {
        return Err(ExactError::new(
            "InvalidRepresentationSource",
            format!("unsupported representation source {source}"),
        ));
    }
    let mode = optional_ordered_lookup(compiler_input, "RepresentationMode")
        .and_then(Value::as_str)
        .unwrap_or("DirectProduct");
    let raw = ordered_lookup(compiler_input, "RepresentationInformation")?;
    if mode == "Induced" {
        output.insert(
            "representation_mode".to_owned(),
            Value::String("Induced".to_owned()),
        );
        let discrete = optional_ordered_lookup(raw, "Discrete").unwrap_or(raw);
        let specifications = tagged_items(discrete)?
            .iter()
            .map(|specification| {
                Ok(json!({
                    "reference_site_index": ordered_lookup(specification, "ReferenceSiteIndex")?,
                    "site_symmetry_operation_indices": decode_tagged_lists(
                        ordered_lookup(specification, "SiteSymmetryOperationIndices")?
                    )?,
                    "site_symmetry_matrices": decode_tagged_lists(
                        ordered_lookup(specification, "SiteSymmetryMatrices")?
                    )?
                }))
            })
            .collect::<ExactResult<Vec<_>>>()?;
        output.insert(
            "induced_specifications".to_owned(),
            Value::Array(specifications),
        );
        if let Some(continuous) = optional_ordered_lookup(raw, "Continuous") {
            output.insert(
                "continuous_site_matrices".to_owned(),
                decode_tagged_lists(ordered_lookup(continuous, "SiteMatrices")?)?,
            );
            output.insert(
                "continuous_parameter".to_owned(),
                Value::String(
                    continuous_parameter
                        .ok_or_else(|| {
                            ExactError::new(
                                "InvalidContinuousSymmetry",
                                "continuous representation has no matching symmetry parameter",
                            )
                        })?
                        .to_owned(),
                ),
            );
        }
        return Ok(());
    }
    if mode != "DirectProduct" {
        return Err(ExactError::new(
            "InvalidRepresentationMode",
            format!("unsupported representation mode {mode}"),
        ));
    }
    let discrete = optional_ordered_lookup(raw, "Discrete").unwrap_or(raw);
    let discrete = decode_tagged_lists(discrete)?;
    let orbit_matrices = array(&discrete)?;
    if orbit_matrices.len() != site_orbits.len() {
        return Err(ExactError::new(
            "InvalidDirectProductData",
            "one discrete matrix list is required for every site orbit",
        ));
    }
    let mut blocks = Vec::new();
    let mut dimensions = Vec::new();
    for (matrices, orbit) in orbit_matrices.iter().zip(site_orbits) {
        let matrices = array(matrices)?;
        let first = matrices
            .first()
            .ok_or_else(|| malformed("discrete matrix list must be nonempty"))?;
        let rows: usize = from_value(required(object(first)?, "rows")?, "matrix rows")?;
        dimensions.push(rows);
        blocks.push(Value::Array(
            matrices
                .iter()
                .map(|matrix| Value::Array(vec![matrix.clone(); orbit.len()]))
                .collect(),
        ));
    }
    output.insert("local_blocks".to_owned(), Value::Array(blocks));
    output.insert("local_dimensions".to_owned(), json!(dimensions));
    if let Some(continuous) = optional_ordered_lookup(raw, "Continuous") {
        output.insert(
            "continuous_site_matrices".to_owned(),
            decode_tagged_lists(ordered_lookup(continuous, "SiteMatrices")?)?,
        );
        output.insert(
            "continuous_parameter".to_owned(),
            Value::String(
                continuous_parameter
                    .ok_or_else(|| {
                        ExactError::new(
                            "InvalidContinuousSymmetry",
                            "continuous representation has no matching symmetry parameter",
                        )
                    })?
                    .to_owned(),
            ),
        );
    }
    Ok(())
}

fn compiler_input_model_result(request: &serde_json::Map<String, Value>) -> ExactResult<Value> {
    let target_shell: usize = from_value(required(request, "target_shell")?, "target_shell")?;
    let mut results = compiler_input_model_results(request, &[target_shell])?;
    Ok(results.remove(0))
}

fn compiler_input_all_model_results(
    request: &serde_json::Map<String, Value>,
) -> ExactResult<Value> {
    let compiler_input = required(request, "compiler_input")?;
    let initial_shells: usize = from_value(
        ordered_lookup(compiler_input, "InitialBondShells")?,
        "InitialBondShells",
    )?;
    if initial_shells == 0 {
        return Err(ExactError::new(
            "InvalidPeriodicBondInput",
            "InitialBondShells must be positive",
        ));
    }
    let targets = (1..=initial_shells).collect::<Vec<_>>();
    Ok(json!({"shell_results": compiler_input_model_results(request, &targets)?}))
}

fn compiler_input_model_results(
    request: &serde_json::Map<String, Value>,
    target_shells: &[usize],
) -> ExactResult<Vec<Value>> {
    let context = context_from_json(required(request, "context")?, 128)?;
    let compiler_input = required(request, "compiler_input")?;
    let bindings = lattice_parameter_bindings(compiler_input)?;
    let lattice = substitute_exact_symbols(ordered_lookup(compiler_input, "Lattice")?, &bindings)?;
    let ResolvedLattice {
        matrix: lattice_matrix,
        encoded: lattice,
        approximations,
    } = resolve_model_lattice(&context, &lattice)?;
    let lattice_rows = quadratic_matrix_rows(&lattice_matrix)?;

    if target_shells.is_empty() || target_shells.contains(&0) {
        return Err(ExactError::new(
            "InvalidPeriodicBondInput",
            "target_shell must be positive",
        ));
    }
    let initial_shells = optional_ordered_lookup(compiler_input, "InitialBondShells").map_or_else(
        || Ok(*target_shells.iter().max().expect("nonempty target shells")),
        |value| from_value(value, "InitialBondShells"),
    )?;
    if target_shells.iter().any(|&target| initial_shells < target) {
        return Err(ExactError::new(
            "InvalidPeriodicBondInput",
            "target_shell exceeds InitialBondShells",
        ));
    }

    let mut compiled = serde_json::Map::new();
    compiled.insert("context".to_owned(), required(request, "context")?.clone());
    propagate_hamiltonian_options(request, &mut compiled);
    let continuous_parameter = compiler_symmetry_arguments(compiler_input, &mut compiled)?;
    let symmetry = explicit_symmetry(&context, &compiled)?;
    let site_orbits = site_orbits_from_seeds(
        &context,
        ordered_lookup(compiler_input, "WyckoffPosition")?,
        &symmetry.spatial,
    )?;
    let encoded_orbits = Value::Array(
        site_orbits
            .iter()
            .map(|orbit| {
                Value::Array(
                    orbit
                        .iter()
                        .map(|site| Value::Array(site.iter().map(element_to_json).collect()))
                        .collect(),
                )
            })
            .collect(),
    );
    let flattened_sites = Value::Array(
        site_orbits
            .iter()
            .flatten()
            .map(|site| {
                Ok(Value::Array(
                    site.iter()
                        .map(encoded_rational_from_element)
                        .collect::<ExactResult<Vec<_>>>()?,
                ))
            })
            .collect::<ExactResult<Vec<_>>>()?,
    );
    compiled.insert("site_orbits".to_owned(), encoded_orbits.clone());
    compiled.insert("sites".to_owned(), flattened_sites);
    compiled.insert("lattice".to_owned(), json!(lattice_rows));
    compiled.insert("basis_lattice".to_owned(), lattice);
    compiled.insert("requested_shells".to_owned(), json!(initial_shells));
    compiler_representation_arguments(
        compiler_input,
        &site_orbits,
        continuous_parameter.as_deref(),
        &mut compiled,
    )?;
    let shell_indices = target_shells
        .iter()
        .map(|target| target - 1)
        .collect::<Vec<_>>();
    let mut results = explicit_block_model_results(&compiled, &shell_indices)?;
    for (result, &target_shell) in results.iter_mut().zip(target_shells) {
        result["compiled_from_raw_input"] = Value::Bool(true);
        result["site_orbits"] = encoded_orbits.clone();
        result["evaluated_lattice"] = matrix_to_json(&lattice_matrix);
        record_lattice_approximation(result, &approximations);
        result["target_shell"] = json!(target_shell);
    }
    Ok(results)
}

#[allow(clippy::too_many_lines)]
fn data_compiler_input_model_result(
    request: &serde_json::Map<String, Value>,
    catalog: &DataCatalog,
) -> ExactResult<Value> {
    let target_shell: usize = from_value(required(request, "target_shell")?, "target_shell")?;
    let mut results = data_compiler_input_model_results(request, catalog, &[target_shell])?;
    Ok(results.remove(0))
}

fn data_compiler_input_all_model_results(
    request: &serde_json::Map<String, Value>,
    catalog: &DataCatalog,
) -> ExactResult<Value> {
    let compiler_input = required(request, "compiler_input")?;
    let initial_shells: usize = from_value(
        ordered_lookup(compiler_input, "InitialBondShells")?,
        "InitialBondShells",
    )?;
    if initial_shells == 0 {
        return Err(ExactError::new(
            "InvalidPeriodicBondInput",
            "InitialBondShells must be positive",
        ));
    }
    let targets = (1..=initial_shells).collect::<Vec<_>>();
    Ok(json!({
        "shell_results": data_compiler_input_model_results(request, catalog, &targets)?
    }))
}

#[allow(clippy::too_many_lines)]
fn data_compiler_input_model_results(
    request: &serde_json::Map<String, Value>,
    catalog: &DataCatalog,
    target_shells: &[usize],
) -> ExactResult<Vec<Value>> {
    let context = context_from_json(required(request, "context")?, 128)?;
    let compiler_input = required(request, "compiler_input")?;
    let data_trace = required(request, "data_trace")?;
    let msg_id = ordered_lookup(data_trace, "msg_stable_id")?
        .as_str()
        .ok_or_else(|| malformed("msg_stable_id must be a string"))?;
    let selections = if let Some(value) = optional_ordered_lookup(data_trace, "wyckoff_selections")
    {
        tagged_items(value)?
            .iter()
            .map(|selection| {
                let source_ordinal = from_value(
                    ordered_lookup(selection, "source_ordinal_1based")?,
                    "source_ordinal_1based",
                )?;
                let letter = ordered_lookup(selection, "letter")?
                    .as_str()
                    .ok_or_else(|| malformed("Wyckoff selection letter must be a string"))?;
                Ok((source_ordinal, letter.to_owned()))
            })
            .collect::<ExactResult<Vec<_>>>()?
    } else {
        vec![(
            from_value(
                ordered_lookup(data_trace, "wyckoff_source_ordinal_1based")?,
                "wyckoff_source_ordinal_1based",
            )?,
            ordered_lookup(data_trace, "wyckoff_letter")?
                .as_str()
                .ok_or_else(|| malformed("wyckoff_letter must be a string"))?
                .to_owned(),
        )]
    };
    let raw_seeds = tagged_items(ordered_lookup(compiler_input, "WyckoffPosition")?)?;
    if selections.is_empty() || selections.len() != raw_seeds.len() {
        return Err(ExactError::new(
            "InvalidWyckoffSeed",
            "Data model requires one ordered Wyckoff seed per Data selection",
        ));
    }
    let mut site_data_records = Vec::with_capacity(selections.len());
    let mut site_orbits = Vec::with_capacity(selections.len());
    for ((source_ordinal, expected_letter), raw_seed) in selections.iter().zip(raw_seeds) {
        let seed = encoded_matrix_from_json(&context, raw_seed)?;
        if seed.rows() != 2 || seed.columns() != 3 {
            return Err(ExactError::new(
                "InvalidWyckoffSeed",
                "each Data model Wyckoff seed must be an exact 2 by 3 matrix",
            ));
        }
        let mut wyckoff_bindings = BTreeMap::new();
        for (column, name) in ["MagneticTB`x", "MagneticTB`y", "MagneticTB`z"]
            .into_iter()
            .enumerate()
        {
            wyckoff_bindings.insert(name.to_owned(), seed.entry(0, column)?.clone());
        }
        let site_data = compile_msg_wyckoff_sites_by_letter(
            catalog,
            msg_id,
            *source_ordinal,
            expected_letter,
            &wyckoff_bindings,
        )
        .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
        if site_data.letter() != expected_letter {
            return Err(ExactError::new(
                "WyckoffIdentityMismatch",
                format!(
                    "Data ordinal {source_ordinal} has letter {}, expected {expected_letter}",
                    site_data.letter()
                ),
            ));
        }
        site_orbits.extend(site_data.sites().site_orbits().iter().cloned());
        site_data_records.push(site_data);
    }

    let lattice_bindings = lattice_parameter_bindings(compiler_input)?;
    let lattice = substitute_exact_symbols(
        ordered_lookup(compiler_input, "Lattice")?,
        &lattice_bindings,
    )?;
    let ResolvedLattice {
        matrix: lattice_matrix,
        encoded: lattice,
        approximations,
    } = resolve_model_lattice(&context, &lattice)?;
    let lattice_rows = quadratic_matrix_rows(&lattice_matrix)?;
    if target_shells.is_empty() || target_shells.contains(&0) {
        return Err(ExactError::new(
            "InvalidPeriodicBondInput",
            "target_shell must be positive",
        ));
    }
    let initial_shells = optional_ordered_lookup(compiler_input, "InitialBondShells").map_or_else(
        || Ok(*target_shells.iter().max().expect("nonempty target shells")),
        |value| from_value(value, "InitialBondShells"),
    )?;
    if target_shells.iter().any(|&target| initial_shells < target) {
        return Err(ExactError::new(
            "InvalidPeriodicBondInput",
            "target_shell exceeds InitialBondShells",
        ));
    }

    let encoded_orbits = Value::Array(
        site_orbits
            .iter()
            .map(|orbit| {
                Value::Array(
                    orbit
                        .iter()
                        .map(|site| Value::Array(site.iter().map(element_to_json).collect()))
                        .collect(),
                )
            })
            .collect(),
    );
    let flattened_sites = Value::Array(
        site_orbits
            .iter()
            .flatten()
            .map(|site| {
                Ok(Value::Array(
                    site.iter()
                        .map(encoded_rational_from_element)
                        .collect::<ExactResult<Vec<_>>>()?,
                ))
            })
            .collect::<ExactResult<Vec<_>>>()?,
    );
    let symmetry = site_data_records[0].symmetry();
    let mut compiled = serde_json::Map::new();
    compiled.insert("context".to_owned(), required(request, "context")?.clone());
    propagate_hamiltonian_options(request, &mut compiled);
    compiled.insert(
        "multiplication_table".to_owned(),
        json!(symmetry.group().multiplication_table()),
    );
    compiled.insert(
        "antiunitary_flags".to_owned(),
        json!(symmetry.antiunitary_flags()),
    );
    compiled.insert(
        "operation_labels".to_owned(),
        Value::Array(
            symmetry
                .operations()
                .iter()
                .map(|operation| Value::String(operation.label().to_owned()))
                .collect(),
        ),
    );
    compiled.insert(
        "spatial_actions".to_owned(),
        Value::Array(
            symmetry
                .operations()
                .iter()
                .map(|operation| {
                    json!({
                        "rotation": matrix_to_json(operation.spatial().rotation()),
                        "translation": operation.spatial().translation().iter()
                            .map(element_to_json).collect::<Vec<_>>()
                    })
                })
                .collect(),
        ),
    );
    compiled.insert("site_orbits".to_owned(), encoded_orbits.clone());
    compiled.insert("sites".to_owned(), flattened_sites);
    compiled.insert("lattice".to_owned(), json!(lattice_rows));
    compiled.insert("basis_lattice".to_owned(), lattice);
    compiled.insert("requested_shells".to_owned(), json!(initial_shells));
    compiler_representation_arguments(compiler_input, &site_orbits, None, &mut compiled)?;
    let shell_indices = target_shells
        .iter()
        .map(|target| target - 1)
        .collect::<Vec<_>>();
    let mut results = explicit_block_model_results(&compiled, &shell_indices)?;
    let encoded_selections = Value::Array(
        selections
            .iter()
            .map(|(source_ordinal, letter)| {
                json!({"source_ordinal": source_ordinal, "letter": letter})
            })
            .collect(),
    );
    for (result, &target_shell) in results.iter_mut().zip(target_shells) {
        result["compiled_from_raw_input"] = Value::Bool(true);
        result["data_loaded"] = Value::Bool(true);
        result["msg_id"] = Value::String(msg_id.to_owned());
        result["wyckoff_selections"] = encoded_selections.clone();
        if selections.len() == 1 {
            result["wyckoff_source_ordinal"] = json!(selections[0].0);
            result["wyckoff_letter"] = Value::String(selections[0].1.clone());
        }
        result["site_orbits"] = encoded_orbits.clone();
        result["evaluated_lattice"] = matrix_to_json(&lattice_matrix);
        record_lattice_approximation(result, &approximations);
        result["target_shell"] = json!(target_shell);
    }
    Ok(results)
}

#[allow(clippy::too_many_lines)]
fn explicit_block_model_result(request: &serde_json::Map<String, Value>) -> ExactResult<Value> {
    let shell_index: usize = from_value(required(request, "shell_index")?, "shell_index")?;
    let mut results = explicit_block_model_results(request, &[shell_index])?;
    Ok(results.remove(0))
}

#[allow(clippy::too_many_lines)]
fn explicit_block_model_results(
    request: &serde_json::Map<String, Value>,
    shell_indices: &[usize],
) -> ExactResult<Vec<Value>> {
    let (hermitian, kernel_method, validation_kind, kernel_method_name, validation_level) =
        hamiltonian_options(request)?;
    let (requested_shells, shells) = bond_shells(request)?;
    if shell_indices.is_empty() || shell_indices.iter().any(|&index| index >= requested_shells) {
        return Err(ExactError::new(
            "InvalidPeriodicBondInput",
            "shell_index is outside requested_shells",
        ));
    }
    let context = context_from_json(required(request, "context")?, 128)?;
    let explicit = explicit_symmetry(&context, request)?;
    let symmetry = explicit.group;
    let spatial = explicit.spatial;
    let antiunitary_flags = explicit.antiunitary_flags;
    let symmetry_spin_rotations = explicit.spin_rotations;
    let operation_labels = explicit.operation_labels;
    let model_lattice = encoded_matrix_from_json(
        &context,
        request
            .get("basis_lattice")
            .unwrap_or(required(request, "lattice")?),
    )?;
    let resolved_spin_actions = if let Some(rotations) = symmetry_spin_rotations.as_ref() {
        rotations.clone()
    } else {
        spatial
            .iter()
            .map(|operation| compile_spatial_spin_action(operation.rotation(), &model_lattice))
            .collect::<Result<Vec<_>, _>>()
            .map_err(|error| ExactError::new(error.tag(), error.to_string()))?
    };
    let serialized_spin_actions = resolved_spin_actions
        .iter()
        .map(matrix_to_json)
        .collect::<Vec<_>>();
    let serialized_spin_rotations = symmetry_spin_rotations
        .as_ref()
        .map(|rotations| rotations.iter().map(matrix_to_json).collect::<Vec<_>>());
    let site_orbits = exact_site_orbits(&context, required(request, "site_orbits")?)?;
    let sites = compile_site_permutations(symmetry.group(), site_orbits.clone(), spatial.clone())
        .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
    let mut induced_reference_indices: Option<Vec<usize>> = None;
    let mut induced_coset_representatives: Option<Vec<Vec<usize>>> = None;
    let mut transported_basis_states: Option<Vec<Vec<Vec<Value>>>> = None;
    let (representation, local_blocks, local_dimensions, spinor_basis_flags) =
        if let Some(induced_value) = request.get("induced_specifications") {
            let raw_specifications = array(induced_value)?;
            if raw_specifications.len() != sites.actions().len() {
                return Err(ExactError::new(
                    "InvalidInducedData",
                    "one induced specification is required for every site orbit",
                ));
            }
            let specifications = raw_specifications
                .iter()
                .enumerate()
                .map(|(orbit_index, value)| {
                    let record = object(value)?;
                    let reference: usize = from_value(
                        required(record, "reference_site_index")?,
                        "reference_site_index",
                    )?;
                    let coset_representatives =
                        if let Some(value) = record.get("coset_representative_indices") {
                            from_value(value, "coset_representative_indices")?
                        } else {
                            automatic_coset_representatives(
                                &sites.actions()[orbit_index],
                                &antiunitary_flags,
                                reference,
                            )
                            .map_err(|error| ExactError::new(error.tag(), error.to_string()))?
                        };
                    Ok(InducedOrbitSpec {
                        subgroup_indices: from_value(
                            required(record, "site_symmetry_operation_indices")?,
                            "site_symmetry_operation_indices",
                        )?,
                        site_symmetry_matrices: exact_matrices(
                            &context,
                            required(record, "site_symmetry_matrices")?,
                        )?,
                        coset_representatives,
                    })
                })
                .collect::<ExactResult<Vec<_>>>()?;
            induced_reference_indices = Some(
                raw_specifications
                    .iter()
                    .map(|value| {
                        from_value(
                            required(object(value)?, "reference_site_index")?,
                            "reference_site_index",
                        )
                    })
                    .collect::<ExactResult<Vec<_>>>()?,
            );
            induced_coset_representatives = Some(
                specifications
                    .iter()
                    .map(|specification| specification.coset_representatives.clone())
                    .collect(),
            );
            let compiled = compile_induced(
                symmetry.group(),
                sites.actions(),
                &specifications,
                &antiunitary_flags,
            )
            .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
            let local_blocks = compiled.representation.local_blocks().to_vec();
            let local_dimensions = compiled.representation.local_dimensions().to_vec();
            (
                compiled.representation,
                local_blocks,
                local_dimensions,
                None,
            )
        } else if let Some(basis_value) = request.get("basis_functions") {
            let basis_items = array(basis_value)?;
            if basis_items.len() != site_orbits.len() {
                return Err(ExactError::new(
                    "InvalidFunctionBasis",
                    "one ordered basis list is required for every site orbit",
                ));
            }
            if request.get("representation_mode").and_then(Value::as_str) == Some("Induced") {
                let specifications = sites
                    .actions()
                    .iter()
                    .zip(basis_items)
                    .map(|(action, basis)| {
                        let coset_representatives =
                            automatic_coset_representatives(action, &antiunitary_flags, 0)
                                .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
                        let compiled = compile_induced_basis_item(
                            &context,
                            basis,
                            &spatial,
                            symmetry_spin_rotations.as_deref(),
                            &antiunitary_flags,
                            action.action_table(),
                            &coset_representatives,
                            &model_lattice,
                        )?;
                        Ok((
                            coset_representatives,
                            compiled.local_blocks().to_vec(),
                            compiled.spinor(),
                            compiled
                                .transported_basis_states()
                                .iter()
                                .map(|states| basis_states_to_json(states))
                                .collect::<ExactResult<Vec<_>>>()?,
                        ))
                    })
                    .collect::<ExactResult<Vec<_>>>()?;
                let mut coset_representatives = Vec::with_capacity(specifications.len());
                let mut local_blocks = Vec::with_capacity(specifications.len());
                let mut local_dimensions = Vec::with_capacity(specifications.len());
                let mut spinor_basis_flags = Vec::with_capacity(specifications.len());
                let mut basis_states = Vec::with_capacity(specifications.len());
                for (representatives, blocks, spinor, states) in specifications {
                    local_dimensions.push(states.first().map_or(0, Vec::len));
                    coset_representatives.push(representatives);
                    local_blocks.push(blocks);
                    spinor_basis_flags.push(spinor);
                    basis_states.push(states);
                }
                if local_dimensions.contains(&0) {
                    return Err(ExactError::new(
                        "InvalidFunctionBasis",
                        "an induced transported basis must be nonempty",
                    ));
                }
                induced_reference_indices = Some(vec![0; local_blocks.len()]);
                induced_coset_representatives = Some(coset_representatives);
                transported_basis_states = Some(basis_states);
                let representation = compile_block_monomial(
                    "Induced",
                    symmetry.group(),
                    sites.actions(),
                    local_blocks.clone(),
                    local_dimensions.clone(),
                    &antiunitary_flags,
                )
                .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
                (
                    representation,
                    local_blocks,
                    local_dimensions,
                    Some(spinor_basis_flags),
                )
            } else {
                let compiled = basis_items
                    .iter()
                    .map(|basis| {
                        compile_basis_item(
                            &context,
                            basis,
                            &spatial,
                            symmetry_spin_rotations.as_deref(),
                            &antiunitary_flags,
                            &model_lattice,
                        )
                    })
                    .collect::<ExactResult<Vec<_>>>()?;
                let dimensions = compiled
                    .iter()
                    .map(|action| action.local_matrices()[0].rows())
                    .collect::<Vec<_>>();
                let spinor_basis_flags = compiled
                    .iter()
                    .map(CompiledBasisAction::spinor)
                    .collect::<Vec<_>>();
                let blocks = compiled
                    .iter()
                    .zip(&site_orbits)
                    .map(|(action, orbit)| {
                        action
                            .local_matrices()
                            .iter()
                            .map(|matrix| vec![matrix.clone(); orbit.len()])
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>();
                let representation = compile_block_monomial(
                    "BasisCatalog",
                    symmetry.group(),
                    sites.actions(),
                    blocks.clone(),
                    dimensions.clone(),
                    &antiunitary_flags,
                )
                .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
                transported_basis_states = Some(
                    basis_items
                        .iter()
                        .zip(&site_orbits)
                        .map(|(basis, orbit)| {
                            let resolved =
                                resolve_basis_item_states(&context, basis, &model_lattice)?;
                            Ok(vec![resolved; orbit.len()])
                        })
                        .collect::<ExactResult<Vec<_>>>()?,
                );
                (representation, blocks, dimensions, Some(spinor_basis_flags))
            }
        } else {
            let blocks = exact_matrix_cube(&context, required(request, "local_blocks")?)?;
            let dimensions: Vec<usize> =
                from_value(required(request, "local_dimensions")?, "local_dimensions")?;
            let representation = compile_block_monomial(
                "ExplicitBlocks",
                symmetry.group(),
                sites.actions(),
                blocks.clone(),
                dimensions.clone(),
                &antiunitary_flags,
            )
            .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
            (representation, blocks, dimensions, None)
        };
    let global_dimensions = sites
        .site_orbits()
        .iter()
        .zip(&local_dimensions)
        .flat_map(|(orbit, &dimension)| vec![dimension; orbit.len()])
        .collect::<Vec<_>>();
    let site_actions = (0..symmetry.group().order())
        .map(|operation| {
            local_blocks
                .iter()
                .flat_map(|orbit| orbit[operation].clone())
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let continuous_generators = if let Some(matrices) = request.get("continuous_site_matrices") {
        let parameter = required(request, "continuous_parameter")?
            .as_str()
            .ok_or_else(|| malformed("continuous_parameter must be a symbol name string"))?;
        Some(continuous_site_generators(&context, matrices, parameter)?)
    } else {
        request
            .get("continuous_site_generators")
            .map(|value| exact_matrices(&context, value))
            .transpose()?
    };
    let cell_translations = sites
        .cell_translations()
        .iter()
        .map(|orbit| {
            orbit
                .iter()
                .map(|operation| {
                    operation
                        .iter()
                        .map(|translation| {
                            translation.iter().map(element_to_json).collect::<Vec<_>>()
                        })
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let mut orbital_layout = Vec::new();
    let mut global_site_index = 0;
    let induced_mode = request.get("representation_mode").and_then(Value::as_str)
        == Some("Induced")
        || request.get("induced_specifications").is_some();
    for (orbit_index, (orbit, &local_dimension)) in
        site_orbits.iter().zip(&local_dimensions).enumerate()
    {
        for (site_in_orbit, position) in orbit.iter().enumerate() {
            let cartesian_position = (0..3)
                .map(|column| {
                    position.iter().enumerate().try_fold(
                        Element::zero(&context)?,
                        |total, (row, coordinate)| {
                            total.add(&coordinate.multiply(model_lattice.entry(row, column)?)?)
                        },
                    )
                })
                .collect::<ExactResult<Vec<_>>>()?;
            for local_orbital_index in 0..local_dimension {
                let mut record = serde_json::Map::new();
                record.insert("orbital_index".to_owned(), json!(orbital_layout.len()));
                record.insert("orbit_index".to_owned(), json!(orbit_index));
                record.insert("site_in_orbit".to_owned(), json!(site_in_orbit));
                record.insert("site_index".to_owned(), json!(global_site_index));
                record.insert("local_orbital_index".to_owned(), json!(local_orbital_index));
                record.insert(
                    "fractional_position".to_owned(),
                    Value::Array(position.iter().map(element_to_json).collect()),
                );
                record.insert(
                    "cartesian_position".to_owned(),
                    Value::Array(cartesian_position.iter().map(element_to_json).collect()),
                );
                if let Some(states) = &transported_basis_states {
                    let state_record = states
                        .get(orbit_index)
                        .and_then(|sites| sites.get(site_in_orbit))
                        .and_then(|basis| basis.get(local_orbital_index))
                        .ok_or_else(|| {
                            ExactError::new(
                                "InvalidFunctionBasis",
                                "transported basis states do not align with orbital layout",
                            )
                        })?;
                    let state_record = object(state_record)?;
                    for field in [
                        "basis_state",
                        "spatial_orbital",
                        "spin_state",
                        "spin_structure",
                    ] {
                        record.insert(field.to_owned(), required(state_record, field)?.clone());
                    }
                }
                if induced_mode {
                    let reference = induced_reference_indices
                        .as_ref()
                        .and_then(|indices| indices.get(orbit_index))
                        .copied()
                        .ok_or_else(|| {
                            ExactError::new(
                                "InvalidInducedData",
                                "reference sites do not align with orbital layout",
                            )
                        })?;
                    let transport = induced_coset_representatives
                        .as_ref()
                        .and_then(|orbits| orbits.get(orbit_index))
                        .and_then(|operations| operations.get(site_in_orbit))
                        .copied()
                        .ok_or_else(|| {
                            ExactError::new(
                                "InvalidInducedData",
                                "coset representatives do not align with orbital layout",
                            )
                        })?;
                    record.insert("reference_site_index".to_owned(), json!(reference));
                    record.insert("transport_operation_index".to_owned(), json!(transport));
                    record.insert(
                        "transport_operation_label".to_owned(),
                        Value::String(operation_labels[transport].clone()),
                    );
                }
                orbital_layout.push(Value::Object(record));
            }
            global_site_index += 1;
        }
    }
    let generator_indices = symmetry
        .group()
        .find_generator_indices(None)
        .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
    let supplied_generators = request
        .get("seitz_generators")
        .or_else(|| request.get("spin_space_generators"));
    let generator_display_records = if let Some(generators) = supplied_generators {
        array(generators)?
            .iter()
            .map(|generator| {
                let generator = object(generator)?;
                Ok(json!({
                    "label": required(generator, "label")?.clone(),
                    "antiunitary": required(generator, "antiunitary")?.clone()
                }))
            })
            .collect::<ExactResult<Vec<_>>>()?
    } else {
        generator_indices
            .iter()
            .map(|&index| {
                json!({
                    "label": &operation_labels[index],
                    "antiunitary": antiunitary_flags[index]
                })
            })
            .collect()
    };
    let common = json!({
        "field_context": context_to_json(&context),
        "model_lattice": request.get("basis_lattice").unwrap_or(required(request, "lattice")?).clone(),
        "evaluated_model_lattice": matrix_to_json(&model_lattice),
        "representation_mode": request.get("representation_mode").and_then(Value::as_str).unwrap_or("DirectProduct"),
        "representation_source": if request.get("basis_functions").is_some() { "BasisFunctions" } else { "Matrices" },
        "representation_method": representation.method(),
        "group_order": symmetry.group().order(),
        "multiplication_table": symmetry.group().multiplication_table(),
        "antiunitary_flags": &antiunitary_flags,
        "operation_labels": operation_labels,
        "generator_indices": generator_indices,
        "generator_display_records": generator_display_records,
        "spatial_actions": spatial.iter().map(|operation| json!({
            "rotation": matrix_to_json(operation.rotation()),
            "translation": operation.translation().iter().map(element_to_json).collect::<Vec<_>>()
        })).collect::<Vec<_>>(),
        "spin_rotations": &serialized_spin_rotations,
        "spin_actions": &serialized_spin_actions,
        "image_site_indices": sites.image_site_indices(),
        "cell_translations": cell_translations,
        "local_dimensions": &local_dimensions,
        "spinor_basis_flags": &spinor_basis_flags,
        "site_dimensions": &global_dimensions,
        "orbital_layout": orbital_layout,
        "representation_matrices": representation.representation_matrices().iter()
            .map(matrix_to_json).collect::<Vec<_>>(),
    });
    let common = object(&common)?.clone();

    shell_indices
        .iter()
        .map(|&shell_index| {
            let bonds = shells[shell_index].bonds();
            let directed = compile_directed_bond_orbits(&symmetry, &sites, bonds)?;
            let constraints = if let Some(generators) = &continuous_generators {
                compile_bond_constraints_with_continuous(
                    &symmetry,
                    bonds,
                    &directed,
                    &global_dimensions,
                    &site_actions,
                    Some(generators),
                    hermitian,
                    kernel_method,
                    validation_kind,
                )?
            } else {
                compile_bond_constraints(
                    &symmetry,
                    bonds,
                    &directed,
                    &global_dimensions,
                    &site_actions,
                    hermitian,
                    kernel_method,
                    validation_kind,
                )?
            };
            let solved = parameter_space_result(&ParameterSpaceInput {
                context: &context,
                symmetry: &symmetry,
                bonds,
                directed: &directed,
                constraints: &constraints,
                local_dimensions: &global_dimensions,
                site_actions: &site_actions,
                representation_matrices: representation.representation_matrices(),
            })?;
            let directed_bond_orbits = directed
                .orbits()
                .iter()
                .map(|orbit| {
                    json!({
                        "representative_bond_index": orbit.representative_bond_index(),
                        "members": orbit.members().iter().map(|member| json!({
                            "bond_index": member.bond_index(),
                            "transporter_operation": member.transporter_operation()
                        })).collect::<Vec<_>>()
                    })
                })
                .collect::<Vec<_>>();
            let constraint_orbits = constraints
                .orbits()
                .iter()
                .enumerate()
                .map(|(orbit_index, orbit)| {
                    json!({
                        "orbit_index": orbit_index,
                        "source_directed_orbit_indices": orbit.source_directed_orbit_indices(),
                        "representative_bond_index": orbit.representative_bond_index(),
                        "dimensions": [orbit.dimensions().0, orbit.dimensions().1],
                        "stabilizer_generator_indices": orbit.stabilizer_generator_indices(),
                        "reverse_representative_operation": orbit.reverse_representative_operation(),
                        "nullity": orbit.kernel().nullity,
                        "residual_verified": orbit.kernel().residual_computed
                            .then_some(orbit.kernel().exact_residual_verified),
                        "independent_verified": orbit.kernel().independent_checked
                            .then_some(orbit.kernel().independent_verified),
                        "rank_nullity_verified": orbit.kernel().rank_nullity_checked
                            .then_some(orbit.kernel().rank_nullity_verified)
                    })
                })
                .collect::<Vec<_>>();
            let mut parameter_order = Vec::with_capacity(constraints.parameter_count());
            for (constraint_orbit_index, orbit) in constraints.orbits().iter().enumerate() {
                for orbit_parameter_index in 0..orbit.kernel().nullity {
                    parameter_order.push(json!({
                        "parameter_index": parameter_order.len(),
                        "constraint_orbit_index": constraint_orbit_index,
                        "orbit_parameter_index": orbit_parameter_index,
                        "representative_bond_index": orbit.representative_bond_index()
                    }));
                }
            }
            let mut result = common.clone();
            result.insert("bond_count".to_owned(), json!(bonds.len()));
            result.insert(
                "bonds".to_owned(),
                Value::Array(bonds.iter().map(periodic_bond_to_json).collect()),
            );
            result.insert(
                "directed_bond_action_table".to_owned(),
                json!(directed.action().action_table()),
            );
            result.insert(
                "reverse_bond_indices".to_owned(),
                json!(directed.reverse_indices()),
            );
            result.insert(
                "directed_bond_orbits".to_owned(),
                Value::Array(directed_bond_orbits),
            );
            result.insert(
                "constraint_orbits".to_owned(),
                Value::Array(constraint_orbits),
            );
            result.insert("parameter_order".to_owned(), Value::Array(parameter_order));
            result.insert("parameter_count".to_owned(), json!(constraints.parameter_count()));
            result.insert("parameter_space".to_owned(), solved);
            result.insert("hermitian".to_owned(), Value::Bool(hermitian));
            result.insert(
                "kernel_method".to_owned(),
                Value::String(kernel_method_name.clone()),
            );
            result.insert(
                "validation_level".to_owned(),
                Value::String(validation_level.clone()),
            );
            let hamiltonian =
                symbolic_hamiltonian_from_result(&context, &result, shell_index + 1)?;
            result.insert(
                "model_identity_sha256".to_owned(),
                hamiltonian["model_identity_sha256"].clone(),
            );
            result.insert("hamiltonian".to_owned(), hamiltonian);
            Ok(Value::Object(result))
        })
        .collect()
}
