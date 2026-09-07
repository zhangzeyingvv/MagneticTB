fn exact_vector_numeric(values: &[Element], field: &str) -> ExactResult<[f64; 3]> {
    if values.len() != 3 {
        return Err(ExactError::new(
            "InvalidCrystalStructureData",
            format!("{field} must have three components"),
        ));
    }
    let mut output = [0.0; 3];
    for (index, value) in values.iter().enumerate() {
        let (real, imaginary) = approximate_element(value)?;
        if !real.is_finite() || !imaginary.is_finite() || imaginary.abs() > 1.0e-12 {
            return Err(ExactError::new(
                "InvalidCrystalStructureData",
                format!("{field} must evaluate to a finite real vector"),
            ));
        }
        output[index] = real;
    }
    Ok(output)
}
fn row_times_matrix(vector: [f64; 3], matrix: [[f64; 3]; 3]) -> [f64; 3] {
    std::array::from_fn(|column| (0..3).map(|row| vector[row] * matrix[row][column]).sum())
}

#[allow(clippy::cast_precision_loss)]
fn crystal_translation_f64(value: i64) -> ExactResult<f64> {
    if value.unsigned_abs() > 9_007_199_254_740_991 {
        return Err(ExactError::new(
            "InvalidCrystalCellRange",
            "CellRange translation exceeds the exact integer range of f64",
        ));
    }
    Ok(value as f64)
}

#[allow(clippy::too_many_lines)]
fn crystal_structure_data(
    request: &serde_json::Map<String, Value>,
    catalog: Option<&DataCatalog>,
) -> ExactResult<Value> {
    let result = object(required(request, "model_result")?)?;
    let compiler_input = required(request, "compiler_input")?;
    let context = context_from_json(required(result, "field_context")?, 128)?;
    let lattice_exact = matrix_from_json(&context, required(result, "evaluated_model_lattice")?)?;
    if lattice_exact.rows() != 3 || lattice_exact.columns() != 3 {
        return Err(ExactError::new(
            "InvalidCrystalStructureData",
            "the current model lattice must be 3 by 3",
        ));
    }
    let mut lattice = [[0.0; 3]; 3];
    for (row, output_row) in lattice.iter_mut().enumerate() {
        for (column, output) in output_row.iter_mut().enumerate() {
            let (real, imaginary) = approximate_element(lattice_exact.entry(row, column)?)?;
            if !real.is_finite() || !imaginary.is_finite() || imaginary.abs() > 1.0e-12 {
                return Err(ExactError::new(
                    "InvalidCrystalStructureData",
                    "the current model lattice must evaluate to finite real values",
                ));
            }
            *output = real;
        }
    }
    let spatial = array(required(result, "spatial_actions")?)?
        .iter()
        .map(|operation| {
            let operation = object(operation)?;
            SpatialOperation::new(
                matrix_from_json(&context, required(operation, "rotation")?)?,
                exact_elements(&context, required(operation, "translation")?)?,
            )
            .map_err(|error| ExactError::new(error.tag(), error.to_string()))
        })
        .collect::<ExactResult<Vec<_>>>()?;
    let antiunitary_flags: Vec<bool> =
        from_value(required(result, "antiunitary_flags")?, "antiunitary_flags")?;
    let spin_rotations = result
        .get("spin_rotations")
        .filter(|value| !value.is_null())
        .map(|value| exact_matrices(&context, value))
        .transpose()?;
    if spatial.is_empty()
        || antiunitary_flags.len() != spatial.len()
        || spin_rotations
            .as_ref()
            .is_some_and(|rotations| rotations.len() != spatial.len())
    {
        return Err(ExactError::new(
            "InvalidCrystalStructureData",
            "the current model symmetry data are not aligned",
        ));
    }
    let cell_range: [[i64; 2]; 3] = from_value(required(request, "cell_range")?, "cell_range")?;
    if cell_range.iter().any(|bounds| bounds[0] > bounds[1]) {
        return Err(ExactError::new(
            "InvalidCrystalCellRange",
            "CellRange bounds must be ordered integers",
        ));
    }
    let minimum_length = lattice
        .iter()
        .map(|row| row.iter().map(|value| value * value).sum::<f64>().sqrt())
        .fold(f64::INFINITY, f64::min);
    if !minimum_length.is_finite() || minimum_length <= 0.0 {
        return Err(ExactError::new(
            "InvalidCrystalStructureData",
            "the current model lattice has an invalid row length",
        ));
    }
    let positive_option = |name: &str, default: f64| -> ExactResult<f64> {
        let value =
            request
                .get(name)
                .filter(|value| !value.is_null())
                .map_or(Ok(default), |value| {
                    value.as_f64().ok_or_else(|| {
                        ExactError::new(
                            "InvalidCrystalStructureOption",
                            format!("{name} must be Automatic or a positive number"),
                        )
                    })
                })?;
        if !value.is_finite() || value <= 0.0 {
            return Err(ExactError::new(
                "InvalidCrystalStructureOption",
                format!("{name} must be Automatic or a positive number"),
            ));
        }
        Ok(value)
    };
    let moment_scale = positive_option("moment_scale", 0.55 * minimum_length)?;
    let atom_radius = positive_option("atom_radius", 0.08 * minimum_length)?;
    let raw_seeds = tagged_items(ordered_lookup(compiler_input, "WyckoffPosition")?)?;
    if raw_seeds.is_empty() {
        return Err(ExactError::new(
            "InvalidCrystalStructureData",
            "the current model has no atomic positions",
        ));
    }
    let compiled_orbits = exact_site_orbits(&context, required(result, "site_orbits")?)?;
    if compiled_orbits.len() != raw_seeds.len() {
        return Err(ExactError::new(
            "InvalidCrystalStructureData",
            "compiled site orbits are not aligned with the ordered Wyckoff inputs",
        ));
    }
    let operation_indices = (0..spatial.len()).collect::<Vec<_>>();
    let mut base_records = Vec::new();
    let data_loaded = result
        .get("data_loaded")
        .and_then(Value::as_bool)
        .unwrap_or(false);
    let data_selections = if data_loaded {
        let selections = array(required(result, "wyckoff_selections")?)?;
        if selections.len() != raw_seeds.len() {
            return Err(ExactError::new(
                "InvalidCrystalStructureData",
                "Data Wyckoff selections are not aligned with compiled site orbits",
            ));
        }
        Some(selections)
    } else {
        None
    };
    for (orbit_index, (raw_seed, compiled_orbit)) in
        raw_seeds.iter().zip(&compiled_orbits).enumerate()
    {
        let seed = encoded_matrix_from_json(&context, raw_seed)?;
        if seed.rows() != 2 || seed.columns() != 3 {
            return Err(ExactError::new(
                "InvalidCrystalStructureData",
                "each Wyckoff seed must be an exact 2 by 3 position/moment matrix",
            ));
        }
        if compiled_orbit.is_empty() {
            return Err(ExactError::new(
                "InvalidCrystalStructureData",
                "compiled site orbits must be nonempty",
            ));
        }
        if let Some(selections) = data_selections {
            let catalog = catalog.ok_or_else(|| {
                ExactError::new(
                    "MissingDataCatalog",
                    "Data crystal structure rendering requires the frozen Rust Data catalog",
                )
            })?;
            let msg_id = required(result, "msg_id")?
                .as_str()
                .ok_or_else(|| malformed("msg_id must be a string"))?;
            let selection = object(&selections[orbit_index])?;
            let source_ordinal: usize =
                from_value(required(selection, "source_ordinal")?, "source_ordinal")?;
            let letter = required(selection, "letter")?
                .as_str()
                .ok_or_else(|| malformed("Wyckoff selection letter must be a string"))?;
            let wyckoff = catalog.wyckoff(msg_id).ok_or_else(|| {
                ExactError::new(
                    "UnknownStableId",
                    format!("no Wyckoff data close over MSG stable ID {msg_id}"),
                )
            })?;
            let entry = wyckoff
                .entries()
                .iter()
                .find(|entry| {
                    entry.source_ordinal() == source_ordinal && entry.letter() == letter
                })
                .ok_or_else(|| {
                    ExactError::new(
                        "UnknownWyckoffOrdinal",
                        format!(
                            "MSG {msg_id} has no Wyckoff source ordinal {source_ordinal} with letter {letter}"
                        ),
                    )
                })?;
            if entry.positions().len() != compiled_orbit.len() {
                return Err(ExactError::new(
                    "InvalidCrystalStructureData",
                    "Data Wyckoff positions are not aligned with the compiled site orbit",
                ));
            }
            let mut bindings = BTreeMap::new();
            for (column, name) in ["x", "y", "z"].into_iter().enumerate() {
                bindings.insert(format!("MagneticTB`{name}"), seed.entry(0, column)?.clone());
            }
            for (column, name) in ["mx", "my", "mz"].into_iter().enumerate() {
                bindings.insert(format!("MagneticTB`{name}"), seed.entry(1, column)?.clone());
            }
            for (equivalent_index, (position, compiled_position)) in
                entry.positions().iter().zip(compiled_orbit).enumerate()
            {
                let evaluated_position = position
                    .coordinates()
                    .iter()
                    .map(|value| {
                        let evaluated = evaluate_data_expression(&context, value, &bindings)
                            .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
                        fractional_part(&evaluated)
                            .map_err(|error| ExactError::new(error.tag(), error.to_string()))
                    })
                    .collect::<ExactResult<Vec<_>>>()?;
                if evaluated_position != *compiled_position {
                    return Err(ExactError::new(
                        "InvalidCrystalStructureData",
                        "Data Wyckoff order differs from the compiled site-orbit order",
                    ));
                }
                let moment = position
                    .moment()
                    .iter()
                    .map(|value| {
                        evaluate_data_expression(&context, value, &bindings)
                            .map_err(|error| ExactError::new(error.tag(), error.to_string()))
                    })
                    .collect::<ExactResult<Vec<_>>>()?;
                base_records.push((
                    orbit_index,
                    equivalent_index,
                    MagneticSite {
                        position: compiled_position.clone(),
                        moment,
                    },
                ));
            }
        } else {
            let magnetic_seed = MagneticSite {
                position: compiled_orbit[0].clone(),
                moment: (0..3)
                    .map(|column| seed.entry(1, column).cloned())
                    .collect::<ExactResult<Vec<_>>>()?,
            };
            let generated = magnetic_orbit(
                &magnetic_seed,
                &operation_indices,
                &spatial,
                spin_rotations.as_deref(),
                &antiunitary_flags,
            )?;
            for (equivalent_index, compiled_position) in compiled_orbit.iter().enumerate() {
                let mut matches = generated
                    .iter()
                    .filter(|site| site.position == *compiled_position);
                let site = matches.next().ok_or_else(|| {
                    ExactError::new(
                        "InvalidCrystalStructureData",
                        "compiled site position has no matching magnetic-site image",
                    )
                })?;
                if matches.next().is_some() {
                    return Err(ExactError::new(
                        "InvalidCrystalStructureData",
                        "compiled site position has multiple magnetic-site images",
                    ));
                }
                base_records.push((orbit_index, equivalent_index, site.clone()));
            }
        }
    }
    let mut translations = Vec::new();
    for first in cell_range[0][0]..=cell_range[0][1] {
        for second in cell_range[1][0]..=cell_range[1][1] {
            for third in cell_range[2][0]..=cell_range[2][1] {
                translations.push([first, second, third]);
            }
        }
    }
    let mut atom_records = Vec::new();
    let mut magnetic_count = 0usize;
    for translation in &translations {
        for (orbit_index, equivalent_index, site) in &base_records {
            let fractional_position = exact_vector_numeric(&site.position, "fractional position")?;
            let fractional_moment = exact_vector_numeric(&site.moment, "fractional moment")?;
            let translation_numeric = [
                crystal_translation_f64(translation[0])?,
                crystal_translation_f64(translation[1])?,
                crystal_translation_f64(translation[2])?,
            ];
            let shifted =
                std::array::from_fn(|axis| fractional_position[axis] + translation_numeric[axis]);
            let cartesian_position = row_times_matrix(shifted, lattice);
            let cartesian_moment = row_times_matrix(fractional_moment, lattice);
            let norm = cartesian_moment
                .iter()
                .map(|value| value * value)
                .sum::<f64>()
                .sqrt();
            let magnetic = norm > 1.0e-12;
            let arrow_end: Option<[f64; 3]> = magnetic.then(|| {
                std::array::from_fn(|axis| {
                    cartesian_position[axis] + moment_scale * cartesian_moment[axis] / norm
                })
            });
            magnetic_count += usize::from(magnetic);
            atom_records.push(json!({
                "OrbitIndex": orbit_index + 1,
                "EquivalentIndex": equivalent_index + 1,
                "FractionalPosition": site.position.iter().map(element_to_json).collect::<Vec<_>>(),
                "FractionalMoment": site.moment.iter().map(element_to_json).collect::<Vec<_>>(),
                "CellTranslation": translation,
                "CartesianPosition": cartesian_position,
                "CartesianMoment": cartesian_moment,
                "Magnetic": magnetic,
                "ArrowEnd": arrow_end
            }));
        }
    }
    let mut cell_edges = Vec::new();
    for translation in &translations {
        for axis in 0..3 {
            for first in 0..=1 {
                for second in 0..=1 {
                    let mut corner = [0.0; 3];
                    let other_axes = (0..3)
                        .filter(|candidate| *candidate != axis)
                        .collect::<Vec<_>>();
                    corner[other_axes[0]] = f64::from(first);
                    corner[other_axes[1]] = f64::from(second);
                    let translation_numeric = [
                        crystal_translation_f64(translation[0])?,
                        crystal_translation_f64(translation[1])?,
                        crystal_translation_f64(translation[2])?,
                    ];
                    let start =
                        std::array::from_fn(|index| corner[index] + translation_numeric[index]);
                    let mut end = start;
                    end[axis] += 1.0;
                    cell_edges.push([
                        row_times_matrix(start, lattice),
                        row_times_matrix(end, lattice),
                    ]);
                }
            }
        }
    }
    Ok(json!({
        "Schema": "MagneticTBCrystalStructureData",
        "SchemaVersion": 1,
        "Lattice": lattice,
        "CellRange": cell_range,
        "Translations": translations,
        "AtomRadius": atom_radius,
        "MomentScale": moment_scale,
        "AtomRecords": atom_records,
        "AtomCount": atom_records.len(),
        "MagneticAtomCount": magnetic_count,
        "CellEdges": cell_edges
    }))
}

fn split_magnetic_orbit(
    full_orbit: &[MagneticSite],
    subgroup_indices: &[usize],
    spatial: &[SpatialOperation],
    spin_rotations: Option<&[ExactMatrix]>,
    antiunitary_flags: &[bool],
) -> ExactResult<Vec<Vec<usize>>> {
    let mut remaining = (0..full_orbit.len()).collect::<Vec<_>>();
    let mut groups = Vec::new();
    while let Some(&seed_index) = remaining.first() {
        let images = magnetic_orbit(
            &full_orbit[seed_index],
            subgroup_indices,
            spatial,
            spin_rotations,
            antiunitary_flags,
        )?;
        let mut indices = images
            .iter()
            .map(|image| {
                let matches = full_orbit
                    .iter()
                    .enumerate()
                    .filter_map(|(index, candidate)| (candidate == image).then_some(index))
                    .collect::<Vec<_>>();
                if matches.len() == 1 {
                    Ok(matches[0])
                } else {
                    Err(ExactError::new(
                        "OrbitSplitFailed",
                        "a retained-subgroup image does not identify exactly one parent site",
                    ))
                }
            })
            .collect::<ExactResult<Vec<_>>>()?;
        indices.sort_unstable();
        indices.dedup();
        if !indices.contains(&seed_index) || indices.iter().any(|index| !remaining.contains(index))
        {
            return Err(ExactError::new(
                "OrbitSplitFailed",
                "retained-subgroup site orbits overlap or omit their seed",
            ));
        }
        remaining.retain(|index| !indices.contains(index));
        groups.push(indices);
    }
    if groups.iter().flatten().count() != full_orbit.len() {
        return Err(ExactError::new(
            "OrbitSplitFailed",
            "retained-subgroup site orbits do not partition the parent orbit",
        ));
    }
    Ok(groups)
}

fn exact_matrix_expression(rows: &[Vec<Element>]) -> ExactResult<Value> {
    if rows.is_empty() || rows[0].is_empty() || rows.iter().any(|row| row.len() != rows[0].len()) {
        return Err(malformed(
            "exact matrix expression must be nonempty and rectangular",
        ));
    }
    Ok(json!({
        "kind": "matrix",
        "rows": rows.len(),
        "columns": rows[0].len(),
        "entries": rows.iter().flatten().map(|value| {
            quadratic_expression_from_element(value)
                .unwrap_or_else(|_| cyclotomic_expression_from_element(value))
        }).collect::<Vec<_>>()
    }))
}

fn selected_symmetry_information(
    compiler_input: &Value,
    subgroup_indices: &[usize],
    spatial: &[SpatialOperation],
    antiunitary_flags: &[bool],
    operation_labels: &[String],
) -> ExactResult<Value> {
    if let Some(source) = optional_ordered_lookup(compiler_input, "SymmetryInformation") {
        if let Ok(items) = tagged_items(source) {
            return Ok(json!({
                "kind": "list",
                "items": subgroup_indices.iter().map(|&index| {
                    items.get(index).cloned().ok_or_else(|| ExactError::new(
                        "InvalidBrokenSymmetryInput",
                        "symmetry-information order is not aligned with the compiled group"
                    ))
                }).collect::<ExactResult<Vec<_>>>()?
            }));
        }
        let source = object(source)?;
        if source.get("kind").and_then(Value::as_str) == Some("ordered_association") {
            let entries = array(required(source, "entries")?)?;
            return Ok(json!({
                "kind": "ordered_association",
                "entries": subgroup_indices.iter().map(|&index| {
                    entries.get(index).cloned().ok_or_else(|| ExactError::new(
                        "InvalidBrokenSymmetryInput",
                        "spin-space information order is not aligned with the compiled group"
                    ))
                }).collect::<ExactResult<Vec<_>>>()?
            }));
        }
    }
    let items = subgroup_indices
        .iter()
        .map(|&index| {
            let operation = spatial.get(index).ok_or_else(|| {
                ExactError::new(
                    "InvalidBrokenSymmetryInput",
                    "subgroup operation index is out of range",
                )
            })?;
            let rotation = (0..3)
                .map(|row| {
                    (0..3)
                        .map(|column| operation.rotation().entry(row, column).cloned())
                        .collect::<ExactResult<Vec<_>>>()
                })
                .collect::<ExactResult<Vec<_>>>()?;
            Ok(json!({
                "kind": "list",
                "items": [
                    operation_labels[index],
                    exact_matrix_expression(&rotation)?,
                    {
                        "kind": "list",
                        "items": operation.translation().iter().map(|value| {
                            quadratic_expression_from_element(value)
                                .unwrap_or_else(|_| cyclotomic_expression_from_element(value))
                        }).collect::<Vec<_>>()
                    },
                    if antiunitary_flags[index] { "T" } else { "F" }
                ]
            }))
        })
        .collect::<ExactResult<Vec<_>>>()?;
    Ok(json!({"kind": "list", "items": items}))
}

#[allow(clippy::too_many_lines)]
fn broken_symmetry_init_rules(request: &serde_json::Map<String, Value>) -> ExactResult<Value> {
    let result = object(required(request, "model_result")?)?;
    let compiler_input = required(request, "compiler_input")?;
    let context = context_from_json(required(result, "field_context")?, 128)?;
    let multiplication_table: Vec<Vec<usize>> = from_value(
        required(result, "multiplication_table")?,
        "multiplication_table",
    )?;
    let group = GroupAlgebra::new(multiplication_table)
        .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
    let labels: Vec<String> =
        from_value(required(result, "operation_labels")?, "operation_labels")?;
    let selector = array(required(request, "selector")?)?;
    let mut generators = Vec::new();
    for item in selector {
        let index =
            if let Some(index) = item.as_u64().and_then(|value| usize::try_from(value).ok()) {
                (index < group.order()).then_some(index)
            } else if let Some(label) = item.as_str() {
                let matches = labels
                    .iter()
                    .enumerate()
                    .filter_map(|(index, candidate)| (candidate == label).then_some(index))
                    .collect::<Vec<_>>();
                (matches.len() == 1).then_some(matches[0])
            } else {
                None
            }
            .ok_or_else(|| {
                ExactError::new(
                    "InvalidBrokenSymmetrySelector",
                    "retained generators must be valid operation indices or unambiguous labels",
                )
            })?;
        if !generators.contains(&index) {
            generators.push(index);
        }
    }
    let subgroup_indices = group
        .generated_subgroup(&generators)
        .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
    let spatial = array(required(result, "spatial_actions")?)?
        .iter()
        .map(|operation| {
            let operation = object(operation)?;
            SpatialOperation::new(
                matrix_from_json(&context, required(operation, "rotation")?)?,
                exact_elements(&context, required(operation, "translation")?)?,
            )
            .map_err(|error| ExactError::new(error.tag(), error.to_string()))
        })
        .collect::<ExactResult<Vec<_>>>()?;
    let antiunitary_flags: Vec<bool> =
        from_value(required(result, "antiunitary_flags")?, "antiunitary_flags")?;
    if spatial.len() != group.order()
        || labels.len() != group.order()
        || antiunitary_flags.len() != group.order()
    {
        return Err(ExactError::new(
            "InvalidBrokenSymmetryInput",
            "compiled group operations and metadata are not aligned",
        ));
    }
    let spin_rotations = result
        .get("spin_rotations")
        .filter(|value| !value.is_null())
        .map(|value| exact_matrices(&context, value))
        .transpose()?;
    if spin_rotations
        .as_ref()
        .is_some_and(|rotations| rotations.len() != group.order())
    {
        return Err(ExactError::new(
            "InvalidBrokenSymmetryInput",
            "spin rotations are not aligned with the compiled group",
        ));
    }
    let raw_seeds = tagged_items(ordered_lookup(compiler_input, "WyckoffPosition")?)?;
    let raw_basis = optional_ordered_lookup(compiler_input, "BasisFunctions").ok_or_else(|| {
        ExactError::new(
            "ExplicitRepresentationInput",
            "brokenSymmetryInitRules cannot regenerate initfromrep representation matrices",
        )
    })?;
    let raw_basis = tagged_items(raw_basis)?;
    if raw_basis.len() != raw_seeds.len() {
        return Err(ExactError::new(
            "InvalidBrokenSymmetryInput",
            "basis specifications are not aligned with Wyckoff seeds",
        ));
    }
    let all_indices = (0..group.order()).collect::<Vec<_>>();
    let mut split_seeds = Vec::new();
    let mut split_basis = Vec::new();
    let mut source_orbit_indices = Vec::new();
    for (orbit_index, raw_seed) in raw_seeds.iter().enumerate() {
        let seed = encoded_matrix_from_json(&context, raw_seed)?;
        if seed.rows() != 2 || seed.columns() != 3 {
            return Err(ExactError::new(
                "InvalidBrokenSymmetryInput",
                "each Wyckoff seed must be an exact 2 by 3 position/moment matrix",
            ));
        }
        let magnetic_seed = MagneticSite {
            position: (0..3)
                .map(|column| seed.entry(0, column).cloned())
                .collect::<ExactResult<Vec<_>>>()?,
            moment: (0..3)
                .map(|column| seed.entry(1, column).cloned())
                .collect::<ExactResult<Vec<_>>>()?,
        };
        let full_orbit = magnetic_orbit(
            &magnetic_seed,
            &all_indices,
            &spatial,
            spin_rotations.as_deref(),
            &antiunitary_flags,
        )?;
        let groups = split_magnetic_orbit(
            &full_orbit,
            &subgroup_indices,
            &spatial,
            spin_rotations.as_deref(),
            &antiunitary_flags,
        )?;
        for group in groups {
            let representative = &full_orbit[group[0]];
            split_seeds.push(exact_matrix_expression(&[
                representative.position.clone(),
                representative.moment.clone(),
            ])?);
            split_basis.push(raw_basis[orbit_index].clone());
            source_orbit_indices.push(orbit_index);
        }
    }
    let symmetry_information = selected_symmetry_information(
        compiler_input,
        &subgroup_indices,
        &spatial,
        &antiunitary_flags,
        &labels,
    )?;
    let canonical_input = json!({
        "kind": "ordered_association",
        "entries": [
            {"key": "Lattice", "value": ordered_lookup(compiler_input, "Lattice")?},
            {"key": "LatticeParameters", "value": ordered_lookup(compiler_input, "LatticeParameters")?},
            {"key": "WyckoffPosition", "value": {"kind": "list", "items": split_seeds}},
            {"key": "SymmetryInformation", "value": symmetry_information},
            {"key": "BasisFunctions", "value": {"kind": "list", "items": split_basis}},
            {"key": "Debug", "value": optional_ordered_lookup(compiler_input, "Debug").cloned().unwrap_or(Value::Bool(false))},
            {"key": "InitialBondShells", "value": optional_ordered_lookup(compiler_input, "InitialBondShells").cloned().unwrap_or_else(|| json!(10))},
            {"key": "GenerateSymmetryGroup", "value": false},
            {"key": "RepresentationMode", "value": optional_ordered_lookup(compiler_input, "RepresentationMode").cloned().unwrap_or_else(|| Value::String("DirectProduct".to_owned()))}
        ]
    });
    Ok(json!({
        "schema": "magnetictb.broken_symmetry_init_rules.v1",
        "version": 1,
        "subgroup_indices": subgroup_indices,
        "source_orbit_indices": source_orbit_indices,
        "compiler_input": canonical_input
    }))
}
