fn exact_numeric_matrix(matrix: &ExactMatrix, field: &str) -> ExactResult<DMatrix<Complex<f64>>> {
    let entries = matrix
        .entries()
        .iter()
        .map(approximate_element)
        .map(|value| value.map(|(real, imaginary)| Complex::new(real, imaginary)))
        .collect::<ExactResult<Vec<_>>>()?;
    if entries
        .iter()
        .any(|value| !value.re.is_finite() || !value.im.is_finite())
    {
        return Err(ExactError::new(
            "NonFiniteBandCorepInput",
            format!("{field} is not finite in f64"),
        ));
    }
    Ok(DMatrix::from_row_slice(
        matrix.rows(),
        matrix.columns(),
        &entries,
    ))
}
fn standard_k_path_result(request: &serde_json::Map<String, Value>) -> ExactResult<Value> {
    let context = context_from_json(required(request, "context")?, 128)?;
    let lattice = encoded_matrix_from_json(&context, required(request, "lattice")?)?;
    let lattice = exact_numeric_matrix(&lattice, "lattice")?;
    if lattice.nrows() != 3 || lattice.ncols() != 3 {
        return Err(ExactError::new(
            "InvalidStandardKPathLattice",
            "lattice must be a 3 by 3 matrix",
        ));
    }
    let mut rows = [[0.0; 3]; 3];
    for row in 0..3 {
        for column in 0..3 {
            let value = lattice[(row, column)];
            if value.im.abs() > 1.0e-12 {
                return Err(ExactError::new(
                    "InvalidStandardKPathLattice",
                    "lattice entries must be real",
                ));
            }
            rows[row][column] = value.re;
        }
    }
    let requested = required(request, "bravais_type")?
        .as_str()
        .ok_or_else(|| malformed("bravais_type must be a string"))?;
    let requested = (requested != "Automatic").then_some(requested);
    let tolerance = required(request, "tolerance")?
        .as_f64()
        .ok_or_else(|| malformed("tolerance must be a real number"))?;
    let data = standard_k_path(rows, requested, tolerance)
        .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
    let path = data
        .path
        .iter()
        .map(|segment| {
            json!([
                [segment.start, segment.end],
                [segment.start_label, segment.end_label]
            ])
        })
        .collect::<Vec<_>>();
    Ok(json!({
        "BravaisType": data.bravais_type,
        "BZType": data.bz_type,
        "Parameters": data.parameters,
        "GramMatrix": data.gram_matrix,
        "ReciprocalCosines": data.reciprocal_cosines,
        "Points": data.points,
        "Sequences": data.sequences,
        "Sequence": data.sequences,
        "Path": path,
        "LabelConvention": "Bradley-CracknellWithSetyawan-CurtaroloFallback"
    }))
}

fn brillouin_zone_data(request: &serde_json::Map<String, Value>) -> ExactResult<Value> {
    let context = context_from_json(required(request, "context")?, 128)?;
    let lattice = encoded_matrix_from_json(&context, required(request, "lattice")?)?;
    let lattice = exact_numeric_matrix(&lattice, "lattice")?;
    if lattice.nrows() != 3 || lattice.ncols() != 3 {
        return Err(ExactError::new(
            "InvalidBrillouinZoneInput",
            "lattice must be a 3 by 3 matrix",
        ));
    }
    let mut rows = [[0.0; 3]; 3];
    for row in 0..3 {
        for column in 0..3 {
            let value = lattice[(row, column)];
            if !value.re.is_finite() || !value.im.is_finite() || value.im.abs() > 1.0e-12 {
                return Err(ExactError::new(
                    "InvalidBrillouinZoneInput",
                    "lattice entries must evaluate to finite real numbers",
                ));
            }
            rows[row][column] = value.re;
        }
    }
    let translation_range = required(request, "translation_range")?
        .as_u64()
        .and_then(|value| usize::try_from(value).ok())
        .ok_or_else(|| {
            ExactError::new(
                "InvalidBrillouinZoneOption",
                "TranslationRange must be a positive integer",
            )
        })?;
    let tolerance = required(request, "tolerance")?.as_f64().ok_or_else(|| {
        ExactError::new(
            "InvalidBrillouinZoneOption",
            "Tolerance must be a positive real number",
        )
    })?;
    let zone = first_brillouin_zone(rows, translation_range, tolerance)
        .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
    let raw_path: Vec<([[f64; 3]; 2], [String; 2])> =
        from_value(required(request, "path")?, "path")?;
    let path = raw_path
        .iter()
        .map(|(points, labels)| (points[0], points[1], labels[0].clone(), labels[1].clone()))
        .collect::<Vec<_>>();
    let displayed =
        fold_band_path_to_first_bz(zone.reciprocal_lattice, &path, translation_range, tolerance)
            .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
    Ok(json!({
        "Schema": "MagneticTBBrillouinZoneData",
        "SchemaVersion": 1,
        "Lattice": rows,
        "ReciprocalLattice": zone.reciprocal_lattice,
        "Vertices": zone.vertices,
        "Facets": zone.facets,
        "FacetCount": zone.facets.len(),
        "TranslationRange": translation_range,
        "Tolerance": tolerance,
        "BravaisType": request.get("bravais_type").cloned().unwrap_or(Value::Null),
        "BZType": request.get("bz_type").cloned().unwrap_or(Value::Null),
        "KPath": raw_path,
        "DisplayedKPath": displayed.iter().map(|segment| json!([
            [segment.start, segment.end],
            [segment.start_label, segment.end_label]
        ])).collect::<Vec<_>>(),
        "CartesianKPath": displayed.iter().map(|segment| [
            segment.cartesian_start,
            segment.cartesian_end
        ]).collect::<Vec<_>>()
    }))
}

fn real_matrix_three(matrix: &ExactMatrix, field: &str) -> ExactResult<Matrix3<f64>> {
    if matrix.rows() != 3 || matrix.columns() != 3 {
        return Err(ExactError::new(
            "InvalidBandCorepSymmetry",
            format!("{field} must be 3 by 3"),
        ));
    }
    let numeric = exact_numeric_matrix(matrix, field)?;
    let mut entries = [0.0; 9];
    for row in 0..3 {
        for column in 0..3 {
            let value = numeric[(row, column)];
            if value.im.abs() > 1.0e-12 {
                return Err(ExactError::new(
                    "InvalidBandCorepSymmetry",
                    format!("{field} must be real"),
                ));
            }
            entries[row * 3 + column] = value.re;
        }
    }
    Ok(Matrix3::from_row_slice(&entries))
}

fn real_vector_three(values: &[Element], field: &str) -> ExactResult<Vector3<f64>> {
    if values.len() != 3 {
        return Err(ExactError::new(
            "InvalidBandCorepInput",
            format!("{field} must have three coordinates"),
        ));
    }
    let mut numeric = [0.0; 3];
    for (index, value) in values.iter().enumerate() {
        let (real, imaginary) = approximate_element(value)?;
        if imaginary.abs() > 1.0e-12 {
            return Err(ExactError::new(
                "InvalidBandCorepInput",
                format!("{field} must be real"),
            ));
        }
        numeric[index] = real;
    }
    Ok(Vector3::from_row_slice(&numeric))
}

fn numeric_json_matrix(value: &Value, field: &str) -> ExactResult<DMatrix<Complex<f64>>> {
    let rows = array(value)?;
    if rows.is_empty() {
        return Err(malformed(format!("{field} must be nonempty")));
    }
    let columns = array(&rows[0])?.len();
    if columns == 0
        || rows
            .iter()
            .any(|row| array(row).map_or(true, |row| row.len() != columns))
    {
        return Err(malformed(format!(
            "{field} must be rectangular and nonempty"
        )));
    }
    let entries = rows
        .iter()
        .flat_map(|row| array(row).expect("validated numeric matrix row"))
        .map(|entry| {
            let record = object(entry)?;
            let real = required(record, "real")?
                .as_f64()
                .ok_or_else(|| malformed(format!("{field} real part must be numeric")))?;
            let imaginary = required(record, "imaginary")?
                .as_f64()
                .ok_or_else(|| malformed(format!("{field} imaginary part must be numeric")))?;
            finite_complex(real, imaginary, field)
                .map(|(real, imaginary)| Complex::new(real, imaginary))
        })
        .collect::<ExactResult<Vec<_>>>()?;
    Ok(DMatrix::from_row_slice(rows.len(), columns, &entries))
}

fn numeric_matrix_json(matrix: &DMatrix<Complex<f64>>) -> Value {
    Value::Array(
        (0..matrix.nrows())
            .map(|row| {
                Value::Array(
                    (0..matrix.ncols())
                        .map(|column| {
                            let value = matrix[(row, column)];
                            json!({"real": value.re, "imaginary": value.im})
                        })
                        .collect(),
                )
            })
            .collect(),
    )
}

fn trace_sign(
    canonical: &DMatrix<Complex<f64>>,
    magnetictb: &DMatrix<Complex<f64>>,
) -> ExactResult<f64> {
    if canonical.nrows() != 2
        || canonical.ncols() != 2
        || magnetictb.nrows() != 2
        || magnetictb.ncols() != 2
    {
        return Err(ExactError::new(
            "InvalidBandCorepSpinRotation",
            "SOC spin rotations must be 2 by 2",
        ));
    }
    let inverse = magnetictb.clone().try_inverse().ok_or_else(|| {
        ExactError::new(
            "InvalidBandCorepSpinRotation",
            "MagneticTB spin rotation is singular",
        )
    })?;
    let relative = canonical * inverse;
    let plus = DMatrix::identity(2, 2);
    let minus = -&plus;
    let distance = |target: &DMatrix<Complex<f64>>| {
        (&relative - target)
            .iter()
            .map(|value| value.norm())
            .fold(0.0, f64::max)
    };
    if distance(&plus) <= 1.0e-8 {
        Ok(1.0)
    } else if distance(&minus) <= 1.0e-8 {
        Ok(-1.0)
    } else {
        Err(ExactError::new(
            "IncompatibleBandCorepSpinConvention",
            "MSGCorep and MagneticTB spin lifts do not differ by a central sign",
        ))
    }
}

fn numeric_hamiltonian_matrix(value: &Value) -> ExactResult<DMatrix<Complex<f64>>> {
    let record = object(value)?;
    numeric_json_matrix(required(record, "matrix")?, "Hamiltonian matrix")
}

fn aligned_band_operations(
    context: &Arc<CyclotomicContext>,
    model: &serde_json::Map<String, Value>,
    request: &serde_json::Map<String, Value>,
) -> ExactResult<(Vec<SpatialOperation>, Vec<bool>)> {
    let model_operations = spatial_operations(context, required(model, "spatial_actions")?)?;
    let model_antiunitary: Vec<bool> =
        from_value(required(model, "antiunitary_flags")?, "antiunitary_flags")?;
    let external_value = required(request, "msgcorep_operations")?;
    let external_operations = spatial_operations(context, external_value)?;
    let external_antiunitary = array(external_value)?
        .iter()
        .map(|operation| {
            required(object(operation)?, "antiunitary")?
                .as_bool()
                .ok_or_else(|| malformed("MSGCorep antiunitary flag must be bool"))
        })
        .collect::<ExactResult<Vec<_>>>()?;
    if model_operations.len() != external_operations.len()
        || model_antiunitary != external_antiunitary
        || model_operations
            .iter()
            .zip(&external_operations)
            .any(|(left, right)| left != right)
    {
        return Err(ExactError::new(
            "MSGCorepOperationOrderMismatch",
            "the prepared model must use the complete MSGCorep operation list in its original order",
        ));
    }
    Ok((model_operations, model_antiunitary))
}

fn band_orbital_positions(
    context: &Arc<CyclotomicContext>,
    model: &serde_json::Map<String, Value>,
) -> ExactResult<Vec<Vector3<f64>>> {
    array(required(model, "orbital_layout")?)?
        .iter()
        .enumerate()
        .map(|(index, orbital)| {
            let orbital = object(orbital)?;
            let position = exact_elements(context, required(orbital, "fractional_position")?)?;
            real_vector_three(&position, &format!("orbital {} position", index + 1))
        })
        .collect()
}

fn band_kpoints(
    context: &Arc<CyclotomicContext>,
    request: &serde_json::Map<String, Value>,
) -> ExactResult<Vec<Vector3<f64>>> {
    let points = array(required(request, "kset")?)?
        .iter()
        .enumerate()
        .map(|(point_index, point)| {
            let point = array(point)?;
            if point.len() != 3 {
                return Err(ExactError::new(
                    "InvalidBandCorepKSet",
                    format!("k point {} must have three coordinates", point_index + 1),
                ));
            }
            Ok(Vector3::new(
                numeric_real(context, &point[0], "k point x")?,
                numeric_real(context, &point[1], "k point y")?,
                numeric_real(context, &point[2], "k point z")?,
            ))
        })
        .collect::<ExactResult<Vec<_>>>()?;
    if points.is_empty() {
        return Err(ExactError::new(
            "InvalidBandCorepKSet",
            "kset must be nonempty",
        ));
    }
    Ok(points)
}

fn band_hamiltonians(
    hamiltonian: &serde_json::Map<String, Value>,
    request: &serde_json::Map<String, Value>,
    reduced_kpoints: &[Vector3<f64>],
) -> ExactResult<Vec<DMatrix<Complex<f64>>>> {
    let parameters = required(request, "parameters")?;
    reduced_kpoints
        .iter()
        .map(|point| {
            let evaluation = serde_json::Map::from_iter([
                ("hamiltonian".to_owned(), Value::Object(hamiltonian.clone())),
                ("parameters".to_owned(), parameters.clone()),
                (
                    "momentum".to_owned(),
                    json!([
                        std::f64::consts::TAU * point[0],
                        std::f64::consts::TAU * point[1],
                        std::f64::consts::TAU * point[2]
                    ]),
                ),
            ]);
            numeric_hamiltonian_matrix(&evaluate_symbolic_hamiltonian(&evaluation)?)
        })
        .collect()
}

fn band_spin_rotations(
    request: &serde_json::Map<String, Value>,
    operation_count: usize,
    soc: bool,
) -> ExactResult<Option<Vec<DMatrix<Complex<f64>>>>> {
    if !soc {
        return Ok(None);
    }
    let result = array(required(request, "msgcorep_spin_rotations")?)?
        .iter()
        .enumerate()
        .map(|(index, matrix)| {
            numeric_json_matrix(matrix, &format!("MSGCorep spin rotation {}", index + 1))
        })
        .collect::<ExactResult<Vec<_>>>()?;
    if result.len() != operation_count {
        return Err(ExactError::new(
            "InvalidBandCorepSpinRotation",
            "MSGCorep spin rotations must align with the operation list",
        ));
    }
    Ok(Some(result))
}

type NumericBandOperations = (Vec<NumericSymmetryOperation>, Vec<DMatrix<Complex<f64>>>);

fn band_numeric_operations(
    context: &Arc<CyclotomicContext>,
    model: &serde_json::Map<String, Value>,
    request: &serde_json::Map<String, Value>,
    model_operations: &[SpatialOperation],
    model_antiunitary: &[bool],
) -> ExactResult<NumericBandOperations> {
    let representations = exact_matrices(context, required(model, "representation_matrices")?)?;
    if representations.len() != model_operations.len() {
        return Err(ExactError::new(
            "InvalidBandCorepRepresentation",
            "representation matrices must align with MSGCorep operations",
        ));
    }
    let lattice = encoded_matrix_from_json(context, required(model, "lattice")?)?;
    let soc = required(request, "soc")?
        .as_bool()
        .ok_or_else(|| malformed("soc must be bool"))?;
    let canonical_spin = band_spin_rotations(request, model_operations.len(), soc)?;
    let mut tb_spin_rotations = Vec::new();
    let operations = model_operations
        .iter()
        .enumerate()
        .map(|(index, operation)| {
            let representation = exact_numeric_matrix(
                &representations[index],
                &format!("representation {}", index + 1),
            )?;
            let tb_spin = if soc {
                Some(exact_numeric_matrix(
                    &compile_spatial_spinor_matrix(operation.rotation(), &lattice)
                        .map_err(|error| ExactError::new(error.tag(), error.to_string()))?,
                    &format!("MagneticTB spin rotation {}", index + 1),
                )?)
            } else {
                None
            };
            let sign = if let (Some(canonical), Some(tb_spin)) = (
                canonical_spin.as_ref().map(|items| &items[index]),
                tb_spin.as_ref(),
            ) {
                trace_sign(canonical, tb_spin)?
            } else {
                1.0
            };
            if let Some(tb_spin) = tb_spin {
                tb_spin_rotations.push(tb_spin);
            }
            Ok(NumericSymmetryOperation {
                rotation: real_matrix_three(operation.rotation(), "spatial rotation")?,
                translation: real_vector_three(operation.translation(), "spatial translation")?,
                antiunitary: model_antiunitary[index],
                representation,
                trace_sign: sign,
            })
        })
        .collect::<ExactResult<Vec<_>>>()?;
    Ok((operations, tb_spin_rotations))
}

fn band_corep_trace(request: &serde_json::Map<String, Value>) -> ExactResult<Value> {
    let hamiltonian = object(required(request, "hamiltonian")?)?;
    let context = context_from_json(required(hamiltonian, "field_context")?, 128)?;
    let model = object(required(request, "model_identity")?)?;
    let model_context = context_from_json(required(model, "field_context")?, 128)?;
    if context != model_context {
        return Err(ExactError::new(
            "BandCorepModelMismatch",
            "Hamiltonian and prepared model use different exact fields",
        ));
    }
    let (model_operations, model_antiunitary) = aligned_band_operations(&context, model, request)?;
    let reduced_kpoints = band_kpoints(&context, request)?;
    let hamiltonians = band_hamiltonians(hamiltonian, request, &reduced_kpoints)?;
    let orbital_positions = band_orbital_positions(&context, model)?;
    let (operations, tb_spin_rotations) = band_numeric_operations(
        &context,
        model,
        request,
        &model_operations,
        &model_antiunitary,
    )?;
    let traces = compute_band_trace(&BandTraceInput {
        hamiltonians: &hamiltonians,
        reduced_kpoints: &reduced_kpoints,
        orbital_positions: &orbital_positions,
        operations: &operations,
    })?;
    Ok(json!({
        "schema": "magnetictb.band_corep_trace.v1",
        "energies": traces.iter().map(|point| &point.energies).collect::<Vec<_>>(),
        "degeneracies": traces.iter().map(|point| &point.degeneracies).collect::<Vec<_>>(),
        "little_symmetry_counts": traces.iter().map(|point| point.little_symmetry_indices.len()).collect::<Vec<_>>(),
        "little_symmetry_indices": traces.iter().map(|point| &point.little_symmetry_indices).collect::<Vec<_>>(),
        "traces": traces.iter().map(|point| {
            point.traces.iter().map(|band| {
                band.iter().map(|value| json!({"real": value.re, "imaginary": value.im})).collect::<Vec<_>>()
            }).collect::<Vec<_>>()
        }).collect::<Vec<_>>(),
        "spin_rotations": tb_spin_rotations.iter().map(numeric_matrix_json).collect::<Vec<_>>()
    }))
}
