use magnetictb_properties::{
    BerryCurvatureOptions, CellMatrix, Complex, CubeSurfaceMesh, DMatrix, FiniteGeometry,
    GaplessSearchData, GaplessSearchOptions, HoppingData, PointChernOptions, PropertiesError,
    PropertiesResult, RefinementMethod, SurfaceOptions, SurfaceSide, Translation,
    WilsonLoopOptions, berry_curvature_from_samples, berry_phase, berry_plaquette_path,
    build_bloch_hamiltonian, build_real_space_hamiltonian, build_slab_hamiltonian,
    cube_surface_mesh, format_wannier90_hr, hermitian_eigenvalues, legacy_berryph, legacy_wloop,
    legacy_z2_path, parse_wannier90_hr, point_chern_number_from_samples, surface_green_function,
    transform_hoppings, wilson_loop,
};
use serde_json::{Map, Value, json};

fn malformed(detail: impl Into<String>) -> PropertiesError {
    PropertiesError::new("MalformedPropertiesInput", detail)
}

fn object(value: &Value) -> PropertiesResult<&Map<String, Value>> {
    value
        .as_object()
        .ok_or_else(|| malformed("expected a JSON object"))
}

fn array<'a>(value: &'a Value, field: &str) -> PropertiesResult<&'a [Value]> {
    value
        .as_array()
        .map(Vec::as_slice)
        .ok_or_else(|| malformed(format!("{field} must be an array")))
}

fn required<'a>(record: &'a Map<String, Value>, field: &str) -> PropertiesResult<&'a Value> {
    record
        .get(field)
        .ok_or_else(|| malformed(format!("missing required field {field}")))
}

fn finite_real(value: &Value, field: &str) -> PropertiesResult<f64> {
    let number = value
        .as_f64()
        .ok_or_else(|| malformed(format!("{field} must be a real JSON number")))?;
    if number.is_finite() {
        Ok(number)
    } else {
        Err(malformed(format!("{field} must be finite")))
    }
}

fn usize_value(value: &Value, field: &str) -> PropertiesResult<usize> {
    let number = value
        .as_u64()
        .and_then(|number| usize::try_from(number).ok())
        .ok_or_else(|| malformed(format!("{field} must be a nonnegative integer")))?;
    Ok(number)
}

fn i64_value(value: &Value, field: &str) -> PropertiesResult<i64> {
    value
        .as_i64()
        .ok_or_else(|| malformed(format!("{field} must be an integer in the i64 range")))
}

fn u64_value(value: &Value, field: &str) -> PropertiesResult<u64> {
    value
        .as_u64()
        .ok_or_else(|| malformed(format!("{field} must be a nonnegative integer")))
}

fn real_vector(value: &Value, field: &str) -> PropertiesResult<Vec<f64>> {
    array(value, field)?
        .iter()
        .enumerate()
        .map(|(index, value)| finite_real(value, &format!("{field}[{index}]")))
        .collect()
}

fn real_vectors(value: &Value, field: &str) -> PropertiesResult<Vec<Vec<f64>>> {
    array(value, field)?
        .iter()
        .enumerate()
        .map(|(index, value)| real_vector(value, &format!("{field}[{index}]")))
        .collect()
}

fn complex_value(value: &Value, field: &str) -> PropertiesResult<Complex<f64>> {
    if value.is_number() {
        return finite_real(value, field).map(|real| Complex::new(real, 0.0));
    }
    let record = object(value)?;
    Ok(Complex::new(
        finite_real(required(record, "real")?, &format!("{field}.real"))?,
        finite_real(
            required(record, "imaginary")?,
            &format!("{field}.imaginary"),
        )?,
    ))
}

fn complex_matrix(value: &Value, field: &str) -> PropertiesResult<DMatrix<Complex<f64>>> {
    let rows = array(value, field)?;
    if rows.is_empty() {
        return Err(malformed(format!("{field} must be nonempty")));
    }
    let columns = array(&rows[0], &format!("{field}[0]"))?.len();
    if columns == 0 {
        return Err(malformed(format!("{field} rows must be nonempty")));
    }
    let mut entries = Vec::with_capacity(rows.len() * columns);
    for (row_index, row) in rows.iter().enumerate() {
        let row = array(row, &format!("{field}[{row_index}]"))?;
        if row.len() != columns {
            return Err(malformed(format!("{field} rows must have equal length")));
        }
        for (column_index, value) in row.iter().enumerate() {
            entries.push(complex_value(
                value,
                &format!("{field}[{row_index}][{column_index}]"),
            )?);
        }
    }
    Ok(DMatrix::from_row_slice(rows.len(), columns, &entries))
}

fn complex_matrices(value: &Value, field: &str) -> PropertiesResult<Vec<DMatrix<Complex<f64>>>> {
    array(value, field)?
        .iter()
        .enumerate()
        .map(|(index, value)| complex_matrix(value, &format!("{field}[{index}]")))
        .collect()
}

fn complex_json(value: Complex<f64>) -> Value {
    json!({"real": value.re, "imaginary": value.im})
}

fn matrix_json(matrix: &DMatrix<Complex<f64>>) -> Value {
    Value::Array(
        (0..matrix.nrows())
            .map(|row| {
                Value::Array(
                    (0..matrix.ncols())
                        .map(|column| complex_json(matrix[(row, column)]))
                        .collect(),
                )
            })
            .collect(),
    )
}

fn translation(value: &Value, field: &str) -> PropertiesResult<Translation> {
    let values = array(value, field)?;
    if values.len() != 3 {
        return Err(malformed(format!("{field} must have length three")));
    }
    Ok([
        i64_value(&values[0], &format!("{field}[0]"))?,
        i64_value(&values[1], &format!("{field}[1]"))?,
        i64_value(&values[2], &format!("{field}[2]"))?,
    ])
}

fn translation_list(value: &Value, field: &str) -> PropertiesResult<Vec<Translation>> {
    array(value, field)?
        .iter()
        .enumerate()
        .map(|(index, value)| translation(value, &format!("{field}[{index}]")))
        .collect()
}

fn real_three(value: &Value, field: &str) -> PropertiesResult<[f64; 3]> {
    let values = real_vector(value, field)?;
    if values.len() != 3 {
        return Err(malformed(format!("{field} must have length three")));
    }
    Ok([values[0], values[1], values[2]])
}

fn real_three_matrix(value: &Value, field: &str) -> PropertiesResult<[[f64; 3]; 3]> {
    let rows = array(value, field)?;
    if rows.len() != 3 {
        return Err(malformed(format!("{field} must have three rows")));
    }
    Ok([
        real_three(&rows[0], &format!("{field}[0]"))?,
        real_three(&rows[1], &format!("{field}[1]"))?,
        real_three(&rows[2], &format!("{field}[2]"))?,
    ])
}

fn cell_matrix(value: &Value, field: &str) -> PropertiesResult<CellMatrix> {
    let rows = array(value, field)?;
    if rows.len() != 3 {
        return Err(malformed(format!("{field} must have three rows")));
    }
    let mut result = [[0; 3]; 3];
    for (row_index, row) in rows.iter().enumerate() {
        let row = array(row, &format!("{field}[{row_index}]"))?;
        if row.len() != 3 {
            return Err(malformed(format!(
                "{field}[{row_index}] must have length three"
            )));
        }
        for (column_index, value) in row.iter().enumerate() {
            result[row_index][column_index] =
                i64_value(value, &format!("{field}[{row_index}][{column_index}]"))?;
        }
    }
    Ok(result)
}

fn hopping_data(value: &Value) -> PropertiesResult<HoppingData> {
    let record = object(value)?;
    let num_wannier = usize_value(required(record, "NumWannier")?, "NumWannier")?;
    let translations = translation_list(required(record, "Translations")?, "Translations")?;
    let degeneracies = array(required(record, "Degeneracies")?, "Degeneracies")?
        .iter()
        .enumerate()
        .map(|(index, value)| u64_value(value, &format!("Degeneracies[{index}]")))
        .collect::<PropertiesResult<Vec<_>>>()?;
    let matrices = complex_matrices(required(record, "HoppingMatrices")?, "HoppingMatrices")?;
    let lattice = record
        .get("Lattice")
        .map(|value| real_three_matrix(value, "Lattice"))
        .transpose()?;
    let wannier_centers = record
        .get("WannierCenters")
        .map(|value| {
            array(value, "WannierCenters")?
                .iter()
                .enumerate()
                .map(|(index, value)| real_three(value, &format!("WannierCenters[{index}]")))
                .collect::<PropertiesResult<Vec<_>>>()
        })
        .transpose()?;
    let cell = record
        .get("CellMatrix")
        .map(|value| cell_matrix(value, "CellMatrix"))
        .transpose()?;
    let representatives = record
        .get("CellRepresentatives")
        .map(|value| translation_list(value, "CellRepresentatives"))
        .transpose()?
        .unwrap_or_default();
    let data = HoppingData {
        num_wannier,
        translations,
        degeneracies,
        matrices,
        lattice,
        wannier_centers,
        source: record
            .get("Source")
            .and_then(Value::as_str)
            .map(str::to_owned),
        cell_matrix: cell,
        cell_representatives: representatives,
    };
    data.validate()?;
    Ok(data)
}

fn hopping_data_json(data: &HoppingData) -> PropertiesResult<Value> {
    let residual = data.hermiticity_residual()?;
    let mut record = Map::from_iter([
        ("Schema".to_owned(), json!("MagneticTBHoppingData")),
        ("SchemaVersion".to_owned(), json!(1)),
        ("NumWannier".to_owned(), json!(data.num_wannier)),
        ("NumTranslations".to_owned(), json!(data.translations.len())),
        ("Translations".to_owned(), json!(data.translations)),
        ("Degeneracies".to_owned(), json!(data.degeneracies)),
        (
            "HoppingMatrices".to_owned(),
            Value::Array(data.matrices.iter().map(matrix_json).collect()),
        ),
        ("HermitianResidual".to_owned(), json!(residual)),
        (
            "Convention".to_owned(),
            json!("H(R)[a,b]=<0,a|H|R,b>; H(k)=Sum_R H(R) Exp[2 Pi I k.R]"),
        ),
    ]);
    if let Some(source) = &data.source {
        record.insert("Source".to_owned(), json!(source));
    }
    if let Some(lattice) = data.lattice {
        record.insert("Lattice".to_owned(), json!(lattice));
    }
    if let Some(centers) = &data.wannier_centers {
        record.insert("WannierCenters".to_owned(), json!(centers));
    }
    if let Some(cell) = data.cell_matrix {
        record.insert("CellMatrix".to_owned(), json!(cell));
        let determinant = i128::from(cell[0][0])
            * (i128::from(cell[1][1]) * i128::from(cell[2][2])
                - i128::from(cell[1][2]) * i128::from(cell[2][1]))
            - i128::from(cell[0][1])
                * (i128::from(cell[1][0]) * i128::from(cell[2][2])
                    - i128::from(cell[1][2]) * i128::from(cell[2][0]))
            + i128::from(cell[0][2])
                * (i128::from(cell[1][0]) * i128::from(cell[2][1])
                    - i128::from(cell[1][1]) * i128::from(cell[2][0]));
        record.insert(
            "CellVolumeFactor".to_owned(),
            json!(determinant.unsigned_abs()),
        );
    }
    if !data.cell_representatives.is_empty() {
        record.insert(
            "CellRepresentatives".to_owned(),
            json!(data.cell_representatives),
        );
    }
    Ok(Value::Object(record))
}

fn optional_tolerance(
    record: &Map<String, Value>,
    field: &str,
    default: f64,
) -> PropertiesResult<f64> {
    record
        .get(field)
        .map_or(Ok(default), |value| finite_real(value, field))
}

fn bool_value(value: &Value, field: &str) -> PropertiesResult<bool> {
    value
        .as_bool()
        .ok_or_else(|| malformed(format!("{field} must be boolean")))
}

fn point_chern_options(record: &Map<String, Value>) -> PropertiesResult<PointChernOptions> {
    let defaults = PointChernOptions::default();
    Ok(PointChernOptions {
        subdivisions: record
            .get("surface_subdivisions")
            .map_or(Ok(defaults.subdivisions), |value| {
                usize_value(value, "surface_subdivisions")
            })?,
        hermitian_tolerance: optional_tolerance(
            record,
            "hermitian_tolerance",
            defaults.hermitian_tolerance,
        )?,
        surface_gap_tolerance: optional_tolerance(
            record,
            "surface_gap_tolerance",
            defaults.surface_gap_tolerance,
        )?,
        center_gap_tolerance: optional_tolerance(
            record,
            "center_gap_tolerance",
            defaults.center_gap_tolerance,
        )?,
        overlap_tolerance: optional_tolerance(
            record,
            "overlap_tolerance",
            defaults.overlap_tolerance,
        )?,
        integer_tolerance: optional_tolerance(
            record,
            "integer_tolerance",
            defaults.integer_tolerance,
        )?,
        require_gapless_center: record
            .get("require_gapless_center")
            .map_or(Ok(defaults.require_gapless_center), |value| {
                bool_value(value, "require_gapless_center")
            })?,
    })
}

fn transformed_input(
    request: &Map<String, Value>,
    tolerance: f64,
) -> PropertiesResult<HoppingData> {
    let data = hopping_data(required(request, "data")?)?;
    request
        .get("cell_matrix")
        .map(|cell| transform_hoppings(&data, cell_matrix(cell, "cell_matrix")?, tolerance))
        .unwrap_or(Ok(data))
}

fn options(record: &Map<String, Value>) -> PropertiesResult<WilsonLoopOptions> {
    let defaults = WilsonLoopOptions::default();
    Ok(WilsonLoopOptions {
        hermitian_tolerance: record
            .get("hermitian_tolerance")
            .map_or(Ok(defaults.hermitian_tolerance), |value| {
                finite_real(value, "hermitian_tolerance")
            })?,
        gap_tolerance: record
            .get("gap_tolerance")
            .map_or(Ok(defaults.gap_tolerance), |value| {
                finite_real(value, "gap_tolerance")
            })?,
        covariance_tolerance: record
            .get("covariance_tolerance")
            .map_or(Ok(defaults.covariance_tolerance), |value| {
                finite_real(value, "covariance_tolerance")
            })?,
        overlap_tolerance: record
            .get("overlap_tolerance")
            .map_or(Ok(defaults.overlap_tolerance), |value| {
                finite_real(value, "overlap_tolerance")
            })?,
    })
}

fn wilson_data_json(data: &magnetictb_properties::WilsonLoopData) -> Value {
    json!({
        "Schema": "MagneticTBWilsonLoop",
        "SchemaVersion": 1,
        "Eigenvalues": data.eigenvalues.iter().map(|&(real, imaginary)| {
            complex_json(Complex::new(real, imaginary))
        }).collect::<Vec<_>>(),
        "PhasesOverPi": data.phases_over_pi,
        "WannierCenters": data.wannier_centers,
        "WilsonMatrix": matrix_json(&data.wilson_matrix),
        "UnitarityResidual": data.unitarity_residual,
        "MinimumLinkSingularValue": data.minimum_link_singular_value,
        "MinimumDirectGap": data.minimum_direct_gap,
        "Path": data.path,
        "PathSubdivisions": data.path_subdivisions,
        "OccupiedBands": data.occupied_bands,
        "ClosureVector": data.closure_vector,
        "ReciprocalCoordinates": data.reciprocal_coordinates,
        "EndpointCovarianceResidual": data.endpoint_covariance_residual
    })
}

fn berry_phase_json(data: &magnetictb_properties::BerryPhaseData) -> Value {
    json!({
        "Schema": "MagneticTBBerryPhase",
        "SchemaVersion": 1,
        "Phase": data.phase,
        "PhaseOverPi": data.phase_over_pi,
        "WilsonDeterminant": complex_json(Complex::new(
            data.wilson_determinant.0,
            data.wilson_determinant.1,
        )),
        "WilsonMatrix": matrix_json(&data.wilson_matrix),
        "UnitarityResidual": data.unitarity_residual,
        "MinimumLinkSingularValue": data.minimum_link_singular_value,
        "MinimumDirectGap": data.minimum_direct_gap,
        "Path": data.path,
        "OccupiedBands": data.occupied_bands,
        "ClosureVector": data.closure_vector,
        "ReciprocalCoordinates": data.reciprocal_coordinates,
        "EndpointCovarianceResidual": data.endpoint_covariance_residual,
        "Convention": "occupied-subspace links follow the forward path; phase=Arg Det[WilsonMatrix] in radians"
    })
}

#[allow(clippy::too_many_lines)]
fn compute_value(payload: &str) -> PropertiesResult<Value> {
    let value: Value = serde_json::from_str(payload)
        .map_err(|error| malformed(format!("invalid JSON: {error}")))?;
    let request = object(&value)?;
    let operation = required(request, "operation")?
        .as_str()
        .ok_or_else(|| malformed("operation must be a string"))?;
    let option_record = request
        .get("options")
        .map(object)
        .transpose()?
        .cloned()
        .unwrap_or_default();
    let wilson_options = options(&option_record)?;
    if operation == "format_wannier90_hr" {
        let data = hopping_data(required(request, "data")?)?;
        let digits = request
            .get("real_digits")
            .map_or(Ok(12), |value| usize_value(value, "real_digits"))?;
        return Ok(json!({"Text": format_wannier90_hr(&data, digits)?}));
    }
    if operation == "parse_wannier90_hr" {
        let text = required(request, "text")?
            .as_str()
            .ok_or_else(|| malformed("text must be a string"))?;
        let source = request
            .get("source")
            .and_then(Value::as_str)
            .map(str::to_owned);
        let precision = optional_tolerance(request, "precision_cutoff", 1.0e-10)?;
        let cutoff = request
            .get("cell_cutoff")
            .filter(|value| !value.is_null())
            .map(|value| i64_value(value, "cell_cutoff"))
            .transpose()?;
        let parsed = parse_wannier90_hr(text, source, precision, cutoff)?;
        let mut result = hopping_data_json(&parsed.data)?;
        let result = result
            .as_object_mut()
            .expect("hopping data JSON is an object");
        result.insert("Schema".to_owned(), json!("Wannier90HRData"));
        result.insert("Header".to_owned(), json!(parsed.header));
        result.insert(
            "CellCutoff".to_owned(),
            cutoff.map_or(Value::String("All".to_owned()), Value::from),
        );
        result.insert("PrecisionCutoff".to_owned(), json!(precision));
        return Ok(Value::Object(result.clone()));
    }
    if operation == "band_eigenvalues" {
        let tolerance = optional_tolerance(request, "hermitian_tolerance", 1.0e-10)?;
        let matrices = complex_matrices(required(request, "matrices")?, "matrices")?;
        let values = matrices
            .iter()
            .map(|matrix| hermitian_eigenvalues(matrix, tolerance))
            .collect::<PropertiesResult<Vec<_>>>()?;
        return Ok(json!({"Eigenvalues": values}));
    }
    if operation == "legacy_z2_path" || operation == "legacy_berryph" {
        let tolerance = optional_tolerance(request, "hermitian_tolerance", 1.0e-10)?;
        let matrices = complex_matrices(required(request, "matrices")?, "matrices")?;
        let occupied = usize_value(required(request, "occupied")?, "occupied")?;
        let values = if operation == "legacy_z2_path" {
            legacy_z2_path(&matrices, occupied, tolerance)?
        } else {
            legacy_berryph(&matrices, occupied, tolerance)?
        };
        return Ok(json!({"Values": values}));
    }
    if operation == "legacy_wloop" {
        let tolerance = optional_tolerance(request, "hermitian_tolerance", 1.0e-10)?;
        let occupied = usize_value(required(request, "occupied")?, "occupied")?;
        let matrix_grid = array(required(request, "matrix_grid")?, "matrix_grid")?
            .iter()
            .enumerate()
            .map(|(index, value)| complex_matrices(value, &format!("matrix_grid[{index}]")))
            .collect::<PropertiesResult<Vec<_>>>()?;
        return Ok(json!({
            "Values": legacy_wloop(&matrix_grid, occupied, tolerance)?
        }));
    }
    if operation == "transform_hoppings" {
        let tolerance = optional_tolerance(request, "hermiticity_tolerance", 1.0e-9)?;
        let data = hopping_data(required(request, "data")?)?;
        let transformed = transform_hoppings(
            &data,
            cell_matrix(required(request, "cell_matrix")?, "cell_matrix")?,
            tolerance,
        )?;
        return hopping_data_json(&transformed);
    }
    if operation == "build_bloch_hamiltonian" {
        let tolerance = optional_tolerance(request, "hermiticity_tolerance", 1.0e-9)?;
        let data = transformed_input(request, tolerance)?;
        let momentum = real_three(required(request, "momentum")?, "momentum")?;
        let result = build_bloch_hamiltonian(&data, momentum, tolerance)?;
        return Ok(json!({
            "Schema": "MagneticTBBlochHamiltonian",
            "SchemaVersion": 1,
            "Hamiltonian": matrix_json(&result.hamiltonian),
            "Dimension": result.hamiltonian.nrows(),
            "CrystalMomentum": result.momentum,
            "CellMatrix": request.get("cell_matrix").cloned().unwrap_or(Value::Null),
            "CellVolumeFactor": result.cell_volume_factor,
            "CellRepresentatives": result.cell_representatives,
            "HermitianResidual": result.hermitian_residual,
            "Convention": "H(k)=Sum_R H(R) Exp[2 Pi I k.R], with k in reciprocal fractional coordinates of the active cell"
        }));
    }
    if operation == "build_real_space_hamiltonian" {
        let tolerance = optional_tolerance(request, "hermiticity_tolerance", 1.0e-9)?;
        let data = hopping_data(required(request, "data")?)?;
        let geometry_value = required(request, "geometry")?;
        let geometry_values = array(geometry_value, "geometry")?;
        let geometry = if geometry_values.len() == 3
            && geometry_values.iter().all(|value| value.as_u64().is_some())
        {
            FiniteGeometry::Rectangular([
                usize_value(&geometry_values[0], "geometry[0]")?,
                usize_value(&geometry_values[1], "geometry[1]")?,
                usize_value(&geometry_values[2], "geometry[2]")?,
            ])
        } else {
            FiniteGeometry::ExplicitCells(translation_list(geometry_value, "geometry")?)
        };
        let boundary_values = array(
            required(request, "boundary_conditions")?,
            "boundary_conditions",
        )?;
        if boundary_values.len() != 3 {
            return Err(malformed("boundary_conditions must have length three"));
        }
        let mut boundary = [false; 3];
        for (index, value) in boundary_values.iter().enumerate() {
            boundary[index] = match value.as_str() {
                Some("Open") => false,
                Some("Periodic") => true,
                _ => return Err(malformed("boundary conditions must be Open or Periodic")),
            };
        }
        let result = build_real_space_hamiltonian(&data, geometry, boundary, tolerance)?;
        let basis_records = result
            .basis_records
            .iter()
            .map(|record| {
                json!({
                    "BasisIndex": record.basis_index,
                    "CellIndex": record.cell_index,
                    "Cell": record.cell,
                    "OrbitalIndex": record.orbital_index,
                    "WannierCenter": record.wannier_center,
                    "FractionalPosition": record.fractional_position,
                    "CartesianPosition": record.cartesian_position,
                })
            })
            .collect::<Vec<_>>();
        let mut output = json!({
            "Schema": "MagneticTBFiniteRealSpaceHamiltonian",
            "SchemaVersion": 2,
            "Hamiltonian": matrix_json(&result.hamiltonian),
            "Shape": result.shape,
            "Cells": result.cells,
            "NumCells": result.cells.len(),
            "OrbitalsPerCell": result.orbitals_per_cell,
            "Dimension": result.hamiltonian.nrows(),
            "BoundaryConditions": result.boundary_conditions.map(|periodic| if periodic {"Periodic"} else {"Open"}),
            "Size": result.size,
            "CellBounds": result.cell_bounds,
            "HermitianResidual": result.hermitian_residual,
            "Convention": "matrix row=(source cell,orbital), column=(target cell,orbital)",
            "PositionDataAvailable": result.position_data_available,
            "BasisRecords": basis_records,
        });
        if result.position_data_available {
            let output = output
                .as_object_mut()
                .ok_or_else(|| malformed("real-space output must be an object"))?;
            output.insert("Lattice".to_owned(), json!(result.lattice));
            output.insert("WannierCenters".to_owned(), json!(result.wannier_centers));
        }
        return Ok(output);
    }
    if operation == "build_slab_hamiltonian" {
        let tolerance = optional_tolerance(request, "hermiticity_tolerance", 1.0e-9)?;
        let data = transformed_input(request, tolerance)?;
        let size_value = array(required(request, "size")?, "size")?;
        if size_value.len() != 3 {
            return Err(malformed("size must have length three"));
        }
        let size = [
            usize_value(&size_value[0], "size[0]")?,
            usize_value(&size_value[1], "size[1]")?,
            usize_value(&size_value[2], "size[2]")?,
        ];
        let directions = array(
            required(request, "periodic_directions")?,
            "periodic_directions",
        )?
        .iter()
        .enumerate()
        .map(|(index, value)| {
            let one_based = usize_value(value, &format!("periodic_directions[{index}]"))?;
            one_based
                .checked_sub(1)
                .ok_or_else(|| malformed("periodic directions are one-based"))
        })
        .collect::<PropertiesResult<Vec<_>>>()?;
        let momentum = real_vector(required(request, "momentum")?, "momentum")?;
        let result = build_slab_hamiltonian(&data, size, &directions, &momentum, tolerance)?;
        return Ok(json!({
            "Schema": "MagneticTBHybridSpaceHamiltonian",
            "SchemaVersion": 1,
            "Hamiltonian": matrix_json(&result.hamiltonian),
            "Cells": result.cells,
            "NumFiniteCells": result.cells.len(),
            "OrbitalsPerCell": result.orbitals_per_cell,
            "Dimension": result.hamiltonian.nrows(),
            "Size": result.size,
            "PeriodicDirections": result.periodic_directions.iter().map(|index| index + 1).collect::<Vec<_>>(),
            "OpenDirections": result.open_directions.iter().map(|index| index + 1).collect::<Vec<_>>(),
            "CrystalMomentum": result.momentum,
            "EmbeddedCrystalMomentum": result.embedded_momentum,
            "CellMatrix": request.get("cell_matrix").cloned().unwrap_or(Value::Null),
            "HermitianResidual": result.hermitian_residual,
            "Convention": "H(R)[a,b]=<0,a|H|R,b>; periodic hopping carries Exp[2 Pi I k.R], while nonperiodic directions use finite open cells"
        }));
    }
    if operation == "surface_green_function" {
        let hermiticity_tolerance = optional_tolerance(request, "hermiticity_tolerance", 1.0e-9)?;
        let data = transformed_input(request, hermiticity_tolerance)?;
        let momentum_values = real_vector(required(request, "momentum")?, "momentum")?;
        if momentum_values.len() != 2 {
            return Err(malformed("surface momentum must have length two"));
        }
        let energy = finite_real(required(request, "energy")?, "energy")?;
        let side = match request
            .get("surface")
            .and_then(Value::as_str)
            .unwrap_or("Positive")
        {
            "Positive" => SurfaceSide::Positive,
            "Negative" => SurfaceSide::Negative,
            _ => return Err(malformed("surface must be Positive or Negative")),
        };
        let surface_options = SurfaceOptions {
            broadening: optional_tolerance(request, "broadening", 1.0e-3)?,
            tolerance: optional_tolerance(request, "tolerance", 1.0e-10)?,
            max_iterations: request
                .get("max_iterations")
                .map_or(Ok(200), |value| usize_value(value, "max_iterations"))?,
            side,
            hermiticity_tolerance,
        };
        let result = surface_green_function(
            &data,
            [momentum_values[0], momentum_values[1]],
            energy,
            surface_options,
        )?;
        return Ok(json!({
            "Schema": "MagneticTBSurfaceGreenFunction",
            "SchemaVersion": 1,
            "GreenFunction": matrix_json(&result.green_function),
            "SpectralWeight": result.spectral_weight,
            "Energy": result.energy,
            "Broadening": surface_options.broadening,
            "SurfaceMomentum": result.momentum,
            "Surface": match surface_options.side { SurfaceSide::Positive => "Positive", SurfaceSide::Negative => "Negative" },
            "PrincipalLayerThickness": result.blocks.thickness,
            "OrbitalsPerCell": result.blocks.orbitals_per_cell,
            "PrincipalLayerDimension": result.blocks.h00.nrows(),
            "H00": matrix_json(&result.blocks.h00),
            "H01": matrix_json(&result.blocks.h01),
            "H10": matrix_json(&result.blocks.h10),
            "Iterations": result.iterations,
            "CouplingResidual": result.coupling_residual,
            "DysonResidual": result.dyson_residual,
            "HoppingHermitianResidual": result.hopping_hermitian_residual,
            "BlockHermitianResidual": result.block_hermitian_residual,
            "Convention": "surface plane=cell axes 1,2; stacking direction=cell axis 3"
        }));
    }
    if operation == "point_chern_mesh" {
        let point = real_three(required(request, "point")?, "point")?;
        let radius = finite_real(required(request, "radius")?, "radius")?;
        let chern_options = point_chern_options(request)?;
        let mesh = cube_surface_mesh(point, radius, chern_options.subdivisions)?;
        return Ok(json!({
            "Points": mesh.points,
            "Triangles": mesh.triangles,
            "SurfaceSubdivisions": chern_options.subdivisions,
        }));
    }
    if operation == "point_chern_number" {
        let point = real_three(required(request, "point")?, "point")?;
        let radius = finite_real(required(request, "radius")?, "radius")?;
        let occupied = usize_value(required(request, "occupied")?, "occupied")?;
        let chern_options = point_chern_options(request)?;
        let points = array(required(request, "points")?, "points")?
            .iter()
            .enumerate()
            .map(|(index, value)| real_three(value, &format!("points[{index}]")))
            .collect::<PropertiesResult<Vec<_>>>()?;
        let triangles = array(required(request, "triangles")?, "triangles")?
            .iter()
            .enumerate()
            .map(|(index, value)| {
                let values = array(value, &format!("triangles[{index}]"))?;
                if values.len() != 3 {
                    return Err(malformed(format!(
                        "triangles[{index}] must have length three"
                    )));
                }
                Ok([
                    usize_value(&values[0], &format!("triangles[{index}][0]"))?,
                    usize_value(&values[1], &format!("triangles[{index}][1]"))?,
                    usize_value(&values[2], &format!("triangles[{index}][2]"))?,
                ])
            })
            .collect::<PropertiesResult<Vec<_>>>()?;
        let mesh = CubeSurfaceMesh { points, triangles };
        let center_matrix = complex_matrix(required(request, "center_matrix")?, "center_matrix")?;
        let surface_matrices = complex_matrices(required(request, "matrices")?, "matrices")?;
        let result = point_chern_number_from_samples(
            &center_matrix,
            &surface_matrices,
            &mesh,
            point,
            radius,
            occupied,
            chern_options,
        )?;
        return Ok(json!({
            "Schema": "MagneticTBPointChernNumber",
            "SchemaVersion": 1,
            "ChernNumber": result.chern_number,
            "RawChernNumber": result.raw_chern_number,
            "QuantizationError": result.quantization_error,
            "Point": result.point,
            "Radius": result.radius,
            "OccupiedBands": result.occupied_bands,
            "HamiltonianDimension": result.hamiltonian_dimension,
            "CenterGap": result.center_gap,
            "RequireGaplessCenter": result.require_gapless_center,
            "MinimumSurfaceGap": result.minimum_surface_gap,
            "MinimumOverlapSingularValue": result.minimum_overlap_singular_value,
            "SurfaceSubdivisions": result.surface_subdivisions,
            "TriangleCount": result.triangle_count,
            "UniqueVertexCount": result.unique_vertex_count,
            "Method": "gauge-invariant occupied-subspace flux on an oriented cube mesh"
        }));
    }
    if operation == "berry_plaquette_path" {
        let point = real_vector(required(request, "point")?, "point")?;
        let directions = array(required(request, "directions")?, "directions")?;
        let step = array(required(request, "step_size")?, "step_size")?;
        if directions.len() != 2 || step.len() != 2 {
            return Err(malformed("directions and step_size must have length two"));
        }
        let curvature_options = BerryCurvatureOptions {
            directions: [
                usize_value(&directions[0], "directions[0]")?,
                usize_value(&directions[1], "directions[1]")?,
            ],
            step_size: [
                finite_real(&step[0], "step_size[0]")?,
                finite_real(&step[1], "step_size[1]")?,
            ],
            wilson: wilson_options,
        };
        return berry_plaquette_path(&point, curvature_options)
            .map(|path| json!({"PlaquettePath": path}));
    }
    let matrices = complex_matrices(required(request, "matrices")?, "matrices")?;
    let centers = real_vectors(required(request, "centers")?, "centers")?;
    let occupied = usize_value(required(request, "occupied")?, "occupied")?;
    match operation {
        "wilson_loop" => {
            let path = real_vectors(required(request, "path")?, "path")?;
            wilson_loop(&matrices, &path, &centers, occupied, wilson_options)
                .map(|data| wilson_data_json(&data))
        }
        "berry_phase" => {
            let path = real_vectors(required(request, "path")?, "path")?;
            berry_phase(&matrices, &path, &centers, occupied, wilson_options)
                .map(|data| berry_phase_json(&data))
        }
        "berry_curvature" => {
            let point = real_vector(required(request, "point")?, "point")?;
            let directions = array(required(request, "directions")?, "directions")?;
            let step = array(required(request, "step_size")?, "step_size")?;
            if directions.len() != 2 || step.len() != 2 {
                return Err(malformed("directions and step_size must have length two"));
            }
            let curvature_options = BerryCurvatureOptions {
                directions: [
                    usize_value(&directions[0], "directions[0]")?,
                    usize_value(&directions[1], "directions[1]")?,
                ],
                step_size: [
                    finite_real(&step[0], "step_size[0]")?,
                    finite_real(&step[1], "step_size[1]")?,
                ],
                wilson: wilson_options,
            };
            let data = berry_curvature_from_samples(
                &matrices,
                &centers,
                occupied,
                &point,
                curvature_options,
            )?;
            Ok(json!({
                "Schema": "MagneticTBBerryCurvature",
                "SchemaVersion": 1,
                "Curvature": data.curvature,
                "Flux": data.flux,
                "Area": data.area,
                "Point": data.point,
                "Directions": data.directions.iter().map(|index| index + 1).collect::<Vec<_>>(),
                "StepSize": data.step_size,
                "PlaquettePath": data.plaquette_path,
                "PhaseData": berry_phase_json(&data.phase_data),
                "Convention": "positive orientation follows Directions[[1]] then Directions[[2]]"
            }))
        }
        _ => Err(PropertiesError::new(
            "UnknownPropertiesOperation",
            format!("unknown properties operation {operation}"),
        )),
    }
}

pub(crate) fn compute(payload: &str) -> PropertiesResult<String> {
    serde_json::to_string(&compute_value(payload)?)
        .map_err(|error| malformed(format!("could not serialize result: {error}")))
}

pub(crate) fn parse_gapless_request(
    payload: &str,
) -> PropertiesResult<(usize, GaplessSearchOptions)> {
    let value: Value = serde_json::from_str(payload)
        .map_err(|error| malformed(format!("invalid JSON: {error}")))?;
    let record = object(&value)?;
    let occupied = usize_value(required(record, "occupied")?, "occupied")?;
    let zone = array(required(record, "brillouin_zone")?, "brillouin_zone")?
        .iter()
        .enumerate()
        .map(|(index, value)| {
            let interval = real_vector(value, &format!("brillouin_zone[{index}]"))?;
            if interval.len() != 2 {
                return Err(malformed(format!(
                    "brillouin_zone[{index}] must have length two"
                )));
            }
            Ok([interval[0], interval[1]])
        })
        .collect::<PropertiesResult<Vec<_>>>()?;
    let grid_size = array(required(record, "grid_size")?, "grid_size")?
        .iter()
        .enumerate()
        .map(|(index, value)| usize_value(value, &format!("grid_size[{index}]")))
        .collect::<PropertiesResult<Vec<_>>>()?;
    let method = match required(record, "refinement_method")?.as_str() {
        Some("PrincipalAxis") => RefinementMethod::PrincipalAxis,
        Some("QuasiNewton") => RefinementMethod::QuasiNewton,
        _ => {
            return Err(malformed(
                "refinement_method must be PrincipalAxis or QuasiNewton",
            ));
        }
    };
    Ok((
        occupied,
        GaplessSearchOptions {
            brillouin_zone: zone,
            grid_size,
            candidate_count: usize_value(required(record, "candidate_count")?, "candidate_count")?,
            gap_tolerance: finite_real(required(record, "gap_tolerance")?, "gap_tolerance")?,
            merge_tolerance: finite_real(required(record, "merge_tolerance")?, "merge_tolerance")?,
            hermitian_tolerance: finite_real(
                required(record, "hermitian_tolerance")?,
                "hermitian_tolerance",
            )?,
            max_iterations: usize_value(required(record, "max_iterations")?, "max_iterations")?,
            refinement_method: method,
        },
    ))
}

pub(crate) fn gapless_result_json(data: &GaplessSearchData) -> PropertiesResult<String> {
    let method = match data.refinement_method {
        RefinementMethod::PrincipalAxis => "PrincipalAxis",
        RefinementMethod::QuasiNewton => "QuasiNewton",
    };
    let records = data
        .point_records
        .iter()
        .map(|record| {
            json!({
                "Point": record.point,
                "Gap": record.gap,
                "Source": record.source,
                "GridSeed": record.grid_seed,
                "GridGap": record.grid_gap,
                "Refined": record.refined,
                "RefinementStatus": record.refinement_status,
                "ObjectiveMinimum": record.objective_minimum,
                "Eigenvalues": record.eigenvalues,
            })
        })
        .collect::<Vec<_>>();
    serde_json::to_string(&json!({
        "Schema": "MagneticTBGaplessPointSearch",
        "SchemaVersion": 1,
        "Points": data.points,
        "PointRecords": records,
        "OccupiedBands": data.occupied_bands,
        "HamiltonianDimension": data.hamiltonian_dimension,
        "BrillouinZone": data.brillouin_zone,
        "GridSize": data.grid_size,
        "GridPointCount": data.grid_point_count,
        "LocalMinimumCount": data.local_minimum_count,
        "SeedCount": data.seed_count,
        "FailedRefinements": data.failed_refinements,
        "GapTolerance": data.gap_tolerance,
        "MergeTolerance": data.merge_tolerance,
        "RefinementMethod": method,
        "MinimumGap": data.minimum_gap,
        "Method": "periodic grid local minima plus configurable gap-squared refinement",
        "CompletenessGuaranteed": false,
    }))
    .map_err(|error| malformed(format!("could not serialize gapless result: {error}")))
}
