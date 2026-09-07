fn symbolic_terms(
    expression: &ExactExpression,
) -> ExactResult<Vec<(QuadraticReal, Option<Value>)>> {
    let empty_bindings = BTreeMap::new();
    if let Ok(scalar) = QuadraticReal::from_expression(expression, &empty_bindings) {
        return Ok(vec![(scalar, None)]);
    }
    if let ExactExpression::Call { head, arguments } = expression {
        if head == "System`Plus" {
            return arguments
                .iter()
                .try_fold(Vec::new(), |mut terms, argument| {
                    terms.extend(symbolic_terms(argument)?);
                    Ok(terms)
                });
        }
        if head == "System`Times" {
            let mut scalar = QuadraticReal::one();
            let mut residual = Vec::new();
            for argument in arguments {
                if let Ok(value) = QuadraticReal::from_expression(argument, &empty_bindings) {
                    scalar = scalar.multiply(&value);
                } else {
                    residual.push(exact_expression_value(argument)?);
                }
            }
            let residual = match residual.len() {
                0 => None,
                1 => residual.pop(),
                _ => Some(json!({
                    "kind": "call",
                    "head": "System`Times",
                    "arguments": residual
                })),
            };
            return Ok(vec![(scalar, residual)]);
        }
    }
    Ok(vec![(
        QuadraticReal::one(),
        Some(exact_expression_value(expression)?),
    )])
}
fn symbolic_cartesian_displacement(
    displacement: &Value,
    symbolic_lattice: &Value,
) -> ExactResult<Value> {
    let displacement = array(displacement)?
        .iter()
        .map(quadratic_from_json)
        .collect::<ExactResult<Vec<_>>>()?;
    let lattice_record = object(symbolic_lattice)?;
    if lattice_record.get("kind").and_then(Value::as_str) != Some("matrix") {
        return Err(malformed("symbolic model lattice must be an exact matrix"));
    }
    let rows = size_field(lattice_record, "rows")?;
    let columns = size_field(lattice_record, "columns")?;
    let entries: Vec<ExactExpression> = from_value(
        required(lattice_record, "entries")?,
        "symbolic model lattice entries",
    )?;
    if rows.checked_mul(columns) != Some(entries.len()) || rows != displacement.len() {
        return Err(ExactError::new(
            "UnsupportedCoordinateMode",
            "symbolic lattice and fractional displacement dimensions do not align",
        ));
    }
    let mut output = Vec::with_capacity(columns);
    for column in 0..columns {
        let mut combined: Vec<(String, QuadraticReal, Option<Value>)> = Vec::new();
        for (row, coordinate) in displacement.iter().enumerate() {
            for (coefficient, residual) in symbolic_terms(&entries[row * columns + column])? {
                let coefficient = coordinate.multiply(&coefficient);
                if coefficient.is_zero() {
                    continue;
                }
                let key = residual
                    .as_ref()
                    .map_or_else(|| "null".to_owned(), Value::to_string);
                if let Some((_, existing, _)) = combined
                    .iter_mut()
                    .find(|(existing_key, _, _)| existing_key == &key)
                {
                    *existing = existing.add(&coefficient);
                } else {
                    combined.push((key, coefficient, residual));
                }
            }
        }
        let mut terms = combined
            .into_iter()
            .filter(|(_, coefficient, _)| !coefficient.is_zero())
            .map(|(_, coefficient, residual)| match residual {
                None => quadratic_expression(&coefficient),
                Some(residual) => {
                    if coefficient == QuadraticReal::one() {
                        residual
                    } else {
                        json!({
                            "kind": "call",
                            "head": "System`Times",
                            "arguments": [quadratic_expression(&coefficient), residual]
                        })
                    }
                }
            })
            .collect::<Vec<_>>();
        output.push(match terms.len() {
            0 => json!({"kind": "integer", "value": "0"}),
            1 => terms.remove(0),
            _ => json!({"kind": "call", "head": "System`Plus", "arguments": terms}),
        });
    }
    Ok(Value::Array(output))
}

fn cartesian_symbolic_hamiltonian(request: &serde_json::Map<String, Value>) -> ExactResult<Value> {
    let mut hamiltonian = required(request, "hamiltonian")?.clone();
    let lattice = required(request, "lattice")?;
    let mut symbolic_lattice = request
        .get("symbolic_lattice")
        .cloned()
        .unwrap_or_else(|| lattice.clone());
    normalize_stable_inexact(&mut symbolic_lattice)?;
    let record = hamiltonian
        .as_object_mut()
        .ok_or_else(|| malformed("hamiltonian must be an object"))?;
    if required(record, "schema")?.as_str() != Some("magnetictb.symbolic_hamiltonian.v1") {
        return Err(ExactError::new(
            "UnsupportedCoordinateMode",
            "only symbolic Hamiltonian v1 can be converted to Cartesian momentum",
        ));
    }
    for term in record
        .get_mut("terms")
        .and_then(Value::as_array_mut)
        .ok_or_else(|| malformed("Hamiltonian terms must be an array"))?
    {
        for coefficient in term
            .get_mut("fourier_coefficients")
            .and_then(Value::as_array_mut)
            .ok_or_else(|| malformed("Fourier coefficients must be an array"))?
        {
            let coefficient = coefficient
                .as_object_mut()
                .ok_or_else(|| malformed("Fourier coefficient must be an object"))?;
            let fractional = required(coefficient, "displacement")?.clone();
            let numerical = cartesian_displacement(&fractional, lattice)?;
            let symbolic = symbolic_cartesian_displacement(&fractional, &symbolic_lattice)?;
            coefficient.insert("displacement".to_owned(), symbolic);
            coefficient.insert("numerical_displacement".to_owned(), numerical);
        }
    }
    for row in record
        .get_mut("matrix_entries")
        .and_then(Value::as_array_mut)
        .ok_or_else(|| malformed("Hamiltonian matrix entries must be an array"))?
    {
        for cell in row
            .as_array_mut()
            .ok_or_else(|| malformed("Hamiltonian matrix row must be an array"))?
        {
            for contribution in cell
                .as_array_mut()
                .ok_or_else(|| malformed("Hamiltonian matrix cell must be an array"))?
            {
                let contribution = contribution
                    .as_object_mut()
                    .ok_or_else(|| malformed("Hamiltonian contribution must be an object"))?;
                let fractional = required(contribution, "displacement")?.clone();
                let numerical = cartesian_displacement(&fractional, lattice)?;
                let symbolic = symbolic_cartesian_displacement(&fractional, &symbolic_lattice)?;
                contribution.insert("displacement".to_owned(), symbolic);
                contribution.insert("numerical_displacement".to_owned(), numerical);
            }
        }
    }
    record.insert(
        "gauge".to_owned(),
        json!({
            "momentum_coordinates": "cartesian_kx_ky_kz",
            "phase_convention": "exp(+i*k_dot_cartesian_displacement)",
            "row_axis": "bra_destination",
            "column_axis": "ket_source",
            "bz_folding": false
        }),
    );
    Ok(hamiltonian)
}

fn vector_to_json(values: &[QuadraticReal]) -> Value {
    Value::Array(values.iter().map(quadratic_to_json).collect())
}

fn parameter_name(shell: usize, number: usize) -> String {
    match shell {
        1 => format!("e{number}"),
        2 => format!("t{number}"),
        3 => format!("r{number}"),
        4 => format!("s{number}"),
        _ => format!("p{shell}n{number}"),
    }
}

fn hamiltonian_options(
    request: &serde_json::Map<String, Value>,
) -> ExactResult<(
    bool,
    ConstraintKernelMethod,
    ConstraintValidationLevel,
    String,
    String,
)> {
    let hermitian = request
        .get("hermitian")
        .map_or(Ok(true), |value| from_value(value, "hermitian"))?;
    let method_name = request
        .get("kernel_method")
        .and_then(Value::as_str)
        .unwrap_or("Iterative");
    let method = match method_name {
        "Iterative" => ConstraintKernelMethod::Iterative,
        "Stacked" => ConstraintKernelMethod::Stacked,
        "Cyclotomic" => ConstraintKernelMethod::Cyclotomic,
        _ => {
            return Err(ExactError::new(
                "InvalidKernelMethod",
                format!("unsupported constraint-kernel method {method_name}"),
            ));
        }
    };
    let validation = request
        .get("validation_level")
        .and_then(Value::as_str)
        .unwrap_or("Basic");
    if !matches!(validation, "None" | "Basic" | "Full") {
        return Err(ExactError::new(
            "InvalidValidationLevel",
            format!("unsupported validation level {validation}"),
        ));
    }
    let validation_kind = match validation {
        "None" => ConstraintValidationLevel::None,
        "Basic" => ConstraintValidationLevel::Basic,
        "Full" => ConstraintValidationLevel::Full,
        _ => unreachable!("validation was checked above"),
    };
    Ok((
        hermitian,
        method,
        validation_kind,
        method_name.to_owned(),
        validation.to_owned(),
    ))
}

fn propagate_hamiltonian_options(
    source: &serde_json::Map<String, Value>,
    target: &mut serde_json::Map<String, Value>,
) {
    for key in ["hermitian", "kernel_method", "validation_level"] {
        if let Some(value) = source.get(key) {
            target.insert(key.to_owned(), value.clone());
        }
    }
}

fn hamiltonian_gauge() -> Value {
    json!({
        "momentum_coordinates": "stable_kx_ky_kz",
        "phase_convention": "exp(+i*k_dot_fractional_displacement)",
        "row_axis": "bra_destination",
        "column_axis": "ket_source",
        "bz_folding": false
    })
}

fn model_identity_payload(
    context: &Arc<CyclotomicContext>,
    result: &serde_json::Map<String, Value>,
) -> ExactResult<Value> {
    Ok(json!({
        "schema": "magnetictb.model_identity.v1",
        "field_context": context_to_json(context),
        "gauge": hamiltonian_gauge(),
        "lattice": required(result, "model_lattice")?.clone(),
        "multiplication_table": required(result, "multiplication_table")?.clone(),
        "antiunitary_flags": required(result, "antiunitary_flags")?.clone(),
        "operation_labels": required(result, "operation_labels")?.clone(),
        "spatial_actions": required(result, "spatial_actions")?.clone(),
        "spin_rotations": required(result, "spin_rotations")?.clone(),
        "image_site_indices": required(result, "image_site_indices")?.clone(),
        "cell_translations": required(result, "cell_translations")?.clone(),
        "local_dimensions": required(result, "local_dimensions")?.clone(),
        "site_dimensions": required(result, "site_dimensions")?.clone(),
        "orbital_layout": required(result, "orbital_layout")?.clone(),
        "representation_matrices": required(result, "representation_matrices")?.clone()
    }))
}

fn canonical_sha256(value: &Value) -> ExactResult<String> {
    let bytes = serde_json::to_vec(value)
        .map_err(|error| malformed(format!("cannot canonicalize model identity: {error}")))?;
    let mut output = String::with_capacity(64);
    for byte in Sha256::digest(bytes) {
        write!(&mut output, "{byte:02x}")
            .map_err(|error| malformed(format!("cannot encode model identity hash: {error}")))?;
    }
    Ok(output)
}

fn hamiltonian_matrix_entries(
    context: &Arc<CyclotomicContext>,
    dimension: usize,
    terms: &[Value],
) -> ExactResult<Value> {
    let mut entries = vec![vec![Vec::<Value>::new(); dimension]; dimension];
    for term in terms {
        let record = object(term)?;
        let parameter = required(record, "parameter")?.clone();
        for coefficient in array(required(record, "fourier_coefficients")?)? {
            let coefficient = object(coefficient)?;
            let displacement = required(coefficient, "displacement")?.clone();
            let numerical_displacement = coefficient.get("numerical_displacement").cloned();
            let matrix = matrix_from_json(context, required(coefficient, "matrix")?)?;
            if matrix.rows() != dimension || matrix.columns() != dimension {
                return Err(ExactError::new(
                    "IncompatibleHamiltonian",
                    "a Fourier coefficient matrix has the wrong Hamiltonian shape",
                ));
            }
            for (row, row_entries) in entries.iter_mut().enumerate() {
                for (column, cell) in row_entries.iter_mut().enumerate() {
                    let value = matrix.entry(row, column)?;
                    if !value.is_zero() {
                        let mut contribution = json!({
                            "parameter": parameter,
                            "displacement": displacement,
                            "coefficient": element_to_json(value)
                        });
                        if let Some(numerical) = &numerical_displacement {
                            contribution
                                .as_object_mut()
                                .expect("Hamiltonian contribution is an object")
                                .insert("numerical_displacement".to_owned(), numerical.clone());
                        }
                        cell.push(contribution);
                    }
                }
            }
        }
    }
    Ok(Value::Array(
        entries
            .into_iter()
            .map(|row| Value::Array(row.into_iter().map(Value::Array).collect()))
            .collect(),
    ))
}

fn symbolic_hamiltonian_from_result(
    context: &Arc<CyclotomicContext>,
    result: &serde_json::Map<String, Value>,
    shell: usize,
) -> ExactResult<Value> {
    let site_dimensions: Vec<usize> =
        from_value(required(result, "site_dimensions")?, "site_dimensions")?;
    let dimension = site_dimensions.iter().try_fold(0usize, |total, value| {
        total
            .checked_add(*value)
            .ok_or_else(|| ExactError::new("DimensionMismatch", "Hamiltonian dimension overflows"))
    })?;
    if dimension == 0 {
        return Err(ExactError::new(
            "DimensionMismatch",
            "Hamiltonian dimension must be positive",
        ));
    }
    let parameter_order = array(required(result, "parameter_order")?)?;
    let parameter_space = object(required(result, "parameter_space")?)?;
    let solutions = array(required(parameter_space, "parameter_basis_solutions")?)?;
    if parameter_order.len() != solutions.len() {
        return Err(malformed(
            "parameter order and parameter-basis solutions do not align",
        ));
    }
    let terms = parameter_order
        .iter()
        .zip(solutions)
        .enumerate()
        .map(|(index, (order, solution))| {
            let order = object(order)?;
            let solution = object(solution)?;
            let parameter_index: usize =
                from_value(required(order, "parameter_index")?, "parameter_index")?;
            if parameter_index != index
                || from_value::<usize>(
                    required(solution, "parameter_index")?,
                    "solution parameter_index",
                )? != index
            {
                return Err(malformed("parameter indices are not in stable order"));
            }
            Ok(json!({
                "parameter": {
                    "name": parameter_name(shell, index + 1),
                    "shell": shell,
                    "parameter_number": index + 1,
                    "core_parameter_index": index
                },
                "fourier_coefficients": required(solution, "fourier_coefficients")?.clone(),
                "gamma_hamiltonian": required(solution, "gamma_hamiltonian")?.clone(),
                "verified": required(solution, "verified")?.clone()
            }))
        })
        .collect::<ExactResult<Vec<_>>>()?;
    let identity_payload = model_identity_payload(context, result)?;
    let identity_sha256 = canonical_sha256(&identity_payload)?;
    let matrix_entries = hamiltonian_matrix_entries(context, dimension, &terms)?;
    Ok(json!({
        "schema": "magnetictb.symbolic_hamiltonian.v1",
        "version": 1,
        "model_identity_sha256": identity_sha256,
        "field_context": context_to_json(context),
        "gauge": hamiltonian_gauge(),
        "shape": [dimension, dimension],
        "shells": [shell],
        "hermitian": required(result, "hermitian")?.clone(),
        "kernel_method": required(result, "kernel_method")?.clone(),
        "validation_level": required(result, "validation_level")?.clone(),
        "parameter_names": terms.iter().map(|term| {
            term["parameter"]["name"].clone()
        }).collect::<Vec<_>>(),
        "terms": terms,
        "matrix_entries": matrix_entries,
        "covariance_verified": required(parameter_space, "full_parameter_space_verified")?.clone()
    }))
}

struct UnconstrainedParameter<'a> {
    context: &'a Arc<CyclotomicContext>,
    dimension: usize,
    row: usize,
    column: usize,
    coefficient: &'a Element,
    name: &'a str,
    shell: usize,
    parameter_index: usize,
    displacement: &'a Value,
}

fn unconstrained_parameter_term(spec: &UnconstrainedParameter<'_>) -> ExactResult<Value> {
    let mut entries = vec![Element::zero(spec.context)?; spec.dimension * spec.dimension];
    entries[spec.row * spec.dimension + spec.column] = spec.coefficient.clone();
    let matrix = ExactMatrix::new(
        Arc::clone(spec.context),
        spec.dimension,
        spec.dimension,
        entries,
    )?;
    let matrix_json = matrix_to_json(&matrix);
    Ok(json!({
        "parameter": {
            "name": spec.name,
            "shell": spec.shell,
            "parameter_number": spec.parameter_index + 1,
            "core_parameter_index": spec.parameter_index
        },
        "fourier_coefficients": [{
            "displacement": spec.displacement,
            "matrix": matrix_json.clone()
        }],
        "gamma_hamiltonian": matrix_json,
        "verified": false
    }))
}

fn unconstrained_symbolic_hamiltonian(
    request: &serde_json::Map<String, Value>,
) -> ExactResult<Value> {
    let result = object(required(request, "model_result")?)?;
    let shell: usize = from_value(required(request, "shell")?, "shell")?;
    if shell == 0 {
        return Err(ExactError::new(
            "InvalidBondShell",
            "the shell index must be a positive integer",
        ));
    }
    let context = context_from_json(required(result, "field_context")?, 128)?;
    let site_dimensions: Vec<usize> =
        from_value(required(result, "site_dimensions")?, "site_dimensions")?;
    if site_dimensions.is_empty() || site_dimensions.contains(&0) {
        return Err(ExactError::new(
            "DimensionMismatch",
            "unconstrained Hamiltonian site dimensions must be positive",
        ));
    }
    let mut offsets = Vec::with_capacity(site_dimensions.len());
    let mut dimension = 0usize;
    for local_dimension in &site_dimensions {
        offsets.push(dimension);
        dimension = dimension.checked_add(*local_dimension).ok_or_else(|| {
            ExactError::new("DimensionMismatch", "Hamiltonian dimension overflows")
        })?;
    }
    let one = Element::one(&context)?;
    let imaginary = Element::root_of_unity(&context, 4, 1)?;
    let bonds = array(required(result, "bonds")?)?;
    let mut terms = Vec::new();
    for (bond_index, bond) in bonds.iter().enumerate() {
        let bond = object(bond)?;
        let row_site: usize = from_value(required(bond, "source_site")?, "source_site")?;
        let column_site: usize = from_value(required(bond, "target_site")?, "target_site")?;
        let row_dimension = *site_dimensions.get(row_site).ok_or_else(|| {
            ExactError::new(
                "InvalidBondSite",
                "unconstrained bond row site is out of range",
            )
        })?;
        let column_dimension = *site_dimensions.get(column_site).ok_or_else(|| {
            ExactError::new(
                "InvalidBondSite",
                "unconstrained bond column site is out of range",
            )
        })?;
        for local_row in 0..row_dimension {
            for local_column in 0..column_dimension {
                let global_row = offsets[row_site] + local_row;
                let global_column = offsets[column_site] + local_column;
                for (prefix, coefficient) in [("tr", &one), ("ti", &imaginary)] {
                    let name = format!(
                        "{prefix}{}{}{}",
                        bond_index + 1,
                        local_row + 1,
                        local_column + 1
                    );
                    terms.push(unconstrained_parameter_term(&UnconstrainedParameter {
                        context: &context,
                        dimension,
                        row: global_row,
                        column: global_column,
                        coefficient,
                        name: &name,
                        shell,
                        parameter_index: terms.len(),
                        displacement: required(bond, "displacement")?,
                    })?);
                }
            }
        }
    }
    let identity_payload = model_identity_payload(&context, result)?;
    let identity_sha256 = canonical_sha256(&identity_payload)?;
    let matrix_entries = hamiltonian_matrix_entries(&context, dimension, &terms)?;
    Ok(json!({
        "schema": "magnetictb.symbolic_hamiltonian.v1",
        "version": 1,
        "model_identity_sha256": identity_sha256,
        "field_context": context_to_json(&context),
        "gauge": hamiltonian_gauge(),
        "shape": [dimension, dimension],
        "shells": [shell],
        "hermitian": false,
        "kernel_method": "Unconstrained",
        "validation_level": "None",
        "parameter_names": terms.iter().map(|term| {
            term["parameter"]["name"].clone()
        }).collect::<Vec<_>>(),
        "terms": terms,
        "matrix_entries": matrix_entries,
        "covariance_verified": false
    }))
}

fn combine_symbolic_hamiltonians(request: &serde_json::Map<String, Value>) -> ExactResult<Value> {
    let left = object(required(request, "left")?)?;
    let right = object(required(request, "right")?)?;
    for record in [left, right] {
        if required(record, "schema")?.as_str() != Some("magnetictb.symbolic_hamiltonian.v1")
            || required(record, "version")?.as_u64() != Some(1)
        {
            return Err(ExactError::new(
                "IncompatibleHamiltonian",
                "only magnetictb.symbolic_hamiltonian.v1 values can be combined",
            ));
        }
    }
    for field in [
        "model_identity_sha256",
        "field_context",
        "gauge",
        "shape",
        "hermitian",
        "kernel_method",
        "validation_level",
    ] {
        if required(left, field)? != required(right, field)? {
            return Err(ExactError::new(
                "IncompatibleHamiltonian",
                format!("Hamiltonian {field} values differ"),
            ));
        }
    }
    let mut shells: Vec<usize> = from_value(required(left, "shells")?, "left shells")?;
    let right_shells: Vec<usize> = from_value(required(right, "shells")?, "right shells")?;
    if right_shells.iter().any(|shell| shells.contains(shell)) {
        return Err(ExactError::new(
            "DuplicateHamiltonianShell",
            "the same shell cannot be added twice",
        ));
    }
    shells.extend(right_shells);
    shells.sort_unstable();
    let mut terms = array(required(left, "terms")?)?.clone();
    terms.extend(array(required(right, "terms")?)?.clone());
    terms.sort_by_key(|term| {
        let parameter = term.get("parameter");
        (
            parameter
                .and_then(|value| value.get("shell"))
                .and_then(Value::as_u64)
                .unwrap_or(u64::MAX),
            parameter
                .and_then(|value| value.get("parameter_number"))
                .and_then(Value::as_u64)
                .unwrap_or(u64::MAX),
        )
    });
    let mut names = Vec::with_capacity(terms.len());
    for term in &terms {
        let name = required(object(required(object(term)?, "parameter")?)?, "name")?
            .as_str()
            .ok_or_else(|| malformed("Hamiltonian parameter name must be a string"))?;
        if names.contains(&name) {
            return Err(ExactError::new(
                "DuplicateHamiltonianParameter",
                format!("Hamiltonian parameter {name} occurs more than once"),
            ));
        }
        names.push(name);
    }
    let context = context_from_json(required(left, "field_context")?, 128)?;
    let shape: Vec<usize> = from_value(required(left, "shape")?, "shape")?;
    if shape.len() != 2 || shape[0] != shape[1] || shape[0] == 0 {
        return Err(ExactError::new(
            "IncompatibleHamiltonian",
            "Hamiltonian shape must be a positive square matrix",
        ));
    }
    let matrix_entries = hamiltonian_matrix_entries(&context, shape[0], &terms)?;
    let verified = required(left, "covariance_verified")?.as_bool() == Some(true)
        && required(right, "covariance_verified")?.as_bool() == Some(true);
    Ok(json!({
        "schema": "magnetictb.symbolic_hamiltonian.v1",
        "version": 1,
        "model_identity_sha256": required(left, "model_identity_sha256")?.clone(),
        "field_context": required(left, "field_context")?.clone(),
        "gauge": required(left, "gauge")?.clone(),
        "shape": shape,
        "shells": shells,
        "hermitian": required(left, "hermitian")?.clone(),
        "kernel_method": required(left, "kernel_method")?.clone(),
        "validation_level": required(left, "validation_level")?.clone(),
        "parameter_names": names,
        "terms": terms,
        "matrix_entries": matrix_entries,
        "covariance_verified": verified
    }))
}

fn complex_add(left: (f64, f64), right: (f64, f64)) -> (f64, f64) {
    (left.0 + right.0, left.1 + right.1)
}

fn complex_multiply(left: (f64, f64), right: (f64, f64)) -> (f64, f64) {
    (
        left.0.mul_add(right.0, -(left.1 * right.1)),
        left.0.mul_add(right.1, left.1 * right.0),
    )
}

fn finite_complex(real: f64, imaginary: f64, field: &str) -> ExactResult<(f64, f64)> {
    if real.is_finite() && imaginary.is_finite() {
        Ok((real, imaginary))
    } else {
        Err(ExactError::new(
            "NonFiniteHamiltonianValue",
            format!("{field} must evaluate to a finite complex number"),
        ))
    }
}

fn approximate_element(value: &Element) -> ExactResult<(f64, f64)> {
    let conductor_u32 = u32::try_from(value.context().conductor).map_err(|_| {
        ExactError::new(
            "NonFiniteHamiltonianValue",
            "cyclotomic conductor is outside the numerical evaluation range",
        )
    })?;
    let conductor = f64::from(conductor_u32);
    let mut result = (0.0, 0.0);
    for (power, coefficient) in value.coefficients().iter().enumerate() {
        let scalar = coefficient.to_f64().ok_or_else(|| {
            ExactError::new(
                "NonFiniteHamiltonianValue",
                "exact cyclotomic coefficient is outside the finite f64 range",
            )
        })?;
        let power_u32 = u32::try_from(power).map_err(|_| {
            ExactError::new(
                "NonFiniteHamiltonianValue",
                "cyclotomic power is outside the numerical evaluation range",
            )
        })?;
        let phase = std::f64::consts::TAU * f64::from(power_u32) / conductor;
        result = complex_add(result, (scalar * phase.cos(), scalar * phase.sin()));
    }
    finite_complex(result.0, result.1, "exact scalar")
}

fn numeric_complex(
    context: &Arc<CyclotomicContext>,
    value: &Value,
    field: &str,
) -> ExactResult<(f64, f64)> {
    if let Some(number) = value.as_f64() {
        return finite_complex(number, 0.0, field);
    }
    if let Some(record) = value.as_object()
        && record.get("kind").and_then(Value::as_str) == Some("approx_complex")
    {
        let real = required(record, "real")?
            .as_f64()
            .ok_or_else(|| malformed(format!("{field} real part must be numeric")))?;
        let imaginary = required(record, "imaginary")?
            .as_f64()
            .ok_or_else(|| malformed(format!("{field} imaginary part must be numeric")))?;
        return finite_complex(real, imaginary, field);
    }
    approximate_element(&encoded_element_from_json(context, value)?)
}

fn numeric_real(context: &Arc<CyclotomicContext>, value: &Value, field: &str) -> ExactResult<f64> {
    if let Some(number) = value.as_f64() {
        return finite_complex(number, 0.0, field).map(|number| number.0);
    }
    let element = encoded_element_from_json(context, value)?;
    if element.conjugate()? != element {
        return Err(ExactError::new(
            "ComplexMomentum",
            format!("{field} must be real"),
        ));
    }
    approximate_element(&element).map(|number| number.0)
}

fn evaluate_hamiltonian_cell(
    context: &Arc<CyclotomicContext>,
    cell: &Value,
    parameters: &BTreeMap<&str, (f64, f64)>,
    momentum: &[f64],
) -> ExactResult<(f64, f64)> {
    let mut sum = (0.0, 0.0);
    for contribution in array(cell)? {
        let contribution = object(contribution)?;
        let parameter_record = object(required(contribution, "parameter")?)?;
        let parameter_name = required(parameter_record, "name")?
            .as_str()
            .ok_or_else(|| malformed("contribution parameter must be a string"))?;
        let parameter = *parameters.get(parameter_name).ok_or_else(|| {
            ExactError::new(
                "InvalidHamiltonianParameters",
                format!("unknown Hamiltonian parameter {parameter_name}"),
            )
        })?;
        let coefficient = approximate_element(&encoded_element_from_json(
            context,
            required(contribution, "coefficient")?,
        )?)?;
        let displacement_value = contribution
            .get("numerical_displacement")
            .unwrap_or(required(contribution, "displacement")?);
        let displacement = array(displacement_value)?
            .iter()
            .map(quadratic_from_json)
            .collect::<ExactResult<Vec<_>>>()?;
        if displacement.len() != momentum.len() {
            return Err(malformed("Hamiltonian displacement dimension is not three"));
        }
        let phase =
            displacement
                .iter()
                .zip(momentum)
                .try_fold(0.0, |total, (coordinate, momentum)| {
                    coordinate
                        .to_f64()
                        .map(|coordinate| total + coordinate * momentum)
                        .map_err(|error| ExactError::new(error.tag(), error.to_string()))
                })?;
        sum = complex_add(
            sum,
            complex_multiply(
                parameter,
                complex_multiply(coefficient, (phase.cos(), phase.sin())),
            ),
        );
    }
    finite_complex(sum.0, sum.1, "Hamiltonian matrix entry")
}

fn evaluate_symbolic_hamiltonian(request: &serde_json::Map<String, Value>) -> ExactResult<Value> {
    let hamiltonian = object(required(request, "hamiltonian")?)?;
    if required(hamiltonian, "schema")?.as_str() != Some("magnetictb.symbolic_hamiltonian.v1") {
        return Err(ExactError::new(
            "InvalidHamiltonianEvaluation",
            "only symbolic Hamiltonian v1 can be evaluated",
        ));
    }
    let context = context_from_json(required(hamiltonian, "field_context")?, 128)?;
    let shape: Vec<usize> = from_value(required(hamiltonian, "shape")?, "shape")?;
    if shape.len() != 2 || shape[0] == 0 || shape[0] != shape[1] {
        return Err(ExactError::new(
            "InvalidHamiltonianEvaluation",
            "Hamiltonian shape must be positive and square",
        ));
    }
    let momentum = array(required(request, "momentum")?)?;
    if momentum.len() != 3 {
        return Err(ExactError::new(
            "InvalidMomentum",
            "momentum must contain exactly kx, ky, kz",
        ));
    }
    let momentum = momentum
        .iter()
        .enumerate()
        .map(|(index, value)| numeric_real(&context, value, ["kx", "ky", "kz"][index]))
        .collect::<ExactResult<Vec<_>>>()?;
    let raw_parameters = object(required(request, "parameters")?)?;
    let names = array(required(hamiltonian, "parameter_names")?)?;
    if raw_parameters.len() != names.len() {
        return Err(ExactError::new(
            "InvalidHamiltonianParameters",
            "parameter mapping must contain every Hamiltonian parameter exactly once",
        ));
    }
    let mut parameters = BTreeMap::new();
    for name in names {
        let name = name
            .as_str()
            .ok_or_else(|| malformed("Hamiltonian parameter name must be a string"))?;
        let value = raw_parameters.get(name).ok_or_else(|| {
            ExactError::new(
                "InvalidHamiltonianParameters",
                format!("missing Hamiltonian parameter {name}"),
            )
        })?;
        parameters.insert(name, numeric_complex(&context, value, name)?);
    }
    let rows = array(required(hamiltonian, "matrix_entries")?)?;
    if rows.len() != shape[0] {
        return Err(malformed(
            "Hamiltonian matrix row count does not match shape",
        ));
    }
    let mut matrix = Vec::with_capacity(shape[0]);
    for row in rows {
        let cells = array(row)?;
        if cells.len() != shape[1] {
            return Err(malformed(
                "Hamiltonian matrix column count does not match shape",
            ));
        }
        let mut numeric_row = Vec::with_capacity(shape[1]);
        for cell in cells {
            let sum = evaluate_hamiltonian_cell(&context, cell, &parameters, &momentum)?;
            numeric_row.push(json!({"real": sum.0, "imaginary": sum.1}));
        }
        matrix.push(Value::Array(numeric_row));
    }
    Ok(json!({
        "schema": "magnetictb.numeric_hamiltonian.v1",
        "shape": shape,
        "momentum": momentum,
        "matrix": matrix
    }))
}

#[allow(clippy::too_many_lines)]
fn hopping_data_from_shell_results(request: &serde_json::Map<String, Value>) -> ExactResult<Value> {
    let shell_results = array(required(request, "shell_results")?)?;
    if shell_results.is_empty() {
        return Err(ExactError::new(
            "InvalidShellSelection",
            "at least one solved shell is required",
        ));
    }
    let parameters = object(required(request, "parameters")?)?;
    let first = object(&shell_results[0])?;
    let context = context_from_json(required(first, "field_context")?, 128)?;
    let site_dimensions: Vec<usize> =
        from_value(required(first, "site_dimensions")?, "site_dimensions")?;
    if site_dimensions.is_empty() || site_dimensions.contains(&0) {
        return Err(malformed("site dimensions must be nonzero"));
    }
    let mut offsets = Vec::with_capacity(site_dimensions.len());
    let mut dimension = 0usize;
    for &local in &site_dimensions {
        offsets.push(dimension);
        dimension = dimension
            .checked_add(local)
            .ok_or_else(|| malformed("orbital dimension overflows usize"))?;
    }
    let identity = required(first, "model_identity_sha256")?;
    let mut translations = BTreeMap::<[i64; 3], DMatrix<Complex<f64>>>::new();
    let mut shells = Vec::with_capacity(shell_results.len());
    let mut used_parameter_names = Vec::new();
    for raw_result in shell_results {
        let result = object(raw_result)?;
        if required(result, "model_identity_sha256")? != identity
            || required(result, "field_context")? != required(first, "field_context")?
            || from_value::<Vec<usize>>(required(result, "site_dimensions")?, "site_dimensions")?
                != site_dimensions
        {
            return Err(ExactError::new(
                "IncompatibleHamiltonian",
                "selected shells do not belong to one prepared model",
            ));
        }
        let shell: usize = from_value(required(result, "target_shell")?, "target_shell")?;
        shells.push(shell);
        let bonds = array(required(result, "bonds")?)?;
        let hamiltonian = object(required(result, "hamiltonian")?)?;
        let names = array(required(hamiltonian, "parameter_names")?)?;
        let parameter_space = object(required(result, "parameter_space")?)?;
        let solutions = array(required(parameter_space, "parameter_basis_solutions")?)?;
        if names.len() != solutions.len() {
            return Err(malformed(
                "parameter names and basis solutions do not align",
            ));
        }
        for (name, solution) in names.iter().zip(solutions) {
            let name = name
                .as_str()
                .ok_or_else(|| malformed("Hamiltonian parameter name must be a string"))?;
            let parameter = parameters.get(name).ok_or_else(|| {
                ExactError::new(
                    "NonnumericHopping",
                    format!("missing numeric value for Hamiltonian parameter {name}"),
                )
            })?;
            let scalar = numeric_complex(&context, parameter, name)?;
            used_parameter_names.push(name.to_owned());
            let solution = object(solution)?;
            for term in array(required(solution, "terms")?)? {
                let term = object(term)?;
                let bond_index: usize = from_value(required(term, "bond_index")?, "bond_index")?;
                let bond = bonds
                    .get(bond_index)
                    .ok_or_else(|| malformed("term bond index is out of range"))?;
                let bond = object(bond)?;
                let translation_values: Vec<i64> =
                    from_value(required(bond, "translation")?, "translation")?;
                if translation_values.len() != 3 {
                    return Err(malformed("bond translation must have length three"));
                }
                let translation = [
                    translation_values[0],
                    translation_values[1],
                    translation_values[2],
                ];
                let row_block: usize = from_value(required(term, "row_block")?, "row_block")?;
                let column_block: usize =
                    from_value(required(term, "column_block")?, "column_block")?;
                if row_block >= site_dimensions.len() || column_block >= site_dimensions.len() {
                    return Err(malformed("term block index is out of range"));
                }
                let matrix = matrix_from_json(&context, required(term, "matrix")?)?;
                if matrix.rows() != site_dimensions[row_block]
                    || matrix.columns() != site_dimensions[column_block]
                {
                    return Err(malformed("term matrix does not match its site blocks"));
                }
                let target = translations
                    .entry(translation)
                    .or_insert_with(|| DMatrix::zeros(dimension, dimension));
                for row in 0..matrix.rows() {
                    for column in 0..matrix.columns() {
                        let coefficient = approximate_element(matrix.entry(row, column)?)?;
                        let value = complex_multiply(scalar, coefficient);
                        target[(offsets[row_block] + row, offsets[column_block] + column)] +=
                            Complex::new(value.0, value.1);
                    }
                }
            }
        }
    }
    used_parameter_names.sort();
    used_parameter_names.dedup();
    let orbital_layout = array(required(first, "orbital_layout")?)?;
    if orbital_layout.len() != dimension {
        return Err(malformed(
            "orbital layout does not match the Hamiltonian dimension",
        ));
    }
    let centers = orbital_layout
        .iter()
        .enumerate()
        .map(|(orbital, value)| {
            let record = object(value)?;
            let position = array(required(record, "fractional_position")?)?;
            if position.len() != 3 {
                return Err(malformed(
                    "orbital fractional position must have length three",
                ));
            }
            position
                .iter()
                .enumerate()
                .map(|(axis, value)| {
                    approximate_element(&element_from_json(&context, value)?).and_then(
                        |(real, imaginary)| {
                            if imaginary.abs() <= 1.0e-12 {
                                Ok(real)
                            } else {
                                Err(ExactError::new(
                                    "NonnumericHopping",
                                    format!("orbital {orbital} coordinate {axis} is not real"),
                                ))
                            }
                        },
                    )
                })
                .collect::<ExactResult<Vec<_>>>()
        })
        .collect::<ExactResult<Vec<_>>>()?;
    let lattice = matrix_from_json(&context, required(first, "evaluated_lattice")?)?;
    if lattice.rows() != 3 || lattice.columns() != 3 {
        return Err(malformed("model lattice must be 3 by 3"));
    }
    let lattice = (0..3)
        .map(|row| {
            (0..3)
                .map(|column| approximate_element(lattice.entry(row, column)?).map(|value| value.0))
                .collect::<ExactResult<Vec<_>>>()
        })
        .collect::<ExactResult<Vec<_>>>()?;
    let hopping_translations = translations.keys().copied().collect::<Vec<_>>();
    let hopping_matrices = translations
        .values()
        .map(|matrix| {
            Value::Array(
                (0..dimension)
                    .map(|row| {
                        Value::Array(
                            (0..dimension)
                                .map(|column| {
                                    let value = matrix[(row, column)];
                                    json!({"real": value.re, "imaginary": value.im})
                                })
                                .collect(),
                        )
                    })
                    .collect(),
            )
        })
        .collect::<Vec<_>>();
    Ok(json!({
        "Schema": "MagneticTBHoppingData",
        "SchemaVersion": 1,
        "Source": "MagneticTBRealSpaceShellCache",
        "NumWannier": dimension,
        "NumTranslations": hopping_translations.len(),
        "Translations": hopping_translations,
        "Degeneracies": vec![1; hopping_matrices.len()],
        "HoppingMatrices": hopping_matrices,
        "Shells": shells,
        "Parameters": used_parameter_names,
        "Lattice": lattice,
        "WannierCenters": centers,
        "Convention": "H(R)[a,b]=<0,a|H|R,b>; H(k)=Sum_R H(R) Exp[2 Pi I k.R]"
    }))
}

#[allow(clippy::too_many_lines)]
fn hopping_data_from_symbolic_matrix(
    request: &serde_json::Map<String, Value>,
) -> ExactResult<Value> {
    const MAX_SAFE_TRANSLATION: f64 = 9_007_199_254_740_991.0;

    let context = context_from_json(required(request, "context")?, 128)?;
    let rows = array(required(request, "matrix")?)?;
    if rows.is_empty() {
        return Err(ExactError::new(
            "InvalidHamiltonian",
            "the hand-written Hamiltonian must be nonempty and square",
        ));
    }
    let dimension = rows.len();
    if rows
        .iter()
        .any(|row| array(row).map_or(true, |row| row.len() != dimension))
    {
        return Err(ExactError::new(
            "InvalidHamiltonian",
            "the hand-written Hamiltonian must be nonempty and square",
        ));
    }
    let centers = array(required(request, "wannier_centers")?)?;
    if centers.len() != dimension {
        return Err(ExactError::new(
            "InvalidWannierCenters",
            "one Wannier center is required per Hamiltonian row",
        ));
    }
    let centers = centers
        .iter()
        .enumerate()
        .map(|(index, value)| {
            let values = array(value)?;
            if values.len() != 3 {
                return Err(ExactError::new(
                    "InvalidWannierCenters",
                    format!("Wannier center {index} must have length three"),
                ));
            }
            Ok([
                numeric_real(&context, &values[0], "Wannier center")?,
                numeric_real(&context, &values[1], "Wannier center")?,
                numeric_real(&context, &values[2], "Wannier center")?,
            ])
        })
        .collect::<ExactResult<Vec<_>>>()?;
    let parameters = object(required(request, "parameters")?)?;
    let tolerance = request
        .get("translation_tolerance")
        .and_then(Value::as_f64)
        .unwrap_or(1.0e-9);
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(ExactError::new(
            "InvalidWannier90Option",
            "TranslationTolerance must be finite and positive",
        ));
    }
    let mut translations = BTreeMap::<[i64; 3], DMatrix<Complex<f64>>>::new();
    for (row_index, row) in rows.iter().enumerate() {
        for (column_index, cell) in array(row)?.iter().enumerate() {
            for contribution in array(cell)? {
                let contribution = object(contribution)?;
                let parameter_name = required(contribution, "parameter")?
                    .as_str()
                    .ok_or_else(|| malformed("contribution parameter must be a string"))?;
                let parameter = parameters.get(parameter_name).ok_or_else(|| {
                    ExactError::new(
                        "NonnumericHopping",
                        format!("missing numeric value for Hamiltonian parameter {parameter_name}"),
                    )
                })?;
                let parameter = numeric_complex(&context, parameter, parameter_name)?;
                let coefficient = approximate_element(&encoded_element_from_json(
                    &context,
                    required(contribution, "coefficient")?,
                )?)?;
                let displacement = array(required(contribution, "displacement")?)?;
                if displacement.len() != 3 {
                    return Err(malformed("Hamiltonian displacement must have length three"));
                }
                let displacement = displacement
                    .iter()
                    .map(quadratic_from_json)
                    .map(|value| {
                        value.and_then(|value| {
                            value
                                .to_f64()
                                .map_err(|error| ExactError::new(error.tag(), error.to_string()))
                        })
                    })
                    .collect::<ExactResult<Vec<_>>>()?;
                let mut translation = [0_i64; 3];
                for axis in 0..3 {
                    let raw =
                        displacement[axis] + centers[row_index][axis] - centers[column_index][axis];
                    let rounded = raw.round();
                    if !raw.is_finite()
                        || (raw - rounded).abs() > tolerance
                        || rounded.abs() > MAX_SAFE_TRANSLATION
                    {
                        return Err(ExactError::new(
                            "NonintegerTranslation",
                            format!(
                                "Hamiltonian element ({},{}) gives noninteger translation coordinate {raw}",
                                row_index + 1,
                                column_index + 1
                            ),
                        ));
                    }
                    translation[axis] = format!("{rounded:.0}")
                        .parse::<i64>()
                        .map_err(|_| malformed("integer translation cannot be represented"))?;
                }
                let value = complex_multiply(parameter, coefficient);
                if value.0.hypot(value.1) <= tolerance {
                    continue;
                }
                let matrix = translations
                    .entry(translation)
                    .or_insert_with(|| DMatrix::zeros(dimension, dimension));
                matrix[(row_index, column_index)] += Complex::new(value.0, value.1);
            }
        }
    }
    if translations.is_empty() {
        translations.insert([0, 0, 0], DMatrix::zeros(dimension, dimension));
    }
    Ok(json!({
        "Schema": "Wannier90HRData",
        "SchemaVersion": 1,
        "Source": "HamiltonianCoefficientExtraction",
        "NumWannier": dimension,
        "NumTranslations": translations.len(),
        "Translations": translations.keys().collect::<Vec<_>>(),
        "Degeneracies": vec![1; translations.len()],
        "HoppingMatrices": translations.values().map(numeric_matrix_json).collect::<Vec<_>>()
    }))
}

fn validate_symbolic_expression_rows(rows: &[Value]) -> ExactResult<()> {
    if rows.is_empty()
        || rows
            .iter()
            .any(|row| array(row).map_or(true, |row| row.len() != rows.len()))
    {
        return Err(ExactError::new(
            "InvalidHamiltonian",
            "Hamiltonian must be a nonempty square matrix",
        ));
    }
    Ok(())
}

fn evaluate_symbolic_expression_rows(
    context: &Arc<CyclotomicContext>,
    rows: &[Value],
    parameters: &serde_json::Map<String, Value>,
    momentum: &[f64],
) -> ExactResult<Vec<Vec<Value>>> {
    if momentum.len() != 3 {
        return Err(ExactError::new(
            "InvalidMomentum",
            "momentum must contain kx, ky, kz",
        ));
    }
    let mut output = Vec::with_capacity(rows.len());
    for row in rows {
        let mut output_row = Vec::with_capacity(rows.len());
        for cell in array(row)? {
            let mut sum = (0.0, 0.0);
            for contribution in array(cell)? {
                let contribution = object(contribution)?;
                let name = required(contribution, "parameter")?
                    .as_str()
                    .ok_or_else(|| malformed("contribution parameter must be a string"))?;
                let parameter = parameters.get(name).ok_or_else(|| {
                    ExactError::new(
                        "InvalidHamiltonianParameters",
                        format!("missing Hamiltonian parameter {name}"),
                    )
                })?;
                let parameter = numeric_complex(context, parameter, name)?;
                let coefficient = approximate_element(&element_from_json(
                    context,
                    required(contribution, "coefficient")?,
                )?)?;
                let displacement = array(required(contribution, "displacement")?)?;
                if displacement.len() != momentum.len() {
                    return Err(malformed("Hamiltonian displacement must have length three"));
                }
                let phase = displacement.iter().zip(momentum.iter()).try_fold(
                    0.0,
                    |sum, (coordinate, momentum)| {
                        quadratic_from_json(coordinate)?.to_f64().map_or_else(
                            |error| Err(ExactError::new(error.tag(), error.to_string())),
                            |coordinate| Ok(sum + coordinate * momentum),
                        )
                    },
                )?;
                sum = complex_add(
                    sum,
                    complex_multiply(
                        parameter,
                        complex_multiply(coefficient, (phase.cos(), phase.sin())),
                    ),
                );
            }
            let sum = finite_complex(sum.0, sum.1, "Hamiltonian matrix entry")?;
            output_row.push(json!({"real": sum.0, "imaginary": sum.1}));
        }
        output.push(output_row);
    }
    Ok(output)
}

fn evaluate_symbolic_expression_matrix(
    request: &serde_json::Map<String, Value>,
) -> ExactResult<Value> {
    let context = context_from_json(required(request, "context")?, 128)?;
    let rows = array(required(request, "matrix")?)?;
    validate_symbolic_expression_rows(rows)?;
    let parameters = object(required(request, "parameters")?)?;
    let momentum = array(required(request, "momentum")?)?;
    let momentum = momentum
        .iter()
        .enumerate()
        .map(|(index, value)| numeric_real(&context, value, ["kx", "ky", "kz"][index]))
        .collect::<ExactResult<Vec<_>>>()?;
    Ok(json!({
        "matrix": evaluate_symbolic_expression_rows(&context, rows, parameters, &momentum)?
    }))
}

fn evaluate_symbolic_expression_matrices(
    request: &serde_json::Map<String, Value>,
) -> ExactResult<Value> {
    let context = context_from_json(required(request, "context")?, 128)?;
    let rows = array(required(request, "matrix")?)?;
    validate_symbolic_expression_rows(rows)?;
    let parameters = object(required(request, "parameters")?)?;
    let momenta = array(required(request, "momenta")?)?;
    if momenta.is_empty() {
        return Err(ExactError::new(
            "InvalidBandPath",
            "momenta must contain at least one point",
        ));
    }
    let mut matrices = Vec::with_capacity(momenta.len());
    for (point_index, momentum) in momenta.iter().enumerate() {
        let momentum = array(momentum)?;
        if momentum.len() != 3 {
            return Err(ExactError::new(
                "InvalidMomentum",
                format!("momenta[{point_index}] must contain kx, ky, kz"),
            ));
        }
        let momentum = momentum
            .iter()
            .enumerate()
            .map(|(index, value)| {
                numeric_real(&context, value, &format!("momenta[{point_index}][{index}]"))
            })
            .collect::<ExactResult<Vec<_>>>()?;
        matrices.push(evaluate_symbolic_expression_rows(
            &context, rows, parameters, &momentum,
        )?);
    }
    Ok(json!({"matrices": matrices}))
}
