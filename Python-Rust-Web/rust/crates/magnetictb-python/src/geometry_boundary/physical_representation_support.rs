fn periodic_bond_to_json(bond: &PeriodicBondRecord) -> Value {
    json!({
        "source_site": bond.source_site(),
        "target_site": bond.target_site(),
        "translation": bond.translation(),
        "displacement": vector_to_json(bond.displacement()),
        "source_endpoint": vector_to_json(bond.source_endpoint()),
        "target_endpoint": vector_to_json(bond.target_endpoint()),
        "squared_distance": quadratic_to_json(bond.squared_distance())
    })
}
fn exact_quadratic_rows(value: &Value, field: &str) -> ExactResult<Vec<Vec<QuadraticReal>>> {
    let encoded: Vec<Vec<ExactExpression>> = from_value(value, field)?;
    let bindings = BTreeMap::new();
    encoded
        .iter()
        .map(|row| {
            row.iter()
                .map(|coordinate| QuadraticReal::from_expression(coordinate, &bindings))
                .collect::<Result<Vec<_>, _>>()
        })
        .collect::<Result<Vec<_>, _>>()
        .map_err(|error| ExactError::new(error.tag(), error.to_string()))
}

fn exact_elements(context: &Arc<CyclotomicContext>, value: &Value) -> ExactResult<Vec<Element>> {
    array(value)?
        .iter()
        .map(|element| encoded_element_from_json(context, element))
        .collect()
}

fn exact_site_orbits(
    context: &Arc<CyclotomicContext>,
    value: &Value,
) -> ExactResult<Vec<Vec<Vec<Element>>>> {
    array(value)?
        .iter()
        .map(|orbit| {
            array(orbit)?
                .iter()
                .map(|site| exact_elements(context, site))
                .collect()
        })
        .collect()
}

fn exact_matrices(
    context: &Arc<CyclotomicContext>,
    value: &Value,
) -> ExactResult<Vec<ExactMatrix>> {
    array(value)?
        .iter()
        .map(|matrix| encoded_matrix_from_json(context, matrix))
        .collect()
}

fn polynomial_expression(
    context: &Arc<CyclotomicContext>,
    value: &Value,
) -> ExactResult<PolynomialExpression> {
    let record = object(value)?;
    match string_field(record, "kind")? {
        "symbol" => {
            let name = string_field(record, "name")?;
            let variable = name.rsplit('`').next().unwrap_or(name);
            match variable {
                "x" => Ok(PolynomialExpression::Variable(0)),
                "y" => Ok(PolynomialExpression::Variable(1)),
                "z" => Ok(PolynomialExpression::Variable(2)),
                _ => Err(ExactError::new(
                    "UnsupportedFunctionBasisExpression",
                    format!("unsupported polynomial variable {name}"),
                )),
            }
        }
        "call" => {
            let arguments = array(required(record, "arguments")?)?;
            match string_field(record, "head")? {
                "System`Plus" => Ok(PolynomialExpression::Add(
                    arguments
                        .iter()
                        .map(|argument| polynomial_expression(context, argument))
                        .collect::<ExactResult<Vec<_>>>()?,
                )),
                "System`Times" => Ok(PolynomialExpression::Multiply(
                    arguments
                        .iter()
                        .map(|argument| polynomial_expression(context, argument))
                        .collect::<ExactResult<Vec<_>>>()?,
                )),
                "System`Power" if arguments.len() == 2 => {
                    let exponent = object(&arguments[1])?;
                    if string_field(exponent, "kind")? != "integer" {
                        return Err(ExactError::new(
                            "UnsupportedFunctionBasisExpression",
                            "polynomial powers require a nonnegative integer exponent",
                        ));
                    }
                    let exponent =
                        string_field(exponent, "value")?
                            .parse::<u8>()
                            .map_err(|_| {
                                ExactError::new(
                                    "UnsupportedFunctionBasisExpression",
                                    "polynomial exponent is outside the supported u8 range",
                                )
                            })?;
                    Ok(PolynomialExpression::Power(
                        Box::new(polynomial_expression(context, &arguments[0])?),
                        exponent,
                    ))
                }
                "System`Complex" => Ok(PolynomialExpression::Scalar(encoded_element_from_json(
                    context, value,
                )?)),
                head => Err(ExactError::new(
                    "UnsupportedFunctionBasisExpression",
                    format!("unsupported exact polynomial head {head}"),
                )),
            }
        }
        _ => Ok(PolynomialExpression::Scalar(encoded_element_from_json(
            context, value,
        )?)),
    }
}

fn exact_function_basis(
    context: &Arc<CyclotomicContext>,
    value: &Value,
) -> ExactResult<Vec<Vec<PolynomialExpression>>> {
    let record = object(value)?;
    if string_field(record, "kind")? != "matrix" {
        return Err(ExactError::new(
            "InvalidFunctionBasis",
            "an explicit function basis must use one tagged matrix per orbit",
        ));
    }
    let rows = size_field(record, "rows")?;
    let columns = size_field(record, "columns")?;
    let entries = array(required(record, "entries")?)?;
    if rows == 0 || !matches!(columns, 1 | 2) || entries.len() != rows.saturating_mul(columns) {
        return Err(ExactError::new(
            "InvalidFunctionBasis",
            "an explicit function basis must be a nonempty scalar or two-component matrix",
        ));
    }
    entries
        .chunks(columns)
        .map(|row| {
            row.iter()
                .map(|entry| polynomial_expression(context, entry))
                .collect()
        })
        .collect()
}

fn compile_basis_item(
    context: &Arc<CyclotomicContext>,
    basis: &Value,
    spatial: &[SpatialOperation],
    spin_rotations: Option<&[ExactMatrix]>,
    antiunitary_flags: &[bool],
    lattice: &ExactMatrix,
) -> ExactResult<CompiledBasisAction> {
    if basis.is_array() {
        let labels: Vec<String> = from_value(basis, "basis_functions orbit")?;
        compile_catalog_basis_action(&labels, spatial, spin_rotations, antiunitary_flags, lattice)
            .map_err(|error| ExactError::new(error.tag(), error.to_string()))
    } else if basis.is_object() {
        let expressions = exact_function_basis(context, basis)?;
        compile_function_basis_action(
            &expressions,
            spatial,
            spin_rotations,
            antiunitary_flags,
            lattice,
        )
        .map_err(|error| ExactError::new(error.tag(), error.to_string()))
    } else {
        Err(ExactError::new(
            "InvalidFunctionBasis",
            "each orbit must use ordered catalog labels or one exact polynomial matrix",
        ))
    }
}

#[allow(clippy::too_many_arguments)]
fn compile_induced_basis_item(
    context: &Arc<CyclotomicContext>,
    basis: &Value,
    spatial: &[SpatialOperation],
    spin_rotations: Option<&[ExactMatrix]>,
    antiunitary_flags: &[bool],
    image_site_indices: &[Vec<usize>],
    coset_representatives: &[usize],
    lattice: &ExactMatrix,
) -> ExactResult<CompiledInducedBasisAction> {
    if basis.is_array() {
        let labels: Vec<String> = from_value(basis, "basis_functions orbit")?;
        compile_catalog_induced_basis_action(
            &labels,
            spatial,
            spin_rotations,
            antiunitary_flags,
            image_site_indices,
            coset_representatives,
            lattice,
        )
        .map_err(|error| ExactError::new(error.tag(), error.to_string()))
    } else if basis.is_object() {
        let expressions = exact_function_basis(context, basis)?;
        compile_function_induced_basis_action(
            &expressions,
            spatial,
            spin_rotations,
            antiunitary_flags,
            image_site_indices,
            coset_representatives,
            lattice,
        )
        .map_err(|error| ExactError::new(error.tag(), error.to_string()))
    } else {
        Err(ExactError::new(
            "InvalidFunctionBasis",
            "each orbit must use ordered catalog labels or one exact polynomial matrix",
        ))
    }
}

fn polynomial_terms_to_json(
    terms: &[magnetictb_symmetry::PolynomialTerm],
    one: &Element,
) -> Value {
    let mut sum = Vec::with_capacity(terms.len());
    for term in terms {
        if term.coefficient().is_zero() {
            continue;
        }
        let mut factors = Vec::new();
        let exponents = term.exponents();
        let has_monomial = exponents.iter().any(|&value| value != 0);
        if !has_monomial || term.coefficient() != one {
            factors.push(cyclotomic_expression_from_element(term.coefficient()));
        }
        for (name, exponent) in ["x", "y", "z"].into_iter().zip(exponents) {
            if *exponent == 0 {
                continue;
            }
            let variable = json!({"kind": "symbol", "name": format!("MagneticTB`{name}")});
            factors.push(if *exponent == 1 {
                variable
            } else {
                json!({
                    "kind": "call",
                    "head": "System`Power",
                    "arguments": [
                        variable,
                        {"kind": "integer", "value": exponent.to_string()}
                    ]
                })
            });
        }
        sum.push(match factors.len() {
            0 => json!({"kind": "integer", "value": "1"}),
            1 => factors.remove(0),
            _ => json!({
                "kind": "call",
                "head": "System`Times",
                "arguments": factors
            }),
        });
    }
    match sum.len() {
        0 => json!({"kind": "integer", "value": "0"}),
        1 => sum.remove(0),
        _ => json!({
            "kind": "call",
            "head": "System`Plus",
            "arguments": sum
        }),
    }
}

fn polynomial_basis_state_to_json(state: &FunctionBasisState) -> ExactResult<Value> {
    let one = Element::one(
        state
            .components()
            .iter()
            .flatten()
            .next()
            .map(|term| term.coefficient().context())
            .ok_or_else(|| {
                ExactError::new(
                    "InvalidFunctionBasis",
                    "transported basis state contains no exact coefficient",
                )
            })?,
    )?;
    let components = state
        .components()
        .iter()
        .map(|component| polynomial_terms_to_json(component, &one))
        .collect::<Vec<_>>();
    Ok(if components.len() == 1 {
        components.into_iter().next().expect("one component")
    } else {
        json!({"kind": "list", "items": components})
    })
}

fn proportional_polynomial_factor(
    reference: &[magnetictb_symmetry::PolynomialTerm],
    component: &[magnetictb_symmetry::PolynomialTerm],
) -> ExactResult<Option<Element>> {
    let reference_term = reference
        .iter()
        .find(|term| !term.coefficient().is_zero())
        .ok_or_else(|| {
        ExactError::new(
            "InvalidFunctionBasis",
            "spinor factorization requires a nonzero reference component",
        )
    })?;
    if component.iter().all(|term| term.coefficient().is_zero()) {
        return Ok(Some(Element::zero(reference_term.coefficient().context())?));
    }
    let Some(component_term) = component
        .iter()
        .find(|term| {
            !term.coefficient().is_zero() && term.exponents() == reference_term.exponents()
        })
    else {
        return Ok(None);
    };
    let factor = component_term
        .coefficient()
        .divide(reference_term.coefficient())?;
    if reference
        .iter()
        .filter(|term| !term.coefficient().is_zero())
        .count()
        != component
            .iter()
            .filter(|term| !term.coefficient().is_zero())
            .count()
    {
        return Ok(None);
    }
    for reference_value in reference
        .iter()
        .filter(|term| !term.coefficient().is_zero())
    {
        let Some(component_value) = component
            .iter()
            .find(|term| {
                !term.coefficient().is_zero()
                    && term.exponents() == reference_value.exponents()
            })
        else {
            return Ok(None);
        };
        if component_value.coefficient() != &reference_value.coefficient().multiply(&factor)? {
            return Ok(None);
        }
    }
    Ok(Some(factor))
}

fn basis_state_record_to_json(state: &FunctionBasisState) -> ExactResult<Value> {
    let basis_state = polynomial_basis_state_to_json(state)?;
    let components = state.components();
    if components.len() == 1 {
        return Ok(json!({
            "basis_state": &basis_state,
            "spatial_orbital": &basis_state,
            "spin_state": Value::Null,
            "spin_structure": "Spinless"
        }));
    }
    if components.len() != 2 {
        return Ok(json!({
            "basis_state": &basis_state,
            "spatial_orbital": Value::Null,
            "spin_state": &basis_state,
            "spin_structure": "GeneralInternalState"
        }));
    }
    let Some(reference) = components
        .iter()
        .find(|component| component.iter().any(|term| !term.coefficient().is_zero()))
    else {
        return Ok(json!({
            "basis_state": &basis_state,
            "spatial_orbital": Value::Null,
            "spin_state": &basis_state,
            "spin_structure": "SpatialSpinor"
        }));
    };
    let factors = components
        .iter()
        .map(|component| proportional_polynomial_factor(reference, component))
        .collect::<ExactResult<Option<Vec<_>>>>()?;
    if let Some(factors) = factors {
        let one = Element::one(reference[0].coefficient().context())?;
        return Ok(json!({
            "basis_state": &basis_state,
            "spatial_orbital": polynomial_terms_to_json(reference, &one),
            "spin_state": {
                "kind": "list",
                "items": factors.iter().map(cyclotomic_expression_from_element).collect::<Vec<_>>()
            },
            "spin_structure": "FactorizedSpinor"
        }));
    }
    Ok(json!({
        "basis_state": &basis_state,
        "spatial_orbital": Value::Null,
        "spin_state": &basis_state,
        "spin_structure": "SpatialSpinor"
    }))
}

fn basis_states_to_json(states: &[FunctionBasisState]) -> ExactResult<Vec<Value>> {
    states.iter().map(basis_state_record_to_json).collect()
}

fn resolve_basis_item_states(
    context: &Arc<CyclotomicContext>,
    basis: &Value,
    lattice: &ExactMatrix,
) -> ExactResult<Vec<Value>> {
    let states = if basis.is_array() {
        let labels: Vec<String> = from_value(basis, "basis_functions orbit")?;
        resolve_catalog_basis_states(&labels, lattice)
    } else if basis.is_object() {
        resolve_function_basis_states(&exact_function_basis(context, basis)?, lattice)
    } else {
        return Err(ExactError::new(
            "InvalidFunctionBasis",
            "each orbit must use ordered catalog labels or one exact polynomial matrix",
        ));
    }
    .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
    basis_states_to_json(&states)
}

fn point_matrix_result(request: &serde_json::Map<String, Value>) -> ExactResult<Value> {
    let context = context_from_json(required(request, "context")?, 128)?;
    let spatial = spatial_operations(&context, required(request, "spatial_actions")?)?;
    let antiunitary_flags: Vec<bool> =
        from_value(required(request, "antiunitary_flags")?, "antiunitary_flags")?;
    let lattice = encoded_matrix_from_json(&context, required(request, "lattice")?)?;
    let basis = required(request, "basis")?;
    let compiled = if basis.is_array() {
        let labels: Vec<String> = from_value(basis, "basis")?;
        compile_catalog_basis_action(&labels, &spatial, None, &antiunitary_flags, &lattice)
    } else {
        compile_function_basis_action(
            &exact_function_basis(&context, basis)?,
            &spatial,
            None,
            &antiunitary_flags,
            &lattice,
        )
    }
    .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
    Ok(json!({
        "matrices": compiled
            .local_matrices()
            .iter()
            .map(matrix_to_json)
            .collect::<Vec<_>>()
    }))
}

struct DualExpression {
    value: Element,
    derivative: Element,
}

fn dual_expression(
    context: &Arc<CyclotomicContext>,
    value: &Value,
    parameter: &str,
) -> ExactResult<DualExpression> {
    if let Some(record) = value.as_object() {
        match record.get("kind").and_then(Value::as_str) {
            Some("symbol") => {
                let name = string_field(record, "name")?;
                if name == parameter {
                    return Ok(DualExpression {
                        value: Element::zero(context)?,
                        derivative: Element::one(context)?,
                    });
                }
                return Err(ExactError::new(
                    "UnsupportedContinuousExpression",
                    format!("unsupported continuous symbol {name}"),
                ));
            }
            Some("call") => return dual_call(context, record, parameter),
            _ => {}
        }
    }
    Ok(DualExpression {
        value: encoded_element_from_json(context, value)?,
        derivative: Element::zero(context)?,
    })
}

fn dual_call(
    context: &Arc<CyclotomicContext>,
    record: &serde_json::Map<String, Value>,
    parameter: &str,
) -> ExactResult<DualExpression> {
    let arguments = array(required(record, "arguments")?)?;
    match string_field(record, "head")? {
        "System`Plus" => arguments.iter().try_fold(
            DualExpression {
                value: Element::zero(context)?,
                derivative: Element::zero(context)?,
            },
            |total, argument| {
                let argument = dual_expression(context, argument, parameter)?;
                Ok(DualExpression {
                    value: total.value.add(&argument.value)?,
                    derivative: total.derivative.add(&argument.derivative)?,
                })
            },
        ),
        "System`Times" => arguments.iter().try_fold(
            DualExpression {
                value: Element::one(context)?,
                derivative: Element::zero(context)?,
            },
            |total, argument| {
                let argument = dual_expression(context, argument, parameter)?;
                Ok(DualExpression {
                    derivative: total
                        .derivative
                        .multiply(&argument.value)?
                        .add(&total.value.multiply(&argument.derivative)?)?,
                    value: total.value.multiply(&argument.value)?,
                })
            },
        ),
        "System`Complex" if arguments.len() == 2 => Ok(DualExpression {
            value: encoded_element_from_json(context, &Value::Object(record.clone()))?,
            derivative: Element::zero(context)?,
        }),
        "System`Sin" | "System`Cos" if arguments.len() == 1 => {
            let argument = dual_expression(context, &arguments[0], parameter)?;
            if !argument.value.is_zero() {
                return Err(ExactError::new(
                    "UnsupportedContinuousExpression",
                    "Sin/Cos arguments must vanish at the continuous identity",
                ));
            }
            if string_field(record, "head")? == "System`Sin" {
                Ok(DualExpression {
                    value: Element::zero(context)?,
                    derivative: argument.derivative,
                })
            } else {
                Ok(DualExpression {
                    value: Element::one(context)?,
                    derivative: Element::zero(context)?,
                })
            }
        }
        "System`Power" if arguments.len() == 2 => {
            let base = object(&arguments[0])?;
            if base.get("kind").and_then(Value::as_str) != Some("symbol")
                || string_field(base, "name")? != "System`E"
            {
                return Err(ExactError::new(
                    "UnsupportedContinuousExpression",
                    "continuous powers must use the exact exponential base E",
                ));
            }
            let exponent = dual_expression(context, &arguments[1], parameter)?;
            if !exponent.value.is_zero() {
                return Err(ExactError::new(
                    "InvalidContinuousRepresentation",
                    "continuous exponential is not identity at parameter zero",
                ));
            }
            Ok(DualExpression {
                value: Element::one(context)?,
                derivative: exponent.derivative,
            })
        }
        head => Err(ExactError::new(
            "UnsupportedContinuousExpression",
            format!("unsupported exact continuous expression head {head}"),
        )),
    }
}

fn continuous_matrix_generator(
    context: &Arc<CyclotomicContext>,
    value: &Value,
    parameter: &str,
) -> ExactResult<ExactMatrix> {
    let record = object(value)?;
    let rows: usize = from_value(required(record, "rows")?, "rows")?;
    let columns: usize = from_value(required(record, "columns")?, "columns")?;
    if rows == 0 || rows != columns {
        return Err(ExactError::new(
            "InvalidContinuousRepresentation",
            "continuous site matrix must be nonempty and square",
        ));
    }
    let duals = array(required(record, "entries")?)?
        .iter()
        .map(|entry| dual_expression(context, entry, parameter))
        .collect::<ExactResult<Vec<_>>>()?;
    let at_identity = ExactMatrix::new(
        Arc::clone(context),
        rows,
        columns,
        duals.iter().map(|entry| entry.value.clone()).collect(),
    )?;
    if !exact_matrix_equal(&at_identity, &ExactMatrix::identity(context, rows)?) {
        return Err(ExactError::new(
            "InvalidContinuousRepresentation",
            "continuous site matrix is not identity at parameter zero",
        ));
    }
    let imaginary = Element::root_of_unity(context, 4, 1)?;
    let generator = ExactMatrix::new(
        Arc::clone(context),
        rows,
        columns,
        duals
            .iter()
            .map(|entry| imaginary.multiply(&entry.derivative))
            .collect::<ExactResult<Vec<_>>>()?,
    )?;
    if !exact_hermitian_matrix(&generator)? {
        return Err(ExactError::new(
            "InvalidContinuousRepresentation",
            "I times the exact derivative must be Hermitian",
        ));
    }
    Ok(generator)
}

fn continuous_site_generators(
    context: &Arc<CyclotomicContext>,
    value: &Value,
    parameter: &str,
) -> ExactResult<Vec<ExactMatrix>> {
    let mut result = Vec::new();
    for orbit in array(value)? {
        for matrix in array(orbit)? {
            result.push(continuous_matrix_generator(context, matrix, parameter)?);
        }
    }
    Ok(result)
}

fn exact_matrix_cube(
    context: &Arc<CyclotomicContext>,
    value: &Value,
) -> ExactResult<Vec<Vec<Vec<ExactMatrix>>>> {
    array(value)?
        .iter()
        .map(|orbit| {
            array(orbit)?
                .iter()
                .map(|operation| exact_matrices(context, operation))
                .collect()
        })
        .collect()
}
