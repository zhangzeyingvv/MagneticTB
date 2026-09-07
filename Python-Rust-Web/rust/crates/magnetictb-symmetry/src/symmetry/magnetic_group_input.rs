/// Compiles the ordered affine operations of one frozen MSG record into the exact
/// Seitz quotient group used by the rest of the Rust core.
pub fn compile_msg_group(group: &MsgGroup) -> SymmetryResult<CompiledSeitzGroup> {
    let context = make_context(1, 128)?;
    compile_msg_group_in_context(group, &context)
}
/// Compiles a frozen MSG record in a caller-selected cyclotomic superfield.
///
/// MSG affine data are rational, so this preserves their exact values while
/// allowing downstream orbital and hopping matrices to share a field such as
/// Q(i) or a larger cyclotomic extension.
pub fn compile_msg_group_in_context(
    group: &MsgGroup,
    context: &Arc<CyclotomicContext>,
) -> SymmetryResult<CompiledSeitzGroup> {
    let operations = group
        .operations()
        .iter()
        .map(|operation| data_operation(context, operation))
        .collect::<SymmetryResult<Vec<_>>>()?;
    compile_seitz_group(operations)
}

fn data_operation(
    context: &Arc<CyclotomicContext>,
    operation: &DataSymmetryOperation,
) -> SymmetryResult<SeitzOperation> {
    let rotation = data_matrix(context, operation.rotation())?;
    let translation = operation
        .translation()
        .iter()
        .map(|value| data_element(context, value))
        .collect::<SymmetryResult<Vec<_>>>()?;
    Ok(SeitzOperation::new(
        operation.stable_id(),
        operation.label(),
        SpatialOperation::new(rotation, translation)?,
        operation.antiunitary(),
    ))
}

fn data_matrix(
    context: &Arc<CyclotomicContext>,
    matrix: &DataExactMatrix,
) -> SymmetryResult<ExactMatrix> {
    ExactMatrix::new(
        Arc::clone(context),
        matrix.rows(),
        matrix.columns(),
        matrix
            .entries()
            .iter()
            .map(|value| data_element(context, value))
            .collect::<SymmetryResult<Vec<_>>>()?,
    )
    .map_err(Into::into)
}

fn data_element(
    context: &Arc<CyclotomicContext>,
    value: &DataExactExpression,
) -> SymmetryResult<Element> {
    let rational = match value {
        DataExactExpression::Integer { value } => Rational::parse(value, "1")?,
        DataExactExpression::Rational {
            numerator,
            denominator,
        } => Rational::parse(numerator, denominator)?,
        DataExactExpression::Symbol { .. }
        | DataExactExpression::RootOfUnity { .. }
        | DataExactExpression::Call { .. } => {
            return Err(SymmetryError::new(
                "UnsupportedDataExactExpression",
                "MSG rotation and translation fields must be exact rational numbers",
            ));
        }
    };
    Element::from_polynomial(Arc::clone(context), &[rational]).map_err(Into::into)
}

/// Evaluates one tagged expression from the frozen Data snapshot in an exact
/// cyclotomic context. Unsupported heads and missing symbols fail explicitly.
pub fn evaluate_data_expression(
    context: &Arc<CyclotomicContext>,
    value: &DataExactExpression,
    bindings: &BTreeMap<String, Element>,
) -> SymmetryResult<Element> {
    match value {
        DataExactExpression::Integer { .. } | DataExactExpression::Rational { .. } => {
            data_element(context, value)
        }
        DataExactExpression::Symbol { name } => {
            let result = bindings.get(name).ok_or_else(|| {
                SymmetryError::new(
                    "MissingExactSymbol",
                    format!("no exact binding was supplied for {name}"),
                )
            })?;
            if result.context() != context {
                return Err(SymmetryError::new(
                    "ConductorMismatch",
                    format!("binding for {name} uses a different exact field"),
                ));
            }
            Ok(result.clone())
        }
        DataExactExpression::RootOfUnity { order, power } => {
            Element::root_of_unity(context, *order, *power).map_err(Into::into)
        }
        DataExactExpression::Call { head, arguments } if head == "System`Plus" => {
            let mut result = Element::zero(context)?;
            for argument in arguments {
                result = result.add(&evaluate_data_expression(context, argument, bindings)?)?;
            }
            Ok(result)
        }
        DataExactExpression::Call { head, arguments } if head == "System`Times" => {
            let mut result = Element::one(context)?;
            for argument in arguments {
                result =
                    result.multiply(&evaluate_data_expression(context, argument, bindings)?)?;
            }
            Ok(result)
        }
        DataExactExpression::Call { head, arguments }
            if head == "System`Power" && arguments.len() == 2 =>
        {
            let exponent = match &arguments[1] {
                DataExactExpression::Integer { value } => value.parse::<i64>().map_err(|_| {
                    SymmetryError::new(
                        "InvalidExactExponent",
                        "integer exponent is outside the supported i64 range",
                    )
                })?,
                _ => {
                    return Err(SymmetryError::new(
                        "InvalidExactExponent",
                        "cyclotomic expression powers require an exact integer exponent",
                    ));
                }
            };
            evaluate_data_expression(context, &arguments[0], bindings)?
                .power(exponent)
                .map_err(Into::into)
        }
        DataExactExpression::Call { head, .. } => Err(SymmetryError::new(
            "UnsupportedExactCall",
            format!("unsupported exact Data expression head {head}"),
        )),
    }
}

fn matrix_vector(matrix: &ExactMatrix, vector: &[Element]) -> SymmetryResult<Vec<Element>> {
    if matrix.columns() != vector.len()
        || vector
            .iter()
            .any(|value| value.context() != matrix.context())
    {
        return Err(SymmetryError::new(
            "SpatialDimensionMismatch",
            "matrix-vector dimensions or exact fields differ",
        ));
    }
    let mut result = Vec::with_capacity(matrix.rows());
    for row in 0..matrix.rows() {
        let mut value = Element::zero(matrix.context())?;
        for (column, coordinate) in vector.iter().enumerate() {
            value = value.add(&matrix.entry(row, column)?.multiply(coordinate)?)?;
        }
        result.push(value);
    }
    Ok(result)
}

fn add_vectors(left: &[Element], right: &[Element]) -> SymmetryResult<Vec<Element>> {
    if left.len() != right.len() {
        return Err(SymmetryError::new(
            "SpatialDimensionMismatch",
            "translation dimensions differ",
        ));
    }
    left.iter()
        .zip(right)
        .map(|(left, right)| left.add(right).map_err(Into::into))
        .collect()
}

fn subtract_vectors(left: &[Element], right: &[Element]) -> SymmetryResult<Vec<Element>> {
    if left.len() != right.len() {
        return Err(SymmetryError::new(
            "SpatialDimensionMismatch",
            "translation dimensions differ",
        ));
    }
    left.iter()
        .zip(right)
        .map(|(left, right)| left.subtract(right).map_err(Into::into))
        .collect()
}
