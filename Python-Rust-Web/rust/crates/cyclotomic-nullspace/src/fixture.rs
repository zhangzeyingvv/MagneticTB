//! Direct port of the canonical data boundary in `FixtureIO.wl`.

use crate::error::{ExactResult, fail};
use crate::exact::{CyclotomicContext, Element, Polynomial, Rational, make_context};
use crate::matrix::{CommonKernelResult, CompiledProblem, ExactMatrix, KernelResult};
use num_bigint::{BigInt, Sign};
use serde_json::{Map, Value, json};
use std::fs;
use std::path::Path;
use std::str::FromStr;
use std::sync::Arc;

pub fn read_json(path: impl AsRef<Path>) -> ExactResult<Value> {
    let text = fs::read_to_string(path)
        .map_err(|error| crate::ExactError::new("MalformedSerialization", error.to_string()))?;
    serde_json::from_str(&text)
        .map_err(|error| crate::ExactError::new("MalformedSerialization", error.to_string()))
}

pub fn write_json(path: impl AsRef<Path>, value: &Value) -> ExactResult<()> {
    let text = serde_json::to_string_pretty(value)
        .map_err(|error| crate::ExactError::new("MalformedSerialization", error.to_string()))?;
    fs::write(path, format!("{text}\n"))
        .map_err(|error| crate::ExactError::new("MalformedSerialization", error.to_string()))
}

pub fn required<'a>(object: &'a Map<String, Value>, key: &str) -> ExactResult<&'a Value> {
    object.get(key).ok_or_else(|| {
        crate::ExactError::new(
            "MalformedSerialization",
            format!("missing JSON field: {key}"),
        )
    })
}

pub fn object(value: &Value) -> ExactResult<&Map<String, Value>> {
    value
        .as_object()
        .ok_or_else(|| crate::ExactError::new("MalformedSerialization", "expected object"))
}

pub fn array(value: &Value) -> ExactResult<&Vec<Value>> {
    value
        .as_array()
        .ok_or_else(|| crate::ExactError::new("MalformedSerialization", "expected array"))
}

pub fn string_field<'a>(object: &'a Map<String, Value>, key: &str) -> ExactResult<&'a str> {
    required(object, key)?.as_str().ok_or_else(|| {
        crate::ExactError::new("MalformedSerialization", format!("{key} must be a string"))
    })
}

pub fn size_value(value: &Value, field: &str) -> ExactResult<usize> {
    let raw = value.as_u64().ok_or_else(|| {
        crate::ExactError::new(
            "MalformedSerialization",
            format!("{field} must be a nonnegative integer"),
        )
    })?;
    usize::try_from(raw).map_err(|_| {
        crate::ExactError::new("MalformedSerialization", format!("{field} exceeds usize"))
    })
}

pub fn size_field(object: &Map<String, Value>, key: &str) -> ExactResult<usize> {
    size_value(required(object, key)?, key)
}

pub fn bool_field(object: &Map<String, Value>, key: &str) -> ExactResult<bool> {
    required(object, key)?.as_bool().ok_or_else(|| {
        crate::ExactError::new("MalformedSerialization", format!("{key} must be boolean"))
    })
}

pub fn rational_from_json(value: &Value) -> ExactResult<Rational> {
    let object = object(value)?;
    let numerator = string_field(object, "numerator")?;
    let denominator = string_field(object, "denominator")?;
    Rational::parse(numerator, denominator)
}

pub fn polynomial_from_json(value: &Value) -> ExactResult<Polynomial> {
    array(value)?.iter().map(rational_from_json).collect()
}

pub fn context_from_json(value: &Value, max_degree: usize) -> ExactResult<Arc<CyclotomicContext>> {
    let object = object(value)?;
    let conductor = size_field(object, "conductor")?;
    let expected_degree = size_field(object, "degree")?;
    if conductor < 1 || expected_degree < 1 {
        return fail("MalformedSerialization", "invalid context dimensions");
    }
    let expected_polynomial = polynomial_from_json(required(object, "cyclotomic_polynomial")?)?;
    let context = make_context(conductor, max_degree.max(expected_degree))?;
    if context.degree != expected_degree || context.cyclotomic_polynomial != expected_polynomial {
        return fail(
            "MalformedSerialization",
            "fixture context does not match generated context",
        );
    }
    Ok(context)
}

pub fn element_from_json(context: &Arc<CyclotomicContext>, value: &Value) -> ExactResult<Element> {
    let object = object(value)?;
    let coefficients = polynomial_from_json(required(object, "coefficients")?)?;
    if coefficients.len() != context.degree {
        return fail(
            "MalformedSerialization",
            "element coefficient count does not match context degree",
        );
    }
    Element::from_polynomial(Arc::clone(context), &coefficients)
}

/// Decodes the tagged exact-expression subset used by the `CoreDataClosure`
/// checkpoint into one exact cyclotomic element.
pub fn encoded_element_from_json(
    context: &Arc<CyclotomicContext>,
    value: &Value,
) -> ExactResult<Element> {
    if let Some(integer) = bare_json_integer(value)? {
        return Element::from_polynomial(Arc::clone(context), &[integer]);
    }
    let record = object(value)?;
    if record.contains_key("coefficients") {
        return element_from_json(context, value);
    }
    match string_field(record, "kind")? {
        "integer" => Element::from_polynomial(
            Arc::clone(context),
            &[Rational::parse(string_field(record, "value")?, "1")?],
        ),
        "rational" => Element::from_polynomial(
            Arc::clone(context),
            &[Rational::parse(
                string_field(record, "numerator")?,
                string_field(record, "denominator")?,
            )?],
        ),
        "root_of_unity" => {
            let order = size_field(record, "order")?;
            let power = required(record, "power")?.as_i64().ok_or_else(|| {
                crate::ExactError::new(
                    "MalformedSerialization",
                    "root_of_unity power must be a signed integer",
                )
            })?;
            Element::root_of_unity(context, order, power)
        }
        "call" => encoded_call_from_json(context, record),
        kind => fail(
            "UnsupportedExactExpression",
            format!("unsupported encoded exact-expression kind {kind}"),
        ),
    }
}

fn encoded_call_from_json(
    context: &Arc<CyclotomicContext>,
    record: &Map<String, Value>,
) -> ExactResult<Element> {
    let head = string_field(record, "head")?;
    let arguments = array(required(record, "arguments")?)?;
    match head {
        "System`Plus" => arguments
            .iter()
            .try_fold(Element::zero(context)?, |total, value| {
                total.add(&encoded_element_from_json(context, value)?)
            }),
        "System`Times" => arguments
            .iter()
            .try_fold(Element::one(context)?, |total, value| {
                total.multiply(&encoded_element_from_json(context, value)?)
            }),
        "System`Complex" if arguments.len() == 2 => {
            let real = encoded_element_from_json(context, &arguments[0])?;
            let imaginary = encoded_element_from_json(context, &arguments[1])?;
            let unit = Element::root_of_unity(context, 4, 1)?;
            real.add(&imaginary.multiply(&unit)?)
        }
        "System`Power" if arguments.len() == 2 => {
            encoded_power_from_json(context, &arguments[0], &arguments[1])
        }
        "System`Cos" if arguments.len() == 1 => {
            encoded_trigonometric_from_json(context, &arguments[0], TrigonometricFunction::Cos)
        }
        "System`Sin" if arguments.len() == 1 => {
            encoded_trigonometric_from_json(context, &arguments[0], TrigonometricFunction::Sin)
        }
        "System`Csc" if arguments.len() == 1 => {
            encoded_trigonometric_from_json(context, &arguments[0], TrigonometricFunction::Sin)?
                .inverse()
        }
        _ => fail(
            "UnsupportedExactExpression",
            format!("unsupported encoded exact-expression head {head}"),
        ),
    }
}

fn encoded_power_from_json(
    context: &Arc<CyclotomicContext>,
    base: &Value,
    exponent: &Value,
) -> ExactResult<Element> {
    if let Some(power) = encoded_integer(exponent)? {
        return encoded_element_from_json(context, base)?.power(power);
    }
    let exponent_record = object(exponent)?;
    if string_field(exponent_record, "kind")? != "rational"
        || string_field(exponent_record, "denominator")? != "2"
        || !matches!(string_field(exponent_record, "numerator")?, "1" | "-1")
    {
        return fail(
            "UnsupportedExactExpression",
            "encoded power exponent must be an integer or canonical +/-1/2",
        );
    }
    let root = exact_rational_square_root(&encoded_element_from_json(context, base)?)?;
    if string_field(exponent_record, "numerator")? == "1" {
        Ok(root)
    } else {
        root.inverse()
    }
}

#[derive(Clone, Copy)]
enum TrigonometricFunction {
    Cos,
    Sin,
}

fn encoded_trigonometric_from_json(
    context: &Arc<CyclotomicContext>,
    argument: &Value,
    function: TrigonometricFunction,
) -> ExactResult<Element> {
    let multiple = rational_pi_multiple(argument)?;
    let numerator = multiple.numerator_string().parse::<i64>().map_err(|_| {
        crate::ExactError::new(
            "UnsupportedExactExpression",
            "trigonometric Pi multiple numerator exceeds i64",
        )
    })?;
    let denominator = multiple
        .denominator_string()
        .parse::<usize>()
        .map_err(|_| {
            crate::ExactError::new(
                "UnsupportedExactExpression",
                "trigonometric Pi multiple denominator exceeds usize",
            )
        })?;
    let order = denominator.checked_mul(2).ok_or_else(|| {
        crate::ExactError::new(
            "UnsupportedExactExpression",
            "trigonometric root-of-unity order overflows usize",
        )
    })?;
    let positive = Element::root_of_unity(context, order, numerator)?;
    let negative = Element::root_of_unity(context, order, -numerator)?;
    let two = Element::from_polynomial(Arc::clone(context), &[Rational::from_i64(2)])?;
    match function {
        TrigonometricFunction::Cos => positive.add(&negative)?.divide(&two),
        TrigonometricFunction::Sin => {
            let imaginary = Element::root_of_unity(context, 4, 1)?;
            positive
                .subtract(&negative)?
                .divide(&two.multiply(&imaginary)?)
        }
    }
}

fn rational_pi_multiple(value: &Value) -> ExactResult<Rational> {
    let record = object(value)?;
    if record.get("kind").and_then(Value::as_str) == Some("symbol")
        && record.get("name").and_then(Value::as_str) == Some("System`Pi")
    {
        return Ok(Rational::one());
    }
    if record.get("kind").and_then(Value::as_str) != Some("call")
        || record.get("head").and_then(Value::as_str) != Some("System`Times")
    {
        return fail(
            "UnsupportedExactExpression",
            "exact trigonometric arguments must be rational multiples of Pi",
        );
    }
    let mut coefficient = Rational::one();
    let mut pi_count = 0;
    for factor in array(required(record, "arguments")?)? {
        let factor_record = object(factor)?;
        if factor_record.get("kind").and_then(Value::as_str) == Some("symbol")
            && factor_record.get("name").and_then(Value::as_str) == Some("System`Pi")
        {
            pi_count += 1;
        } else {
            coefficient = coefficient.multiply(&encoded_rational_literal(factor)?);
        }
    }
    if pi_count != 1 {
        return fail(
            "UnsupportedExactExpression",
            "exact trigonometric arguments must contain Pi exactly once",
        );
    }
    Ok(coefficient)
}

fn encoded_rational_literal(value: &Value) -> ExactResult<Rational> {
    if let Some(integer) = bare_json_integer(value)? {
        return Ok(integer);
    }
    let record = object(value)?;
    match string_field(record, "kind")? {
        "integer" => Rational::parse(string_field(record, "value")?, "1"),
        "rational" => Rational::parse(
            string_field(record, "numerator")?,
            string_field(record, "denominator")?,
        ),
        _ => fail(
            "UnsupportedExactExpression",
            "expected an exact rational literal",
        ),
    }
}

fn exact_rational_square_root(value: &Element) -> ExactResult<Element> {
    if value
        .coefficients()
        .iter()
        .skip(1)
        .any(|coefficient| !coefficient.is_zero())
    {
        return fail(
            "UnsupportedExactExpression",
            "fractional power is not a rational perfect square or a supported canonical radical",
        );
    }
    let rational = &value.coefficients()[0];
    let numerator = BigInt::from_str(&rational.numerator_string()).map_err(|_| {
        crate::ExactError::new("MalformedSerialization", "invalid rational numerator")
    })?;
    let denominator = BigInt::from_str(&rational.denominator_string()).map_err(|_| {
        crate::ExactError::new("MalformedSerialization", "invalid rational denominator")
    })?;
    if numerator.sign() == Sign::Minus {
        return fail(
            "UnsupportedExactExpression",
            "fractional power of a negative rational is not supported",
        );
    }

    let context = value.context();
    let root_sum = |order| -> ExactResult<Element> {
        Element::root_of_unity(context, order, 1)?.add(&Element::root_of_unity(context, order, -1)?)
    };
    for square_free_part in [1_i64, 2, 3, 5, 6] {
        let quotient = Rational::new(
            numerator.clone(),
            &denominator * BigInt::from(square_free_part),
        )?;
        let quotient_numerator = BigInt::from_str(&quotient.numerator_string()).map_err(|_| {
            crate::ExactError::new("MalformedSerialization", "invalid rational numerator")
        })?;
        let quotient_denominator =
            BigInt::from_str(&quotient.denominator_string()).map_err(|_| {
                crate::ExactError::new("MalformedSerialization", "invalid rational denominator")
            })?;
        let (Some(numerator_root), Some(denominator_root)) = (
            exact_nonnegative_integer_square_root(&quotient_numerator),
            exact_nonnegative_integer_square_root(&quotient_denominator),
        ) else {
            continue;
        };
        let coefficient = Rational::new(numerator_root, denominator_root)?;
        let radical = match square_free_part {
            1 => Element::one(context)?,
            2 => root_sum(8)?,
            3 => root_sum(12)?,
            5 => Element::one(context)?.add(&root_sum(5)?.scale(&Rational::from_i64(2))?)?,
            6 => root_sum(8)?.multiply(&root_sum(12)?)?,
            _ => unreachable!("the supported radical list is fixed"),
        };
        return radical.scale(&coefficient);
    }
    fail(
        "UnsupportedExactExpression",
        "rational square root is outside the supported canonical radicals 2, 3, 5, and 6",
    )
}

fn exact_nonnegative_integer_square_root(value: &BigInt) -> Option<BigInt> {
    if value.sign() == Sign::Minus {
        return None;
    }
    if value == &BigInt::from(0) {
        return Some(BigInt::from(0));
    }
    let mut estimate = value.clone();
    let two = BigInt::from(2);
    loop {
        let next = (&estimate + value / &estimate) / &two;
        if next >= estimate {
            return (&estimate * &estimate == *value).then_some(estimate);
        }
        estimate = next;
    }
}

fn encoded_integer(value: &Value) -> ExactResult<Option<i64>> {
    if let Some(integer) = bare_json_integer(value)? {
        return integer
            .numerator_string()
            .parse::<i64>()
            .map(Some)
            .map_err(|_| {
                crate::ExactError::new(
                    "MalformedSerialization",
                    "encoded integer is outside i64 range",
                )
            });
    }
    let record = object(value)?;
    if string_field(record, "kind")? != "integer" {
        return Ok(None);
    }
    string_field(record, "value")?
        .parse::<i64>()
        .map(Some)
        .map_err(|_| {
            crate::ExactError::new(
                "MalformedSerialization",
                "encoded integer is outside i64 range",
            )
        })
}

/// Parses only a canonical bare JSON integer token and never converts through
/// `f64`.  The `serde_json/arbitrary_precision` feature preserves the original
/// number token, so decimal and exponent spellings remain distinguishable and
/// are rejected even when they denote a mathematical integer.
fn bare_json_integer(value: &Value) -> ExactResult<Option<Rational>> {
    if !value.is_number() {
        return Ok(None);
    }
    let token = value.to_string();
    Rational::parse(&token, "1").map(Some).map_err(|_| {
        crate::ExactError::new(
            "MalformedSerialization",
            format!("exact bare JSON number must use integer syntax, received {token}"),
        )
    })
}

pub fn encoded_matrix_from_json(
    context: &Arc<CyclotomicContext>,
    value: &Value,
) -> ExactResult<ExactMatrix> {
    let record = object(value)?;
    if !record.contains_key("kind") {
        return matrix_from_json(context, value);
    }
    if string_field(record, "kind")? != "matrix" {
        return fail("MalformedSerialization", "expected encoded exact matrix");
    }
    let rows = size_field(record, "rows")?;
    let columns = size_field(record, "columns")?;
    let entries = array(required(record, "entries")?)?;
    if entries.len()
        != rows.checked_mul(columns).ok_or_else(|| {
            crate::ExactError::new("MalformedSerialization", "matrix shape overflows")
        })?
    {
        return fail(
            "MalformedSerialization",
            "encoded matrix entry count does not match shape",
        );
    }
    ExactMatrix::new(
        Arc::clone(context),
        rows,
        columns,
        entries
            .iter()
            .map(|entry| encoded_element_from_json(context, entry))
            .collect::<ExactResult<Vec<_>>>()?,
    )
}

pub fn matrix_from_json(
    context: &Arc<CyclotomicContext>,
    value: &Value,
) -> ExactResult<ExactMatrix> {
    let object = object(value)?;
    let rows = size_field(object, "rows")?;
    let columns = size_field(object, "columns")?;
    let entries_value = array(required(object, "entries")?)?;
    let expected = rows.checked_mul(columns).ok_or_else(|| {
        crate::ExactError::new("MalformedSerialization", "matrix shape overflows")
    })?;
    if entries_value.len() != expected {
        return fail(
            "MalformedSerialization",
            "fixture matrix entry count does not match shape",
        );
    }
    let entries = entries_value
        .iter()
        .map(|entry| element_from_json(context, entry))
        .collect::<ExactResult<Vec<_>>>()?;
    ExactMatrix::new(Arc::clone(context), rows, columns, entries)
}

pub fn problem_from_json(value: &Value) -> ExactResult<CompiledProblem> {
    let object = object(value)?;
    let context = context_from_json(required(object, "context")?, 128)?;
    let coordinate_dimension = size_field(object, "coordinate_dimension")?;
    let raw_constraints = array(required(object, "constraints")?)?;
    let mut constraints = Vec::with_capacity(raw_constraints.len());
    for raw_constraint in raw_constraints {
        let matrix = matrix_from_json(&context, raw_constraint)?;
        if matrix.columns() != coordinate_dimension {
            return fail(
                "DimensionMismatch",
                "constraint columns do not match coordinate dimension",
            );
        }
        constraints.push(matrix);
    }
    Ok(CompiledProblem {
        context,
        coordinate_dimension,
        constraints,
    })
}

pub fn sizes_from_json(value: &Value) -> ExactResult<Vec<usize>> {
    array(value)?
        .iter()
        .map(|item| size_value(item, "array item"))
        .collect()
}

#[must_use]
pub fn rational_to_json(value: &Rational) -> Value {
    json!({
        "numerator": value.numerator_string(),
        "denominator": value.denominator_string()
    })
}

#[must_use]
pub fn polynomial_to_json(polynomial: &[Rational]) -> Value {
    Value::Array(polynomial.iter().map(rational_to_json).collect())
}

#[must_use]
pub fn context_to_json(context: &CyclotomicContext) -> Value {
    json!({
        "conductor": context.conductor,
        "degree": context.degree,
        "cyclotomic_polynomial": polynomial_to_json(&context.cyclotomic_polynomial)
    })
}

#[must_use]
pub fn element_to_json(element: &Element) -> Value {
    json!({"coefficients": polynomial_to_json(element.coefficients())})
}

#[must_use]
pub fn matrix_to_json(matrix: &ExactMatrix) -> Value {
    json!({
        "rows": matrix.rows(),
        "columns": matrix.columns(),
        "entries": matrix.entries().iter().map(element_to_json).collect::<Vec<_>>()
    })
}

#[must_use]
pub fn kernel_result_to_json(result: &KernelResult) -> Value {
    json!({
        "reduced_matrix": matrix_to_json(&result.reduced_matrix),
        "nullspace_rows": matrix_to_json(&result.nullspace_rows),
        "basis_matrix": matrix_to_json(&result.basis_matrix),
        "pivot_columns": result.pivot_columns,
        "free_columns": result.free_columns,
        "rank": result.rank,
        "nullity": result.nullity,
        "exact_residual_verified": result.exact_residual_verified
    })
}

#[must_use]
pub fn common_result_to_json(result: &CommonKernelResult) -> Value {
    json!({
        "iteration_nullities": result.iteration_nullities,
        "basis_matrix": matrix_to_json(&result.basis_matrix),
        "nullspace_rows": matrix_to_json(&result.nullspace_rows),
        "rank": result.rank,
        "nullity": result.nullity,
        "exact_residual_verified": result.exact_residual_verified
    })
}
