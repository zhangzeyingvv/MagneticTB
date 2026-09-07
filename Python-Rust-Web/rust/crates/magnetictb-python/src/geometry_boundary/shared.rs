use std::collections::BTreeMap;
use std::sync::Arc;

use cyclotomic_nullspace::fixture::{
    array, context_from_json, context_to_json, element_from_json, element_to_json,
    encoded_element_from_json, encoded_matrix_from_json, matrix_from_json, matrix_to_json, object,
    size_field, string_field,
};
use cyclotomic_nullspace::{
    CyclotomicContext, Element, ExactError, ExactMatrix, ExactResult, Rational,
};
use magnetictb_abstract_group::GroupAlgebra;
use magnetictb_crystal_geometry::{
    BondSearchOptions, PeriodicBondRecord, PeriodicBondSearchResult, PeriodicBondShell,
    QuadraticReal, SitePermutationData, compile_msg_wyckoff_sites,
    compile_msg_wyckoff_sites_by_letter, compile_site_permutations, find_periodic_bond_shells,
    periodic_bond_shells_in_box,
};
use magnetictb_data::{DataCatalog, ExactExpression};
use magnetictb_linear_algebra::{
    ConstraintKernelMethod, ConstraintValidationLevel, add, exact_hermitian_matrix,
    exact_matrix_equal,
};
use magnetictb_properties::{first_brillouin_zone, fold_band_path_to_first_bz, standard_k_path};
use magnetictb_representation::{
    InducedOrbitSpec, automatic_coset_representatives, compile_block_monomial,
    compile_direct_product, compile_induced,
};
use magnetictb_symmetry::{
    CompiledBasisAction, CompiledInducedBasisAction, CompiledSeitzGroup, FunctionBasisState,
    PolynomialExpression, SeitzOperation, SpatialOperation, SpinSpaceOperation,
    compile_catalog_basis_action, compile_catalog_induced_basis_action,
    compile_function_basis_action, compile_function_induced_basis_action,
    compile_msg_group_in_context, compile_seitz_group, compile_seitz_group_with_table,
    compile_spatial_spin_action, compile_spatial_spinor_matrix, compile_spin_space_group,
    evaluate_data_expression, fractional_part, generate_seitz_group, generate_spin_space_group,
    resolve_catalog_basis_states, resolve_function_basis_states,
};
use magnetictb_tight_binding::{
    BandTraceInput, BondConstraintData, DirectedBondOrbitData, NumericSymmetryOperation,
    compile_bond_constraints, compile_bond_constraints_with_continuous,
    compile_directed_bond_orbits, compile_hermitian_bond_constraints, compute_band_trace,
    fourier_coefficients, reconstruct_solved_bond_model, verify_hamiltonian_symmetry,
};
use nalgebra::{Complex, DMatrix, Matrix3, Vector3};
use serde::de::DeserializeOwned;
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::fmt::Write;

fn malformed(detail: impl Into<String>) -> ExactError {
    ExactError::new("MalformedSerialization", detail)
}
fn stable_machine_real(decimal: &str) -> ExactResult<Value> {
    let value = decimal.parse::<f64>().map_err(|_| {
        ExactError::new(
            "InvalidInexactInput",
            format!("cannot parse finite machine real {decimal}"),
        )
    })?;
    if !value.is_finite() {
        return Err(ExactError::new(
            "InvalidInexactInput",
            "machine-real input must be finite",
        ));
    }
    if value == 0.0 {
        return Ok(json!({"kind": "integer", "value": "0"}));
    }

    // Preserve the complete decimal spelling emitted by Python instead of
    // snapping every machine real to a 0.001 grid. Mathematica's NumericQ
    // lattice boundary distinguishes these inputs; the previous rounding
    // silently changed both the metric and model identity. The reduced
    // rational is only an internal carrier through the exact compiler. Bond
    // topology is still selected by the dedicated f64 geometry path.
    let (mantissa, exponent_text) = decimal
        .split_once(['e', 'E'])
        .map_or((decimal, "0"), |parts| parts);
    let exponent = exponent_text
        .parse::<i32>()
        .map_err(|_| ExactError::new("InvalidInexactInput", "machine-real exponent is invalid"))?;
    if !(-4096..=4096).contains(&exponent) {
        return Err(ExactError::new(
            "InvalidInexactInput",
            "machine-real exponent is outside the finite binary64 range",
        ));
    }
    let (negative, unsigned) = mantissa
        .strip_prefix('-')
        .map_or((false, mantissa), |value| (true, value));
    let unsigned = unsigned.strip_prefix('+').unwrap_or(unsigned);
    let (integer, fraction) = unsigned
        .split_once('.')
        .map_or((unsigned, ""), |parts| parts);
    if integer.is_empty()
        || !integer.bytes().all(|digit| digit.is_ascii_digit())
        || !fraction.bytes().all(|digit| digit.is_ascii_digit())
    {
        return Err(ExactError::new(
            "InvalidInexactInput",
            "machine-real mantissa is invalid",
        ));
    }
    let mut digits = format!("{integer}{fraction}");
    let first_nonzero = digits
        .find(|digit: char| digit != '0')
        .unwrap_or(digits.len());
    digits.drain(..first_nonzero);
    if digits.is_empty() {
        return Ok(json!({"kind": "integer", "value": "0"}));
    }
    let decimal_places = i32::try_from(fraction.len())
        .map_err(|_| malformed("machine-real mantissa is too long"))?;
    let scale = decimal_places - exponent;
    let (numerator, denominator) = if scale <= 0 {
        let zeros =
            usize::try_from(-scale).map_err(|_| malformed("machine-real exponent is too large"))?;
        digits.push_str(&"0".repeat(zeros));
        (digits, "1".to_owned())
    } else {
        let zeros =
            usize::try_from(scale).map_err(|_| malformed("machine-real exponent is too small"))?;
        (digits, format!("1{}", "0".repeat(zeros)))
    };
    let signed_numerator = if negative {
        format!("-{numerator}")
    } else {
        numerator
    };
    let rational = Rational::parse(&signed_numerator, &denominator)?;
    if rational.denominator_string() == "1" {
        Ok(json!({
            "kind": "integer",
            "value": rational.numerator_string()
        }))
    } else {
        Ok(json!({
            "kind": "rational",
            "numerator": rational.numerator_string(),
            "denominator": rational.denominator_string()
        }))
    }
}

fn normalize_stable_inexact(value: &mut Value) -> ExactResult<()> {
    match value {
        Value::Array(items) => {
            for item in items {
                normalize_stable_inexact(item)?;
            }
        }
        Value::Object(record) => {
            match record.get("kind").and_then(Value::as_str) {
                Some("machine_real") => {
                    let decimal = record
                        .get("decimal")
                        .and_then(Value::as_str)
                        .ok_or_else(|| malformed("machine_real requires decimal text"))?;
                    *value = stable_machine_real(decimal)?;
                    return Ok(());
                }
                Some("machine_complex") => {
                    let real = record
                        .get("real_decimal")
                        .and_then(Value::as_str)
                        .ok_or_else(|| malformed("machine_complex requires real_decimal"))?;
                    let imaginary = record
                        .get("imaginary_decimal")
                        .and_then(Value::as_str)
                        .ok_or_else(|| malformed("machine_complex requires imaginary_decimal"))?;
                    let real = stable_machine_real(real)?;
                    let imaginary = stable_machine_real(imaginary)?;
                    *value = if imaginary == json!({"kind": "integer", "value": "0"}) {
                        real
                    } else {
                        json!({
                            "kind": "call",
                            "head": "System`Complex",
                            "arguments": [real, imaginary]
                        })
                    };
                    return Ok(());
                }
                _ => {}
            }
            for item in record.values_mut() {
                normalize_stable_inexact(item)?;
            }
        }
        _ => {}
    }
    Ok(())
}

fn from_value<T: DeserializeOwned>(value: &Value, field: &str) -> ExactResult<T> {
    serde_json::from_value(value.clone())
        .map_err(|error| malformed(format!("invalid {field}: {error}")))
}

fn required<'a>(record: &'a serde_json::Map<String, Value>, field: &str) -> ExactResult<&'a Value> {
    record
        .get(field)
        .ok_or_else(|| malformed(format!("missing field {field}")))
}

fn tagged_items(value: &Value) -> ExactResult<&Vec<Value>> {
    if let Some(items) = value.as_array() {
        return Ok(items);
    }
    let record = object(value)?;
    if record.get("kind").and_then(Value::as_str) != Some("list") {
        return Err(malformed("expected a tagged list"));
    }
    array(required(record, "items")?)
}

fn ordered_lookup<'a>(value: &'a Value, key: &str) -> ExactResult<&'a Value> {
    let record = object(value)?;
    if record.get("kind").and_then(Value::as_str) != Some("ordered_association") {
        return Err(malformed("expected an ordered association"));
    }
    for entry in array(required(record, "entries")?)? {
        let entry = object(entry)?;
        if required(entry, "key")?.as_str() == Some(key) {
            return required(entry, "value");
        }
    }
    Err(malformed(format!("missing ordered-association key {key}")))
}

fn optional_ordered_lookup<'a>(value: &'a Value, key: &str) -> Option<&'a Value> {
    ordered_lookup(value, key).ok()
}

fn decode_tagged_lists(value: &Value) -> ExactResult<Value> {
    if let Some(record) = value.as_object()
        && record.get("kind").and_then(Value::as_str) == Some("list")
    {
        return Ok(Value::Array(
            array(required(record, "items")?)?
                .iter()
                .map(decode_tagged_lists)
                .collect::<ExactResult<Vec<_>>>()?,
        ));
    }
    match value {
        Value::Array(items) => Ok(Value::Array(
            items
                .iter()
                .map(decode_tagged_lists)
                .collect::<ExactResult<Vec<_>>>()?,
        )),
        _ => Ok(value.clone()),
    }
}

fn lattice_parameter_bindings(compiler_input: &Value) -> ExactResult<BTreeMap<String, Value>> {
    let mut result = BTreeMap::new();
    for rule in tagged_items(ordered_lookup(compiler_input, "LatticeParameters")?)? {
        let record = object(rule)?;
        if string_field(record, "kind")? != "call" || string_field(record, "head")? != "System`Rule"
        {
            return Err(malformed(
                "LatticeParameters must contain exact Rule expressions",
            ));
        }
        let arguments = array(required(record, "arguments")?)?;
        if arguments.len() != 2 {
            return Err(malformed(
                "a lattice-parameter Rule must have two arguments",
            ));
        }
        let symbol = object(&arguments[0])?;
        if string_field(symbol, "kind")? != "symbol" {
            return Err(malformed("a lattice-parameter Rule must bind a symbol"));
        }
        let name = string_field(symbol, "name")?.to_owned();
        if result.insert(name, arguments[1].clone()).is_some() {
            return Err(malformed("duplicate lattice-parameter binding"));
        }
    }
    Ok(result)
}

fn substitute_exact_symbols(
    value: &Value,
    bindings: &BTreeMap<String, Value>,
) -> ExactResult<Value> {
    if let Some(record) = value.as_object()
        && record.get("kind").and_then(Value::as_str) == Some("symbol")
    {
        let name = string_field(record, "name")?;
        if !bindings.contains_key(name) && matches!(name, "System`Pi" | "System`E") {
            // Constants in a literal lattice need no user parameter rule.
            // Their supported exact/numerical treatment belongs to the
            // lattice resolver, just like constants inside a bound value.
            return Ok(value.clone());
        }
        return bindings.get(name).cloned().ok_or_else(|| {
            ExactError::new(
                "MissingExactSymbol",
                format!("no exact lattice-parameter binding was supplied for {name}"),
            )
        });
    }
    match value {
        Value::Array(items) => Ok(Value::Array(
            items
                .iter()
                .map(|item| substitute_exact_symbols(item, bindings))
                .collect::<ExactResult<Vec<_>>>()?,
        )),
        Value::Object(record) => Ok(Value::Object(
            record
                .iter()
                .map(|(key, item)| Ok((key.clone(), substitute_exact_symbols(item, bindings)?)))
                .collect::<ExactResult<serde_json::Map<_, _>>>()?,
        )),
        _ => Ok(value.clone()),
    }
}

fn quadratic_to_json(value: &QuadraticReal) -> Value {
    json!({
        "basis": ["1", "sqrt(2)", "sqrt(3)", "sqrt(6)"],
        "coefficients": value.coefficients().iter().map(|coefficient| json!({
            "numerator": coefficient.numerator_string(),
            "denominator": coefficient.denominator_string()
        })).collect::<Vec<_>>()
    })
}

fn quadratic_from_json(value: &Value) -> ExactResult<QuadraticReal> {
    let record = object(value)?;
    let basis = array(required(record, "basis")?)?;
    if basis.as_slice()
        != [
            Value::String("1".to_owned()),
            Value::String("sqrt(2)".to_owned()),
            Value::String("sqrt(3)".to_owned()),
            Value::String("sqrt(6)".to_owned()),
        ]
    {
        return Err(malformed("unsupported quadratic-real basis"));
    }
    let coefficients = array(required(record, "coefficients")?)?;
    if coefficients.len() != 4 {
        return Err(malformed(
            "quadratic-real value must contain four coefficients",
        ));
    }
    let coefficients = coefficients
        .iter()
        .map(|coefficient| {
            let coefficient = object(coefficient)?;
            Rational::parse(
                string_field(coefficient, "numerator")?,
                string_field(coefficient, "denominator")?,
            )
        })
        .collect::<ExactResult<Vec<_>>>()?;
    let basis = [
        QuadraticReal::one(),
        QuadraticReal::sqrt_two(),
        QuadraticReal::sqrt_three(),
        QuadraticReal::sqrt_two().multiply(&QuadraticReal::sqrt_three()),
    ];
    Ok(basis
        .iter()
        .zip(coefficients)
        .fold(QuadraticReal::zero(), |sum, (basis, coefficient)| {
            sum.add(&basis.scale(&coefficient))
        }))
}

fn cartesian_displacement(displacement: &Value, lattice: &Value) -> ExactResult<Value> {
    let displacement = array(displacement)?
        .iter()
        .map(quadratic_from_json)
        .collect::<ExactResult<Vec<_>>>()?;
    let lattice_record = object(lattice)?;
    if lattice_record.get("kind").and_then(Value::as_str) != Some("matrix") {
        return Err(malformed("model lattice must be an exact matrix"));
    }
    let rows = size_field(lattice_record, "rows")?;
    let columns = size_field(lattice_record, "columns")?;
    let entries: Vec<ExactExpression> = from_value(
        required(lattice_record, "entries")?,
        "model lattice entries",
    )?;
    if rows.checked_mul(columns) != Some(entries.len()) {
        return Err(malformed(
            "model lattice dimensions do not match its entries",
        ));
    }
    let bindings = BTreeMap::new();
    let lattice = entries
        .chunks(columns)
        .map(|row| {
            row.iter()
                .map(|entry| {
                    QuadraticReal::from_expression(entry, &bindings)
                        .map_err(|error| ExactError::new(error.tag(), error.to_string()))
                })
                .collect::<ExactResult<Vec<_>>>()
        })
        .collect::<ExactResult<Vec<_>>>()?;
    if lattice.len() != displacement.len()
        || lattice.is_empty()
        || lattice.iter().any(|row| row.len() != lattice[0].len())
    {
        return Err(ExactError::new(
            "UnsupportedCoordinateMode",
            "model lattice and fractional displacement dimensions do not align",
        ));
    }
    let mut cartesian = vec![QuadraticReal::zero(); lattice[0].len()];
    for (coordinate, row) in displacement.iter().zip(&lattice) {
        for (result, lattice_value) in cartesian.iter_mut().zip(row) {
            *result = result.add(&coordinate.multiply(lattice_value));
        }
    }
    Ok(vector_to_json(&cartesian))
}

fn exact_expression_value(expression: &ExactExpression) -> ExactResult<Value> {
    serde_json::to_value(expression)
        .map_err(|error| malformed(format!("cannot encode exact expression: {error}")))
}

fn quadratic_expression(value: &QuadraticReal) -> Value {
    let radical = |radicand: &str| {
        json!({
            "kind": "call",
            "head": "System`Power",
            "arguments": [
                {"kind": "integer", "value": radicand},
                {"kind": "rational", "numerator": "1", "denominator": "2"}
            ]
        })
    };
    let mut terms = Vec::new();
    for (index, coefficient) in value.coefficients().iter().enumerate() {
        if coefficient.is_zero() {
            continue;
        }
        let encoded = encoded_rational_value(coefficient);
        if index == 0 {
            terms.push(encoded);
            continue;
        }
        let mut factors = Vec::new();
        if coefficient.numerator_string() != "1" || coefficient.denominator_string() != "1" {
            factors.push(encoded);
        }
        if index & 1 != 0 {
            factors.push(radical("2"));
        }
        if index & 2 != 0 {
            factors.push(radical("3"));
        }
        terms.push(if factors.len() == 1 {
            factors.remove(0)
        } else {
            json!({"kind": "call", "head": "System`Times", "arguments": factors})
        });
    }
    match terms.len() {
        0 => json!({"kind": "integer", "value": "0"}),
        1 => terms.remove(0),
        _ => json!({"kind": "call", "head": "System`Plus", "arguments": terms}),
    }
}
