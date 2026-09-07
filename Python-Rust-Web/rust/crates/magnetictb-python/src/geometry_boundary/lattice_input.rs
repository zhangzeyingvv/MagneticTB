// Owner-approved lattice policy: retain supported exact entries; otherwise
// evaluate the complete substituted entry numerically and carry its seven
// significant decimal digits as a rational. This is not an exact-algebra or
// representation fallback. All subsequent symmetry checks remain authoritative.
struct ResolvedLattice {
    matrix: ExactMatrix,
    encoded: Value,
    approximations: Vec<Value>,
}

fn lattice_rational_literal(value: &Value) -> ExactResult<Option<Rational>> {
    if value.is_number() {
        return Rational::parse(&value.to_string(), "1").map(Some);
    }
    let record = object(value)?;
    match record.get("kind").and_then(Value::as_str) {
        Some("integer") => Rational::parse(string_field(record, "value")?, "1").map(Some),
        Some("rational") => Rational::parse(
            string_field(record, "numerator")?,
            string_field(record, "denominator")?,
        )
        .map(Some),
        _ => Ok(None),
    }
}

// Python's AST retains nested multiplication. Canonicalize rational factors
// before the exact attempt so 7*Pi/18 and (7/18)*Pi take the same path, and
// equivalent spellings of crystallographic angles keep their exact radicals.
fn normalize_lattice_expression(value: &Value) -> ExactResult<Value> {
    let Some(record) = value.as_object() else {
        return Ok(value.clone());
    };
    if record.get("kind").and_then(Value::as_str) != Some("call") {
        return Ok(value.clone());
    }
    let head = string_field(record, "head")?;
    let arguments = array(required(record, "arguments")?)?
        .iter()
        .map(normalize_lattice_expression)
        .collect::<ExactResult<Vec<_>>>()?;
    if head != "System`Times" {
        return Ok(json!({"kind": "call", "head": head, "arguments": arguments}));
    }
    let mut factors = Vec::new();
    for argument in arguments {
        if argument.get("head").and_then(Value::as_str) == Some("System`Times") {
            factors.extend(
                array(required(object(&argument)?, "arguments")?)?
                    .iter()
                    .cloned(),
            );
        } else {
            factors.push(argument);
        }
    }
    let mut coefficient = Rational::one();
    let mut remaining = Vec::new();
    for factor in factors {
        if let Some(rational) = lattice_rational_literal(&factor)? {
            coefficient = coefficient.multiply(&rational);
        } else {
            remaining.push(factor);
        }
    }
    if remaining.is_empty() {
        return Ok(encoded_rational_value(&coefficient));
    }
    if coefficient != Rational::one() {
        remaining.insert(0, encoded_rational_value(&coefficient));
    }
    if remaining.len() == 1 {
        return Ok(remaining.remove(0));
    }
    Ok(json!({"kind": "call", "head": head, "arguments": remaining}))
}

fn finite_lattice_number(value: Complex<f64>) -> ExactResult<Complex<f64>> {
    if value.re.is_finite() && value.im.is_finite() {
        Ok(value)
    } else {
        Err(ExactError::new(
            "InvalidLattice",
            "lattice expression must evaluate to a finite number",
        ))
    }
}

// Numerical evaluation is confined to constant lattice ASTs after parameter
// substitution. Rust/num-complex supplies the numerical elementary functions.
fn numeric_lattice_expression(value: &Value) -> ExactResult<Complex<f64>> {
    if let Some(rational) = lattice_rational_literal(value)? {
        return rational
            .to_f64()
            .map(|real| Complex::new(real, 0.0))
            .ok_or_else(|| {
                ExactError::new(
                    "InvalidLattice",
                    "lattice scalar is outside the finite f64 range",
                )
            });
    }
    let record = object(value)?;
    let result = match string_field(record, "kind")? {
        "symbol" => match string_field(record, "name")? {
            "System`Pi" => Complex::new(std::f64::consts::PI, 0.0),
            "System`E" => Complex::new(std::f64::consts::E, 0.0),
            name => {
                return Err(ExactError::new(
                    "InvalidLattice",
                    format!("unresolved lattice symbol {name}"),
                ));
            }
        },
        "root_of_unity" => {
            let order = size_field(record, "order")?;
            let order =
                i64::try_from(order).map_err(|_| malformed("lattice root order exceeds i64"))?;
            if order <= 0 {
                return Err(ExactError::new(
                    "InvalidRootOfUnity",
                    "root order must be positive",
                ));
            }
            let power = required(record, "power")?
                .as_i64()
                .ok_or_else(|| malformed("root power must be a signed integer"))?;
            let mut power = power.rem_euclid(order);
            if power > order / 2 {
                power -= order;
            }
            let fraction = Rational::parse(&power.to_string(), &order.to_string())?
                .to_f64()
                .ok_or_else(|| malformed("invalid root phase"))?;
            Complex::from_polar(1.0, std::f64::consts::TAU * fraction)
        }
        "call" => {
            let head = string_field(record, "head")?;
            let arguments = array(required(record, "arguments")?)?
                .iter()
                .map(numeric_lattice_expression)
                .collect::<ExactResult<Vec<_>>>()?;
            match (head, arguments.as_slice()) {
                ("System`Plus", args) => args.iter().copied().sum(),
                ("System`Times", args) => args.iter().copied().product(),
                ("System`Sin", [argument]) => argument.sin(),
                ("System`Cos", [argument]) => argument.cos(),
                ("System`Csc", [argument]) => Complex::new(1.0, 0.0) / argument.sin(),
                ("System`Exp", [argument]) => argument.exp(),
                ("System`Power", [base, exponent])
                    if base.im == 0.0
                        && exponent.im == 0.0
                        && (base.re >= 0.0 || exponent.re.fract() == 0.0) =>
                {
                    Complex::new(base.re.powf(exponent.re), 0.0)
                }
                ("System`Power", [base, exponent]) => base.powc(*exponent),
                ("System`Complex", [real, imaginary]) if real.im == 0.0 && imaginary.im == 0.0 => {
                    Complex::new(real.re, imaginary.re)
                }
                _ => {
                    return Err(ExactError::new(
                        "UnsupportedExactExpression",
                        format!("unsupported numeric lattice call {head}"),
                    ));
                }
            }
        }
        kind => {
            return Err(ExactError::new(
                "UnsupportedExactExpression",
                format!("unsupported numeric lattice kind {kind}"),
            ));
        }
    };
    finite_lattice_number(result)
}

fn rounded_lattice_rational(value: Complex<f64>) -> ExactResult<Value> {
    let value = finite_lattice_number(value)?;
    if value.im != 0.0 {
        return Err(ExactError::new(
            "InvalidLattice",
            "lattice expression must be real",
        ));
    }
    // Scientific notation uses one leading digit plus six fractional digits.
    // Decimal-to-Q is exact after this single, explicit rounding step.
    stable_machine_real(&format!("{:.6e}", value.re))
}

fn resolve_model_lattice(
    context: &Arc<CyclotomicContext>,
    value: &Value,
) -> ExactResult<ResolvedLattice> {
    let record = object(value)?;
    if record.get("kind").and_then(Value::as_str) != Some("matrix") {
        return Err(malformed("model lattice requires an encoded matrix"));
    }
    if size_field(record, "rows")? != 3 || size_field(record, "columns")? != 3 {
        return Err(ExactError::new(
            "InvalidLattice",
            "the evaluated lattice must be a 3 by 3 matrix",
        ));
    }
    let entries = array(required(record, "entries")?)?;
    if entries.len() != 9 {
        return Err(malformed("lattice entry count does not match its shape"));
    }
    let mut resolved = Vec::with_capacity(9);
    let mut approximations = Vec::new();
    for (index, source) in entries.iter().enumerate() {
        let expression = normalize_lattice_expression(source)?;
        let exact = encoded_element_from_json(context, &expression).and_then(|element| {
            if element.conjugate()? != element {
                return Err(ExactError::new(
                    "InvalidLattice",
                    "lattice expression must be real",
                ));
            }
            rational_quadratic_coordinates(&element)?;
            Ok(element)
        });
        match exact {
            Ok(_) => resolved.push(expression),
            Err(error)
                if matches!(
                    error.tag(),
                    "UnsupportedExactExpression"
                        | "ConductorMismatch"
                        | "UnsupportedBondExactExpression"
                ) =>
            {
                let rational = rounded_lattice_rational(numeric_lattice_expression(&expression)?)?;
                approximations.push(json!({
                    "row_index": index / 3,
                    "column_index": index % 3,
                    "source_expression": source,
                    "rational_value": &rational
                }));
                resolved.push(rational);
            }
            Err(error) => return Err(error),
        }
    }
    let encoded = json!({"kind": "matrix", "rows": 3, "columns": 3, "entries": resolved});
    Ok(ResolvedLattice {
        matrix: encoded_matrix_from_json(context, &encoded)?,
        encoded,
        approximations,
    })
}

fn record_lattice_approximation(result: &mut Value, entries: &[Value]) {
    if !entries.is_empty() {
        result["lattice_approximation"] = json!({
            "method": "numeric_then_rational",
            "significant_digits": 7,
            "entries": entries
        });
    }
}

#[cfg(test)]
mod lattice_input_tests {
    use super::*;

    fn number(value: i64) -> Value {
        json!({"kind": "integer", "value": value.to_string()})
    }

    fn rational(numerator: i64, denominator: i64) -> Value {
        encoded_rational_value(
            &Rational::parse(&numerator.to_string(), &denominator.to_string()).unwrap(),
        )
    }

    fn call(head: &str, arguments: Vec<Value>) -> Value {
        let mut value = json!({"kind": "call", "head": format!("System`{head}")});
        value["arguments"] = Value::Array(arguments);
        value
    }

    fn lattice(last_row: [Value; 3]) -> Value {
        let mut value = json!({"kind": "matrix", "rows": 3, "columns": 3});
        let mut entries = vec![
            number(1),
            number(0),
            number(0),
            number(0),
            number(2),
            number(0),
        ];
        entries.extend(last_row);
        value["entries"] = Value::Array(entries);
        value
    }

    #[test]
    fn seven_digits_are_significant_digits_across_scales() {
        for (input, expected) in [
            (1.234_567_89, rational(1_234_568, 1_000_000)),
            (-0.342_020_143_325, rational(-3_420_201, 10_000_000)),
            (0.000_001_234_567_89, rational(1_234_568, 1_000_000_000_000)),
            (12_345_678.9, number(12_345_680)),
            (9.999_999_9, number(10)),
            (0.0, number(0)),
        ] {
            assert_eq!(
                rounded_lattice_rational(Complex::new(input, 0.0)).unwrap(),
                expected
            );
        }
        for input in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(
                rounded_lattice_rational(Complex::new(input, 0.0))
                    .unwrap_err()
                    .tag(),
                "InvalidLattice"
            );
        }
    }

    #[test]
    fn supported_hexagonal_sqrt_three_and_nested_pi_stay_exact() {
        let context = cyclotomic_nullspace::make_context(24, 128).unwrap();
        let sqrt_three = call("Power", vec![number(3), rational(1, 2)]);
        let resolved = resolve_model_lattice(
            &context,
            &lattice([number(1), number(0), sqrt_three.clone()]),
        )
        .unwrap();
        assert!(resolved.approximations.is_empty());
        assert_eq!(
            resolved.matrix.entry(2, 2).unwrap(),
            &encoded_element_from_json(&context, &sqrt_three).unwrap()
        );

        let pi = json!({"kind": "symbol", "name": "System`Pi"});
        let angle = call(
            "Times",
            vec![call("Times", vec![number(2), pi]), rational(1, 6)],
        );
        let resolved = resolve_model_lattice(
            &context,
            &lattice([
                call("Cos", vec![angle.clone()]),
                number(0),
                call("Sin", vec![angle]),
            ]),
        )
        .unwrap();
        assert!(resolved.approximations.is_empty());
        assert_eq!(
            resolved.matrix.entry(2, 0).unwrap(),
            &encoded_element_from_json(&context, &rational(1, 2)).unwrap()
        );
        let expected = call("Times", vec![rational(1, 2), sqrt_three]);
        assert_eq!(
            resolved.matrix.entry(2, 2).unwrap(),
            &encoded_element_from_json(&context, &expected).unwrap()
        );
    }

    #[test]
    fn general_angle_is_evaluated_before_rounding_the_complete_entry() {
        let context = cyclotomic_nullspace::make_context(24, 128).unwrap();
        let angle = call(
            "Times",
            vec![
                call(
                    "Times",
                    vec![number(7), json!({"kind": "symbol", "name": "System`Pi"})],
                ),
                rational(1, 18),
            ],
        );
        let input = lattice([
            call("Times", vec![number(3), call("Cos", vec![angle.clone()])]),
            number(0),
            call("Times", vec![number(3), call("Sin", vec![angle])]),
        ]);
        let resolved = resolve_model_lattice(&context, &input).unwrap();
        assert_eq!(resolved.approximations.len(), 2);
        assert_eq!(resolved.encoded["entries"][6], rational(51_303, 50_000));
        assert_eq!(resolved.encoded["entries"][8], rational(1_409_539, 500_000));
        assert_eq!(
            resolved.approximations[0]["source_expression"],
            input["entries"][6]
        );
        assert_eq!(resolved.approximations[0]["row_index"], 2);
        assert_eq!(resolved.approximations[0]["column_index"], 0);
    }

    #[test]
    fn lattice_conversion_does_not_relax_exact_matrix_decoding() {
        let context = cyclotomic_nullspace::make_context(24, 128).unwrap();
        let sqrt_seven = call("Power", vec![number(7), rational(1, 2)]);
        let input = lattice([number(0), number(0), sqrt_seven.clone()]);
        let resolved = resolve_model_lattice(&context, &input).unwrap();
        assert_eq!(
            resolved.encoded["entries"][8],
            rational(2_645_751, 1_000_000)
        );
        assert!(encoded_matrix_from_json(&context, &input).is_err());
        for (entry, tag) in [
            (
                call("Power", vec![number(-7), rational(1, 2)]),
                "InvalidLattice",
            ),
            (
                json!({"kind": "symbol", "name": "MagneticTB`unbound"}),
                "InvalidLattice",
            ),
            (call("Power", vec![number(0), number(-1)]), "DivisionByZero"),
            (
                call("Sin", vec![number(1), number(2)]),
                "UnsupportedExactExpression",
            ),
        ] {
            let error = resolve_model_lattice(&context, &lattice([number(0), number(0), entry]))
                .err()
                .unwrap();
            assert_eq!(error.tag(), tag);
        }
    }
}
