fn evaluate(
    request: &serde_json::Map<String, Value>,
    catalog: Option<&DataCatalog>,
) -> ExactResult<Value> {
    let operation = required(request, "operation")?
        .as_str()
        .ok_or_else(|| malformed("operation must be a string"))?;
    match operation {
        "periodic_bond_shells_in_box" => box_result(request),
        "find_periodic_bond_shells" => adaptive_result(request),
        "compile_data_directed_bond_orbits" => directed_orbit_result(
            request,
            catalog.ok_or_else(|| {
                ExactError::new(
                    "MissingDataCatalog",
                    "this geometry operation requires the frozen Rust Data catalog",
                )
            })?,
        ),
        "compile_data_scalar_model_closure" => explicit_scalar_model_result(
            request,
            catalog.ok_or_else(|| {
                ExactError::new(
                    "MissingDataCatalog",
                    "this model operation requires the frozen Rust Data catalog",
                )
            })?,
        ),
        "compile_explicit_block_model_closure" => explicit_block_model_result(request),
        "compile_model_input_closure" => compiler_input_model_result(request),
        "compile_model_input_closures" => compiler_input_all_model_results(request),
        "compile_data_model_input_closure" => data_compiler_input_model_result(
            request,
            catalog.ok_or_else(|| {
                ExactError::new(
                    "MissingDataCatalog",
                    "this model-input operation requires the frozen Rust Data catalog",
                )
            })?,
        ),
        "compile_data_model_input_closures" => data_compiler_input_all_model_results(
            request,
            catalog.ok_or_else(|| {
                ExactError::new(
                    "MissingDataCatalog",
                    "this model-input operation requires the frozen Rust Data catalog",
                )
            })?,
        ),
        "combine_symbolic_hamiltonians" => combine_symbolic_hamiltonians(request),
        "unconstrained_symbolic_hamiltonian" => unconstrained_symbolic_hamiltonian(request),
        "cartesian_symbolic_hamiltonian" => cartesian_symbolic_hamiltonian(request),
        "symham_ii" => convention_ii_matrix(request),
        "broken_symmetry_init_rules" => broken_symmetry_init_rules(request),
        "crystal_structure_data" => crystal_structure_data(request, catalog),
        "evaluate_symbolic_hamiltonian" => evaluate_symbolic_hamiltonian(request),
        "hopping_data_from_shell_results" => hopping_data_from_shell_results(request),
        "hopping_data_from_symbolic_matrix" => hopping_data_from_symbolic_matrix(request),
        "evaluate_symbolic_expression_matrix" => evaluate_symbolic_expression_matrix(request),
        "evaluate_symbolic_expression_matrices" => evaluate_symbolic_expression_matrices(request),
        "standard_k_path" => standard_k_path_result(request),
        "brillouin_zone_data" => brillouin_zone_data(request),
        "point_matrix" => point_matrix_result(request),
        "band_corep_trace" => band_corep_trace(request),
        _ => Err(ExactError::new(
            "UnsupportedOperation",
            format!("unsupported geometry operation: {operation}"),
        )),
    }
}

fn parse_request(payload: &str) -> ExactResult<Value> {
    let mut value: Value = serde_json::from_str(payload)
        .map_err(|error| malformed(format!("invalid JSON: {error}")))?;
    normalize_stable_inexact(&mut value)?;
    Ok(value)
}

fn serialize_result(value: &Value) -> ExactResult<String> {
    serde_json::to_string(value)
        .map_err(|error| malformed(format!("result serialization failed: {error}")))
}

pub fn compute(payload: &str) -> ExactResult<String> {
    let value = parse_request(payload)?;
    let request = value
        .as_object()
        .ok_or_else(|| malformed("request must be an object"))?;
    serialize_result(&evaluate(request, None)?)
}

pub fn compute_with_catalog(payload: &str, catalog: &DataCatalog) -> ExactResult<String> {
    let value = parse_request(payload)?;
    let request = value
        .as_object()
        .ok_or_else(|| malformed("request must be an object"))?;
    serialize_result(&evaluate(request, Some(catalog))?)
}

#[cfg(test)]
mod tests {
    use cyclotomic_nullspace::fixture::{context_to_json, element_to_json, matrix_to_json};
    use cyclotomic_nullspace::{Element, ExactMatrix, make_context};
    use serde_json::{Value, json};

    use super::{
        cartesian_displacement, combine_symbolic_hamiltonians, convention_ii_matrix,
        evaluate_symbolic_hamiltonian, normalize_stable_inexact, stable_machine_real,
        symbolic_cartesian_displacement, vector_to_json,
    };
    use magnetictb_crystal_geometry::QuadraticReal;

    #[test]
    fn cartesian_momentum_uses_row_vector_lattice_displacements() {
        let one = QuadraticReal::one();
        let two = one.add(&one);
        let three = two.add(&one);
        let displacement = vector_to_json(&[one.clone(), one]);
        let lattice = json!({
            "kind": "matrix",
            "rows": 2,
            "columns": 2,
            "entries": [
                {"kind": "integer", "value": "2"},
                {"kind": "integer", "value": "0"},
                {"kind": "integer", "value": "0"},
                {"kind": "integer", "value": "3"}
            ]
        });
        assert_eq!(
            cartesian_displacement(&displacement, &lattice).expect("Cartesian displacement"),
            vector_to_json(&[two, three])
        );

        let half = QuadraticReal::from_rational(
            cyclotomic_nullspace::Rational::parse("1", "2").expect("half"),
        );
        let symbolic_displacement = vector_to_json(&[half, QuadraticReal::one()]);
        let symbolic_lattice = json!({
            "kind": "matrix",
            "rows": 2,
            "columns": 2,
            "entries": [
                {"kind": "symbol", "name": "MagneticTB`a"},
                {"kind": "integer", "value": "0"},
                {"kind": "integer", "value": "0"},
                {"kind": "integer", "value": "3"}
            ]
        });
        assert_eq!(
            symbolic_cartesian_displacement(&symbolic_displacement, &symbolic_lattice)
                .expect("symbolic Cartesian displacement"),
            json!([
                {
                    "kind": "call",
                    "head": "System`Times",
                    "arguments": [
                        {"kind": "rational", "numerator": "1", "denominator": "2"},
                        {"kind": "symbol", "name": "MagneticTB`a"}
                    ]
                },
                {"kind": "integer", "value": "3"}
            ])
        );
    }

    #[test]
    fn stable_machine_numbers_preserve_the_full_decimal_metric() {
        assert_eq!(
            stable_machine_real("1.2344").expect("finite real"),
            json!({"kind": "rational", "numerator": "1543", "denominator": "1250"})
        );
        assert_eq!(
            stable_machine_real("-0.0004").expect("finite real"),
            json!({"kind": "rational", "numerator": "-1", "denominator": "2500"})
        );
        assert_eq!(
            stable_machine_real("1.25").expect("finite real"),
            json!({"kind": "rational", "numerator": "5", "denominator": "4"})
        );
        assert_eq!(
            stable_machine_real("1e-6").expect("finite real"),
            json!({"kind": "rational", "numerator": "1", "denominator": "1000000"})
        );
        assert_eq!(
            stable_machine_real("nan").expect_err("NaN must fail").tag(),
            "InvalidInexactInput"
        );

        let mut complex = json!({
            "kind": "machine_complex",
            "real_decimal": "0.5",
            "imaginary_decimal": "-0.25"
        });
        normalize_stable_inexact(&mut complex).expect("finite complex");
        assert_eq!(
            complex,
            json!({
                "kind": "call",
                "head": "System`Complex",
                "arguments": [
                    {"kind": "rational", "numerator": "1", "denominator": "2"},
                    {"kind": "rational", "numerator": "-1", "denominator": "4"}
                ]
            })
        );
    }

    fn zero_displacement() -> Value {
        Value::Array(
            (0..3)
                .map(|_| {
                    json!({
                        "basis": ["1", "sqrt(2)", "sqrt(3)", "sqrt(6)"],
                        "coefficients": (0..4).map(|_| json!({
                            "numerator": "0", "denominator": "1"
                        })).collect::<Vec<_>>()
                    })
                })
                .collect(),
        )
    }

    fn hamiltonian(shell: usize, name: &str, identity: &str) -> Value {
        let context = make_context(1, 8).expect("rational context");
        let one = Element::one(&context).expect("one");
        let matrix =
            ExactMatrix::new(context.clone(), 1, 1, vec![one.clone()]).expect("one by one matrix");
        let parameter = json!({
            "name": name,
            "shell": shell,
            "parameter_number": 1,
            "core_parameter_index": 0
        });
        let displacement = zero_displacement();
        json!({
            "schema": "magnetictb.symbolic_hamiltonian.v1",
            "version": 1,
            "model_identity_sha256": identity,
            "field_context": context_to_json(&context),
            "gauge": {
                "momentum_coordinates": "stable_kx_ky_kz",
                "phase_convention": "exp(+i*k_dot_fractional_displacement)",
                "row_axis": "bra_destination",
                "column_axis": "ket_source",
                "bz_folding": false
            },
            "shape": [1, 1],
            "shells": [shell],
            "hermitian": true,
            "kernel_method": "Iterative",
            "validation_level": "Basic",
            "parameter_names": [name],
            "terms": [{
                "parameter": parameter.clone(),
                "fourier_coefficients": [{
                    "displacement": displacement.clone(),
                    "matrix": matrix_to_json(&matrix)
                }],
                "gamma_hamiltonian": matrix_to_json(&matrix),
                "verified": true
            }],
            "matrix_entries": [[[{
                "parameter": parameter,
                "displacement": displacement,
                "coefficient": element_to_json(&one)
            }]]],
            "covariance_verified": true
        })
    }

    #[test]
    fn rust_combines_distinct_shells_and_rejects_duplicates_or_mismatches() {
        let left = hamiltonian(1, "e1", "same-model");
        let right = hamiltonian(2, "t1", "same-model");
        let request = json!({"left": left, "right": right});
        let combined =
            combine_symbolic_hamiltonians(request.as_object().expect("combination request"))
                .expect("compatible shells");
        assert_eq!(combined["shells"], json!([1, 2]));
        assert_eq!(combined["parameter_names"], json!(["e1", "t1"]));
        assert_eq!(
            combined["matrix_entries"][0][0].as_array().map(Vec::len),
            Some(2)
        );

        let duplicate = json!({
            "left": hamiltonian(1, "e1", "same-model"),
            "right": hamiltonian(1, "e2", "same-model")
        });
        assert_eq!(
            combine_symbolic_hamiltonians(duplicate.as_object().expect("duplicate request"))
                .expect_err("duplicate shell")
                .tag(),
            "DuplicateHamiltonianShell"
        );

        let mismatch = json!({
            "left": hamiltonian(1, "e1", "model-a"),
            "right": hamiltonian(2, "t1", "model-b")
        });
        assert_eq!(
            combine_symbolic_hamiltonians(mismatch.as_object().expect("mismatch request"))
                .expect_err("different model")
                .tag(),
            "IncompatibleHamiltonian"
        );
    }

    #[test]
    fn rust_evaluates_the_canonical_symbolic_hamiltonian() {
        let request = json!({
            "hamiltonian": hamiltonian(1, "e1", "same-model"),
            "parameters": {"e1": {"kind": "rational", "numerator": "3", "denominator": "2"}},
            "momentum": [0, 0, 0]
        });
        let result = evaluate_symbolic_hamiltonian(
            request.as_object().expect("Hamiltonian evaluation request"),
        )
        .expect("numeric Hamiltonian");
        assert_eq!(result["shape"], json!([1, 1]));
        assert_eq!(result["matrix"][0][0]["real"], json!(1.5));
        assert_eq!(result["matrix"][0][0]["imaginary"], json!(0.0));
    }

    #[test]
    fn convention_ii_applies_row_minus_column_wannier_centers() {
        let context = make_context(1, 8).expect("rational context");
        let one = Element::one(&context).expect("one");
        let zero = zero_displacement();
        let request = json!({
            "context": context_to_json(&context),
            "matrix_entries": [
                [[], [{
                    "parameter": {"name": "t1"},
                    "displacement": zero,
                    "coefficient": element_to_json(&one)
                }]],
                [[], []]
            ],
            "centers": [
                [
                    {"kind": "integer", "value": "0"},
                    {"kind": "integer", "value": "0"},
                    {"kind": "integer", "value": "0"}
                ],
                [
                    {"kind": "rational", "numerator": "1", "denominator": "2"},
                    {"kind": "integer", "value": "0"},
                    {"kind": "integer", "value": "0"}
                ]
            ]
        });
        let result = convention_ii_matrix(request.as_object().expect("symhamII request"))
            .expect("convention-II matrix");
        assert_eq!(
            result["matrix_entries"][0][1][0]["displacement"],
            json!([
                {"kind": "rational", "numerator": "-1", "denominator": "2"},
                super::quadratic_to_json(&QuadraticReal::zero()),
                super::quadratic_to_json(&QuadraticReal::zero())
            ])
        );
    }
}
