struct ScalarConstraintResult {
    constraints: Value,
    solved: Option<Value>,
}

fn scalar_constraint_result(
    request: &serde_json::Map<String, Value>,
    context: &Arc<CyclotomicContext>,
    symmetry: &CompiledSeitzGroup,
    sites: &SitePermutationData,
    bonds: &[PeriodicBondRecord],
    directed: &DirectedBondOrbitData,
) -> ExactResult<ScalarConstraintResult> {
    let site_count = sites.site_orbits().iter().map(Vec::len).sum::<usize>();
    let identity = ExactMatrix::identity(context, 1)?;
    let site_actions = (0..symmetry.group().order())
        .map(|_| vec![identity.clone(); site_count])
        .collect::<Vec<_>>();
    let constraints = compile_hermitian_bond_constraints(
        symmetry,
        bonds,
        directed,
        &vec![1; site_count],
        &site_actions,
    )?;
    let constraint_json = json!({
        "parameter_count": constraints.parameter_count(),
        "orbits": constraints.orbits().iter().map(|orbit| json!({
            "source_directed_orbit_indices": orbit.source_directed_orbit_indices(),
            "representative_bond_index": orbit.representative_bond_index(),
            "dimensions": [orbit.dimensions().0, orbit.dimensions().1],
            "stabilizer_generator_indices": orbit.stabilizer_generator_indices(),
            "reverse_representative_operation": orbit.reverse_representative_operation(),
            "nullity": orbit.kernel().nullity,
            "constraint_blocks": orbit.kernel().constraint_blocks.iter()
                .map(matrix_to_json).collect::<Vec<_>>(),
            "basis_matrix": matrix_to_json(&orbit.kernel().basis_matrix),
            "exact_residual_verified": orbit.kernel().exact_residual_verified
        })).collect::<Vec<_>>()
    });
    let solved = request
        .get("solve_and_verify")
        .is_some_and(|value| value == &Value::Bool(true))
        .then(|| {
            scalar_solved_result(
                context,
                symmetry,
                sites,
                bonds,
                directed,
                &constraints,
                &site_actions,
            )
        })
        .transpose()?;
    Ok(ScalarConstraintResult {
        constraints: constraint_json,
        solved,
    })
}

fn scalar_solved_result(
    context: &Arc<CyclotomicContext>,
    symmetry: &CompiledSeitzGroup,
    sites: &SitePermutationData,
    bonds: &[PeriodicBondRecord],
    directed: &DirectedBondOrbitData,
    constraints: &BondConstraintData,
    site_actions: &[Vec<ExactMatrix>],
) -> ExactResult<Value> {
    let site_count = sites.site_orbits().iter().map(Vec::len).sum::<usize>();
    let identity = ExactMatrix::identity(context, 1)?;
    let local_representations = vec![vec![identity; symmetry.group().order()]];
    let representation = compile_direct_product(
        symmetry.group(),
        sites.actions(),
        &local_representations,
        &symmetry.antiunitary_flags(),
    )
    .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
    parameter_space_result(&ParameterSpaceInput {
        context,
        symmetry,
        bonds,
        directed,
        constraints,
        local_dimensions: &vec![1; site_count],
        site_actions,
        representation_matrices: representation.representation_matrices(),
    })
}

struct ParameterSpaceInput<'a> {
    context: &'a Arc<CyclotomicContext>,
    symmetry: &'a CompiledSeitzGroup,
    bonds: &'a [PeriodicBondRecord],
    directed: &'a DirectedBondOrbitData,
    constraints: &'a BondConstraintData,
    local_dimensions: &'a [usize],
    site_actions: &'a [Vec<ExactMatrix>],
    representation_matrices: &'a [ExactMatrix],
}

fn parameter_space_result(input: &ParameterSpaceInput<'_>) -> ExactResult<Value> {
    let zero = Element::zero(input.context)?;
    let one = Element::one(input.context)?;
    let total_dimension = input.local_dimensions.iter().sum();
    let parameter_count = input.constraints.parameter_count();
    let mut parameter_solutions = Vec::with_capacity(parameter_count);
    let mut full_parameter_space_verified = true;
    for selected in 0..parameter_count {
        let parameters_by_orbit =
            selected_parameter_vector(input.constraints, &zero, &one, selected);
        let solved = reconstruct_solved_bond_model(
            input.symmetry,
            input.bonds,
            input.directed,
            input.constraints,
            input.local_dimensions,
            input.site_actions,
            &parameters_by_orbit,
        )?;
        let coefficients = fourier_coefficients(&solved)?;
        let gamma = coefficients.iter().try_fold(
            ExactMatrix::zero(input.context, total_dimension, total_dimension)?,
            |total, coefficient| add(&total, coefficient.matrix()),
        )?;
        let verification =
            verify_hamiltonian_symmetry(&solved, input.symmetry, input.representation_matrices)?;
        full_parameter_space_verified &= verification.verified();
        parameter_solutions.push(solution_to_json(
            selected,
            &solved,
            &coefficients,
            &gamma,
            &verification,
        ));
    }
    Ok(json!({
        "parameter_count": parameter_count,
        "representation_matrices": input.representation_matrices.iter()
            .map(matrix_to_json).collect::<Vec<_>>(),
        "parameter_basis_solutions": parameter_solutions,
        "full_parameter_space_verified": full_parameter_space_verified
    }))
}

fn selected_parameter_vector(
    constraints: &BondConstraintData,
    zero: &Element,
    one: &Element,
    selected: usize,
) -> Vec<Vec<Element>> {
    let mut result = constraints
        .orbits()
        .iter()
        .map(|orbit| vec![zero.clone(); orbit.kernel().nullity])
        .collect::<Vec<_>>();
    let mut offset = 0;
    for (orbit, parameters) in constraints.orbits().iter().zip(&mut result) {
        if selected >= offset && selected < offset + orbit.kernel().nullity {
            parameters[selected - offset] = one.clone();
            break;
        }
        offset += orbit.kernel().nullity;
    }
    result
}

fn solution_to_json(
    selected: usize,
    solved: &magnetictb_tight_binding::SolvedBondModel,
    coefficients: &[magnetictb_tight_binding::FourierCoefficient],
    gamma: &ExactMatrix,
    verification: &magnetictb_tight_binding::HamiltonianSymmetryVerification,
) -> Value {
    json!({
        "parameter_index": selected,
        "representative_hoppings": solved.representative_hoppings().iter()
            .map(matrix_to_json).collect::<Vec<_>>(),
        "terms": solved.terms().iter().map(|term| json!({
            "bond_index": term.bond_index(),
            "orbit_index": term.orbit_index(),
            "row_block": term.row_block(),
            "column_block": term.column_block(),
            "displacement": vector_to_json(term.displacement()),
            "matrix": matrix_to_json(term.matrix())
        })).collect::<Vec<_>>(),
        "fourier_coefficients": coefficients.iter().map(|coefficient| json!({
            "displacement": vector_to_json(coefficient.displacement()),
            "matrix": matrix_to_json(coefficient.matrix())
        })).collect::<Vec<_>>(),
        "gamma_hamiltonian": matrix_to_json(gamma),
        "gamma_semilinear_by_operation": verification.gamma_semilinear_by_operation(),
        "gamma_unitary_commutator_by_operation":
            verification.gamma_unitary_commutator_by_operation(),
        "fourier_covariance_by_operation": verification.fourier_covariance_by_operation(),
        "verified": verification.verified()
    })
}

fn directed_orbit_result(
    request: &serde_json::Map<String, Value>,
    catalog: &DataCatalog,
) -> ExactResult<Value> {
    let (requested_shells, shells) = bond_shells(request)?;
    let shell_index: usize = from_value(required(request, "shell_index")?, "shell_index")?;
    if shell_index >= requested_shells {
        return Err(ExactError::new(
            "InvalidPeriodicBondInput",
            "shell_index is outside requested_shells",
        ));
    }
    let context = context_from_json(required(request, "context")?, 128)?;
    let encoded_bindings = required(request, "bindings")?
        .as_object()
        .ok_or_else(|| malformed("bindings must be an object"))?;
    let bindings = encoded_bindings
        .iter()
        .map(|(name, value)| {
            element_from_json(&context, value).map(|element| (name.clone(), element))
        })
        .collect::<ExactResult<BTreeMap<_, _>>>()?;
    let msg_id = required(request, "msg_id")?
        .as_str()
        .ok_or_else(|| malformed("msg_id must be a string"))?;
    let source_ordinal: usize = from_value(required(request, "source_ordinal")?, "source_ordinal")?;
    let site_data = compile_msg_wyckoff_sites(catalog, msg_id, source_ordinal, &bindings)
        .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
    let compiled = compile_directed_bond_orbits(
        site_data.symmetry(),
        site_data.sites(),
        shells[shell_index].bonds(),
    )?;
    let mut bond_to_orbit = vec![usize::MAX; shells[shell_index].bonds().len()];
    for (orbit_index, orbit) in compiled.orbits().iter().enumerate() {
        for member in orbit.members() {
            bond_to_orbit[member.bond_index()] = orbit_index;
        }
    }
    let reverse_orbit_indices = compiled
        .orbits()
        .iter()
        .map(|orbit| bond_to_orbit[compiled.reverse_indices()[orbit.representative_bond_index()]])
        .collect::<Vec<_>>();
    let mut result = json!({
        "action_table": compiled.action().action_table(),
        "reverse_indices": compiled.reverse_indices(),
        "reverse_orbit_indices": reverse_orbit_indices,
        "orbits": compiled.orbits().iter().map(|orbit| json!({
            "representative_bond_index": orbit.representative_bond_index(),
            "members": orbit.members().iter().map(|member| json!({
                "bond_index": member.bond_index(),
                "transporter_operation": member.transporter_operation()
            })).collect::<Vec<_>>()
        })).collect::<Vec<_>>()
    });
    if request
        .get("scalar_constraints")
        .is_some_and(|value| value == &Value::Bool(true))
    {
        let scalar = scalar_constraint_result(
            request,
            &context,
            site_data.symmetry(),
            site_data.sites(),
            shells[shell_index].bonds(),
            &compiled,
        )?;
        result["scalar_constraints"] = scalar.constraints;
        if let Some(solved) = scalar.solved {
            result["scalar_solved_model"] = solved;
        }
    }
    Ok(result)
}

fn explicit_scalar_model_result(
    request: &serde_json::Map<String, Value>,
    catalog: &DataCatalog,
) -> ExactResult<Value> {
    let (requested_shells, shells) = bond_shells(request)?;
    let shell_index: usize = from_value(required(request, "shell_index")?, "shell_index")?;
    if shell_index >= requested_shells {
        return Err(ExactError::new(
            "InvalidPeriodicBondInput",
            "shell_index is outside requested_shells",
        ));
    }
    let context = context_from_json(required(request, "context")?, 128)?;
    let msg_id = required(request, "msg_id")?
        .as_str()
        .ok_or_else(|| malformed("msg_id must be a string"))?;
    let source = catalog
        .msg(msg_id)
        .ok_or_else(|| ExactError::new("UnknownStableId", format!("unknown MSG ID {msg_id}")))?;
    let symmetry = compile_msg_group_in_context(source, &context)
        .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
    let site_orbits = exact_site_orbits(&context, required(request, "site_orbits")?)?;
    let sites = compile_site_permutations(
        symmetry.group(),
        site_orbits,
        symmetry
            .operations()
            .iter()
            .map(|operation| operation.spatial().clone())
            .collect(),
    )
    .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
    let bonds = shells[shell_index].bonds();
    let directed = compile_directed_bond_orbits(&symmetry, &sites, bonds)?;
    let scalar = scalar_constraint_result(request, &context, &symmetry, &sites, bonds, &directed)?;
    Ok(json!({
        "msg_id": msg_id,
        "group_order": symmetry.group().order(),
        "antiunitary_flags": symmetry.antiunitary_flags(),
        "image_site_indices": sites.image_site_indices(),
        "bond_count": bonds.len(),
        "scalar_constraints": scalar.constraints,
        "scalar_solved_model": scalar.solved
    }))
}

fn encoded_rational_from_element(value: &Element) -> ExactResult<Value> {
    if value
        .coefficients()
        .iter()
        .skip(1)
        .any(|coefficient| !coefficient.is_zero())
    {
        return Err(ExactError::new(
            "NonrationalSiteCoordinate",
            "periodic site coordinates must be exact rational values",
        ));
    }
    Ok(encoded_rational_value(&value.coefficients()[0]))
}

fn encoded_rational_value(coefficient: &Rational) -> Value {
    if coefficient.denominator_string() == "1" {
        json!({
            "kind": "integer",
            "value": coefficient.numerator_string()
        })
    } else {
        json!({
            "kind": "rational",
            "numerator": coefficient.numerator_string(),
            "denominator": coefficient.denominator_string()
        })
    }
}

fn rational_quadratic_coordinates(value: &Element) -> ExactResult<Vec<Rational>> {
    let context = value.context();
    let one = Element::one(context)?;
    let sqrt_two =
        Element::root_of_unity(context, 8, 1)?.add(&Element::root_of_unity(context, 8, -1)?)?;
    let sqrt_three =
        Element::root_of_unity(context, 12, 1)?.add(&Element::root_of_unity(context, 12, -1)?)?;
    let sqrt_six = sqrt_two.multiply(&sqrt_three)?;
    let basis = [one, sqrt_two, sqrt_three, sqrt_six];
    let mut augmented = (0..context.degree)
        .map(|row| {
            let mut entries = basis
                .iter()
                .map(|element| element.coefficients()[row].clone())
                .collect::<Vec<_>>();
            entries.push(value.coefficients()[row].clone());
            entries
        })
        .collect::<Vec<_>>();
    let mut pivot_row = 0;
    let mut pivot_rows = Vec::new();
    for column in 0..4 {
        let Some(row) = (pivot_row..augmented.len()).find(|&row| !augmented[row][column].is_zero())
        else {
            return Err(ExactError::new(
                "UnsupportedBondExactExpression",
                "quadratic radical basis is singular in the selected cyclotomic field",
            ));
        };
        augmented.swap(pivot_row, row);
        let pivot = augmented[pivot_row][column].clone();
        for entry in &mut augmented[pivot_row][column..=4] {
            *entry = entry.divide(&pivot)?;
        }
        for row in 0..augmented.len() {
            if row == pivot_row || augmented[row][column].is_zero() {
                continue;
            }
            let factor = augmented[row][column].clone();
            let pivot_values = augmented[pivot_row][column..=4].to_vec();
            for (entry, pivot_value) in augmented[row][column..=4].iter_mut().zip(pivot_values) {
                *entry = entry.subtract(&factor.multiply(&pivot_value));
            }
        }
        pivot_rows.push(pivot_row);
        pivot_row += 1;
    }
    if augmented
        .iter()
        .any(|row| row[..4].iter().all(Rational::is_zero) && !row[4].is_zero())
    {
        return Err(ExactError::new(
            "UnsupportedBondExactExpression",
            "exact lattice value lies outside Q(sqrt(2), sqrt(3))",
        ));
    }
    Ok(pivot_rows
        .into_iter()
        .map(|row| augmented[row][4].clone())
        .collect())
}

fn quadratic_expression_from_element(value: &Element) -> ExactResult<Value> {
    let coordinates = rational_quadratic_coordinates(value)?;
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
    for (index, coordinate) in coordinates.iter().enumerate() {
        if coordinate.is_zero() {
            continue;
        }
        let coefficient = encoded_rational_value(coordinate);
        if index == 0 {
            terms.push(coefficient);
            continue;
        }
        let mut factors = Vec::new();
        if coordinate.numerator_string() != "1" || coordinate.denominator_string() != "1" {
            factors.push(coefficient);
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
            json!({
                "kind": "call",
                "head": "System`Times",
                "arguments": factors
            })
        });
    }
    Ok(match terms.len() {
        0 => json!({"kind": "integer", "value": "0"}),
        1 => terms.remove(0),
        _ => json!({
            "kind": "call",
            "head": "System`Plus",
            "arguments": terms
        }),
    })
}

fn cyclotomic_expression_from_element(value: &Element) -> Value {
    let mut terms = Vec::new();
    for (power, coefficient) in value.coefficients().iter().enumerate() {
        if coefficient.is_zero() {
            continue;
        }
        let coefficient_value = encoded_rational_value(coefficient);
        if power == 0 {
            terms.push(coefficient_value);
            continue;
        }
        let root = json!({
            "kind": "root_of_unity",
            "order": value.context().conductor,
            "power": power
        });
        terms.push(
            if coefficient.numerator_string() == "1" && coefficient.denominator_string() == "1" {
                root
            } else {
                json!({
                    "kind": "call",
                    "head": "System`Times",
                    "arguments": [coefficient_value, root]
                })
            },
        );
    }
    match terms.len() {
        0 => json!({"kind": "integer", "value": "0"}),
        1 => terms.remove(0),
        _ => json!({
            "kind": "call",
            "head": "System`Plus",
            "arguments": terms
        }),
    }
}

fn quadratic_from_element(value: &Element) -> ExactResult<QuadraticReal> {
    let coefficients = rational_quadratic_coordinates(value)?;
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

fn phase_coordinate_expression(value: &Value) -> ExactResult<Value> {
    if value
        .as_object()
        .is_some_and(|record| record.contains_key("basis"))
    {
        return Ok(quadratic_expression(&quadratic_from_json(value)?));
    }
    let record = object(value)?;
    if record.contains_key("kind") {
        Ok(value.clone())
    } else {
        Err(malformed("phase coordinate must be exact"))
    }
}

fn shifted_phase_coordinate(displacement: &Value, shift: &Element) -> ExactResult<Value> {
    if shift.is_zero() {
        return Ok(displacement.clone());
    }
    if let Ok(displacement) = quadratic_from_json(displacement)
        && let Ok(shift) = quadratic_from_element(shift)
    {
        return Ok(quadratic_to_json(&displacement.add(&shift)));
    }
    let displacement = phase_coordinate_expression(displacement)?;
    let shift = quadratic_expression_from_element(shift)
        .unwrap_or_else(|_| cyclotomic_expression_from_element(shift));
    if displacement == json!({"kind": "integer", "value": "0"}) {
        return Ok(shift);
    }
    Ok(json!({
        "kind": "call",
        "head": "System`Plus",
        "arguments": [displacement, shift]
    }))
}

fn convention_ii_matrix(request: &serde_json::Map<String, Value>) -> ExactResult<Value> {
    let context = context_from_json(required(request, "context")?, 128)?;
    let rows = array(required(request, "matrix_entries")?)?;
    if rows.is_empty() || rows.iter().any(|row| array(row).is_err()) {
        return Err(ExactError::new(
            "InvalidConventionIIInput",
            "symhamII requires a nonempty rectangular Hamiltonian matrix",
        ));
    }
    let dimension = rows.len();
    if rows
        .iter()
        .any(|row| array(row).is_ok_and(|row| row.len() != dimension))
    {
        return Err(ExactError::new(
            "InvalidConventionIIInput",
            "symhamII requires a square Hamiltonian matrix",
        ));
    }
    let centers = array(required(request, "centers")?)?;
    if centers.len() != dimension {
        return Err(ExactError::new(
            "InvalidWannierCenters",
            format!(
                "symhamII requires one Wannier center per Hamiltonian row; expected {dimension}, received {}",
                centers.len()
            ),
        ));
    }
    let centers = centers
        .iter()
        .map(|center| {
            let center = array(center)?;
            if center.len() != 3 {
                return Err(ExactError::new(
                    "InvalidWannierCenters",
                    "each symhamII Wannier center must have three fractional coordinates",
                ));
            }
            center
                .iter()
                .map(|coordinate| encoded_element_from_json(&context, coordinate))
                .collect::<ExactResult<Vec<_>>>()
        })
        .collect::<ExactResult<Vec<_>>>()?;
    let transformed = rows
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            array(row)?
                .iter()
                .enumerate()
                .map(|(column_index, cell)| {
                    array(cell)?
                        .iter()
                        .map(|contribution| {
                            let mut transformed = object(contribution)?.clone();
                            let displacement = array(required(&transformed, "displacement")?)?;
                            if displacement.len() != 3 {
                                return Err(ExactError::new(
                                    "InvalidConventionIIInput",
                                    "every symhamII Fourier displacement must have three coordinates",
                                ));
                            }
                            let shifted = displacement
                                .iter()
                                .enumerate()
                                .map(|(axis, coordinate)| {
                                    let shift = centers[row_index][axis]
                                        .subtract(&centers[column_index][axis])?;
                                    shifted_phase_coordinate(coordinate, &shift)
                                })
                                .collect::<ExactResult<Vec<_>>>()?;
                            transformed.insert("displacement".to_owned(), Value::Array(shifted));
                            Ok(Value::Object(transformed))
                        })
                        .collect::<ExactResult<Vec<_>>>()
                        .map(Value::Array)
                })
                .collect::<ExactResult<Vec<_>>>()
                .map(Value::Array)
        })
        .collect::<ExactResult<Vec<_>>>()?;
    Ok(json!({
        "schema": "magnetictb.symham_ii_matrix.v1",
        "version": 1,
        "field_context": context_to_json(&context),
        "gauge": {
            "momentum_coordinates": "stable_kx_ky_kz",
            "phase_convention": "exp(+i*k_dot_fractional_displacement_with_orbital_centers)",
            "row_axis": "bra_destination",
            "column_axis": "ket_source",
            "bz_folding": false
        },
        "matrix_entries": transformed
    }))
}
