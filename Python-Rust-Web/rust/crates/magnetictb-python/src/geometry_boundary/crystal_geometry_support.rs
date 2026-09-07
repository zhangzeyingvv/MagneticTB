fn spatial_operations(
    context: &Arc<CyclotomicContext>,
    value: &Value,
) -> ExactResult<Vec<SpatialOperation>> {
    array(value)?
        .iter()
        .map(|operation| {
            let record = object(operation)?;
            SpatialOperation::new(
                encoded_matrix_from_json(context, required(record, "rotation")?)?,
                exact_elements(context, required(record, "translation")?)?,
            )
            .map_err(|error| ExactError::new(error.tag(), error.to_string()))
        })
        .collect()
}

struct BondInputs {
    sites: Vec<Vec<QuadraticReal>>,
    lattice: Vec<Vec<QuadraticReal>>,
    requested_shells: usize,
}

fn bond_inputs(request: &serde_json::Map<String, Value>) -> ExactResult<BondInputs> {
    let sites = exact_quadratic_rows(required(request, "sites")?, "sites")?;
    let lattice = exact_quadratic_rows(required(request, "lattice")?, "lattice")?;
    let requested_shells: usize =
        from_value(required(request, "requested_shells")?, "requested_shells")?;
    if requested_shells == 0 {
        return Err(ExactError::new(
            "InvalidPeriodicBondInput",
            "requested_shells must be positive",
        ));
    }
    Ok(BondInputs {
        sites,
        lattice,
        requested_shells,
    })
}

fn adaptive_options(request: &serde_json::Map<String, Value>) -> ExactResult<BondSearchOptions> {
    let mut options = BondSearchOptions::default();
    if let Some(value) = request.get("minimum_neighbors_per_site") {
        options.minimum_neighbors_per_site = from_value(value, "minimum_neighbors_per_site")?;
    }
    if let Some(value) = request.get("initial_radius") {
        options.initial_radius = Some(from_value(value, "initial_radius")?);
    }
    if let Some(value) = request.get("growth_factor") {
        options.growth_factor = from_value(value, "growth_factor")?;
    }
    if let Some(value) = request.get("maximum_iterations") {
        options.maximum_iterations = from_value(value, "maximum_iterations")?;
    }
    if let Some(value) = request.get("maximum_radius") {
        options.maximum_radius = Some(from_value(value, "maximum_radius")?);
    }
    if let Some(value) = request.get("distance_tolerance") {
        options.distance_tolerance = from_value(value, "distance_tolerance")?;
    }
    Ok(options)
}

fn bond_shells(
    request: &serde_json::Map<String, Value>,
) -> ExactResult<(usize, Vec<PeriodicBondShell>)> {
    let inputs = bond_inputs(request)?;
    let requested_shells = inputs.requested_shells;
    let shells = if let Some(value) = request.get("translation_bound") {
        let bound: i64 = from_value(value, "translation_bound")?;
        periodic_bond_shells_in_box(&inputs.sites, &inputs.lattice, bound)
            .map_err(|error| ExactError::new(error.tag(), error.to_string()))?
    } else {
        find_periodic_bond_shells(
            &inputs.sites,
            &inputs.lattice,
            requested_shells,
            &adaptive_options(request)?,
        )
        .map_err(|error| ExactError::new(error.tag(), error.to_string()))?
        .shells()
        .to_vec()
    };
    if shells.len() < requested_shells {
        return Err(ExactError::new(
            "InsufficientPeriodicBondShells",
            format!(
                "bond search contains {} shells but {requested_shells} were requested",
                shells.len()
            ),
        ));
    }
    Ok((requested_shells, shells))
}

fn box_result(request: &serde_json::Map<String, Value>) -> ExactResult<Value> {
    let inputs = bond_inputs(request)?;
    let requested_shells = inputs.requested_shells;
    let translation_bound: i64 =
        from_value(required(request, "translation_bound")?, "translation_bound")?;
    let shells = periodic_bond_shells_in_box(&inputs.sites, &inputs.lattice, translation_bound)
        .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
    if shells.len() < requested_shells {
        return Err(ExactError::new(
            "InsufficientPeriodicBondShells",
            format!(
                "translation box contains {} shells but {requested_shells} were requested",
                shells.len()
            ),
        ));
    }
    let enumerated_translation_count = shells
        .iter()
        .map(|shell| shell.bonds().len())
        .sum::<usize>();
    Ok(json!({
        "translation_bound": translation_bound,
        "box_complete": true,
        "global_search_complete": false,
        "enumerated_translation_count": enumerated_translation_count,
        "shells": shells.iter().take(requested_shells).map(|shell| json!({
            "squared_distance": quadratic_to_json(shell.squared_distance()),
            "bonds": shell.bonds().iter().map(periodic_bond_to_json).collect::<Vec<_>>()
        })).collect::<Vec<_>>()
    }))
}

fn adaptive_result(request: &serde_json::Map<String, Value>) -> ExactResult<Value> {
    if request.contains_key("translation_bound") {
        return Err(ExactError::new(
            "InvalidPeriodicBondInput",
            "translation_bound belongs to the explicit in-box diagnostic endpoint",
        ));
    }
    let inputs = bond_inputs(request)?;
    let requested_shells = inputs.requested_shells;
    let result: PeriodicBondSearchResult = find_periodic_bond_shells(
        &inputs.sites,
        &inputs.lattice,
        requested_shells,
        &adaptive_options(request)?,
    )
    .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
    Ok(json!({
        "adaptive_search_complete": true,
        "radius": result.radius(),
        "requested_shells": result.requested_shells(),
        "minimum_neighbors_per_site": result.minimum_neighbors_per_site(),
        "neighbor_counts_per_site": result.neighbor_counts_per_site(),
        "candidate_bond_count": result.candidate_bond_count(),
        "enumerated_translation_count": result.enumerated_translation_count(),
        "iterations": result.iterations(),
        "shells": result.shells().iter().take(requested_shells).map(|shell| json!({
            "numeric_distance": shell.numeric_distance(),
            "squared_distance": quadratic_to_json(shell.squared_distance()),
            "bonds": shell.bonds().iter().map(periodic_bond_to_json).collect::<Vec<_>>()
        })).collect::<Vec<_>>()
    }))
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct MagneticSite {
    position: Vec<Element>,
    moment: Vec<Element>,
}

fn determinant_three_exact(matrix: &ExactMatrix) -> ExactResult<Element> {
    if matrix.rows() != 3 || matrix.columns() != 3 {
        return Err(ExactError::new(
            "InvalidBrokenSymmetryInput",
            "magnetic-site rotations must be exact 3 by 3 matrices",
        ));
    }
    let positive = matrix
        .entry(0, 0)?
        .multiply(matrix.entry(1, 1)?)?
        .multiply(matrix.entry(2, 2)?)?
        .add(
            &matrix
                .entry(0, 1)?
                .multiply(matrix.entry(1, 2)?)?
                .multiply(matrix.entry(2, 0)?)?,
        )?
        .add(
            &matrix
                .entry(0, 2)?
                .multiply(matrix.entry(1, 0)?)?
                .multiply(matrix.entry(2, 1)?)?,
        )?;
    let negative = matrix
        .entry(0, 2)?
        .multiply(matrix.entry(1, 1)?)?
        .multiply(matrix.entry(2, 0)?)?
        .add(
            &matrix
                .entry(0, 1)?
                .multiply(matrix.entry(1, 0)?)?
                .multiply(matrix.entry(2, 2)?)?,
        )?
        .add(
            &matrix
                .entry(0, 0)?
                .multiply(matrix.entry(1, 2)?)?
                .multiply(matrix.entry(2, 1)?)?,
        )?;
    positive.subtract(&negative)
}

fn matrix_vector_exact(matrix: &ExactMatrix, vector: &[Element]) -> ExactResult<Vec<Element>> {
    if matrix.columns() != vector.len() {
        return Err(ExactError::new(
            "InvalidBrokenSymmetryInput",
            "magnetic-site vector dimension does not match its operation",
        ));
    }
    (0..matrix.rows())
        .map(|row| {
            vector
                .iter()
                .enumerate()
                .try_fold(Element::zero(matrix.context())?, |sum, (column, value)| {
                    sum.add(&matrix.entry(row, column)?.multiply(value)?)
                })
        })
        .collect()
}

fn transform_magnetic_site(
    site: &MagneticSite,
    spatial: &SpatialOperation,
    spin_rotation: Option<&ExactMatrix>,
    antiunitary: bool,
) -> ExactResult<MagneticSite> {
    let mut position = matrix_vector_exact(spatial.rotation(), &site.position)?;
    if position.len() != spatial.translation().len() {
        return Err(ExactError::new(
            "InvalidBrokenSymmetryInput",
            "magnetic-site position and translation dimensions differ",
        ));
    }
    for (coordinate, translation) in position.iter_mut().zip(spatial.translation()) {
        *coordinate = fractional_part(&coordinate.add(translation)?)
            .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
    }
    let mut moment = if let Some(spin_rotation) = spin_rotation {
        matrix_vector_exact(spin_rotation, &site.moment)?
    } else {
        let determinant = determinant_three_exact(spatial.rotation())?;
        matrix_vector_exact(spatial.rotation(), &site.moment)?
            .into_iter()
            .map(|value| value.multiply(&determinant))
            .collect::<ExactResult<Vec<_>>>()?
    };
    if antiunitary {
        moment = moment
            .into_iter()
            .map(|value| value.negate())
            .collect::<ExactResult<Vec<_>>>()?;
    }
    Ok(MagneticSite { position, moment })
}

fn magnetic_orbit(
    seed: &MagneticSite,
    operation_indices: &[usize],
    spatial: &[SpatialOperation],
    spin_rotations: Option<&[ExactMatrix]>,
    antiunitary_flags: &[bool],
) -> ExactResult<Vec<MagneticSite>> {
    let mut orbit = Vec::new();
    for &index in operation_indices {
        let image = transform_magnetic_site(
            seed,
            spatial.get(index).ok_or_else(|| {
                ExactError::new(
                    "InvalidBrokenSymmetryInput",
                    "subgroup operation index is out of range",
                )
            })?,
            spin_rotations.and_then(|rotations| rotations.get(index)),
            *antiunitary_flags.get(index).ok_or_else(|| {
                ExactError::new(
                    "InvalidBrokenSymmetryInput",
                    "antiunitary metadata is not aligned with operations",
                )
            })?,
        )?;
        if !orbit.contains(&image) {
            orbit.push(image);
        }
    }
    Ok(orbit)
}
