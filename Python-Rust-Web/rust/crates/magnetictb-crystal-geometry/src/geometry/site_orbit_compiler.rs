pub fn compile_msg_wyckoff_sites(
    catalog: &DataCatalog,
    msg_id: &str,
    source_ordinal: usize,
    bindings: &BTreeMap<String, Element>,
) -> GeometryResult<MsgWyckoffSiteData> {
    compile_msg_wyckoff_site_selection(catalog, msg_id, source_ordinal, None, bindings)
}
/// Compiles a frozen Wyckoff entry using both source ordinal and letter.
/// Some Mathematica source records share the same first numeric field, so the
/// letter is part of the stable entry identity for an unambiguous Data lookup.
pub fn compile_msg_wyckoff_sites_by_letter(
    catalog: &DataCatalog,
    msg_id: &str,
    source_ordinal: usize,
    letter: &str,
    bindings: &BTreeMap<String, Element>,
) -> GeometryResult<MsgWyckoffSiteData> {
    compile_msg_wyckoff_site_selection(catalog, msg_id, source_ordinal, Some(letter), bindings)
}

fn compile_msg_wyckoff_site_selection(
    catalog: &DataCatalog,
    msg_id: &str,
    source_ordinal: usize,
    letter: Option<&str>,
    bindings: &BTreeMap<String, Element>,
) -> GeometryResult<MsgWyckoffSiteData> {
    let msg = catalog.msg(msg_id).ok_or_else(|| {
        GeometryError::new("UnknownStableId", format!("unknown MSG stable ID {msg_id}"))
    })?;
    let wyckoff = catalog.wyckoff(msg_id).ok_or_else(|| {
        GeometryError::new(
            "UnknownStableId",
            format!("no Wyckoff data closes over MSG stable ID {msg_id}"),
        )
    })?;
    let entry = wyckoff
        .entries()
        .iter()
        .find(|entry| {
            entry.source_ordinal() == source_ordinal
                && letter.is_none_or(|letter| entry.letter() == letter)
        })
        .ok_or_else(|| {
            GeometryError::new(
                "UnknownWyckoffOrdinal",
                format!(
                    "MSG {msg_id} has no Wyckoff source ordinal {source_ordinal}{}",
                    letter.map_or_else(String::new, |letter| format!(" with letter {letter}"))
                ),
            )
        })?;
    if entry.multiplicity() != entry.positions().len() {
        return Err(GeometryError::new(
            "InvalidWyckoffMultiplicity",
            "Wyckoff multiplicity differs from the ordered position count",
        ));
    }
    let context = if let Some(first) = bindings.values().next() {
        if bindings
            .values()
            .any(|value| value.context() != first.context())
        {
            return Err(GeometryError::new(
                "ConductorMismatch",
                "all Wyckoff bindings must use one exact cyclotomic field",
            ));
        }
        Arc::clone(first.context())
    } else {
        make_context(1, 128)?
    };
    let symmetry = compile_msg_group_in_context(msg, &context)?;
    let context = symmetry.operations()[0].spatial().rotation().context();
    let orbit = entry
        .positions()
        .iter()
        .map(|position| {
            position
                .coordinates()
                .iter()
                .map(|coordinate| {
                    let value = evaluate_data_expression(context, coordinate, bindings)?;
                    fractional_part(&value)
                })
                .collect::<Result<Vec<_>, _>>()
                .map_err(Into::into)
        })
        .collect::<GeometryResult<Vec<_>>>()?;
    let sites = compile_site_permutations(
        symmetry.group(),
        vec![orbit],
        symmetry
            .operations()
            .iter()
            .map(|operation| operation.spatial().clone())
            .collect(),
    )?;
    Ok(MsgWyckoffSiteData {
        source_ordinal,
        letter: entry.letter().to_owned(),
        symmetry,
        sites,
    })
}

fn validate_inputs(
    group: &GroupAlgebra,
    site_orbits: &[Vec<Vec<Element>>],
    spatial_actions: &[SpatialOperation],
) -> GeometryResult<usize> {
    let dimension = site_orbits
        .first()
        .and_then(|orbit| orbit.first())
        .map_or(0, Vec::len);
    if dimension == 0
        || site_orbits.iter().any(Vec::is_empty)
        || spatial_actions.len() != group.order()
    {
        return Err(GeometryError::new(
            "InvalidSiteOrbitData",
            "site orbits must be nonempty and actions must align with the group",
        ));
    }
    let context = &site_orbits[0][0][0].context();
    if site_orbits
        .iter()
        .flatten()
        .any(|site| site.len() != dimension || site.iter().any(|value| value.context() != *context))
        || spatial_actions.iter().any(|operation| {
            operation.rotation().rows() != dimension
                || operation.rotation().context() != *context
                || operation
                    .translation()
                    .iter()
                    .any(|value| value.context() != *context)
        })
    {
        return Err(GeometryError::new(
            "InvalidSiteOrbitData",
            "all coordinates and spatial actions must share one exact dimension and field",
        ));
    }
    Ok(dimension)
}

fn affine_image(operation: &SpatialOperation, vector: &[Element]) -> GeometryResult<Vec<Element>> {
    let mut result = Vec::with_capacity(vector.len());
    for row in 0..operation.rotation().rows() {
        let mut value = operation.translation()[row].clone();
        for (column, coordinate) in vector.iter().enumerate() {
            value = value.add(
                &operation
                    .rotation()
                    .entry(row, column)?
                    .multiply(coordinate)?,
            )?;
        }
        result.push(value);
    }
    Ok(result)
}

fn subtract_vectors(left: &[Element], right: &[Element]) -> GeometryResult<Vec<Element>> {
    if left.len() != right.len() {
        return Err(GeometryError::new(
            "DimensionMismatch",
            "coordinate vectors have different dimensions",
        ));
    }
    left.iter()
        .zip(right)
        .map(|(left, right)| left.subtract(right).map_err(Into::into))
        .collect()
}
