/// Expands discrete spin-space generators, retaining distinct spin elements
/// that share a spatial action, then compiles their exact ordered group table.
pub fn generate_spin_space_group(
    generators: &[SpinSpaceOperation],
) -> SymmetryResult<CompiledSpinSpaceGroup> {
    let context = generators
        .first()
        .ok_or_else(|| SymmetryError::new("EmptySpinSpaceGenerators", "generator list is empty"))?
        .spatial()
        .rotation()
        .context();
    let zero = Element::zero(context)?;
    let identity = SpinSpaceOperation::new(
        SpatialOperation::new(
            ExactMatrix::identity(context, 3)?,
            vec![zero.clone(), zero.clone(), zero],
        )?,
        ExactMatrix::identity(context, 3)?,
        false,
    )?;
    let mut elements = generate_group_by(
        generators,
        &identity,
        SpinSpaceOperation::product,
        |left, right| left.equivalent_mod_lattice(right).unwrap_or(false),
    )?;
    elements.sort_by(|left, right| {
        let left_identity = left.equivalent_mod_lattice(&identity).unwrap_or(false);
        let right_identity = right.equivalent_mod_lattice(&identity).unwrap_or(false);
        match (left_identity, right_identity) {
            (true, false) => Ordering::Less,
            (false, true) => Ordering::Greater,
            _ => compare_spin_space(left, right),
        }
    });
    compile_spin_space_group(elements)
}

/// Compiles an already complete ordered discrete spin-space group without
/// changing its caller-supplied element order.
pub fn compile_spin_space_group(
    elements: Vec<SpinSpaceOperation>,
) -> SymmetryResult<CompiledSpinSpaceGroup> {
    if elements.is_empty() {
        return Err(SymmetryError::new(
            "EmptySpinSpaceGroup",
            "at least one complete ordered spin-space element is required",
        ));
    }
    let mut table = vec![vec![0; elements.len()]; elements.len()];
    for left in 0..elements.len() {
        for right in 0..elements.len() {
            let product = elements[left].product(&elements[right])?;
            let matches = elements
                .iter()
                .enumerate()
                .filter_map(|(index, candidate)| {
                    product
                        .equivalent_mod_lattice(candidate)
                        .ok()
                        .is_some_and(|equal| equal)
                        .then_some(index)
                })
                .collect::<Vec<_>>();
            if matches.len() != 1 {
                return Err(SymmetryError::new(
                    "SpinSpaceClosureFailure",
                    format!(
                        "spin-space product ({left}, {right}) has {} representatives",
                        matches.len()
                    ),
                ));
            }
            table[left][right] = matches[0];
        }
    }
    let seitz_operations = elements
        .iter()
        .enumerate()
        .map(|(index, element)| {
            SeitzOperation::new(
                format!("spin-space-{index}"),
                format!("SpinSpace{}", index + 1),
                element.spatial.clone(),
                element.antiunitary,
            )
        })
        .collect();
    let seitz_group = compile_seitz_group_with_table(seitz_operations, table)?;
    Ok(CompiledSpinSpaceGroup {
        elements,
        seitz_group,
    })
}
