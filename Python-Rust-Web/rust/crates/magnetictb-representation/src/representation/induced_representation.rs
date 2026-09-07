pub fn automatic_coset_representatives(
    action: &GroupAction,
    antiunitary_flags: &[bool],
    reference: usize,
) -> RepresentationResult<Vec<usize>> {
    if antiunitary_flags.len() != action.group().order() || reference >= action.object_count() {
        return Err(RepresentationError::new(
            "InvalidCosetRepresentativeData",
            "site action, flags, and reference index do not align",
        ));
    }
    (0..action.object_count())
        .map(|target| {
            let transporters = action.transporters(reference, target)?;
            transporters
                .iter()
                .find(|&&operation| !antiunitary_flags[operation])
                .or_else(|| transporters.first())
                .copied()
                .ok_or_else(|| {
                    RepresentationError::new(
                        "MissingCosetRepresentative",
                        format!("no operation maps reference {reference} to target {target}"),
                    )
                })
        })
        .collect()
}

pub fn compile_induced(
    group: &GroupAlgebra,
    actions: &[GroupAction],
    specifications: &[InducedOrbitSpec],
    antiunitary_flags: &[bool],
) -> RepresentationResult<CompiledInducedRepresentation> {
    if actions.is_empty() || actions.len() != specifications.len() {
        return Err(RepresentationError::new(
            "InvalidInducedData",
            "one induced specification is required for every orbit",
        ));
    }
    let mut orbit_data = Vec::with_capacity(actions.len());
    for (action, specification) in actions.iter().zip(specifications) {
        orbit_data.push(compile_induced_orbit(
            group,
            action,
            specification,
            antiunitary_flags,
        )?);
    }
    let blocks = orbit_data
        .iter()
        .map(|orbit| orbit.local_blocks.clone())
        .collect();
    let dimensions = orbit_data
        .iter()
        .map(|orbit| orbit.local_dimension)
        .collect();
    let representation = assemble(
        "Induced",
        group,
        actions,
        blocks,
        dimensions,
        antiunitary_flags,
        // Stable Mathematica validates the ordered site-symmetry matrices as exact
        // unitaries and validates the Schreier induction data, but deliberately does
        // not impose the ordinary corepresentation law on the assembled matrices.
        // That distinction is required for double-valued/projective local data, for
        // example D(C2)^2 = -I while the spatial operation satisfies C2^2 = E.
        false,
    )?;
    Ok(CompiledInducedRepresentation {
        representation,
        site_symmetry_data: orbit_data,
    })
}

#[allow(clippy::too_many_lines)]
fn compile_induced_orbit(
    group: &GroupAlgebra,
    action: &GroupAction,
    specification: &InducedOrbitSpec,
    antiunitary_flags: &[bool],
) -> RepresentationResult<InducedOrbitData> {
    if action.group() != group
        || antiunitary_flags.len() != group.order()
        || specification.coset_representatives.len() != action.object_count()
    {
        return Err(RepresentationError::new(
            "InvalidInducedData",
            "group, action, flags, and coset representatives do not align",
        ));
    }
    let mut subgroup = Vec::new();
    for &operation in &specification.subgroup_indices {
        if !subgroup.contains(&operation) {
            subgroup.push(operation);
        }
    }
    group.find_generator_indices(Some(&subgroup))?;
    if subgroup.is_empty()
        || specification.site_symmetry_matrices.len() != subgroup.len()
        || specification
            .coset_representatives
            .iter()
            .any(|&index| index >= group.order())
    {
        return Err(RepresentationError::new(
            "InvalidInducedSubgroup",
            "subgroup matrices or coset representatives are invalid",
        ));
    }
    let first = &specification.site_symmetry_matrices[0];
    if first.rows() == 0
        || first.rows() != first.columns()
        || specification.site_symmetry_matrices.iter().any(|matrix| {
            matrix.context() != first.context()
                || matrix.rows() != first.rows()
                || matrix.columns() != first.columns()
        })
    {
        return Err(RepresentationError::new(
            "NonunitaryLocalRepresentation",
            "site-symmetry matrices must be same-size exact unitaries",
        ));
    }
    for matrix in &specification.site_symmetry_matrices {
        if !exact_unitary_matrix(matrix)? {
            return Err(RepresentationError::new(
                "NonunitaryLocalRepresentation",
                "site-symmetry matrices must be exact unitaries",
            ));
        }
    }

    let mut records = Vec::with_capacity(group.order());
    let mut images = Vec::with_capacity(group.order());
    let mut blocks = Vec::with_capacity(group.order());
    for operation in 0..group.order() {
        let mut operation_records = Vec::with_capacity(action.object_count());
        let mut operation_images = Vec::with_capacity(action.object_count());
        let mut operation_blocks = Vec::with_capacity(action.object_count());
        for source in 0..action.object_count() {
            let (target, subgroup_element) = group.schreier_decomposition(
                &subgroup,
                &specification.coset_representatives,
                operation,
                source,
            )?;
            let position = subgroup
                .iter()
                .position(|&element| element == subgroup_element)
                .ok_or_else(|| {
                    RepresentationError::new(
                        "InvalidSchreierDecomposition",
                        "Schreier subgroup element is absent from local data",
                    )
                })?;
            let matrix = &specification.site_symmetry_matrices[position];
            operation_blocks.push(
                if antiunitary_flags[specification.coset_representatives[target]] {
                    conjugate(matrix)?
                } else {
                    matrix.clone()
                },
            );
            operation_images.push(target);
            operation_records.push((target, subgroup_element));
        }
        if operation_images != action.action_table()[operation] {
            return Err(RepresentationError::new(
                "InducedActionMismatch",
                "Schreier images differ from the supplied ordered site action",
            ));
        }
        records.push(operation_records);
        images.push(operation_images);
        blocks.push(operation_blocks);
    }
    Ok(InducedOrbitData {
        subgroup_indices: subgroup,
        coset_representatives: specification.coset_representatives.clone(),
        schreier_records: records,
        image_site_indices: images,
        local_blocks: blocks,
        local_dimension: first.rows(),
    })
}

pub(crate) fn assemble(
    method: &str,
    group: &GroupAlgebra,
    actions: &[GroupAction],
    local_blocks: Vec<Vec<Vec<ExactMatrix>>>,
    local_dimensions: Vec<usize>,
    antiunitary_flags: &[bool],
    verify_law: bool,
) -> RepresentationResult<CompiledRepresentation> {
    if actions.len() != local_blocks.len()
        || actions.len() != local_dimensions.len()
        || antiunitary_flags.len() != group.order()
    {
        return Err(RepresentationError::new(
            "InvalidBlockMonomialData",
            "orbit actions, blocks, dimensions, and flags do not align",
        ));
    }
    let mut orbit_matrices = Vec::with_capacity(actions.len());
    for ((action, blocks), &dimension) in actions.iter().zip(&local_blocks).zip(&local_dimensions) {
        if dimension == 0 || action.group() != group || blocks.len() != group.order() {
            return Err(RepresentationError::new(
                "InvalidBlockMonomialData",
                "an orbit action or local dimension is invalid",
            ));
        }
        let mut matrices = Vec::with_capacity(group.order());
        for (operation, operation_blocks) in blocks.iter().enumerate() {
            if operation_blocks.len() != action.object_count() {
                return Err(RepresentationError::new(
                    "InvalidBlockMonomialData",
                    "local blocks must be operation-by-source aligned",
                ));
            }
            matrices.push(assemble_orbit_matrix(
                &action.action_table()[operation],
                operation_blocks,
                dimension,
            )?);
        }
        orbit_matrices.push(matrices);
    }

    let mut representation_matrices = Vec::with_capacity(group.order());
    for operation in 0..group.order() {
        let diagonal: Vec<&ExactMatrix> = orbit_matrices
            .iter()
            .map(|orbit| &orbit[operation])
            .collect();
        representation_matrices.push(block_diagonal(&diagonal)?);
    }
    if verify_law {
        verify_corepresentation(group, &representation_matrices, antiunitary_flags)?;
    }

    let block_dimensions: Vec<usize> = actions
        .iter()
        .zip(&local_dimensions)
        .flat_map(|(action, &dimension)| vec![dimension; action.object_count()])
        .collect();
    let dimension = block_dimensions.iter().sum();
    Ok(CompiledRepresentation {
        method: method.to_owned(),
        group: group.clone(),
        antiunitary_flags: antiunitary_flags.to_vec(),
        local_dimensions,
        local_blocks,
        orbit_representation_matrices: orbit_matrices,
        representation_matrices,
        block_dimensions,
        dimension,
    })
}

fn assemble_orbit_matrix(
    images: &[usize],
    blocks: &[ExactMatrix],
    local_dimension: usize,
) -> RepresentationResult<ExactMatrix> {
    let Some(first) = blocks.first() else {
        return Err(RepresentationError::new(
            "InvalidBlockMonomialData",
            "an orbit must contain at least one local block",
        ));
    };
    if images.len() != blocks.len()
        || blocks.iter().any(|block| {
            block.context() != first.context()
                || block.rows() != local_dimension
                || block.columns() != local_dimension
        })
    {
        return Err(RepresentationError::new(
            "InvalidBlockMonomialData",
            "local blocks must be aligned exact unitaries",
        ));
    }
    for block in blocks {
        if !exact_unitary_matrix(block)? {
            return Err(RepresentationError::new(
                "InvalidBlockMonomialData",
                "local blocks must be exact unitaries",
            ));
        }
    }
    let site_count = images.len();
    let total_dimension = site_count.checked_mul(local_dimension).ok_or_else(|| {
        RepresentationError::new("DimensionOverflow", "orbit dimension overflows")
    })?;
    let mut entries = vec![Element::zero(first.context())?; total_dimension * total_dimension];
    for (source, (&target, block)) in images.iter().zip(blocks).enumerate() {
        if target >= site_count {
            return Err(RepresentationError::new(
                "SiteIndexOutOfRange",
                "an image site is outside its orbit",
            ));
        }
        for row in 0..local_dimension {
            for column in 0..local_dimension {
                let global_row = target * local_dimension + row;
                let global_column = source * local_dimension + column;
                entries[global_row * total_dimension + global_column] =
                    block.entry(row, column)?.clone();
            }
        }
    }
    ExactMatrix::new(
        Arc::clone(first.context()),
        total_dimension,
        total_dimension,
        entries,
    )
    .map_err(Into::into)
}

fn block_diagonal(blocks: &[&ExactMatrix]) -> RepresentationResult<ExactMatrix> {
    let Some(first) = blocks.first() else {
        return Err(RepresentationError::new(
            "InvalidBlockMonomialData",
            "at least one orbit block is required",
        ));
    };
    let dimension: usize = blocks.iter().map(|block| block.rows()).sum();
    if blocks.iter().any(|block| {
        block.context() != first.context() || block.rows() == 0 || block.rows() != block.columns()
    }) {
        return Err(RepresentationError::new(
            "InvalidBlockMonomialData",
            "full representation blocks must be square in one exact field",
        ));
    }
    let mut entries = vec![Element::zero(first.context())?; dimension * dimension];
    let mut offset = 0;
    for block in blocks {
        for row in 0..block.rows() {
            for column in 0..block.columns() {
                entries[(offset + row) * dimension + offset + column] =
                    block.entry(row, column)?.clone();
            }
        }
        offset += block.rows();
    }
    ExactMatrix::new(Arc::clone(first.context()), dimension, dimension, entries).map_err(Into::into)
}
