pub fn compile_direct_product(
    group: &GroupAlgebra,
    actions: &[GroupAction],
    local_representations: &[Vec<ExactMatrix>],
    antiunitary_flags: &[bool],
) -> RepresentationResult<CompiledRepresentation> {
    if actions.is_empty()
        || actions.len() != local_representations.len()
        || antiunitary_flags.len() != group.order()
    {
        return Err(RepresentationError::new(
            "InvalidDirectProductData",
            "one action and local representation are required for every orbit",
        ));
    }
    let mut dimensions = Vec::with_capacity(actions.len());
    let mut blocks = Vec::with_capacity(actions.len());
    for (action, local) in actions.iter().zip(local_representations) {
        if action.group() != group || local.len() != group.order() || local.is_empty() {
            return Err(RepresentationError::new(
                "InvalidDirectProductData",
                "local representation or action does not align with the group",
            ));
        }
        verify_corepresentation(group, local, antiunitary_flags)?;
        dimensions.push(local[0].rows());
        blocks.push(
            local
                .iter()
                .map(|matrix| vec![matrix.clone(); action.object_count()])
                .collect(),
        );
    }
    assemble(
        "DirectProduct",
        group,
        actions,
        blocks,
        dimensions,
        antiunitary_flags,
        true,
    )
}
/// Assembles already-compiled orbit-local blocks into the exact full block-monomial
/// representation used by physical site actions.
pub fn compile_block_monomial(
    method: &str,
    group: &GroupAlgebra,
    actions: &[GroupAction],
    local_blocks: Vec<Vec<Vec<ExactMatrix>>>,
    local_dimensions: Vec<usize>,
    antiunitary_flags: &[bool],
) -> RepresentationResult<CompiledRepresentation> {
    assemble(
        method,
        group,
        actions,
        local_blocks,
        local_dimensions,
        antiunitary_flags,
        false,
    )
}
