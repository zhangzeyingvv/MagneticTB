#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CompiledRepresentation {
    method: String,
    group: GroupAlgebra,
    antiunitary_flags: Vec<bool>,
    local_dimensions: Vec<usize>,
    local_blocks: Vec<Vec<Vec<ExactMatrix>>>,
    orbit_representation_matrices: Vec<Vec<ExactMatrix>>,
    representation_matrices: Vec<ExactMatrix>,
    block_dimensions: Vec<usize>,
    dimension: usize,
}
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct InducedOrbitSpec {
    pub subgroup_indices: Vec<usize>,
    pub site_symmetry_matrices: Vec<ExactMatrix>,
    pub coset_representatives: Vec<usize>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct InducedOrbitData {
    pub subgroup_indices: Vec<usize>,
    pub coset_representatives: Vec<usize>,
    pub schreier_records: Vec<Vec<(usize, usize)>>,
    pub image_site_indices: Vec<Vec<usize>>,
    pub local_blocks: Vec<Vec<ExactMatrix>>,
    pub local_dimension: usize,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CompiledInducedRepresentation {
    pub representation: CompiledRepresentation,
    pub site_symmetry_data: Vec<InducedOrbitData>,
}

impl CompiledRepresentation {
    #[must_use]
    pub fn method(&self) -> &str {
        &self.method
    }

    #[must_use]
    pub const fn group(&self) -> &GroupAlgebra {
        &self.group
    }

    #[must_use]
    pub fn antiunitary_flags(&self) -> &[bool] {
        &self.antiunitary_flags
    }

    #[must_use]
    pub fn local_dimensions(&self) -> &[usize] {
        &self.local_dimensions
    }

    #[must_use]
    pub fn local_blocks(&self) -> &[Vec<Vec<ExactMatrix>>] {
        &self.local_blocks
    }

    #[must_use]
    pub fn orbit_representation_matrices(&self) -> &[Vec<ExactMatrix>] {
        &self.orbit_representation_matrices
    }

    #[must_use]
    pub fn representation_matrices(&self) -> &[ExactMatrix] {
        &self.representation_matrices
    }

    #[must_use]
    pub fn block_dimensions(&self) -> &[usize] {
        &self.block_dimensions
    }

    #[must_use]
    pub const fn dimension(&self) -> usize {
        self.dimension
    }
}

pub fn verify_corepresentation(
    group: &GroupAlgebra,
    matrices: &[ExactMatrix],
    antiunitary_flags: &[bool],
) -> RepresentationResult<()> {
    if matrices.len() != group.order()
        || antiunitary_flags.len() != group.order()
        || matrices.is_empty()
    {
        return Err(RepresentationError::new(
            "RepresentationDimensionMismatch",
            "group order, matrices, and antiunitary flags must align",
        ));
    }
    let first = &matrices[0];
    if first.rows() == 0
        || first.rows() != first.columns()
        || matrices.iter().any(|matrix| {
            matrix.context() != first.context()
                || matrix.rows() != first.rows()
                || matrix.columns() != first.columns()
        })
    {
        return Err(RepresentationError::new(
            "NonunitaryRepresentation",
            "representation matrices must be same-size exact unitaries",
        ));
    }
    for matrix in matrices {
        if !exact_unitary_matrix(matrix)? {
            return Err(RepresentationError::new(
                "NonunitaryRepresentation",
                "representation matrices must be exact unitaries",
            ));
        }
    }
    for left in 0..group.order() {
        for right in 0..group.order() {
            let product = group.multiplication_table()[left][right];
            if antiunitary_flags[product] != (antiunitary_flags[left] ^ antiunitary_flags[right]) {
                return Err(RepresentationError::new(
                    "AntiunitaryParityViolation",
                    "antiunitary flags are not a group homomorphism",
                ));
            }
            let right_matrix = if antiunitary_flags[left] {
                conjugate(&matrices[right])?
            } else {
                matrices[right].clone()
            };
            let composed = multiply(&matrices[left], &right_matrix)?;
            if !exact_matrix_equal(&composed, &matrices[product]) {
                return Err(RepresentationError::new(
                    "RepresentationLawViolation",
                    format!("corepresentation law fails for ({left}, {right})"),
                ));
            }
        }
    }
    Ok(())
}
