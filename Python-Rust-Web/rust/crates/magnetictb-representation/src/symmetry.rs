use magnetictb_abstract_group::{GroupAction, GroupAlgebra};

use crate::{RepresentationError, RepresentationResult};

/// Builds the Cayley table of an explicitly ordered list of zero-based
/// permutations.  This is table compilation, not group enumeration.
pub fn ordered_permutation_group(
    permutations: &[Vec<usize>],
) -> RepresentationResult<GroupAlgebra> {
    let degree = permutations.first().map_or(0, Vec::len);
    if degree == 0 || permutations.is_empty() {
        return Err(RepresentationError::new(
            "InvalidPermutationData",
            "the ordered permutation list must be nonempty",
        ));
    }
    let expected: Vec<usize> = (0..degree).collect();
    for permutation in permutations {
        let mut sorted = permutation.clone();
        sorted.sort_unstable();
        if sorted != expected {
            return Err(RepresentationError::new(
                "InvalidPermutationData",
                "each ordered element must be a permutation of the same objects",
            ));
        }
    }

    let mut table = Vec::with_capacity(permutations.len());
    for left in permutations {
        let mut row = Vec::with_capacity(permutations.len());
        for right in permutations {
            let product: Vec<usize> = right.iter().map(|&source| left[source]).collect();
            let index = permutations
                .iter()
                .position(|candidate| candidate == &product)
                .ok_or_else(|| {
                    RepresentationError::new(
                        "PermutationClosureFailure",
                        "an ordered permutation product is absent from the table",
                    )
                })?;
            row.push(index);
        }
        table.push(row);
    }
    GroupAlgebra::new(table).map_err(Into::into)
}

pub fn permutation_action(
    group: &GroupAlgebra,
    permutations: &[Vec<usize>],
) -> RepresentationResult<GroupAction> {
    GroupAction::compile(group, permutations.to_vec()).map_err(Into::into)
}
