#![forbid(unsafe_code)]
#![allow(clippy::missing_errors_doc)]

mod error;
mod physical;
mod representation;
mod symmetry;

pub use error::{RepresentationError, RepresentationResult};
pub use physical::{SiteActionData, SiteActionRecord};
pub use representation::{
    CompiledInducedRepresentation, CompiledRepresentation, InducedOrbitData, InducedOrbitSpec,
    automatic_coset_representatives, compile_block_monomial, compile_direct_product,
    compile_induced, verify_corepresentation,
};
pub use symmetry::{ordered_permutation_group, permutation_action};

#[cfg(test)]
mod tests {
    use cyclotomic_nullspace::{Element, ExactMatrix, make_context, multiply};
    use magnetictb_abstract_group::{GroupAction, GroupAlgebra};

    use super::{
        InducedOrbitSpec, SiteActionData, automatic_coset_representatives, compile_direct_product,
        compile_induced, ordered_permutation_group, permutation_action, verify_corepresentation,
    };

    fn s3() -> (Vec<Vec<usize>>, Vec<bool>) {
        (
            vec![
                vec![0, 1, 2],
                vec![1, 2, 0],
                vec![2, 0, 1],
                vec![0, 2, 1],
                vec![2, 1, 0],
                vec![1, 0, 2],
            ],
            vec![false, false, false, true, true, true],
        )
    }

    #[test]
    fn s3_direct_and_induced_site_action_close_exactly() {
        let (permutations, flags) = s3();
        let group = ordered_permutation_group(&permutations).expect("ordered S3");
        assert_eq!(
            group.multiplication_table(),
            [
                vec![0, 1, 2, 3, 4, 5],
                vec![1, 2, 0, 5, 3, 4],
                vec![2, 0, 1, 4, 5, 3],
                vec![3, 4, 5, 0, 1, 2],
                vec![4, 5, 3, 2, 0, 1],
                vec![5, 3, 4, 1, 2, 0],
            ]
        );
        let first = permutation_action(&group, &permutations).expect("first orbit action");
        let second = GroupAction::compile(&group, vec![vec![0]; 6]).expect("fixed orbit action");
        let actions = vec![first, second];
        let context = make_context(1, 8).expect("rational context");
        let one = ExactMatrix::identity(&context, 1).expect("one");
        let local = vec![vec![one.clone(); 6], vec![one.clone(); 6]];

        let direct = compile_direct_product(&group, &actions, &local, &flags)
            .expect("direct product representation");
        assert_eq!(direct.dimension(), 4);
        assert_eq!(direct.local_dimensions(), [1, 1]);

        let mut specifications = Vec::new();
        for action in &actions {
            let subgroup = action.stabilizer(0).expect("site stabilizer");
            let representatives =
                automatic_coset_representatives(action, &flags, 0).expect("representatives");
            specifications.push(InducedOrbitSpec {
                site_symmetry_matrices: vec![one.clone(); subgroup.len()],
                subgroup_indices: subgroup,
                coset_representatives: representatives,
            });
        }
        let induced = compile_induced(&group, &actions, &specifications, &flags)
            .expect("induced representation");
        assert_eq!(
            induced.representation.representation_matrices(),
            direct.representation_matrices()
        );

        let translations = actions
            .iter()
            .map(|action| vec![vec![vec![0_i64; 3]; action.object_count()]; group.order()])
            .collect();
        let site_action = SiteActionData::new(
            "Induced",
            &group,
            actions,
            translations,
            induced.representation.local_blocks().to_vec(),
            induced.representation.local_dimensions().to_vec(),
            flags,
        )
        .expect("site action data");
        let record = site_action.site_action(0, 1, 4).expect("site action");
        assert_eq!(record.image_site_index, 1);
        assert!(record.antiunitary);
        assert_eq!(record.cell_translation, [0, 0, 0]);
    }

    #[test]
    fn invalid_ordered_permutations_fail_explicitly() {
        let error =
            ordered_permutation_group(&[vec![0, 0]]).expect_err("non-permutation must fail");
        assert_eq!(error.tag(), "InvalidPermutationData");
    }

    #[test]
    fn explicit_induced_accepts_stable_double_valued_local_data() {
        // Ordered operations: E, tau, C2z, C2z tau.  The site action sees the
        // half translation, while C2z fixes both sites.
        let group = GroupAlgebra::new(vec![
            vec![0, 1, 2, 3],
            vec![1, 0, 3, 2],
            vec![2, 3, 0, 1],
            vec![3, 2, 1, 0],
        ])
        .expect("ordered V4 group");
        let action =
            GroupAction::compile(&group, vec![vec![0, 1], vec![1, 0], vec![0, 1], vec![1, 0]])
                .expect("two-site action");
        let context = make_context(4, 8).expect("Gaussian exact field");
        let zero = Element::zero(&context).expect("zero");
        let one = Element::one(&context).expect("one");
        let imaginary = Element::root_of_unity(&context, 4, 1).expect("I");
        let minus_imaginary = imaginary.negate().expect("-I");
        let local_identity = ExactMatrix::identity(&context, 2).expect("I2");
        let local_c2 = ExactMatrix::new(
            context.clone(),
            2,
            2,
            vec![
                minus_imaginary.clone(),
                zero.clone(),
                zero.clone(),
                imaginary.clone(),
            ],
        )
        .expect("double-valued C2z");
        let compiled = compile_induced(
            &group,
            std::slice::from_ref(&action),
            &[InducedOrbitSpec {
                subgroup_indices: vec![0, 2],
                site_symmetry_matrices: vec![local_identity, local_c2],
                coset_representatives: vec![0, 1],
            }],
            &[false; 4],
        )
        .expect("stable double-valued induced representation");

        assert_eq!(compiled.representation.dimension(), 4);
        let c2 = &compiled.representation.representation_matrices()[2];
        let c2_squared = multiply(c2, c2).expect("C2 square");
        let negative_identity = ExactMatrix::new(
            context.clone(),
            4,
            4,
            (0..16)
                .map(|offset| {
                    if offset / 4 == offset % 4 {
                        one.negate().expect("-1")
                    } else {
                        zero.clone()
                    }
                })
                .collect(),
        )
        .expect("-I4");
        assert_eq!(c2_squared, negative_identity);

        // This is intentionally projective relative to the spatial V4 table.
        // DirectProduct and callers that request an ordinary corepresentation
        // continue to use this strict verifier.
        let error = verify_corepresentation(
            &group,
            compiled.representation.representation_matrices(),
            &[false; 4],
        )
        .expect_err("ordinary representation law must distinguish the projective case");
        assert_eq!(error.tag(), "RepresentationLawViolation");
    }
}
