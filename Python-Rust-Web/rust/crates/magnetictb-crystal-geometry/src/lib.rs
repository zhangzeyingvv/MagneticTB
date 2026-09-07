#![forbid(unsafe_code)]
#![allow(clippy::missing_errors_doc)]

mod bond;

pub use bond::{
    BondSearchOptions, PeriodicBondRecord, PeriodicBondSearchResult, PeriodicBondShell,
    QuadraticReal, find_periodic_bond_shells, periodic_bond_shells_in_box,
};

use std::collections::BTreeMap;
use std::error::Error;
use std::fmt::{Display, Formatter};
use std::sync::Arc;

use cyclotomic_nullspace::{Element, ExactError, make_context};
use magnetictb_abstract_group::{GroupAction, GroupAlgebra, GroupError};
use magnetictb_data::DataCatalog;
use magnetictb_linear_algebra::integer_vector;
pub use magnetictb_symmetry::SpatialOperation;
use magnetictb_symmetry::{
    CompiledSeitzGroup, SymmetryError, compile_msg_group_in_context, evaluate_data_expression,
    fractional_part,
};

include!("geometry/site_permutation_compiler.rs");
include!("geometry/site_orbit_compiler.rs");

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;
    use std::sync::Arc;

    use cyclotomic_nullspace::{Element, ExactMatrix, Rational, make_context};
    use magnetictb_data::DataCatalog;
    use magnetictb_representation::ordered_permutation_group;

    use super::{SpatialOperation, compile_msg_wyckoff_sites, compile_site_permutations};

    fn rational(
        context: &Arc<cyclotomic_nullspace::CyclotomicContext>,
        numerator: i64,
        denominator: i64,
    ) -> Element {
        Element::from_polynomial(
            Arc::clone(context),
            &[
                Rational::parse(&numerator.to_string(), &denominator.to_string())
                    .expect("rational"),
            ],
        )
        .expect("element")
    }

    #[test]
    fn exact_s3_site_images_match_the_stable_representation_fixture() {
        let permutations = vec![
            vec![0, 1, 2],
            vec![1, 2, 0],
            vec![2, 0, 1],
            vec![0, 2, 1],
            vec![2, 1, 0],
            vec![1, 0, 2],
        ];
        let group = ordered_permutation_group(&permutations).expect("S3");
        let context = make_context(1, 8).expect("rational context");
        let zero = rational(&context, 0, 1);
        let quarter = rational(&context, 1, 4);
        let first_orbit = vec![
            vec![quarter.clone(), zero.clone(), zero.clone()],
            vec![zero.clone(), quarter.clone(), zero.clone()],
            vec![zero.clone(), zero.clone(), quarter],
        ];
        let second_orbit = vec![vec![zero.clone(), zero.clone(), zero.clone()]];
        let operations = permutations
            .iter()
            .map(|permutation| {
                let mut entries = vec![zero.clone(); 9];
                for (source, &target) in permutation.iter().enumerate() {
                    entries[target * 3 + source] = rational(&context, 1, 1);
                }
                SpatialOperation::new(
                    ExactMatrix::new(Arc::clone(&context), 3, 3, entries).expect("rotation"),
                    vec![zero.clone(); 3],
                )
                .expect("operation")
            })
            .collect();
        let data = compile_site_permutations(&group, vec![first_orbit, second_orbit], operations)
            .expect("site permutations");
        assert_eq!(
            data.image_site_indices(),
            [
                permutations,
                vec![vec![0], vec![0], vec![0], vec![0], vec![0], vec![0]],
            ]
        );
        assert_eq!(
            data.site_symmetry_operations(0, 0).expect("stabilizer"),
            [0, 3]
        );
        assert!(
            data.cell_translations()
                .iter()
                .flatten()
                .flatten()
                .flatten()
                .all(Element::is_zero)
        );
    }

    #[test]
    fn invalid_spatial_dimension_fails_explicitly() {
        let context = make_context(1, 8).expect("context");
        let rotation = ExactMatrix::identity(&context, 2).expect("rotation");
        let error = SpatialOperation::new(rotation, vec![Element::zero(&context).expect("zero")])
            .expect_err("dimension mismatch must fail");
        assert_eq!(error.tag(), "InvalidSpatialOperation");
    }

    #[test]
    fn frozen_type_iv_wyckoff_orbit_has_exact_site_action() {
        let catalog = DataCatalog::load_embedded().expect("frozen Data");
        let context = make_context(1, 8).expect("rational context");
        let bindings = BTreeMap::from([
            ("MagneticTB`x".to_owned(), rational(&context, 1, 7)),
            ("MagneticTB`y".to_owned(), rational(&context, 2, 11)),
            ("MagneticTB`z".to_owned(), rational(&context, 3, 13)),
        ]);
        let compiled = compile_msg_wyckoff_sites(&catalog, "msg-0003", 2, &bindings)
            .expect("Data to exact site action");
        assert_eq!(compiled.letter(), "a");
        assert_eq!(
            compiled.sites().image_site_indices(),
            [vec![vec![0, 1], vec![1, 0]]]
        );
        assert_eq!(
            compiled.sites().cell_translations()[0][1],
            [
                vec![
                    rational(&context, 0, 1),
                    rational(&context, 0, 1),
                    rational(&context, 0, 1),
                ],
                vec![
                    rational(&context, 0, 1),
                    rational(&context, 0, 1),
                    rational(&context, 1, 1),
                ],
            ]
        );
    }
}
