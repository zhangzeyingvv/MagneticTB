#![forbid(unsafe_code)]
#![allow(clippy::missing_errors_doc)]

mod basis;

pub use basis::{
    CompiledBasisAction, CompiledInducedBasisAction, FunctionBasisState, PolynomialExpression,
    PolynomialTerm, compile_catalog_basis_action, compile_catalog_induced_basis_action,
    compile_function_basis_action, compile_function_induced_basis_action,
    compile_spatial_spin_action, compile_spatial_spinor_matrix, resolve_catalog_basis_states,
    resolve_function_basis_states, transform_catalog_basis_states, transform_function_basis_states,
};

use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::error::Error;
use std::fmt::{Display, Formatter};
use std::str::FromStr;
use std::sync::Arc;

use cyclotomic_nullspace::{
    CyclotomicContext, Element, ExactError, ExactMatrix, Rational, make_context, multiply, rref,
};
use magnetictb_abstract_group::{GroupAlgebra, GroupError, generate_group_by};
use magnetictb_data::{
    ExactExpression as DataExactExpression, ExactMatrix as DataExactMatrix, MsgGroup,
    SymmetryOperation as DataSymmetryOperation,
};
use magnetictb_linear_algebra::{exact_matrix_equal, integer_vector};
use num_bigint::BigInt;

include!("symmetry/symmetry_algebra.rs");
include!("symmetry/spin_space_group_input.rs");
include!("symmetry/magnetic_group_input.rs");

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use cyclotomic_nullspace::{Element, ExactMatrix, Rational, make_context};
    use magnetictb_data::DataCatalog;

    use super::{
        SeitzOperation, SpatialOperation, SpinSpaceOperation, compile_msg_group,
        compile_seitz_group, generate_spin_space_group,
    };

    fn integer(context: &Arc<cyclotomic_nullspace::CyclotomicContext>, value: i64) -> Element {
        Element::from_polynomial(Arc::clone(context), &[Rational::from_i64(value)])
            .expect("integer")
    }

    #[test]
    fn ordered_c2_seitz_group_records_lattice_shifts_exactly() {
        let context = make_context(1, 8).expect("rational context");
        let identity = ExactMatrix::identity(&context, 1).expect("identity");
        let operations = vec![
            SeitzOperation::new(
                "c2.e",
                "E",
                SpatialOperation::new(identity.clone(), vec![integer(&context, 0)]).expect("E"),
                false,
            ),
            SeitzOperation::new(
                "c2.t",
                "t(1/2)'",
                SpatialOperation::new(
                    identity,
                    vec![
                        Element::from_polynomial(
                            Arc::clone(&context),
                            &[Rational::parse("1", "2").expect("half")],
                        )
                        .expect("half"),
                    ],
                )
                .expect("translation"),
                true,
            ),
        ];
        let compiled = compile_seitz_group(operations).expect("C2 Seitz group");
        assert_eq!(
            compiled.group().multiplication_table(),
            [vec![0, 1], vec![1, 0]]
        );
        assert_eq!(compiled.antiunitary_flags(), [false, true]);
        assert_eq!(
            compiled.product_records()[1][1].lattice_translation[0],
            integer(&context, 1)
        );
    }

    #[test]
    fn nonclosed_seitz_list_fails_explicitly() {
        let context = make_context(1, 8).expect("context");
        let operation = SeitzOperation::new(
            "half",
            "t(1/2)",
            SpatialOperation::new(
                ExactMatrix::identity(&context, 1).expect("identity"),
                vec![
                    Element::from_polynomial(
                        Arc::clone(&context),
                        &[Rational::parse("1", "2").expect("half")],
                    )
                    .expect("half"),
                ],
            )
            .expect("operation"),
            false,
        );
        let error = compile_seitz_group(vec![operation]).expect_err("no identity representative");
        assert_eq!(error.tag(), "SeitzClosureFailure");
    }

    #[test]
    fn spin_space_generators_expand_with_distinct_shared_spatial_actions() {
        let context = make_context(1, 8).expect("context");
        let identity = ExactMatrix::identity(&context, 3).expect("identity");
        let c4 = ExactMatrix::new(
            Arc::clone(&context),
            3,
            3,
            [0, -1, 0, 1, 0, 0, 0, 0, 1]
                .into_iter()
                .map(|value| integer(&context, value))
                .collect(),
        )
        .expect("C4");
        let zero = vec![integer(&context, 0); 3];
        let generators = vec![
            SpinSpaceOperation::new(
                SpatialOperation::new(identity.clone(), zero.clone()).expect("space"),
                c4,
                false,
            )
            .expect("C4 spin"),
            SpinSpaceOperation::new(
                SpatialOperation::new(
                    identity.clone(),
                    vec![
                        Element::from_polynomial(
                            Arc::clone(&context),
                            &[Rational::parse("1", "2").expect("half")],
                        )
                        .expect("half"),
                        integer(&context, 0),
                        integer(&context, 0),
                    ],
                )
                .expect("half translation"),
                identity,
                true,
            )
            .expect("half T"),
        ];
        let compiled = generate_spin_space_group(&generators).expect("SSG closure");
        assert_eq!(compiled.elements().len(), 8);
        assert_eq!(compiled.seitz_group().group().order(), 8);
        assert_eq!(
            compiled
                .elements()
                .iter()
                .filter(|element| element.antiunitary())
                .count(),
            4
        );
    }

    #[test]
    fn frozen_type_iv_msg_compiles_with_order_and_lattice_shift() {
        let catalog = DataCatalog::load_embedded().expect("frozen data catalog");
        let source = catalog.msg("msg-0003").expect("stable type-IV MSG");
        let compiled = compile_msg_group(source).expect("exact MSG quotient group");
        assert_eq!(
            compiled.group().multiplication_table(),
            [vec![0, 1], vec![1, 0]]
        );
        assert_eq!(compiled.antiunitary_flags(), [false, true]);
        assert_eq!(
            compiled
                .operations()
                .iter()
                .map(SeitzOperation::stable_id)
                .collect::<Vec<_>>(),
            ["msg-0003.op-000", "msg-0003.op-001"]
        );
        assert_eq!(
            compiled.product_records()[1][1].lattice_translation[2],
            integer(compiled.operations()[0].spatial().rotation().context(), 1)
        );
    }
}
