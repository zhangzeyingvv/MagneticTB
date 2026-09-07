use std::collections::BTreeMap;
use std::sync::Arc;

use cyclotomic_nullspace::{CyclotomicContext, Element, ExactMatrix, Rational, multiply, rref};
use magnetictb_linear_algebra::exact_matrix_equal;

use crate::{SpatialOperation, SymmetryError, SymmetryResult};

type Exponents = [u8; 3];
type Polynomial = BTreeMap<Exponents, Element>;
type BasisFunction = Vec<Polynomial>;

include!("basis/function_basis_representation.rs");
include!("basis/basis_catalog.rs");
include!("basis/spin_representation.rs");

#[cfg(test)]
mod tests {
    use super::*;
    use cyclotomic_nullspace::make_context;

    fn matrix(context: &Arc<CyclotomicContext>, rows: [[i64; 3]; 3]) -> ExactMatrix {
        ExactMatrix::new(
            Arc::clone(context),
            3,
            3,
            rows.into_iter()
                .flatten()
                .map(|value| scalar(context, value).expect("integer"))
                .collect(),
        )
        .expect("matrix")
    }

    #[test]
    fn p_orbit_and_c4_spinor_are_exact() {
        let context = make_context(24, 128).expect("context");
        let identity = matrix(&context, [[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
        let c4 = matrix(&context, [[0, -1, 0], [1, 0, 0], [0, 0, 1]]);
        let operations = [
            SpatialOperation::new(identity.clone(), vec![scalar(&context, 0).unwrap(); 3]).unwrap(),
            SpatialOperation::new(identity.clone(), vec![scalar(&context, 0).unwrap(); 3]).unwrap(),
        ];
        let orbital = compile_catalog_basis_action(
            &["px".to_owned(), "py".to_owned()],
            &[
                SpatialOperation::new(identity.clone(), vec![scalar(&context, 0).unwrap(); 3])
                    .unwrap(),
                SpatialOperation::new(c4.clone(), vec![scalar(&context, 0).unwrap(); 3]).unwrap(),
            ],
            None,
            &[false, false],
            &identity,
        )
        .expect("p action");
        let p_c4 = ExactMatrix::new(
            Arc::clone(&context),
            2,
            2,
            vec![
                scalar(&context, 0).unwrap(),
                scalar(&context, -1).unwrap(),
                scalar(&context, 1).unwrap(),
                scalar(&context, 0).unwrap(),
            ],
        )
        .unwrap();
        assert_eq!(orbital.local_matrices()[1], p_c4);

        let spin = compile_catalog_basis_action(
            &["sup".to_owned(), "sdn".to_owned()],
            &operations,
            Some(&[identity.clone(), c4]),
            &[false, false],
            &identity,
        )
        .expect("spin action");
        assert_eq!(spin.local_matrices()[1].rows(), 2);
        assert!(
            magnetictb_linear_algebra::exact_unitary_matrix(&spin.local_matrices()[1]).unwrap()
        );
    }

    #[test]
    fn explicit_two_component_polynomial_basis_uses_the_exact_core() {
        let context = make_context(24, 128).expect("context");
        let identity = ExactMatrix::identity(&context, 3).expect("identity");
        let zero = PolynomialExpression::Scalar(Element::zero(&context).expect("zero"));
        let imaginary = Element::root_of_unity(&context, 4, 1).expect("i");
        let plus = PolynomialExpression::Add(vec![
            PolynomialExpression::Variable(0),
            PolynomialExpression::Multiply(vec![
                PolynomialExpression::Scalar(imaginary.clone()),
                PolynomialExpression::Variable(1),
            ]),
        ]);
        let minus = PolynomialExpression::Add(vec![
            PolynomialExpression::Variable(0),
            PolynomialExpression::Multiply(vec![
                PolynomialExpression::Scalar(imaginary.negate().expect("-i")),
                PolynomialExpression::Variable(1),
            ]),
        ]);
        let action = compile_function_basis_action(
            &[vec![plus, zero.clone()], vec![zero, minus]],
            &[SpatialOperation::new(
                identity.clone(),
                vec![scalar(&context, 0).expect("zero"); 3],
            )
            .expect("operation")],
            Some(std::slice::from_ref(&identity)),
            &[false],
            &identity,
        )
        .expect("explicit basis action");
        assert_eq!(
            action.local_matrices()[0],
            ExactMatrix::identity(&context, 2).expect("basis identity")
        );
    }

    #[test]
    fn double_valued_induced_spinor_uses_transported_local_blocks() {
        let context = make_context(24, 128).expect("context");
        let identity = matrix(&context, [[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
        let c2z = matrix(&context, [[-1, 0, 0], [0, -1, 0], [0, 0, 1]]);
        let zero = vec![scalar(&context, 0).expect("zero"); 3];
        let half_translation = vec![
            rational(&context, 1, 2).expect("half"),
            scalar(&context, 0).expect("zero"),
            scalar(&context, 0).expect("zero"),
        ];
        let operations = vec![
            SpatialOperation::new(identity.clone(), zero.clone()).expect("E"),
            SpatialOperation::new(identity.clone(), half_translation.clone()).expect("tau"),
            SpatialOperation::new(c2z.clone(), zero).expect("C2z"),
            SpatialOperation::new(c2z, half_translation).expect("C2z tau"),
        ];
        let action_table = vec![vec![0, 1], vec![1, 0], vec![0, 1], vec![1, 0]];
        let compiled = compile_catalog_induced_basis_action(
            &["sup".to_owned(), "sdn".to_owned()],
            &operations,
            None,
            &[false; 4],
            &action_table,
            &[0, 1],
            &identity,
        )
        .expect("double-valued induced basis");

        assert!(compiled.spinor());
        assert_eq!(compiled.transported_basis_states().len(), 2);
        assert!(
            compiled
                .transported_basis_states()
                .iter()
                .all(|basis| basis.len() == 2)
        );
        assert_eq!(compiled.local_blocks().len(), 4);
        assert!(compiled.local_blocks().iter().all(|row| row.len() == 2));
        let c2_block = &compiled.local_blocks()[2][0];
        let c2_squared = multiply(c2_block, c2_block).expect("C2 square");
        let minus_identity = ExactMatrix::new(
            Arc::clone(&context),
            2,
            2,
            vec![
                scalar(&context, -1).expect("minus one"),
                scalar(&context, 0).expect("zero"),
                scalar(&context, 0).expect("zero"),
                scalar(&context, -1).expect("minus one"),
            ],
        )
        .expect("minus identity");
        assert_eq!(c2_squared, minus_identity);
    }
}
