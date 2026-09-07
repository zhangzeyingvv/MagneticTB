use std::collections::BTreeMap;
use std::sync::Arc;

use cyclotomic_nullspace::{Element, ExactError, ExactMatrix, ExactResult, multiply};
use magnetictb_crystal_geometry::{PeriodicBondRecord, QuadraticReal};
use magnetictb_linear_algebra::{add, conjugate_transpose, exact_matrix_equal, subtract};
use magnetictb_symmetry::CompiledSeitzGroup;

use crate::bond_orbit::integer_rotation;
use crate::{
    BondConstraintData, DirectedBondOrbitData, propagate_hopping, reconstruct_rectangular_hopping,
};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SolvedBondTerm {
    bond_index: usize,
    orbit_index: usize,
    row_block: usize,
    column_block: usize,
    displacement: Vec<QuadraticReal>,
    matrix: ExactMatrix,
}

impl SolvedBondTerm {
    #[must_use]
    pub const fn bond_index(&self) -> usize {
        self.bond_index
    }

    #[must_use]
    pub const fn orbit_index(&self) -> usize {
        self.orbit_index
    }

    #[must_use]
    pub const fn row_block(&self) -> usize {
        self.row_block
    }

    #[must_use]
    pub const fn column_block(&self) -> usize {
        self.column_block
    }

    #[must_use]
    pub fn displacement(&self) -> &[QuadraticReal] {
        &self.displacement
    }

    #[must_use]
    pub const fn matrix(&self) -> &ExactMatrix {
        &self.matrix
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SolvedBondModel {
    block_dimensions: Vec<usize>,
    parameters_by_orbit: Vec<Vec<Element>>,
    representative_hoppings: Vec<ExactMatrix>,
    terms: Vec<SolvedBondTerm>,
}

impl SolvedBondModel {
    #[must_use]
    pub fn block_dimensions(&self) -> &[usize] {
        &self.block_dimensions
    }

    #[must_use]
    pub fn parameters_by_orbit(&self) -> &[Vec<Element>] {
        &self.parameters_by_orbit
    }

    #[must_use]
    pub fn representative_hoppings(&self) -> &[ExactMatrix] {
        &self.representative_hoppings
    }

    #[must_use]
    pub fn terms(&self) -> &[SolvedBondTerm] {
        &self.terms
    }

    #[must_use]
    pub fn parameter_count(&self) -> usize {
        self.parameters_by_orbit.iter().map(Vec::len).sum()
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FourierCoefficient {
    displacement: Vec<QuadraticReal>,
    matrix: ExactMatrix,
}

impl FourierCoefficient {
    #[must_use]
    pub fn displacement(&self) -> &[QuadraticReal] {
        &self.displacement
    }

    #[must_use]
    pub const fn matrix(&self) -> &ExactMatrix {
        &self.matrix
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct HamiltonianSymmetryVerification {
    semilinear_gamma: Vec<bool>,
    unitary_gamma_commutator: Vec<Option<bool>>,
    fourier_covariance: Vec<bool>,
}

impl HamiltonianSymmetryVerification {
    #[must_use]
    pub fn gamma_semilinear_by_operation(&self) -> &[bool] {
        &self.semilinear_gamma
    }

    #[must_use]
    pub fn gamma_unitary_commutator_by_operation(&self) -> &[Option<bool>] {
        &self.unitary_gamma_commutator
    }

    #[must_use]
    pub fn fourier_covariance_by_operation(&self) -> &[bool] {
        &self.fourier_covariance
    }

    #[must_use]
    pub fn verified(&self) -> bool {
        self.semilinear_gamma.iter().all(|&value| value)
            && self
                .unitary_gamma_commutator
                .iter()
                .flatten()
                .all(|&value| value)
            && self.fourier_covariance.iter().all(|&value| value)
    }
}

#[allow(clippy::too_many_arguments)]
pub fn reconstruct_solved_bond_model(
    symmetry: &CompiledSeitzGroup,
    bonds: &[PeriodicBondRecord],
    directed: &DirectedBondOrbitData,
    constraints: &BondConstraintData,
    local_dimensions: &[usize],
    site_actions: &[Vec<ExactMatrix>],
    parameters_by_orbit: &[Vec<Element>],
) -> ExactResult<SolvedBondModel> {
    validate_reconstruction(
        symmetry,
        bonds,
        directed,
        constraints,
        local_dimensions,
        site_actions,
        parameters_by_orbit,
    )?;
    let mut representative_hoppings = Vec::with_capacity(constraints.orbits().len());
    let mut terms = Vec::with_capacity(bonds.len());
    for (constraint_index, (constraint, parameters)) in constraints
        .orbits()
        .iter()
        .zip(parameters_by_orbit)
        .enumerate()
    {
        let (rows, columns) = constraint.dimensions();
        let representative = reconstruct_rectangular_hopping(
            &constraint.kernel().basis_matrix,
            parameters,
            rows,
            columns,
        )?;
        let directed_index = constraint.source_directed_orbit_indices()[0];
        let orbit = &directed.orbits()[directed_index];
        let representative_bond = &bonds[constraint.representative_bond_index()];
        for member in orbit.members() {
            let operation = member.transporter_operation();
            let hopping = propagate_hopping(
                &representative,
                &site_actions[operation][representative_bond.source_site()],
                &site_actions[operation][representative_bond.target_site()],
                symmetry.operations()[operation].antiunitary(),
            )?;
            terms.push(term_from_bond(
                bonds,
                member.bond_index(),
                constraint_index,
                hopping.clone(),
            ));
            if constraint.source_directed_orbit_indices().len() == 2 {
                terms.push(term_from_bond(
                    bonds,
                    directed.reverse_indices()[member.bond_index()],
                    constraint_index,
                    conjugate_transpose(&hopping)?,
                ));
            }
        }
        representative_hoppings.push(representative);
    }
    terms.sort_by_key(SolvedBondTerm::bond_index);
    if terms.len() != bonds.len()
        || terms
            .iter()
            .enumerate()
            .any(|(index, term)| term.bond_index != index)
    {
        return Err(ExactError::new(
            "InvalidSolvedBondModel",
            "reconstructed terms must cover every selected directed bond exactly once",
        ));
    }
    Ok(SolvedBondModel {
        block_dimensions: local_dimensions.to_vec(),
        parameters_by_orbit: parameters_by_orbit.to_vec(),
        representative_hoppings,
        terms,
    })
}

fn term_from_bond(
    bonds: &[PeriodicBondRecord],
    bond_index: usize,
    orbit_index: usize,
    matrix: ExactMatrix,
) -> SolvedBondTerm {
    let bond = &bonds[bond_index];
    SolvedBondTerm {
        bond_index,
        orbit_index,
        row_block: bond.source_site(),
        column_block: bond.target_site(),
        displacement: bond.displacement().to_vec(),
        matrix,
    }
}

#[allow(clippy::too_many_arguments)]
fn validate_reconstruction(
    symmetry: &CompiledSeitzGroup,
    bonds: &[PeriodicBondRecord],
    directed: &DirectedBondOrbitData,
    constraints: &BondConstraintData,
    local_dimensions: &[usize],
    site_actions: &[Vec<ExactMatrix>],
    parameters_by_orbit: &[Vec<Element>],
) -> ExactResult<()> {
    if bonds.is_empty()
        || local_dimensions.is_empty()
        || local_dimensions.contains(&0)
        || directed.action().object_count() != bonds.len()
        || parameters_by_orbit.len() != constraints.orbits().len()
        || site_actions.len() != symmetry.group().order()
        || site_actions
            .iter()
            .any(|operation| operation.len() != local_dimensions.len())
    {
        return Err(ExactError::new(
            "InvalidSolvedBondModel",
            "bond, constraint, parameter, and site-action layers do not align",
        ));
    }
    for (constraint, parameters) in constraints.orbits().iter().zip(parameters_by_orbit) {
        if parameters.len() != constraint.kernel().nullity
            || constraint.source_directed_orbit_indices().is_empty()
            || constraint.source_directed_orbit_indices().len() > 2
            || constraint
                .source_directed_orbit_indices()
                .iter()
                .any(|&index| index >= directed.orbits().len())
            || directed.orbits()[constraint.source_directed_orbit_indices()[0]]
                .representative_bond_index()
                != constraint.representative_bond_index()
        {
            return Err(ExactError::new(
                "InvalidSolvedBondModel",
                "constraint orbit does not match the directed-bond partition",
            ));
        }
    }
    Ok(())
}

pub fn fourier_coefficients(model: &SolvedBondModel) -> ExactResult<Vec<FourierCoefficient>> {
    let Some(first) = model.terms.first() else {
        return Err(ExactError::new(
            "InvalidSolvedBondModel",
            "a solved bond model must have at least one term",
        ));
    };
    let context = first.matrix.context();
    let total_dimension: usize = model.block_dimensions.iter().sum();
    let offsets = block_offsets(&model.block_dimensions);
    let mut coefficients = BTreeMap::<Vec<QuadraticReal>, ExactMatrix>::new();
    for term in &model.terms {
        let global = embed_block(
            context,
            total_dimension,
            offsets[term.row_block],
            offsets[term.column_block],
            &term.matrix,
        )?;
        if let Some(matrix) = coefficients.get_mut(&term.displacement) {
            *matrix = add(matrix, &global)?;
        } else {
            coefficients.insert(term.displacement.clone(), global);
        }
    }
    Ok(coefficients
        .into_iter()
        .map(|(displacement, matrix)| FourierCoefficient {
            displacement,
            matrix,
        })
        .collect())
}

pub fn verify_hamiltonian_symmetry(
    model: &SolvedBondModel,
    symmetry: &CompiledSeitzGroup,
    representation_matrices: &[ExactMatrix],
) -> ExactResult<HamiltonianSymmetryVerification> {
    if representation_matrices.len() != symmetry.group().order()
        || representation_matrices.is_empty()
    {
        return Err(ExactError::new(
            "InvalidHamiltonianSymmetryInput",
            "one full representation matrix is required per ordered operation",
        ));
    }
    let coefficients = fourier_coefficients(model)?;
    let coefficient_map = coefficients
        .iter()
        .map(|coefficient| (coefficient.displacement.clone(), &coefficient.matrix))
        .collect::<BTreeMap<_, _>>();
    let gamma = coefficients.iter().try_fold(
        ExactMatrix::zero(
            representation_matrices[0].context(),
            representation_matrices[0].rows(),
            representation_matrices[0].columns(),
        )?,
        |total, coefficient| add(&total, &coefficient.matrix),
    )?;
    let mut gamma_semilinear = Vec::with_capacity(symmetry.group().order());
    let mut gamma_commutators = Vec::with_capacity(symmetry.group().order());
    let mut covariance = Vec::with_capacity(symmetry.group().order());
    for (operation, representation) in representation_matrices.iter().enumerate() {
        let antiunitary = symmetry.operations()[operation].antiunitary();
        let transformed_gamma =
            propagate_hopping(&gamma, representation, representation, antiunitary)?;
        gamma_semilinear.push(exact_matrix_equal(&transformed_gamma, &gamma));
        if antiunitary {
            gamma_commutators.push(None);
        } else {
            let commutator = subtract(
                &multiply(representation, &gamma)?,
                &multiply(&gamma, representation)?,
            )?;
            gamma_commutators.push(Some(commutator.entries().iter().all(Element::is_zero)));
        }
        let rotation = integer_rotation(symmetry, operation)?;
        let mut operation_verified = true;
        for coefficient in &coefficients {
            let image_displacement = rotate_displacement(&rotation, &coefficient.displacement)?;
            let Some(expected) = coefficient_map.get(&image_displacement) else {
                operation_verified = false;
                break;
            };
            let transformed = propagate_hopping(
                &coefficient.matrix,
                representation,
                representation,
                antiunitary,
            )?;
            if !exact_matrix_equal(&transformed, expected) {
                operation_verified = false;
                break;
            }
        }
        covariance.push(operation_verified);
    }
    Ok(HamiltonianSymmetryVerification {
        semilinear_gamma: gamma_semilinear,
        unitary_gamma_commutator: gamma_commutators,
        fourier_covariance: covariance,
    })
}

fn block_offsets(dimensions: &[usize]) -> Vec<usize> {
    let mut total = 0;
    dimensions
        .iter()
        .map(|&dimension| {
            let offset = total;
            total += dimension;
            offset
        })
        .collect()
}

fn embed_block(
    context: &Arc<cyclotomic_nullspace::CyclotomicContext>,
    total: usize,
    row_offset: usize,
    column_offset: usize,
    block: &ExactMatrix,
) -> ExactResult<ExactMatrix> {
    if block.context() != context
        || row_offset + block.rows() > total
        || column_offset + block.columns() > total
    {
        return Err(ExactError::new(
            "InvalidSolvedBondModel",
            "term block cannot be embedded in the global Hamiltonian",
        ));
    }
    let mut entries = vec![Element::zero(context)?; total * total];
    for row in 0..block.rows() {
        for column in 0..block.columns() {
            entries[(row_offset + row) * total + column_offset + column] =
                block.entry(row, column)?.clone();
        }
    }
    ExactMatrix::new(Arc::clone(context), total, total, entries)
}

fn rotate_displacement(
    rotation: &[Vec<i64>],
    displacement: &[QuadraticReal],
) -> ExactResult<Vec<QuadraticReal>> {
    if rotation.len() != displacement.len()
        || rotation.iter().any(|row| row.len() != displacement.len())
    {
        return Err(ExactError::new(
            "DimensionMismatch",
            "spatial rotation and bond displacement dimensions differ",
        ));
    }
    Ok(rotation
        .iter()
        .map(|row| {
            row.iter().zip(displacement).fold(
                QuadraticReal::zero(),
                |total, (&coefficient, value)| {
                    total.add(&value.scale(&cyclotomic_nullspace::Rational::from_i64(coefficient)))
                },
            )
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use cyclotomic_nullspace::{Element, ExactMatrix, Rational, make_context};
    use magnetictb_crystal_geometry::{
        QuadraticReal, compile_site_permutations, periodic_bond_shells_in_box,
    };
    use magnetictb_representation::compile_direct_product;
    use magnetictb_symmetry::{SeitzOperation, SpatialOperation, compile_seitz_group};

    use crate::{
        compile_directed_bond_orbits, compile_hermitian_bond_constraints,
        reconstruct_solved_bond_model, verify_hamiltonian_symmetry,
    };

    fn element(
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

    fn real(numerator: i64, denominator: i64) -> QuadraticReal {
        QuadraticReal::from_rational(
            Rational::parse(&numerator.to_string(), &denominator.to_string()).expect("rational"),
        )
    }

    #[test]
    fn half_translation_model_reconstructs_and_is_exactly_covariant() {
        let context = make_context(4, 8).expect("context containing i");
        let identity = ExactMatrix::identity(&context, 1).expect("identity");
        let symmetry = compile_seitz_group(vec![
            SeitzOperation::new(
                "e",
                "E",
                SpatialOperation::new(identity.clone(), vec![element(&context, 0, 1)])
                    .expect("identity"),
                false,
            ),
            SeitzOperation::new(
                "t",
                "t(1/2)'",
                SpatialOperation::new(identity.clone(), vec![element(&context, 1, 2)])
                    .expect("translation"),
                true,
            ),
        ])
        .expect("symmetry");
        let sites = compile_site_permutations(
            symmetry.group(),
            vec![vec![
                vec![element(&context, 0, 1)],
                vec![element(&context, 1, 2)],
            ]],
            symmetry
                .operations()
                .iter()
                .map(|operation| operation.spatial().clone())
                .collect(),
        )
        .expect("sites");
        let shells = periodic_bond_shells_in_box(
            &[vec![real(0, 1)], vec![real(1, 2)]],
            &[vec![real(1, 1)]],
            1,
        )
        .expect("shells");
        let bonds = shells[1].bonds();
        let directed =
            compile_directed_bond_orbits(&symmetry, &sites, bonds).expect("directed bonds");
        let site_actions = vec![
            vec![identity.clone(), identity.clone()],
            vec![identity.clone(), identity.clone()],
        ];
        let constraints =
            compile_hermitian_bond_constraints(&symmetry, bonds, &directed, &[1, 1], &site_actions)
                .expect("constraints");
        let model = reconstruct_solved_bond_model(
            &symmetry,
            bonds,
            &directed,
            &constraints,
            &[1, 1],
            &site_actions,
            &[vec![element(&context, 2, 1), element(&context, 3, 1)]],
        )
        .expect("solved model");
        assert_eq!(model.parameter_count(), 2);
        assert_eq!(model.terms().len(), 4);
        let representation = compile_direct_product(
            symmetry.group(),
            sites.actions(),
            &[vec![identity; 2]],
            &symmetry.antiunitary_flags(),
        )
        .expect("representation");
        let verification = verify_hamiltonian_symmetry(
            &model,
            &symmetry,
            representation.representation_matrices(),
        )
        .expect("symmetry verification");
        assert!(verification.verified());
        assert_eq!(
            verification.gamma_unitary_commutator_by_operation(),
            [Some(true), None]
        );
        assert_eq!(verification.fourier_covariance_by_operation(), [true, true]);
    }
}
