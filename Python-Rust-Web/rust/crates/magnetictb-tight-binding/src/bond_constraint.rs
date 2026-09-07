use std::collections::BTreeSet;

use cyclotomic_nullspace::{ExactError, ExactMatrix, ExactResult};
use magnetictb_crystal_geometry::PeriodicBondRecord;
use magnetictb_linear_algebra::{
    ConstraintKernelMethod, ConstraintKernelResult, ConstraintTarget, ConstraintValidationLevel,
    RectangularConstraintOperation, compile_rectangular_constraint_blocks,
    solve_constraint_kernel_with_options,
};
use magnetictb_symmetry::CompiledSeitzGroup;

use crate::DirectedBondOrbitData;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BondConstraintOrbit {
    source_directed_orbit_indices: Vec<usize>,
    representative_bond_index: usize,
    dimensions: (usize, usize),
    stabilizer_generator_indices: Vec<usize>,
    reverse_representative_operation: Option<usize>,
    kernel: ConstraintKernelResult,
}

impl BondConstraintOrbit {
    #[must_use]
    pub fn source_directed_orbit_indices(&self) -> &[usize] {
        &self.source_directed_orbit_indices
    }

    #[must_use]
    pub const fn representative_bond_index(&self) -> usize {
        self.representative_bond_index
    }

    #[must_use]
    pub const fn dimensions(&self) -> (usize, usize) {
        self.dimensions
    }

    #[must_use]
    pub fn stabilizer_generator_indices(&self) -> &[usize] {
        &self.stabilizer_generator_indices
    }

    #[must_use]
    pub const fn reverse_representative_operation(&self) -> Option<usize> {
        self.reverse_representative_operation
    }

    #[must_use]
    pub const fn kernel(&self) -> &ConstraintKernelResult {
        &self.kernel
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BondConstraintData {
    orbits: Vec<BondConstraintOrbit>,
    parameter_count: usize,
}

impl BondConstraintData {
    #[must_use]
    pub fn orbits(&self) -> &[BondConstraintOrbit] {
        &self.orbits
    }

    #[must_use]
    pub const fn parameter_count(&self) -> usize {
        self.parameter_count
    }
}

pub fn compile_hermitian_bond_constraints(
    symmetry: &CompiledSeitzGroup,
    bonds: &[PeriodicBondRecord],
    directed: &DirectedBondOrbitData,
    local_dimensions: &[usize],
    site_actions: &[Vec<ExactMatrix>],
) -> ExactResult<BondConstraintData> {
    compile_bond_constraints(
        symmetry,
        bonds,
        directed,
        local_dimensions,
        site_actions,
        true,
        ConstraintKernelMethod::Iterative,
        ConstraintValidationLevel::Basic,
    )
}

#[allow(clippy::too_many_arguments)]
pub fn compile_bond_constraints(
    symmetry: &CompiledSeitzGroup,
    bonds: &[PeriodicBondRecord],
    directed: &DirectedBondOrbitData,
    local_dimensions: &[usize],
    site_actions: &[Vec<ExactMatrix>],
    hermitian: bool,
    method: ConstraintKernelMethod,
    validation_level: ConstraintValidationLevel,
) -> ExactResult<BondConstraintData> {
    compile_bond_constraints_with_continuous(
        symmetry,
        bonds,
        directed,
        local_dimensions,
        site_actions,
        None,
        hermitian,
        method,
        validation_level,
    )
}

pub fn compile_hermitian_bond_constraints_with_continuous(
    symmetry: &CompiledSeitzGroup,
    bonds: &[PeriodicBondRecord],
    directed: &DirectedBondOrbitData,
    local_dimensions: &[usize],
    site_actions: &[Vec<ExactMatrix>],
    continuous_site_generators: Option<&[ExactMatrix]>,
) -> ExactResult<BondConstraintData> {
    compile_bond_constraints_with_continuous(
        symmetry,
        bonds,
        directed,
        local_dimensions,
        site_actions,
        continuous_site_generators,
        true,
        ConstraintKernelMethod::Iterative,
        ConstraintValidationLevel::Basic,
    )
}

#[allow(clippy::too_many_arguments)]
pub fn compile_bond_constraints_with_continuous(
    symmetry: &CompiledSeitzGroup,
    bonds: &[PeriodicBondRecord],
    directed: &DirectedBondOrbitData,
    local_dimensions: &[usize],
    site_actions: &[Vec<ExactMatrix>],
    continuous_site_generators: Option<&[ExactMatrix]>,
    hermitian: bool,
    method: ConstraintKernelMethod,
    validation_level: ConstraintValidationLevel,
) -> ExactResult<BondConstraintData> {
    validate_inputs(
        symmetry,
        bonds,
        directed,
        local_dimensions,
        site_actions,
        continuous_site_generators,
    )?;
    let bond_to_orbit = bond_to_orbit(directed, bonds.len())?;
    let mut unseen = (0..directed.orbits().len()).collect::<BTreeSet<_>>();
    let mut orbits = Vec::new();
    while let Some(orbit_index) = unseen.first().copied() {
        let directed_orbit = &directed.orbits()[orbit_index];
        let representative = directed_orbit.representative_bond_index();
        let reverse_bond = directed.reverse_indices()[representative];
        let reverse_orbit = bond_to_orbit[reverse_bond];
        let bond = &bonds[representative];
        let dimensions = (
            local_dimensions[bond.source_site()],
            local_dimensions[bond.target_site()],
        );
        let stabilizer = directed
            .action()
            .stabilizer(representative)
            .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
        let generators = symmetry
            .group()
            .find_generator_indices(Some(&stabilizer))
            .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
        let mut operations = generators
            .iter()
            .map(|&operation| {
                constraint_operation(
                    symmetry,
                    site_actions,
                    bond,
                    operation,
                    ConstraintTarget::Same,
                )
            })
            .collect::<Vec<_>>();
        let reverse_representative_operation = if hermitian {
            append_reverse_constraint(
                symmetry,
                directed,
                site_actions,
                bond,
                representative,
                reverse_bond,
                orbit_index == reverse_orbit,
                dimensions,
                &mut operations,
            )?
        } else {
            None
        };
        let continuous_pair = continuous_site_generators.map(|generators| {
            (
                &generators[bond.source_site()],
                &generators[bond.target_site()],
            )
        });
        let blocks = compile_rectangular_constraint_blocks(
            dimensions.0,
            dimensions.1,
            &operations,
            continuous_pair,
        )?;
        let coordinate_dimension = 2 * dimensions.0 * dimensions.1;
        let kernel = solve_constraint_kernel_with_options(
            site_actions[0][bond.source_site()].context(),
            &blocks,
            coordinate_dimension,
            method,
            validation_level,
        )?;
        let source_directed_orbit_indices = if !hermitian || reverse_orbit == orbit_index {
            vec![orbit_index]
        } else {
            vec![orbit_index, reverse_orbit]
        };
        for index in &source_directed_orbit_indices {
            unseen.remove(index);
        }
        orbits.push(BondConstraintOrbit {
            source_directed_orbit_indices,
            representative_bond_index: representative,
            dimensions,
            stabilizer_generator_indices: generators,
            reverse_representative_operation,
            kernel,
        });
    }
    let parameter_count = orbits.iter().map(|orbit| orbit.kernel.nullity).sum();
    Ok(BondConstraintData {
        orbits,
        parameter_count,
    })
}

fn constraint_operation(
    symmetry: &CompiledSeitzGroup,
    site_actions: &[Vec<ExactMatrix>],
    bond: &PeriodicBondRecord,
    operation: usize,
    target: ConstraintTarget,
) -> RectangularConstraintOperation {
    RectangularConstraintOperation {
        left: site_actions[operation][bond.source_site()].clone(),
        right: site_actions[operation][bond.target_site()].clone(),
        antiunitary: symmetry.operations()[operation].antiunitary(),
        target,
    }
}

#[allow(clippy::too_many_arguments)]
fn append_reverse_constraint(
    symmetry: &CompiledSeitzGroup,
    directed: &DirectedBondOrbitData,
    site_actions: &[Vec<ExactMatrix>],
    bond: &PeriodicBondRecord,
    representative: usize,
    reverse_bond: usize,
    self_reverse_orbit: bool,
    dimensions: (usize, usize),
    operations: &mut Vec<RectangularConstraintOperation>,
) -> ExactResult<Option<usize>> {
    if !self_reverse_orbit {
        return Ok(None);
    }
    if reverse_bond == representative {
        operations.push(RectangularConstraintOperation {
            left: ExactMatrix::identity(
                site_actions[0][bond.source_site()].context(),
                dimensions.0,
            )?,
            right: ExactMatrix::identity(
                site_actions[0][bond.target_site()].context(),
                dimensions.1,
            )?,
            antiunitary: false,
            target: ConstraintTarget::Reverse,
        });
        return Ok(None);
    }
    let operation = directed
        .action()
        .transporters(representative, reverse_bond)
        .map_err(|error| ExactError::new(error.tag(), error.to_string()))?
        .into_iter()
        .next()
        .ok_or_else(|| {
            ExactError::new(
                "MissingReverseBondTransporter",
                "self-reverse directed orbit has no reverse transporter",
            )
        })?;
    operations.push(constraint_operation(
        symmetry,
        site_actions,
        bond,
        operation,
        ConstraintTarget::Reverse,
    ));
    Ok(Some(operation))
}

fn validate_inputs(
    symmetry: &CompiledSeitzGroup,
    bonds: &[PeriodicBondRecord],
    directed: &DirectedBondOrbitData,
    local_dimensions: &[usize],
    site_actions: &[Vec<ExactMatrix>],
    continuous_site_generators: Option<&[ExactMatrix]>,
) -> ExactResult<()> {
    if bonds.is_empty()
        || local_dimensions.is_empty()
        || local_dimensions.contains(&0)
        || directed.action().object_count() != bonds.len()
        || site_actions.len() != symmetry.group().order()
        || site_actions
            .iter()
            .any(|operation| operation.len() != local_dimensions.len())
        || bonds.iter().any(|bond| {
            bond.source_site() >= local_dimensions.len()
                || bond.target_site() >= local_dimensions.len()
        })
        || continuous_site_generators
            .is_some_and(|generators| generators.len() != local_dimensions.len())
    {
        return Err(ExactError::new(
            "InvalidBondConstraintInput",
            "bond, orbit, dimension, and site-action data do not align",
        ));
    }
    let context = site_actions[0][0].context();
    for operation in site_actions {
        for (site, matrix) in operation.iter().enumerate() {
            if matrix.context() != context
                || matrix.rows() != local_dimensions[site]
                || matrix.columns() != local_dimensions[site]
            {
                return Err(ExactError::new(
                    "InvalidBondConstraintInput",
                    "site representation matrices have incompatible exact shape",
                ));
            }
        }
    }
    if let Some(generators) = continuous_site_generators {
        for (site, generator) in generators.iter().enumerate() {
            if generator.context() != context
                || generator.rows() != local_dimensions[site]
                || generator.columns() != local_dimensions[site]
            {
                return Err(ExactError::new(
                    "InvalidBondConstraintInput",
                    "continuous site generators have incompatible exact shape",
                ));
            }
        }
    }
    Ok(())
}

fn bond_to_orbit(directed: &DirectedBondOrbitData, bond_count: usize) -> ExactResult<Vec<usize>> {
    let mut result = vec![usize::MAX; bond_count];
    for (orbit_index, orbit) in directed.orbits().iter().enumerate() {
        for member in orbit.members() {
            if member.bond_index() >= bond_count || result[member.bond_index()] != usize::MAX {
                return Err(ExactError::new(
                    "InvalidBondConstraintInput",
                    "directed bond orbits do not partition the selected shell",
                ));
            }
            result[member.bond_index()] = orbit_index;
        }
    }
    if result.contains(&usize::MAX) {
        return Err(ExactError::new(
            "InvalidBondConstraintInput",
            "directed bond orbits do not cover the selected shell",
        ));
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use cyclotomic_nullspace::{Element, ExactMatrix, Rational, make_context};
    use magnetictb_crystal_geometry::{
        QuadraticReal, compile_site_permutations, periodic_bond_shells_in_box,
    };
    use magnetictb_linear_algebra::{ConstraintKernelMethod, ConstraintValidationLevel};
    use magnetictb_symmetry::{SeitzOperation, SpatialOperation, compile_seitz_group};

    use crate::{
        compile_bond_constraints, compile_directed_bond_orbits, compile_hermitian_bond_constraints,
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
    fn trivial_onsite_hermiticity_has_one_real_parameter() {
        let context = make_context(1, 8).expect("context");
        let symmetry = compile_seitz_group(vec![SeitzOperation::new(
            "e",
            "E",
            SpatialOperation::new(
                ExactMatrix::identity(&context, 1).expect("rotation"),
                vec![element(&context, 0, 1)],
            )
            .expect("operation"),
            false,
        )])
        .expect("symmetry");
        let sites = compile_site_permutations(
            symmetry.group(),
            vec![vec![vec![element(&context, 0, 1)]]],
            vec![symmetry.operations()[0].spatial().clone()],
        )
        .expect("sites");
        let shells = periodic_bond_shells_in_box(&[vec![real(0, 1)]], &[vec![real(1, 1)]], 0)
            .expect("shells");
        let directed =
            compile_directed_bond_orbits(&symmetry, &sites, shells[0].bonds()).expect("directed");
        let constraints = compile_hermitian_bond_constraints(
            &symmetry,
            shells[0].bonds(),
            &directed,
            &[1],
            &[vec![ExactMatrix::identity(&context, 1).expect("identity")]],
        )
        .expect("constraints");
        assert_eq!(constraints.parameter_count(), 1);
        assert_eq!(constraints.orbits()[0].kernel().nullity, 1);

        let nonhermitian = compile_bond_constraints(
            &symmetry,
            shells[0].bonds(),
            &directed,
            &[1],
            &[vec![ExactMatrix::identity(&context, 1).expect("identity")]],
            false,
            ConstraintKernelMethod::Iterative,
            ConstraintValidationLevel::Basic,
        )
        .expect("non-Hermitian constraints");
        assert_eq!(nonhermitian.parameter_count(), 2);
        assert_eq!(nonhermitian.orbits().len(), 1);
        assert_eq!(
            nonhermitian.orbits()[0].source_directed_orbit_indices(),
            &[0]
        );
    }

    #[test]
    fn paired_half_translation_bond_orbits_have_two_real_parameters() {
        let context = make_context(1, 8).expect("context");
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
        let directed =
            compile_directed_bond_orbits(&symmetry, &sites, shells[1].bonds()).expect("directed");
        let site_actions = vec![
            vec![identity.clone(), identity.clone()],
            vec![identity.clone(), identity],
        ];
        let constraints = compile_hermitian_bond_constraints(
            &symmetry,
            shells[1].bonds(),
            &directed,
            &[1, 1],
            &site_actions,
        )
        .expect("constraints");
        assert_eq!(constraints.parameter_count(), 2);
        assert_eq!(
            constraints.orbits()[0].source_directed_orbit_indices(),
            [0, 1]
        );
        assert_eq!(constraints.orbits()[0].kernel().basis_matrix.rows(), 2);
        assert_eq!(constraints.orbits()[0].kernel().basis_matrix.columns(), 2);
    }
}
