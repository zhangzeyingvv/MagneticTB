use cyclotomic_nullspace::{Element, ExactError, ExactResult};
use magnetictb_abstract_group::GroupAction;
use magnetictb_crystal_geometry::{PeriodicBondRecord, SitePermutationData};
use magnetictb_symmetry::CompiledSeitzGroup;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DirectedBondOrbitMember {
    bond_index: usize,
    transporter_operation: usize,
}

impl DirectedBondOrbitMember {
    #[must_use]
    pub const fn bond_index(&self) -> usize {
        self.bond_index
    }

    #[must_use]
    pub const fn transporter_operation(&self) -> usize {
        self.transporter_operation
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DirectedBondOrbit {
    representative_bond_index: usize,
    members: Vec<DirectedBondOrbitMember>,
}

impl DirectedBondOrbit {
    #[must_use]
    pub const fn representative_bond_index(&self) -> usize {
        self.representative_bond_index
    }

    #[must_use]
    pub fn members(&self) -> &[DirectedBondOrbitMember] {
        &self.members
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DirectedBondOrbitData {
    action: GroupAction,
    reverse_indices: Vec<usize>,
    orbits: Vec<DirectedBondOrbit>,
}

impl DirectedBondOrbitData {
    #[must_use]
    pub const fn action(&self) -> &GroupAction {
        &self.action
    }

    #[must_use]
    pub fn reverse_indices(&self) -> &[usize] {
        &self.reverse_indices
    }

    #[must_use]
    pub fn orbits(&self) -> &[DirectedBondOrbit] {
        &self.orbits
    }
}

pub fn compile_directed_bond_orbits(
    symmetry: &CompiledSeitzGroup,
    site_data: &SitePermutationData,
    bonds: &[PeriodicBondRecord],
) -> ExactResult<DirectedBondOrbitData> {
    if bonds.is_empty()
        || symmetry.group().order() != site_data.spatial_actions().len()
        || symmetry.group().order() != site_data.image_records()[0].len()
    {
        return Err(ExactError::new(
            "InvalidDirectedBondInput",
            "nonempty bonds and aligned symmetry/site actions are required",
        ));
    }
    let action_table = compile_bond_action_table(symmetry, site_data, bonds)?;
    let action = GroupAction::compile(symmetry.group(), action_table)
        .map_err(|error| ExactError::new(error.tag(), error.to_string()))?;
    let reverse_indices = compile_reverse_indices(bonds)?;
    let orbits = compile_orbits(&action)?;
    Ok(DirectedBondOrbitData {
        action,
        reverse_indices,
        orbits,
    })
}

fn compile_bond_action_table(
    symmetry: &CompiledSeitzGroup,
    site_data: &SitePermutationData,
    bonds: &[PeriodicBondRecord],
) -> ExactResult<Vec<Vec<usize>>> {
    let offsets = site_offsets(site_data);
    let mut action_table = Vec::with_capacity(symmetry.group().order());
    for operation_index in 0..symmetry.group().order() {
        let rotation = integer_rotation(symmetry, operation_index)?;
        let mut row = Vec::with_capacity(bonds.len());
        for bond in bonds {
            let (source_image, source_cell) =
                global_site_image(site_data, &offsets, operation_index, bond.source_site())?;
            let (target_image, target_cell) =
                global_site_image(site_data, &offsets, operation_index, bond.target_site())?;
            let rotated_translation = multiply_integer_matrix(&rotation, bond.translation())?;
            let image_translation = target_cell
                .iter()
                .zip(&rotated_translation)
                .zip(&source_cell)
                .map(|((&target, &rotated), &source)| {
                    target
                        .checked_add(rotated)
                        .and_then(|value| value.checked_sub(source))
                        .ok_or_else(|| {
                            ExactError::new(
                                "IntegerOverflow",
                                "transformed bond translation overflows i64",
                            )
                        })
                })
                .collect::<ExactResult<Vec<_>>>()?;
            let image_index = bonds
                .iter()
                .position(|candidate| {
                    candidate.source_site() == source_image
                        && candidate.target_site() == target_image
                        && candidate.translation() == image_translation
                })
                .ok_or_else(|| {
                    ExactError::new(
                        "BondImageOutsideShell",
                        format!(
                            "operation {operation_index} maps a bond outside the selected shell"
                        ),
                    )
                })?;
            row.push(image_index);
        }
        action_table.push(row);
    }
    Ok(action_table)
}

fn compile_reverse_indices(bonds: &[PeriodicBondRecord]) -> ExactResult<Vec<usize>> {
    bonds
        .iter()
        .map(|bond| {
            let reverse_translation = bond
                .translation()
                .iter()
                .map(|&value| {
                    value.checked_neg().ok_or_else(|| {
                        ExactError::new("IntegerOverflow", "bond translation negation overflows")
                    })
                })
                .collect::<ExactResult<Vec<_>>>()?;
            bonds
                .iter()
                .position(|candidate| {
                    candidate.source_site() == bond.target_site()
                        && candidate.target_site() == bond.source_site()
                        && candidate.translation() == reverse_translation
                })
                .ok_or_else(|| {
                    ExactError::new(
                        "ReverseBondMissing",
                        "the selected directed shell is not closed under reversal",
                    )
                })
        })
        .collect()
}

fn compile_orbits(action: &GroupAction) -> ExactResult<Vec<DirectedBondOrbit>> {
    action
        .orbits()
        .map_err(|error| ExactError::new(error.tag(), error.to_string()))?
        .into_iter()
        .map(|indices| {
            let representative_bond_index = indices[0];
            let members = indices
                .into_iter()
                .map(|bond_index| {
                    let transporter_operation = action
                        .transporters(representative_bond_index, bond_index)
                        .map_err(|error| ExactError::new(error.tag(), error.to_string()))?
                        .into_iter()
                        .next()
                        .ok_or_else(|| {
                            ExactError::new(
                                "MissingBondTransporter",
                                "bond orbit member has no ordered transporter",
                            )
                        })?;
                    Ok(DirectedBondOrbitMember {
                        bond_index,
                        transporter_operation,
                    })
                })
                .collect::<ExactResult<Vec<_>>>()?;
            Ok(DirectedBondOrbit {
                representative_bond_index,
                members,
            })
        })
        .collect()
}

fn site_offsets(site_data: &SitePermutationData) -> Vec<usize> {
    let mut total = 0;
    site_data
        .site_orbits()
        .iter()
        .map(|orbit| {
            let offset = total;
            total += orbit.len();
            offset
        })
        .collect()
}

fn global_site_image(
    site_data: &SitePermutationData,
    offsets: &[usize],
    operation: usize,
    global_site: usize,
) -> ExactResult<(usize, Vec<i64>)> {
    let (orbit_index, local_site) = offsets
        .iter()
        .enumerate()
        .find_map(|(orbit_index, &offset)| {
            let length = site_data.site_orbits()[orbit_index].len();
            (global_site >= offset && global_site < offset + length)
                .then_some((orbit_index, global_site - offset))
        })
        .ok_or_else(|| ExactError::new("SiteIndexOutOfRange", "bond site is out of range"))?;
    let record = &site_data.image_records()[orbit_index][operation][local_site];
    Ok((
        offsets[orbit_index] + record.target_site,
        record
            .cell_translation
            .iter()
            .map(element_to_i64)
            .collect::<ExactResult<Vec<_>>>()?,
    ))
}

pub(crate) fn integer_rotation(
    symmetry: &CompiledSeitzGroup,
    operation: usize,
) -> ExactResult<Vec<Vec<i64>>> {
    let matrix = symmetry.operations()[operation].spatial().rotation();
    (0..matrix.rows())
        .map(|row| {
            (0..matrix.columns())
                .map(|column| matrix.entry(row, column).and_then(element_to_i64))
                .collect()
        })
        .collect()
}

fn element_to_i64(value: &Element) -> ExactResult<i64> {
    if value.coefficients()[1..]
        .iter()
        .any(|coefficient| !coefficient.is_zero())
        || !value.coefficients()[0].is_integer()
    {
        return Err(ExactError::new(
            "NonintegerSpatialAction",
            "bond cell actions require exact integer rotation and translation data",
        ));
    }
    value.coefficients()[0]
        .numerator_string()
        .parse::<i64>()
        .map_err(|_| ExactError::new("IntegerOverflow", "exact integer is outside i64 range"))
}

fn multiply_integer_matrix(matrix: &[Vec<i64>], vector: &[i64]) -> ExactResult<Vec<i64>> {
    if matrix.len() != vector.len() || matrix.iter().any(|row| row.len() != vector.len()) {
        return Err(ExactError::new(
            "DimensionMismatch",
            "bond translation and spatial rotation dimensions differ",
        ));
    }
    matrix
        .iter()
        .map(|row| {
            row.iter()
                .zip(vector)
                .try_fold(0_i64, |total, (&left, &right)| {
                    left.checked_mul(right)
                        .and_then(|product| total.checked_add(product))
                        .ok_or_else(|| {
                            ExactError::new(
                                "IntegerOverflow",
                                "integer matrix-vector product overflows",
                            )
                        })
                })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use cyclotomic_nullspace::{Element, ExactMatrix, Rational, make_context};
    use magnetictb_crystal_geometry::{
        QuadraticReal, compile_site_permutations, periodic_bond_shells_in_box,
    };
    use magnetictb_symmetry::{SeitzOperation, SpatialOperation, compile_seitz_group};

    use super::compile_directed_bond_orbits;

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
    fn half_translation_compiles_directed_orbits_and_reversals() {
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
                SpatialOperation::new(identity, vec![element(&context, 1, 2)])
                    .expect("translation"),
                true,
            ),
        ])
        .expect("Seitz group");
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
        .expect("site action");
        let shells = periodic_bond_shells_in_box(
            &[vec![real(0, 1)], vec![real(1, 2)]],
            &[vec![real(1, 1)]],
            1,
        )
        .expect("shells");
        let compiled = compile_directed_bond_orbits(&symmetry, &sites, shells[1].bonds())
            .expect("directed bond orbits");
        assert_eq!(
            compiled.action().action_table(),
            [vec![0, 1, 2, 3], vec![2, 3, 0, 1]]
        );
        assert_eq!(compiled.reverse_indices(), [3, 2, 1, 0]);
        assert_eq!(
            compiled
                .orbits()
                .iter()
                .map(|orbit| {
                    orbit
                        .members()
                        .iter()
                        .map(super::DirectedBondOrbitMember::bond_index)
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>(),
            [vec![0, 2], vec![1, 3]]
        );
    }
}
