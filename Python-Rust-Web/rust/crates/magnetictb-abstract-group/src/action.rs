use crate::{GroupAlgebra, GroupError};
use std::collections::BTreeSet;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct GroupAction {
    group: GroupAlgebra,
    action_table: Vec<Vec<usize>>,
    object_count: usize,
}

impl GroupAction {
    /// Tests whether a table is a valid action of `group`.
    #[must_use]
    pub fn is_valid_table(group: &GroupAlgebra, action_table: &[Vec<usize>]) -> bool {
        Self::compile(group, action_table.to_vec()).is_ok()
    }

    /// Compiles and validates an operation-by-source action table.
    ///
    /// # Errors
    ///
    /// Returns an explicit tagged error when the table shape, permutation rows,
    /// identity action, or action law is invalid.
    pub fn compile(
        group: &GroupAlgebra,
        action_table: Vec<Vec<usize>>,
    ) -> Result<Self, GroupError> {
        let object_count = action_table.first().map_or(0, Vec::len);
        if object_count == 0
            || action_table.len() != group.order()
            || action_table.iter().any(|row| row.len() != object_count)
        {
            return Err(GroupError::new(
                "action_shape_mismatch",
                "compile_group_action",
                "expected a nonempty operation-by-source action table",
            ));
        }

        let expected: Vec<usize> = (0..object_count).collect();
        let permutation_rows = action_table.iter().all(|row| {
            let mut sorted = row.clone();
            sorted.sort_unstable();
            sorted == expected
        });
        if !permutation_rows {
            return Err(GroupError::new(
                "action_row_not_permutation",
                "compile_group_action",
                "every action row must be a permutation of the objects",
            ));
        }
        if action_table[group.identity()] != expected {
            return Err(GroupError::new(
                "action_identity_violation",
                "compile_group_action",
                "the identity operation must fix every object",
            ));
        }

        let law_holds = (0..group.order()).all(|left| {
            (0..group.order()).all(|right| {
                (0..object_count).all(|source| {
                    let product = group.multiplication_table()[left][right];
                    action_table[product][source] == action_table[left][action_table[right][source]]
                })
            })
        });
        if !law_holds {
            return Err(GroupError::new(
                "action_law_violation",
                "compile_group_action",
                "the action table does not respect group multiplication",
            ));
        }

        Ok(Self {
            group: group.clone(),
            action_table,
            object_count,
        })
    }

    #[must_use]
    pub const fn group(&self) -> &GroupAlgebra {
        &self.group
    }

    #[must_use]
    pub fn action_table(&self) -> &[Vec<usize>] {
        &self.action_table
    }

    #[must_use]
    pub const fn object_count(&self) -> usize {
        self.object_count
    }

    #[must_use]
    pub fn is_faithful(&self) -> bool {
        self.action_table
            .iter()
            .map(Vec::as_slice)
            .collect::<BTreeSet<_>>()
            .len()
            == self.group.order()
    }

    /// Returns the sorted orbit of an object.
    ///
    /// # Errors
    ///
    /// Returns `object_index_out_of_range` when `source` is invalid.
    pub fn orbit(&self, source: usize) -> Result<Vec<usize>, GroupError> {
        if source >= self.object_count {
            return Err(GroupError::new(
                "object_index_out_of_range",
                "action_orbit",
                "the action source is outside the object list",
            ));
        }
        Ok(self
            .action_table
            .iter()
            .map(|row| row[source])
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect())
    }

    /// Applies an operation to an object.
    ///
    /// # Errors
    ///
    /// Returns a tagged range error when either index is invalid.
    pub fn image(&self, operation: usize, source: usize) -> Result<usize, GroupError> {
        if operation >= self.group.order() {
            return Err(GroupError::new(
                "action_operation_out_of_range",
                "group_action",
                "the action operation is outside the ordered group",
            ));
        }
        if source >= self.object_count {
            return Err(GroupError::new(
                "object_index_out_of_range",
                "group_action",
                "the action source is outside the object list",
            ));
        }
        Ok(self.action_table[operation][source])
    }

    /// Returns all object orbits in deterministic source order.
    ///
    /// # Errors
    ///
    /// Propagates object-index validation errors from orbit construction.
    pub fn orbits(&self) -> Result<Vec<Vec<usize>>, GroupError> {
        let mut unseen: BTreeSet<usize> = (0..self.object_count).collect();
        let mut orbits = Vec::new();
        while let Some(source) = unseen.first().copied() {
            let orbit = self.orbit(source)?;
            for element in &orbit {
                unseen.remove(element);
            }
            orbits.push(orbit);
        }
        Ok(orbits)
    }

    /// Returns the ordered operation indices that fix an object.
    ///
    /// # Errors
    ///
    /// Returns `object_index_out_of_range` when `source` is invalid.
    pub fn stabilizer(&self, source: usize) -> Result<Vec<usize>, GroupError> {
        if source >= self.object_count {
            return Err(GroupError::new(
                "object_index_out_of_range",
                "stabilizer",
                "the stabilizer source is outside the object list",
            ));
        }
        Ok(self
            .action_table
            .iter()
            .enumerate()
            .filter_map(|(operation, row)| (row[source] == source).then_some(operation))
            .collect())
    }

    /// Returns all operation indices that send `source` to `target`.
    ///
    /// A valid pair with no transporter returns an empty vector.
    ///
    /// # Errors
    ///
    /// Returns a tagged source or target range error for invalid indices.
    pub fn transporters(&self, source: usize, target: usize) -> Result<Vec<usize>, GroupError> {
        if source >= self.object_count {
            return Err(GroupError::new(
                "transporter_source_out_of_range",
                "transporter",
                "the transporter source is outside the object list",
            ));
        }
        if target >= self.object_count {
            return Err(GroupError::new(
                "transporter_target_out_of_range",
                "transporter",
                "the transporter target is outside the object list",
            ));
        }
        Ok(self
            .action_table
            .iter()
            .enumerate()
            .filter_map(|(operation, row)| (row[source] == target).then_some(operation))
            .collect())
    }
}
