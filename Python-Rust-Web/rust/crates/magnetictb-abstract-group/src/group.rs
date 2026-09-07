use crate::GroupError;
use std::collections::{BTreeSet, VecDeque};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct GroupAlgebra {
    multiplication_table: Vec<Vec<usize>>,
    identity: usize,
    inverse_indices: Vec<usize>,
}

impl GroupAlgebra {
    /// Tests whether a table defines a valid ordered finite group.
    #[must_use]
    pub fn is_valid_table(multiplication_table: &[Vec<usize>]) -> bool {
        Self::new(multiplication_table.to_vec()).is_ok()
    }

    /// Validates and constructs an ordered finite group from a Cayley table.
    ///
    /// # Errors
    ///
    /// Returns an explicit tagged error when the table is malformed or violates
    /// a finite-group axiom.
    pub fn new(multiplication_table: Vec<Vec<usize>>) -> Result<Self, GroupError> {
        let order = multiplication_table.len();
        if order == 0 || multiplication_table.iter().any(|row| row.len() != order) {
            return Err(GroupError::new(
                "table_not_square",
                "create_group_algebra",
                "expected a nonempty square multiplication table",
            ));
        }

        if multiplication_table
            .iter()
            .flatten()
            .any(|&entry| entry >= order)
        {
            return Err(GroupError::new(
                "table_value_out_of_range",
                "create_group_algebra",
                "a multiplication-table entry is outside the ordered group",
            ));
        }

        let expected: Vec<usize> = (0..order).collect();
        let rows_are_latin = multiplication_table.iter().all(|row| {
            let mut sorted = row.clone();
            sorted.sort_unstable();
            sorted == expected
        });
        let columns_are_latin = (0..order).all(|column| {
            let mut values: Vec<usize> =
                multiplication_table.iter().map(|row| row[column]).collect();
            values.sort_unstable();
            values == expected
        });
        if !rows_are_latin || !columns_are_latin {
            return Err(GroupError::new(
                "non_group_non_latin_table",
                "create_group_algebra",
                "every row and column must be a permutation of all group indices",
            ));
        }

        let identity_candidates: Vec<usize> = (0..order)
            .filter(|&candidate| {
                multiplication_table[candidate] == expected
                    && (0..order).all(|element| multiplication_table[element][candidate] == element)
            })
            .collect();

        let associative = (0..order).all(|left| {
            (0..order).all(|right| {
                (0..order).all(|third| {
                    multiplication_table[multiplication_table[left][right]][third]
                        == multiplication_table[left][multiplication_table[right][third]]
                })
            })
        });
        if !associative {
            let tag = if identity_candidates.is_empty() {
                "no_identity_nonassociative_table"
            } else {
                "table_not_associative"
            };
            return Err(GroupError::new(
                tag,
                "create_group_algebra",
                "the multiplication table is not associative",
            ));
        }

        if identity_candidates.len() != 1 {
            return Err(GroupError::new(
                "identity_not_unique",
                "create_group_algebra",
                "the multiplication table must have exactly one two-sided identity",
            ));
        }
        let identity = identity_candidates[0];

        let inverse_indices: Result<Vec<usize>, GroupError> = (0..order)
            .map(|element| {
                (0..order)
                    .find(|&candidate| {
                        multiplication_table[element][candidate] == identity
                            && multiplication_table[candidate][element] == identity
                    })
                    .ok_or_else(|| {
                        GroupError::new(
                            "inverse_not_found",
                            "create_group_algebra",
                            format!("group element {element} has no two-sided inverse"),
                        )
                    })
            })
            .collect();

        Ok(Self {
            multiplication_table,
            identity,
            inverse_indices: inverse_indices?,
        })
    }

    #[must_use]
    pub fn order(&self) -> usize {
        self.multiplication_table.len()
    }

    #[must_use]
    pub const fn identity(&self) -> usize {
        self.identity
    }

    #[must_use]
    pub fn inverse_indices(&self) -> &[usize] {
        &self.inverse_indices
    }

    #[must_use]
    pub fn multiplication_table(&self) -> &[Vec<usize>] {
        &self.multiplication_table
    }

    /// Computes `left * right` using the frozen zero-based table convention.
    ///
    /// # Errors
    ///
    /// Returns `element_index_out_of_range` when either index is invalid.
    pub fn product(&self, left: usize, right: usize) -> Result<usize, GroupError> {
        if left >= self.order() || right >= self.order() {
            return Err(GroupError::new(
                "element_index_out_of_range",
                "product_index",
                format!(
                    "left={left} and right={right} must both be below {}",
                    self.order()
                ),
            ));
        }
        Ok(self.multiplication_table[left][right])
    }

    /// Returns the sorted subgroup generated by the supplied element indices.
    ///
    /// # Errors
    ///
    /// Returns `generator_index_out_of_range` for an invalid generator.
    pub fn generated_subgroup(&self, generators: &[usize]) -> Result<Vec<usize>, GroupError> {
        if generators
            .iter()
            .any(|&generator| generator >= self.order())
        {
            return Err(GroupError::new(
                "generator_index_out_of_range",
                "generated_subgroup",
                "every generator index must refer to an ordered group element",
            ));
        }

        let mut unique_generators = Vec::new();
        for &generator in generators {
            if !unique_generators.contains(&generator) {
                unique_generators.push(generator);
            }
        }

        let mut seen = BTreeSet::from([self.identity]);
        let mut queue = VecDeque::from([self.identity]);
        while let Some(current) = queue.pop_front() {
            for &generator in &unique_generators {
                for product in [
                    self.multiplication_table[current][generator],
                    self.multiplication_table[generator][current],
                ] {
                    if seen.insert(product) {
                        queue.push_back(product);
                    }
                }
            }
        }
        Ok(seen.into_iter().collect())
    }

    /// Finds deterministic greedy generators for the group or a subgroup.
    ///
    /// # Errors
    ///
    /// Returns `subgroup_not_closed` when the requested indices are not a subgroup.
    pub fn find_generator_indices(
        &self,
        subgroup: Option<&[usize]>,
    ) -> Result<Vec<usize>, GroupError> {
        let indices = match subgroup {
            Some(values) => {
                self.validate_subgroup(values, "subgroup_not_closed", "find_generators")?
            }
            None => (0..self.order()).collect(),
        };

        if indices == [self.identity] {
            return Ok(vec![self.identity]);
        }

        let mut generators = Vec::new();
        let mut generated = vec![self.identity];
        for &candidate in &indices {
            if generated.contains(&candidate) {
                continue;
            }
            let mut candidates = generators.clone();
            candidates.push(candidate);
            let closure = self.generated_subgroup(&candidates)?;
            if !closure.iter().all(|element| indices.contains(element)) {
                return Err(GroupError::new(
                    "subgroup_not_closed",
                    "find_generators",
                    "candidate generators leave the requested subgroup",
                ));
            }
            if closure.len() > generated.len() {
                generators.push(candidate);
                generated = closure;
            }
            if generated == indices {
                break;
            }
        }

        if generated == indices {
            Ok(generators)
        } else {
            Err(GroupError::new(
                "subgroup_not_closed",
                "find_generators",
                "the requested subgroup cannot be generated as a closed subgroup",
            ))
        }
    }

    /// Returns `representative * subgroup` in normalized subgroup order.
    ///
    /// # Errors
    ///
    /// Returns a tagged error for an invalid representative or non-subgroup.
    pub fn right_coset(
        &self,
        subgroup: &[usize],
        representative: usize,
    ) -> Result<Vec<usize>, GroupError> {
        if representative >= self.order() {
            return Err(GroupError::new(
                "coset_representative_out_of_range",
                "right_coset",
                "the coset representative is outside the ordered group",
            ));
        }
        let normalized =
            self.validate_subgroup(subgroup, "coset_subgroup_not_closed", "right_coset")?;
        Ok(normalized
            .iter()
            .map(|&element| self.multiplication_table[representative][element])
            .collect())
    }

    /// Returns `subgroup * representative` in normalized subgroup order.
    ///
    /// # Errors
    ///
    /// Returns a tagged error for an invalid representative or non-subgroup.
    pub fn left_coset(
        &self,
        subgroup: &[usize],
        representative: usize,
    ) -> Result<Vec<usize>, GroupError> {
        if representative >= self.order() {
            return Err(GroupError::new(
                "coset_representative_out_of_range",
                "left_coset",
                "the coset representative is outside the ordered group",
            ));
        }
        let normalized =
            self.validate_subgroup(subgroup, "coset_subgroup_not_closed", "left_coset")?;
        Ok(normalized
            .iter()
            .map(|&element| self.multiplication_table[element][representative])
            .collect())
    }

    /// Solves `operation * rep[source] = rep[target] * h` exactly.
    ///
    /// # Errors
    ///
    /// Returns a tagged error for invalid inputs or when the decomposition is not
    /// unique.
    pub fn schreier_decomposition(
        &self,
        subgroup: &[usize],
        representatives: &[usize],
        operation: usize,
        source: usize,
    ) -> Result<(usize, usize), GroupError> {
        self.validate_subgroup(
            subgroup,
            "schreier_invalid_subgroup",
            "schreier_decomposition",
        )?;
        if representatives
            .iter()
            .any(|&representative| representative >= self.order())
        {
            return Err(GroupError::new(
                "schreier_representative_out_of_range",
                "schreier_decomposition",
                "a Schreier representative is outside the ordered group",
            ));
        }
        if operation >= self.order() {
            return Err(GroupError::new(
                "schreier_operation_out_of_range",
                "schreier_decomposition",
                "the Schreier operation is outside the ordered group",
            ));
        }
        if source >= representatives.len() {
            return Err(GroupError::new(
                "schreier_source_out_of_range",
                "schreier_decomposition",
                "the Schreier source is outside the representative list",
            ));
        }

        let left_product = self.multiplication_table[operation][representatives[source]];
        let mut matches = Vec::new();
        for (target, &representative) in representatives.iter().enumerate() {
            for &subgroup_element in subgroup {
                if subgroup_element < self.order()
                    && self.multiplication_table[representative][subgroup_element] == left_product
                {
                    matches.push((target, subgroup_element));
                }
            }
        }
        if matches.len() == 1 {
            Ok(matches[0])
        } else {
            Err(GroupError::new(
                "schreier_nonunique_decomposition",
                "schreier_decomposition",
                format!(
                    "expected one Schreier decomposition, found {}",
                    matches.len()
                ),
            ))
        }
    }

    fn validate_subgroup(
        &self,
        subgroup: &[usize],
        error_tag: &'static str,
        error_stage: &'static str,
    ) -> Result<Vec<usize>, GroupError> {
        let mut indices = subgroup.to_vec();
        indices.sort_unstable();
        indices.dedup();
        let valid = !indices.is_empty()
            && indices.contains(&self.identity)
            && indices.iter().all(|&element| element < self.order())
            && indices.iter().all(|&left| {
                indices
                    .iter()
                    .all(|&right| indices.contains(&self.multiplication_table[left][right]))
            });
        if !valid {
            return Err(GroupError::new(
                error_tag,
                error_stage,
                "expected a closed subgroup containing the identity",
            ));
        }
        Ok(indices)
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct OrderedElement {
    stable_id: String,
    label: String,
    antiunitary: bool,
}

impl OrderedElement {
    #[must_use]
    pub fn new(stable_id: impl Into<String>, label: impl Into<String>, antiunitary: bool) -> Self {
        Self {
            stable_id: stable_id.into(),
            label: label.into(),
            antiunitary,
        }
    }

    #[must_use]
    pub fn stable_id(&self) -> &str {
        &self.stable_id
    }

    #[must_use]
    pub fn label(&self) -> &str {
        &self.label
    }

    #[must_use]
    pub const fn antiunitary(&self) -> bool {
        self.antiunitary
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct OrderedFiniteGroup {
    algebra: GroupAlgebra,
    elements: Vec<OrderedElement>,
}

impl OrderedFiniteGroup {
    /// Attaches stable ordered metadata to a validated group algebra.
    ///
    /// # Errors
    ///
    /// Returns `element_metadata_mismatch` unless metadata aligns one-to-one with
    /// the table order and has unique nonempty stable identifiers.
    pub fn new(algebra: GroupAlgebra, elements: Vec<OrderedElement>) -> Result<Self, GroupError> {
        let stable_ids: BTreeSet<&str> = elements.iter().map(OrderedElement::stable_id).collect();
        if elements.len() != algebra.order()
            || stable_ids.len() != elements.len()
            || elements
                .iter()
                .any(|element| element.stable_id().is_empty())
        {
            return Err(GroupError::new(
                "element_metadata_mismatch",
                "ordered_group",
                "ordered element metadata must align one-to-one with the Cayley table",
            ));
        }
        Ok(Self { algebra, elements })
    }

    #[must_use]
    pub const fn algebra(&self) -> &GroupAlgebra {
        &self.algebra
    }

    #[must_use]
    pub fn elements(&self) -> &[OrderedElement] {
        &self.elements
    }
}
