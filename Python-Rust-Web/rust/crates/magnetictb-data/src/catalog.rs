use std::collections::{HashMap, HashSet};

use serde::{Deserialize, Serialize};

use crate::DataError;
use crate::exact::{EncodedValue, ExactExpression, ExactMatrix, OrderedAssociation};

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct BravaisLattice {
    stable_id: String,
    conventional_vectors: ExactMatrix,
    primitive_vectors: ExactMatrix,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct SymmetryOperation {
    stable_id: String,
    label: String,
    rotation: ExactMatrix,
    translation: Vec<ExactExpression>,
    antiunitary: bool,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct MsgGroup {
    stable_id: String,
    source_index: usize,
    bns_key: Vec<usize>,
    og_key: Vec<usize>,
    symbol: String,
    classification: String,
    bravais_lattice_id: String,
    operations: Vec<SymmetryOperation>,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct ClassificationMapEntry {
    source_key: Vec<usize>,
    msg_id: String,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct ClassificationMaps {
    bns: Vec<ClassificationMapEntry>,
    og: Vec<ClassificationMapEntry>,
    og_number: Vec<ClassificationMapEntry>,
    gray: Vec<ClassificationMapEntry>,
    type_i: Vec<ClassificationMapEntry>,
    type_iii: Vec<ClassificationMapEntry>,
    type_iv: Vec<ClassificationMapEntry>,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct WyckoffPosition {
    coordinates: Vec<ExactExpression>,
    moment: Vec<ExactExpression>,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct WyckoffEntry {
    source_ordinal: usize,
    multiplicity: usize,
    letter: String,
    positions: Vec<WyckoffPosition>,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct WyckoffGroup {
    msg_id: String,
    bns_key: Vec<usize>,
    entries: Vec<WyckoffEntry>,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct SubperiodicGroup {
    stable_id: String,
    og_key: Vec<usize>,
    symbol: String,
    bravais_lattice_id: String,
    operations: Vec<SymmetryOperation>,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct SubperiodicDataset {
    basic_vectors: OrderedAssociation,
    gray_selectors: OrderedAssociation,
    group_to_bravais: OrderedAssociation,
    groups: Vec<SubperiodicGroup>,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct DataCounts {
    pub bravais_lattices: usize,
    pub msg_groups: usize,
    pub msg_operations: usize,
    pub gray: usize,
    pub type_i: usize,
    pub type_iii: usize,
    pub type_iv: usize,
    pub wyckoff_groups: usize,
    pub wyckoff_entries: usize,
    pub wyckoff_positions: usize,
    pub rod_groups: usize,
    pub rod_operations: usize,
    pub layer_groups: usize,
    pub layer_operations: usize,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct DataResourceManifest {
    pub schema: String,
    pub schema_version: usize,
    pub stable_commit: String,
    pub source_commit: String,
    pub bundle_file_sha256: String,
    pub logical_bundle_sha256: String,
    pub source_schema_sha256: String,
    pub counts: DataCounts,
    pub sections: HashMap<String, ResourceSection>,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct ResourceSection {
    pub file: String,
    pub sha256: String,
}

#[derive(Debug)]
pub struct DataCatalog {
    manifest: DataResourceManifest,
    bravais_lattices: Vec<BravaisLattice>,
    msg_groups: Vec<MsgGroup>,
    classification_maps: ClassificationMaps,
    wyckoff_groups: Vec<WyckoffGroup>,
    rod: SubperiodicDataset,
    layer: SubperiodicDataset,
    bravais_index: HashMap<String, usize>,
    msg_index: HashMap<String, usize>,
    wyckoff_index: HashMap<String, usize>,
    rod_index: HashMap<String, usize>,
    layer_index: HashMap<String, usize>,
}

impl BravaisLattice {
    #[must_use]
    pub fn stable_id(&self) -> &str {
        &self.stable_id
    }

    #[must_use]
    pub const fn conventional_vectors(&self) -> &ExactMatrix {
        &self.conventional_vectors
    }

    #[must_use]
    pub const fn primitive_vectors(&self) -> &ExactMatrix {
        &self.primitive_vectors
    }
}

impl SymmetryOperation {
    #[must_use]
    pub fn stable_id(&self) -> &str {
        &self.stable_id
    }

    #[must_use]
    pub fn label(&self) -> &str {
        &self.label
    }

    #[must_use]
    pub const fn rotation(&self) -> &ExactMatrix {
        &self.rotation
    }

    #[must_use]
    pub fn translation(&self) -> &[ExactExpression] {
        &self.translation
    }

    #[must_use]
    pub const fn antiunitary(&self) -> bool {
        self.antiunitary
    }
}

impl MsgGroup {
    #[must_use]
    pub fn stable_id(&self) -> &str {
        &self.stable_id
    }

    #[must_use]
    pub const fn source_index(&self) -> usize {
        self.source_index
    }

    #[must_use]
    pub fn bns_key(&self) -> &[usize] {
        &self.bns_key
    }

    #[must_use]
    pub fn og_key(&self) -> &[usize] {
        &self.og_key
    }

    #[must_use]
    pub fn symbol(&self) -> &str {
        &self.symbol
    }

    #[must_use]
    pub fn classification(&self) -> &str {
        &self.classification
    }

    #[must_use]
    pub fn bravais_lattice_id(&self) -> &str {
        &self.bravais_lattice_id
    }

    #[must_use]
    pub fn operations(&self) -> &[SymmetryOperation] {
        &self.operations
    }
}

impl WyckoffGroup {
    #[must_use]
    pub fn msg_id(&self) -> &str {
        &self.msg_id
    }

    #[must_use]
    pub fn bns_key(&self) -> &[usize] {
        &self.bns_key
    }

    #[must_use]
    pub fn entries(&self) -> &[WyckoffEntry] {
        &self.entries
    }
}

impl WyckoffEntry {
    #[must_use]
    pub const fn source_ordinal(&self) -> usize {
        self.source_ordinal
    }

    #[must_use]
    pub const fn multiplicity(&self) -> usize {
        self.multiplicity
    }

    #[must_use]
    pub fn letter(&self) -> &str {
        &self.letter
    }

    #[must_use]
    pub fn positions(&self) -> &[WyckoffPosition] {
        &self.positions
    }
}

impl WyckoffPosition {
    #[must_use]
    pub fn coordinates(&self) -> &[ExactExpression] {
        &self.coordinates
    }

    #[must_use]
    pub fn moment(&self) -> &[ExactExpression] {
        &self.moment
    }
}

impl SubperiodicGroup {
    #[must_use]
    pub fn stable_id(&self) -> &str {
        &self.stable_id
    }

    #[must_use]
    pub fn og_key(&self) -> &[usize] {
        &self.og_key
    }

    #[must_use]
    pub fn symbol(&self) -> &str {
        &self.symbol
    }

    #[must_use]
    pub fn bravais_lattice_id(&self) -> &str {
        &self.bravais_lattice_id
    }

    #[must_use]
    pub fn operations(&self) -> &[SymmetryOperation] {
        &self.operations
    }
}

impl SubperiodicDataset {
    #[must_use]
    pub const fn basic_vectors(&self) -> &OrderedAssociation {
        &self.basic_vectors
    }

    #[must_use]
    pub const fn gray_selectors(&self) -> &OrderedAssociation {
        &self.gray_selectors
    }

    #[must_use]
    pub const fn group_to_bravais(&self) -> &OrderedAssociation {
        &self.group_to_bravais
    }

    #[must_use]
    pub fn groups(&self) -> &[SubperiodicGroup] {
        &self.groups
    }
}

impl DataCatalog {
    pub(crate) fn from_parts(
        manifest: DataResourceManifest,
        bravais_lattices: Vec<BravaisLattice>,
        msg_groups: Vec<MsgGroup>,
        classification_maps: ClassificationMaps,
        wyckoff_groups: Vec<WyckoffGroup>,
        rod: SubperiodicDataset,
        layer: SubperiodicDataset,
    ) -> Result<Self, DataError> {
        let bravais_index = unique_index(
            bravais_lattices
                .iter()
                .map(|lattice| lattice.stable_id.as_str()),
            "Bravais lattice",
        )?;
        let msg_index = unique_index(
            msg_groups.iter().map(|group| group.stable_id.as_str()),
            "MSG",
        )?;
        let wyckoff_index = unique_index(
            wyckoff_groups.iter().map(|group| group.msg_id.as_str()),
            "Wyckoff MSG",
        )?;
        let rod_index = unique_index(
            rod.groups.iter().map(|group| group.stable_id.as_str()),
            "rod group",
        )?;
        let layer_index = unique_index(
            layer.groups.iter().map(|group| group.stable_id.as_str()),
            "layer group",
        )?;
        let catalog = Self {
            manifest,
            bravais_lattices,
            msg_groups,
            classification_maps,
            wyckoff_groups,
            rod,
            layer,
            bravais_index,
            msg_index,
            wyckoff_index,
            rod_index,
            layer_index,
        };
        catalog.validate()?;
        Ok(catalog)
    }

    #[must_use]
    pub const fn manifest(&self) -> &DataResourceManifest {
        &self.manifest
    }

    #[must_use]
    pub const fn counts(&self) -> DataCounts {
        self.manifest.counts
    }

    #[must_use]
    pub fn bravais_lattices(&self) -> &[BravaisLattice] {
        &self.bravais_lattices
    }

    #[must_use]
    pub fn msg_groups(&self) -> &[MsgGroup] {
        &self.msg_groups
    }

    #[must_use]
    pub const fn classification_maps(&self) -> &ClassificationMaps {
        &self.classification_maps
    }

    #[must_use]
    pub fn wyckoff_groups(&self) -> &[WyckoffGroup] {
        &self.wyckoff_groups
    }

    #[must_use]
    pub const fn rod(&self) -> &SubperiodicDataset {
        &self.rod
    }

    #[must_use]
    pub const fn layer(&self) -> &SubperiodicDataset {
        &self.layer
    }

    #[must_use]
    pub fn bravais(&self, stable_id: &str) -> Option<&BravaisLattice> {
        self.bravais_index
            .get(stable_id)
            .map(|&index| &self.bravais_lattices[index])
    }

    #[must_use]
    pub fn msg(&self, stable_id: &str) -> Option<&MsgGroup> {
        self.msg_index
            .get(stable_id)
            .map(|&index| &self.msg_groups[index])
    }

    /// Resolves the one-based source index used by Mathematica's private
    /// `MSGOP` association and accepted directly by the public `msgop`.
    #[must_use]
    pub fn msg_by_source_index(&self, source_index: usize) -> Option<&MsgGroup> {
        source_index
            .checked_sub(1)
            .and_then(|index| self.msg_groups.get(index))
            .filter(|group| group.source_index == source_index)
    }

    #[must_use]
    pub fn msg_by_bns(&self, bns_key: &[usize]) -> Option<&MsgGroup> {
        self.msg_by_classification_key("bns", bns_key)
    }

    /// Resolves an exact key in one of Mathematica's frozen MSG dictionaries.
    ///
    /// Unlike [`Self::msg_by_classification_index`], this preserves the actual
    /// Mathematica Association key: gray and OG-number keys have length one,
    /// BNS/type keys have length two, and OG keys have length three.
    #[must_use]
    pub fn msg_by_classification_key(
        &self,
        classification: &str,
        source_key: &[usize],
    ) -> Option<&MsgGroup> {
        let entries = match classification {
            "bns" => &self.classification_maps.bns,
            "og" => &self.classification_maps.og,
            "og_number" => &self.classification_maps.og_number,
            "gray" => &self.classification_maps.gray,
            "typeI" => &self.classification_maps.type_i,
            "typeIII" => &self.classification_maps.type_iii,
            "typeIV" => &self.classification_maps.type_iv,
            _ => return None,
        };
        entries
            .iter()
            .find(|entry| entry.source_key == source_key)
            .and_then(|entry| self.msg(&entry.msg_id))
    }

    /// Resolves a stable Mathematica-style magnetic-group selector by its
    /// one-based position in the frozen classification list.
    #[must_use]
    pub fn msg_by_classification_index(
        &self,
        classification: &str,
        source_index: usize,
    ) -> Option<&MsgGroup> {
        let entries = match classification {
            "gray" => &self.classification_maps.gray,
            "typeI" => &self.classification_maps.type_i,
            "typeIII" => &self.classification_maps.type_iii,
            "typeIV" => &self.classification_maps.type_iv,
            _ => return None,
        };
        source_index
            .checked_sub(1)
            .and_then(|index| entries.get(index))
            .and_then(|entry| self.msg(&entry.msg_id))
    }

    #[must_use]
    pub fn wyckoff(&self, msg_id: &str) -> Option<&WyckoffGroup> {
        self.wyckoff_index
            .get(msg_id)
            .map(|&index| &self.wyckoff_groups[index])
    }

    /// Resolves one stable Wyckoff letter and optional one-based source ordinal.
    ///
    /// # Errors
    ///
    /// Returns a tagged data error when the MSG or selection is unknown, or
    /// when omitting the ordinal would leave the letter ambiguous.
    pub fn wyckoff_entry(
        &self,
        msg_id: &str,
        letter: &str,
        source_ordinal: Option<usize>,
    ) -> Result<&WyckoffEntry, DataError> {
        let group = self.wyckoff(msg_id).ok_or_else(|| {
            DataError::new(
                "UnknownStableId",
                format!("unknown Wyckoff MSG stable ID {msg_id}"),
            )
        })?;
        let matching = group
            .entries()
            .iter()
            .filter(|entry| {
                entry.letter() == letter
                    && source_ordinal.is_none_or(|ordinal| entry.source_ordinal() == ordinal)
            })
            .collect::<Vec<_>>();
        match matching.as_slice() {
            [entry] => Ok(*entry),
            [] => Err(DataError::new(
                "UnknownWyckoffSelection",
                match source_ordinal {
                    Some(ordinal) => format!(
                        "MSG {msg_id} has no Wyckoff entry with letter {letter} and source ordinal {ordinal}"
                    ),
                    None => format!("MSG {msg_id} has no Wyckoff entry with letter {letter}"),
                },
            )),
            _ => Err(DataError::new(
                "AmbiguousWyckoffSelection",
                format!(
                    "MSG {msg_id} has multiple Wyckoff entries with letter {letter}; supply source_ordinal"
                ),
            )),
        }
    }

    #[must_use]
    pub fn rod_group(&self, stable_id: &str) -> Option<&SubperiodicGroup> {
        self.rod_index
            .get(stable_id)
            .map(|&index| &self.rod.groups[index])
    }

    #[must_use]
    pub fn layer_group(&self, stable_id: &str) -> Option<&SubperiodicGroup> {
        self.layer_index
            .get(stable_id)
            .map(|&index| &self.layer.groups[index])
    }

    /// Resolves the exact three-integer OG key accepted by stable `mlgop` and
    /// `mrgop`.  The lookup remains in the Rust-owned ordered database.
    #[must_use]
    pub fn subperiodic_group_by_og(
        &self,
        kind: &str,
        og_key: &[usize],
    ) -> Option<&SubperiodicGroup> {
        let groups = match kind {
            "rod" => self.rod.groups(),
            "layer" => self.layer.groups(),
            _ => return None,
        };
        groups.iter().find(|group| group.og_key() == og_key)
    }

    /// Returns the three-integer OG key stored at a one-based stable gray
    /// selector.  This mirrors `grayrod[n]` and `graylayer[n]` without moving
    /// Association traversal into Python.
    #[must_use]
    pub fn subperiodic_gray_key(&self, kind: &str, source_key: usize) -> Option<Vec<usize>> {
        let association = match kind {
            "rod" => self.rod.gray_selectors(),
            "layer" => self.layer.gray_selectors(),
            _ => return None,
        };
        association.entries().iter().find_map(|entry| {
            let EncodedValue::Integer(entry_key) = entry.key() else {
                return None;
            };
            if usize::try_from(*entry_key).ok() != Some(source_key) {
                return None;
            }
            let EncodedValue::List(list) = entry.value() else {
                return None;
            };
            list.items()
                .iter()
                .map(|value| match value {
                    EncodedValue::Integer(integer) => usize::try_from(*integer).ok(),
                    _ => None,
                })
                .collect()
        })
    }

    #[must_use]
    pub fn subperiodic_basic_vectors(
        &self,
        kind: &str,
        bravais_lattice_id: &str,
    ) -> Option<&EncodedValue> {
        let association = match kind {
            "rod" => self.rod.basic_vectors(),
            "layer" => self.layer.basic_vectors(),
            _ => return None,
        };
        association.entries().iter().find_map(|entry| {
            matches!(entry.key(), EncodedValue::String(key) if key == bravais_lattice_id)
                .then_some(entry.value())
        })
    }

    fn validate(&self) -> Result<(), DataError> {
        self.validate_manifest()?;
        self.validate_msg_closure()?;
        self.validate_wyckoff_closure()?;
        validate_subperiodic(&self.rod, "rod")?;
        validate_subperiodic(&self.layer, "layer")?;
        Ok(())
    }

    fn validate_manifest(&self) -> Result<(), DataError> {
        if self.manifest.schema != "magnetictb.data.resources" || self.manifest.schema_version != 1
        {
            return Err(DataError::new(
                "MalformedResource",
                "unsupported data resource manifest schema",
            ));
        }
        let actual = DataCounts {
            bravais_lattices: self.bravais_lattices.len(),
            msg_groups: self.msg_groups.len(),
            msg_operations: self.msg_groups.iter().map(|g| g.operations.len()).sum(),
            gray: self.classification_maps.gray.len(),
            type_i: self.classification_maps.type_i.len(),
            type_iii: self.classification_maps.type_iii.len(),
            type_iv: self.classification_maps.type_iv.len(),
            wyckoff_groups: self.wyckoff_groups.len(),
            wyckoff_entries: self.wyckoff_groups.iter().map(|g| g.entries.len()).sum(),
            wyckoff_positions: self
                .wyckoff_groups
                .iter()
                .flat_map(|g| &g.entries)
                .map(|entry| entry.positions.len())
                .sum(),
            rod_groups: self.rod.groups.len(),
            rod_operations: self.rod.groups.iter().map(|g| g.operations.len()).sum(),
            layer_groups: self.layer.groups.len(),
            layer_operations: self.layer.groups.iter().map(|g| g.operations.len()).sum(),
        };
        if actual != self.manifest.counts {
            return Err(DataError::new(
                "CountMismatch",
                format!(
                    "resource counts {actual:?} differ from manifest {:?}",
                    self.manifest.counts
                ),
            ));
        }
        if self.classification_maps.bns.len() != self.msg_groups.len()
            || self.classification_maps.og.len() != self.msg_groups.len()
            || self.classification_maps.og_number.len() != self.msg_groups.len()
        {
            return Err(DataError::new(
                "CountMismatch",
                "classification maps do not cover every MSG",
            ));
        }
        Ok(())
    }

    fn validate_msg_closure(&self) -> Result<(), DataError> {
        let msg_ids: HashSet<&str> = self.msg_index.keys().map(String::as_str).collect();
        for (index, group) in self.msg_groups.iter().enumerate() {
            if group.source_index != index + 1 {
                return Err(DataError::new(
                    "OrderingMismatch",
                    format!(
                        "{} has source_index {}",
                        group.stable_id, group.source_index
                    ),
                ));
            }
            if !self.bravais_index.contains_key(&group.bravais_lattice_id) {
                return Err(DataError::new(
                    "DanglingReference",
                    format!(
                        "{} refers to unknown Bravais lattice {}",
                        group.stable_id, group.bravais_lattice_id
                    ),
                ));
            }
            validate_operations(&group.operations, &group.stable_id)?;
        }
        for entry in self
            .classification_maps
            .bns
            .iter()
            .chain(&self.classification_maps.og)
            .chain(&self.classification_maps.og_number)
            .chain(&self.classification_maps.gray)
            .chain(&self.classification_maps.type_i)
            .chain(&self.classification_maps.type_iii)
            .chain(&self.classification_maps.type_iv)
        {
            if !msg_ids.contains(entry.msg_id.as_str()) {
                return Err(DataError::new(
                    "DanglingReference",
                    format!("classification map refers to unknown MSG {}", entry.msg_id),
                ));
            }
        }
        Ok(())
    }

    fn validate_wyckoff_closure(&self) -> Result<(), DataError> {
        for group in &self.wyckoff_groups {
            let msg = self.msg(&group.msg_id).ok_or_else(|| {
                DataError::new(
                    "DanglingReference",
                    format!("Wyckoff data refers to unknown MSG {}", group.msg_id),
                )
            })?;
            if group.bns_key != msg.bns_key {
                return Err(DataError::new(
                    "DanglingReference",
                    format!("Wyckoff BNS key does not match {}", group.msg_id),
                ));
            }
            for entry in &group.entries {
                if entry.source_ordinal == 0 {
                    return Err(DataError::new(
                        "OrderingMismatch",
                        format!("{} has an invalid Wyckoff source ordinal", group.msg_id),
                    ));
                }
                if entry.positions.len() != entry.multiplicity {
                    return Err(DataError::new(
                        "CountMismatch",
                        format!(
                            "{} Wyckoff {} multiplicity mismatch",
                            group.msg_id, entry.letter
                        ),
                    ));
                }
                for position in &entry.positions {
                    if position.coordinates.len() != 3 || position.moment.len() != 3 {
                        return Err(DataError::new(
                            "MalformedResource",
                            format!("{} Wyckoff vectors must be length three", group.msg_id),
                        ));
                    }
                }
            }
        }
        Ok(())
    }
}

fn unique_index<'a>(
    ids: impl IntoIterator<Item = &'a str>,
    kind: &str,
) -> Result<HashMap<String, usize>, DataError> {
    let mut index = HashMap::new();
    for (position, stable_id) in ids.into_iter().enumerate() {
        if index.insert(stable_id.to_owned(), position).is_some() {
            return Err(DataError::new(
                "DuplicateStableId",
                format!("duplicate {kind} stable ID {stable_id}"),
            ));
        }
    }
    Ok(index)
}

fn validate_operations(operations: &[SymmetryOperation], owner: &str) -> Result<(), DataError> {
    let mut ids = HashSet::new();
    for operation in operations {
        if !ids.insert(operation.stable_id.as_str()) {
            return Err(DataError::new(
                "DuplicateStableId",
                format!("duplicate operation ID {} in {owner}", operation.stable_id),
            ));
        }
        operation
            .rotation
            .validate(&format!("{} rotation", operation.stable_id))?;
        if operation.rotation.rows() != 3
            || operation.rotation.columns() != 3
            || operation.translation.len() != 3
        {
            return Err(DataError::new(
                "MalformedResource",
                format!(
                    "{} is not a three-dimensional operation",
                    operation.stable_id
                ),
            ));
        }
    }
    Ok(())
}

fn validate_subperiodic(dataset: &SubperiodicDataset, kind: &str) -> Result<(), DataError> {
    for group in &dataset.groups {
        if group.og_key.len() != 3 || group.bravais_lattice_id.is_empty() {
            return Err(DataError::new(
                "MalformedResource",
                format!("{} has invalid {kind} group metadata", group.stable_id),
            ));
        }
        validate_operations(&group.operations, &group.stable_id)?;
    }
    Ok(())
}
