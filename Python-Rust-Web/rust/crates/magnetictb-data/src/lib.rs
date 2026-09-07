#![forbid(unsafe_code)]

mod catalog;
mod error;
mod exact;

pub use catalog::{
    BravaisLattice, ClassificationMapEntry, ClassificationMaps, DataCatalog, DataCounts,
    DataResourceManifest, MsgGroup, ResourceSection, SubperiodicDataset, SubperiodicGroup,
    SymmetryOperation, WyckoffEntry, WyckoffGroup, WyckoffPosition,
};
pub use error::DataError;
pub use exact::{
    EncodedList, EncodedValue, ExactExpression, ExactMatrix, OrderedAssociation,
    OrderedAssociationEntry,
};

const MANIFEST: &str = include_str!("../resources/manifest.json");
const BRAVAIS_LATTICES: &str = include_str!("../resources/bravais_lattices.json");
const MSG_GROUPS: &str = include_str!("../resources/msg_groups.json");
const CLASSIFICATION_MAPS: &str = include_str!("../resources/classification_maps.json");
const WYCKOFF_GROUPS: &str = include_str!("../resources/wyckoff_groups.json");
const ROD: &str = include_str!("../resources/rod.json");
const LAYER: &str = include_str!("../resources/layer.json");

impl DataCatalog {
    /// Loads the self-contained data snapshot imported from the frozen Mathematica baseline.
    ///
    /// # Errors
    ///
    /// Returns an explicit tagged error if any embedded resource is malformed or its
    /// count, order, stable IDs, or cross-references do not close.
    pub fn load_embedded() -> Result<Self, DataError> {
        Self::from_parts(
            parse(MANIFEST, "manifest")?,
            parse(BRAVAIS_LATTICES, "Bravais lattices")?,
            parse(MSG_GROUPS, "MSG groups")?,
            parse(CLASSIFICATION_MAPS, "classification maps")?,
            parse(WYCKOFF_GROUPS, "Wyckoff groups")?,
            parse(ROD, "rod groups")?,
            parse(LAYER, "layer groups")?,
        )
    }
}

fn parse<T>(source: &str, name: &str) -> Result<T, DataError>
where
    T: serde::de::DeserializeOwned,
{
    serde_json::from_str(source).map_err(|error| {
        DataError::new(
            "MalformedResource",
            format!("cannot decode embedded {name}: {error}"),
        )
    })
}

#[cfg(test)]
mod tests {
    use super::{DataCatalog, MsgGroup};

    #[test]
    fn embedded_catalog_closes_with_frozen_counts_and_order() {
        let catalog = DataCatalog::load_embedded().expect("frozen Data must load");
        let counts = catalog.counts();
        assert_eq!(counts.bravais_lattices, 15);
        assert_eq!(counts.msg_groups, 1_651);
        assert_eq!(counts.msg_operations, 25_033);
        assert_eq!(counts.wyckoff_groups, 1_651);
        assert_eq!(counts.wyckoff_entries, 14_212);
        assert_eq!(counts.wyckoff_positions, 85_509);
        assert_eq!(counts.rod_groups, 394);
        assert_eq!(counts.rod_operations, 4_594);
        assert_eq!(counts.layer_groups, 528);
        assert_eq!(counts.layer_operations, 4_890);

        assert_eq!(catalog.bravais_lattices()[0].stable_id(), "CubicP");
        assert_eq!(catalog.bravais_lattices()[14].stable_id(), "TriclinicP");
        assert_eq!(catalog.msg_groups()[0].stable_id(), "msg-0001");
        assert_eq!(catalog.msg_groups()[1_650].stable_id(), "msg-1651");
        assert_eq!(catalog.msg("msg-0001").map(MsgGroup::source_index), Some(1));
        assert!(catalog.wyckoff("msg-1651").is_some());
        assert!(catalog.rod_group("rod-og-1.1.1").is_some());
        assert!(catalog.layer_group("layer-og-1.1.1").is_some());
    }

    #[test]
    fn unknown_stable_ids_do_not_autocorrect() {
        let catalog = DataCatalog::load_embedded().expect("frozen Data must load");
        assert!(catalog.msg("MSG-0001").is_none());
        assert!(catalog.wyckoff("msg-9999").is_none());
        assert!(catalog.rod_group("rod-og-0.0.0").is_none());
        assert!(catalog.layer_group("layer-og-0.0.0").is_none());
    }

    #[test]
    fn high_level_bns_and_wyckoff_resolution_is_exact() {
        let catalog = DataCatalog::load_embedded().expect("frozen Data must load");
        let group = catalog.msg_by_bns(&[1, 3]).expect("BNS 1.3 exists");
        assert_eq!(group.stable_id(), "msg-0003");
        assert!(catalog.msg_by_bns(&[0, 0]).is_none());

        let entry = catalog
            .wyckoff_entry(group.stable_id(), "a", None)
            .expect("letter a is unambiguous");
        assert_eq!(entry.source_ordinal(), 2);
        assert_eq!(
            catalog
                .wyckoff_entry(group.stable_id(), "missing", None)
                .expect_err("unknown letter must fail")
                .tag(),
            "UnknownWyckoffSelection"
        );
    }

    #[test]
    fn mathematica_style_classification_selectors_are_one_based_and_exact() {
        let catalog = DataCatalog::load_embedded().expect("frozen Data must load");
        for (classification, key, expected) in [
            ("gray", &[1][..], "msg-0002"),
            ("typeI", &[1, 1][..], "msg-0001"),
            ("typeIII", &[2, 6][..], "msg-0006"),
            ("typeIV", &[1, 3][..], "msg-0003"),
            ("bns", &[1, 3][..], "msg-0003"),
            ("og", &[1, 3, 3][..], "msg-0003"),
            ("og_number", &[3][..], "msg-0003"),
        ] {
            assert_eq!(
                catalog
                    .msg_by_classification_key(classification, key)
                    .map(MsgGroup::stable_id),
                Some(expected)
            );
        }
        assert_eq!(
            catalog.msg_by_source_index(3).map(MsgGroup::stable_id),
            Some("msg-0003")
        );
        assert!(catalog.msg_by_source_index(0).is_none());
        assert!(catalog.msg_by_classification_index("gray", 0).is_none());
        assert!(
            catalog
                .msg_by_classification_index("typeIII", usize::MAX)
                .is_none()
        );
        assert!(catalog.msg_by_classification_index("unknown", 1).is_none());
    }
}
