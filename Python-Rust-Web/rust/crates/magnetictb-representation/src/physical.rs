use cyclotomic_nullspace::ExactMatrix;
use magnetictb_abstract_group::{GroupAction, GroupAlgebra};

use crate::{
    CompiledRepresentation, RepresentationError, RepresentationResult, representation::assemble,
};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SiteActionData {
    method: String,
    group: GroupAlgebra,
    actions: Vec<GroupAction>,
    cell_translations: Vec<Vec<Vec<Vec<i64>>>>,
    local_blocks: Vec<Vec<Vec<ExactMatrix>>>,
    local_dimensions: Vec<usize>,
    antiunitary_flags: Vec<bool>,
    coordinate_dimension: usize,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SiteActionRecord {
    pub orbit_index: usize,
    pub source_site_index: usize,
    pub image_site_index: usize,
    pub cell_translation: Vec<i64>,
    pub matrix: ExactMatrix,
    pub antiunitary: bool,
}

impl SiteActionData {
    pub fn new(
        method: impl Into<String>,
        group: &GroupAlgebra,
        actions: Vec<GroupAction>,
        cell_translations: Vec<Vec<Vec<Vec<i64>>>>,
        local_blocks: Vec<Vec<Vec<ExactMatrix>>>,
        local_dimensions: Vec<usize>,
        antiunitary_flags: Vec<bool>,
    ) -> RepresentationResult<Self> {
        if actions.is_empty()
            || actions.len() != cell_translations.len()
            || actions.len() != local_blocks.len()
            || actions.len() != local_dimensions.len()
            || antiunitary_flags.len() != group.order()
        {
            return Err(RepresentationError::new(
                "InvalidSiteActionData",
                "site-action orbit layers do not align",
            ));
        }
        let coordinate_dimension = cell_translations
            .first()
            .and_then(|orbit| orbit.first())
            .and_then(|operation| operation.first())
            .map_or(0, Vec::len);
        if coordinate_dimension == 0 {
            return Err(RepresentationError::new(
                "InvalidSiteActionData",
                "cell translations must have a positive coordinate dimension",
            ));
        }
        for orbit in 0..actions.len() {
            let action = &actions[orbit];
            if action.group() != group
                || cell_translations[orbit].len() != group.order()
                || local_blocks[orbit].len() != group.order()
            {
                return Err(RepresentationError::new(
                    "InvalidSiteActionData",
                    "an orbit does not align with the ordered group",
                ));
            }
            for operation in 0..group.order() {
                if cell_translations[orbit][operation].len() != action.object_count()
                    || local_blocks[orbit][operation].len() != action.object_count()
                    || cell_translations[orbit][operation]
                        .iter()
                        .any(|translation| translation.len() != coordinate_dimension)
                {
                    return Err(RepresentationError::new(
                        "InvalidSiteActionData",
                        "operation-by-source translations and local blocks do not align",
                    ));
                }
            }
        }
        let result = Self {
            method: method.into(),
            group: group.clone(),
            actions,
            cell_translations,
            local_blocks,
            local_dimensions,
            antiunitary_flags,
            coordinate_dimension,
        };
        result.assemble_full()?;
        Ok(result)
    }

    #[must_use]
    pub fn method(&self) -> &str {
        &self.method
    }

    #[must_use]
    pub const fn coordinate_dimension(&self) -> usize {
        self.coordinate_dimension
    }

    pub fn site_action(
        &self,
        orbit: usize,
        source: usize,
        operation: usize,
    ) -> RepresentationResult<SiteActionRecord> {
        let action = self.actions.get(orbit).ok_or_else(|| {
            RepresentationError::new("SiteActionIndexOutOfRange", "orbit index is out of range")
        })?;
        if source >= action.object_count() || operation >= self.group.order() {
            return Err(RepresentationError::new(
                "SiteActionIndexOutOfRange",
                "source or operation index is out of range",
            ));
        }
        Ok(SiteActionRecord {
            orbit_index: orbit,
            source_site_index: source,
            image_site_index: action.image(operation, source)?,
            cell_translation: self.cell_translations[orbit][operation][source].clone(),
            matrix: self.local_blocks[orbit][operation][source].clone(),
            antiunitary: self.antiunitary_flags[operation],
        })
    }

    pub fn assemble_full(&self) -> RepresentationResult<CompiledRepresentation> {
        assemble(
            &self.method,
            &self.group,
            &self.actions,
            self.local_blocks.clone(),
            self.local_dimensions.clone(),
            &self.antiunitary_flags,
            false,
        )
    }
}
