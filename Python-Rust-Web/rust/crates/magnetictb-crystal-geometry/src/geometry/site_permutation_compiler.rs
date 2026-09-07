#[derive(Clone, Debug, Eq, PartialEq)]
pub struct GeometryError {
    tag: String,
    detail: String,
}
impl GeometryError {
    #[must_use]
    pub fn new(tag: impl Into<String>, detail: impl Into<String>) -> Self {
        Self {
            tag: tag.into(),
            detail: detail.into(),
        }
    }

    #[must_use]
    pub fn tag(&self) -> &str {
        &self.tag
    }
}

impl Display for GeometryError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        write!(formatter, "{}: {}", self.tag, self.detail)
    }
}

impl Error for GeometryError {}

impl From<ExactError> for GeometryError {
    fn from(error: ExactError) -> Self {
        Self::new(error.tag(), error.to_string())
    }
}

impl From<GroupError> for GeometryError {
    fn from(error: GroupError) -> Self {
        Self::new(error.tag(), error.to_string())
    }
}

impl From<SymmetryError> for GeometryError {
    fn from(error: SymmetryError) -> Self {
        Self::new(error.tag(), error.to_string())
    }
}

pub type GeometryResult<T> = Result<T, GeometryError>;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SiteImageRecord {
    pub source_site: usize,
    pub target_site: usize,
    pub cell_translation: Vec<Element>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SitePermutationData {
    site_orbits: Vec<Vec<Vec<Element>>>,
    spatial_actions: Vec<SpatialOperation>,
    image_records: Vec<Vec<Vec<SiteImageRecord>>>,
    actions: Vec<GroupAction>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MsgWyckoffSiteData {
    source_ordinal: usize,
    letter: String,
    symmetry: CompiledSeitzGroup,
    sites: SitePermutationData,
}

impl MsgWyckoffSiteData {
    #[must_use]
    pub const fn source_ordinal(&self) -> usize {
        self.source_ordinal
    }

    #[must_use]
    pub fn letter(&self) -> &str {
        &self.letter
    }

    #[must_use]
    pub const fn symmetry(&self) -> &CompiledSeitzGroup {
        &self.symmetry
    }

    #[must_use]
    pub const fn sites(&self) -> &SitePermutationData {
        &self.sites
    }
}

impl SitePermutationData {
    #[must_use]
    pub fn site_orbits(&self) -> &[Vec<Vec<Element>>] {
        &self.site_orbits
    }

    #[must_use]
    pub fn spatial_actions(&self) -> &[SpatialOperation] {
        &self.spatial_actions
    }

    #[must_use]
    pub fn image_records(&self) -> &[Vec<Vec<SiteImageRecord>>] {
        &self.image_records
    }

    #[must_use]
    pub fn actions(&self) -> &[GroupAction] {
        &self.actions
    }

    #[must_use]
    pub fn image_site_indices(&self) -> Vec<Vec<Vec<usize>>> {
        self.image_records
            .iter()
            .map(|orbit| {
                orbit
                    .iter()
                    .map(|operation| operation.iter().map(|record| record.target_site).collect())
                    .collect()
            })
            .collect()
    }

    #[must_use]
    pub fn cell_translations(&self) -> Vec<Vec<Vec<Vec<Element>>>> {
        self.image_records
            .iter()
            .map(|orbit| {
                orbit
                    .iter()
                    .map(|operation| {
                        operation
                            .iter()
                            .map(|record| record.cell_translation.clone())
                            .collect()
                    })
                    .collect()
            })
            .collect()
    }

    pub fn site_symmetry_operations(
        &self,
        orbit: usize,
        site: usize,
    ) -> GeometryResult<Vec<usize>> {
        let action = self.actions.get(orbit).ok_or_else(|| {
            GeometryError::new("SiteIndexOutOfRange", "orbit index is out of range")
        })?;
        action.stabilizer(site).map_err(Into::into)
    }
}

pub fn compile_site_permutations(
    group: &GroupAlgebra,
    site_orbits: Vec<Vec<Vec<Element>>>,
    spatial_actions: Vec<SpatialOperation>,
) -> GeometryResult<SitePermutationData> {
    let dimension = validate_inputs(group, &site_orbits, &spatial_actions)?;
    let mut all_records = Vec::with_capacity(site_orbits.len());
    let mut actions = Vec::with_capacity(site_orbits.len());
    for (orbit_index, orbit) in site_orbits.iter().enumerate() {
        let mut orbit_records = Vec::with_capacity(spatial_actions.len());
        let mut action_table = Vec::with_capacity(spatial_actions.len());
        for (operation_index, operation) in spatial_actions.iter().enumerate() {
            let mut operation_records = Vec::with_capacity(orbit.len());
            let mut images = Vec::with_capacity(orbit.len());
            for source in 0..orbit.len() {
                let transformed = affine_image(operation, &orbit[source])?;
                let mut matches = Vec::new();
                for (target, target_site) in orbit.iter().enumerate() {
                    let difference = subtract_vectors(&transformed, target_site)?;
                    if integer_vector(&difference)? {
                        matches.push((target, difference));
                    }
                }
                if matches.len() != 1 {
                    return Err(GeometryError::new(
                        "SiteImageNotUnique",
                        format!(
                            "operation {operation_index} maps orbit {orbit_index} source {source} to {} targets modulo the lattice",
                            matches.len()
                        ),
                    ));
                }
                let (target, cell_translation) = matches.remove(0);
                images.push(target);
                operation_records.push(SiteImageRecord {
                    source_site: source,
                    target_site: target,
                    cell_translation,
                });
            }
            let mut sorted = images.clone();
            sorted.sort_unstable();
            if sorted != (0..orbit.len()).collect::<Vec<_>>() {
                return Err(GeometryError::new(
                    "SiteActionNotBijective",
                    format!("operation {operation_index} is not bijective on orbit {orbit_index}"),
                ));
            }
            action_table.push(images);
            orbit_records.push(operation_records);
        }
        actions.push(GroupAction::compile(group, action_table)?);
        all_records.push(orbit_records);
    }
    debug_assert!(dimension > 0);
    Ok(SitePermutationData {
        site_orbits,
        spatial_actions,
        image_records: all_records,
        actions,
    })
}
