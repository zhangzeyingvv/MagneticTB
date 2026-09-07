#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SymmetryError {
    tag: String,
    detail: String,
}

impl SymmetryError {
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

impl Display for SymmetryError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        write!(formatter, "{}: {}", self.tag, self.detail)
    }
}

impl Error for SymmetryError {}

impl From<ExactError> for SymmetryError {
    fn from(error: ExactError) -> Self {
        Self::new(error.tag(), error.to_string())
    }
}

impl From<GroupError> for SymmetryError {
    fn from(error: GroupError) -> Self {
        Self::new(error.tag(), error.to_string())
    }
}

pub type SymmetryResult<T> = Result<T, SymmetryError>;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SpatialOperation {
    rotation: ExactMatrix,
    translation: Vec<Element>,
}

impl SpatialOperation {
    pub fn new(rotation: ExactMatrix, translation: Vec<Element>) -> SymmetryResult<Self> {
        if rotation.rows() == 0
            || rotation.rows() != rotation.columns()
            || translation.len() != rotation.rows()
            || translation
                .iter()
                .any(|value| value.context() != rotation.context())
            || rref(&rotation)?.rank != rotation.rows()
        {
            return Err(SymmetryError::new(
                "InvalidSpatialOperation",
                "rotation must be invertible and translation must have the same exact dimension",
            ));
        }
        Ok(Self {
            rotation,
            translation,
        })
    }

    #[must_use]
    pub const fn rotation(&self) -> &ExactMatrix {
        &self.rotation
    }

    #[must_use]
    pub fn translation(&self) -> &[Element] {
        &self.translation
    }

    pub fn product(&self, right: &Self) -> SymmetryResult<Self> {
        if self.rotation.context() != right.rotation.context()
            || self.rotation.rows() != right.rotation.rows()
        {
            return Err(SymmetryError::new(
                "SpatialDimensionMismatch",
                "spatial operations use different exact dimensions or fields",
            ));
        }
        let rotation = multiply(&self.rotation, &right.rotation)?;
        let rotated_translation = matrix_vector(&self.rotation, &right.translation)?;
        let translation = add_vectors(&rotated_translation, &self.translation)?;
        Self::new(rotation, translation)
    }

    pub fn translation_difference(&self, right: &Self) -> SymmetryResult<Vec<Element>> {
        if !exact_matrix_equal(&self.rotation, &right.rotation) {
            return Err(SymmetryError::new(
                "RotationMismatch",
                "spatial rotations differ",
            ));
        }
        subtract_vectors(&self.translation, &right.translation)
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SeitzOperation {
    stable_id: String,
    label: String,
    spatial: SpatialOperation,
    antiunitary: bool,
}

impl SeitzOperation {
    #[must_use]
    pub fn new(
        stable_id: impl Into<String>,
        label: impl Into<String>,
        spatial: SpatialOperation,
        antiunitary: bool,
    ) -> Self {
        Self {
            stable_id: stable_id.into(),
            label: label.into(),
            spatial,
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
    pub const fn spatial(&self) -> &SpatialOperation {
        &self.spatial
    }

    #[must_use]
    pub const fn antiunitary(&self) -> bool {
        self.antiunitary
    }

    pub fn product(&self, right: &Self) -> SymmetryResult<Self> {
        Ok(Self::new(
            "",
            "",
            self.spatial.product(&right.spatial)?,
            self.antiunitary ^ right.antiunitary,
        ))
    }

    pub fn equivalent_mod_lattice(&self, right: &Self) -> SymmetryResult<bool> {
        if self.antiunitary != right.antiunitary
            || !exact_matrix_equal(self.spatial.rotation(), right.spatial.rotation())
        {
            return Ok(false);
        }
        Ok(integer_vector(
            &self.spatial.translation_difference(&right.spatial)?,
        )?)
    }

    pub fn lattice_difference(&self, right: &Self) -> SymmetryResult<Vec<Element>> {
        if self.antiunitary != right.antiunitary {
            return Err(SymmetryError::new(
                "AntiunitaryMismatch",
                "Seitz antiunitary parity differs",
            ));
        }
        self.spatial.translation_difference(&right.spatial)
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SeitzProductRecord {
    pub operation_index: usize,
    pub lattice_translation: Vec<Element>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CompiledSeitzGroup {
    operations: Vec<SeitzOperation>,
    group: GroupAlgebra,
    product_records: Vec<Vec<SeitzProductRecord>>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SpinSpaceOperation {
    spatial: SpatialOperation,
    spin_rotation: ExactMatrix,
    antiunitary: bool,
}

impl SpinSpaceOperation {
    pub fn new(
        spatial: SpatialOperation,
        spin_rotation: ExactMatrix,
        antiunitary: bool,
    ) -> SymmetryResult<Self> {
        if spin_rotation.rows() != 3
            || spin_rotation.columns() != 3
            || spin_rotation.context() != spatial.rotation().context()
            || rref(&spin_rotation)?.rank != 3
        {
            return Err(SymmetryError::new(
                "InvalidSpinSpaceOperation",
                "spin rotation must be an invertible exact 3 by 3 matrix in the spatial field",
            ));
        }
        Ok(Self {
            spatial,
            spin_rotation,
            antiunitary,
        })
    }

    #[must_use]
    pub const fn spatial(&self) -> &SpatialOperation {
        &self.spatial
    }

    #[must_use]
    pub const fn spin_rotation(&self) -> &ExactMatrix {
        &self.spin_rotation
    }

    #[must_use]
    pub const fn antiunitary(&self) -> bool {
        self.antiunitary
    }

    pub fn product(&self, right: &Self) -> SymmetryResult<Self> {
        Self::new(
            normalize_spatial_mod_lattice(&self.spatial.product(&right.spatial)?)?,
            multiply(&self.spin_rotation, &right.spin_rotation)?,
            self.antiunitary ^ right.antiunitary,
        )
    }

    pub fn equivalent_mod_lattice(&self, right: &Self) -> SymmetryResult<bool> {
        if self.antiunitary != right.antiunitary
            || !exact_matrix_equal(&self.spin_rotation, &right.spin_rotation)
            || !exact_matrix_equal(self.spatial.rotation(), right.spatial.rotation())
        {
            return Ok(false);
        }
        Ok(integer_vector(
            &self.spatial.translation_difference(&right.spatial)?,
        )?)
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CompiledSpinSpaceGroup {
    elements: Vec<SpinSpaceOperation>,
    seitz_group: CompiledSeitzGroup,
}

impl CompiledSpinSpaceGroup {
    #[must_use]
    pub fn elements(&self) -> &[SpinSpaceOperation] {
        &self.elements
    }

    #[must_use]
    pub const fn seitz_group(&self) -> &CompiledSeitzGroup {
        &self.seitz_group
    }

    #[must_use]
    pub fn spin_rotations(&self) -> Vec<ExactMatrix> {
        self.elements
            .iter()
            .map(|element| element.spin_rotation.clone())
            .collect()
    }
}

impl CompiledSeitzGroup {
    #[must_use]
    pub fn operations(&self) -> &[SeitzOperation] {
        &self.operations
    }

    #[must_use]
    pub const fn group(&self) -> &GroupAlgebra {
        &self.group
    }

    #[must_use]
    pub fn product_records(&self) -> &[Vec<SeitzProductRecord>] {
        &self.product_records
    }

    #[must_use]
    pub fn antiunitary_flags(&self) -> Vec<bool> {
        self.operations
            .iter()
            .map(SeitzOperation::antiunitary)
            .collect()
    }
}

pub fn compile_seitz_group(operations: Vec<SeitzOperation>) -> SymmetryResult<CompiledSeitzGroup> {
    if operations.is_empty() {
        return Err(SymmetryError::new(
            "EmptySeitzGroup",
            "at least one ordered Seitz operation is required",
        ));
    }
    for left in 0..operations.len() {
        for right in left + 1..operations.len() {
            if operations[left].equivalent_mod_lattice(&operations[right])? {
                return Err(SymmetryError::new(
                    "DuplicateSeitzOperation",
                    format!("operations {left} and {right} are lattice-equivalent"),
                ));
            }
        }
    }
    let mut product_records = Vec::with_capacity(operations.len());
    let mut table = Vec::with_capacity(operations.len());
    for (left_index, left) in operations.iter().enumerate() {
        let mut record_row = Vec::with_capacity(operations.len());
        let mut table_row = Vec::with_capacity(operations.len());
        for (right_index, right) in operations.iter().enumerate() {
            let product = left.product(right)?;
            let mut matches = Vec::new();
            for (index, candidate) in operations.iter().enumerate() {
                if product.equivalent_mod_lattice(candidate)? {
                    matches.push(index);
                }
            }
            if matches.len() != 1 {
                return Err(SymmetryError::new(
                    "SeitzClosureFailure",
                    format!(
                        "product ({left_index}, {right_index}) has {} ordered representatives",
                        matches.len()
                    ),
                ));
            }
            let operation_index = matches[0];
            let lattice_translation = product.lattice_difference(&operations[operation_index])?;
            if !integer_vector(&lattice_translation)? {
                return Err(SymmetryError::new(
                    "NonintegerLatticeTranslation",
                    "Seitz product difference is not an exact integer vector",
                ));
            }
            table_row.push(operation_index);
            record_row.push(SeitzProductRecord {
                operation_index,
                lattice_translation,
            });
        }
        table.push(table_row);
        product_records.push(record_row);
    }
    let group = GroupAlgebra::new(table)?;
    Ok(CompiledSeitzGroup {
        operations,
        group,
        product_records,
    })
}

/// Compiles ordered affine operations against an explicit finite-group table.
///
/// This form supports spin-space groups where distinct ordered elements can
/// share the same spatial Seitz action while their spin actions distinguish
/// the group product.
pub fn compile_seitz_group_with_table(
    operations: Vec<SeitzOperation>,
    multiplication_table: Vec<Vec<usize>>,
) -> SymmetryResult<CompiledSeitzGroup> {
    let group = GroupAlgebra::new(multiplication_table)?;
    if operations.len() != group.order() || operations.is_empty() {
        return Err(SymmetryError::new(
            "InvalidSeitzGroupTable",
            "operation count must match the explicit finite-group table",
        ));
    }
    let mut product_records = Vec::with_capacity(group.order());
    for left in 0..group.order() {
        let mut row = Vec::with_capacity(group.order());
        for right in 0..group.order() {
            let product = operations[left].product(&operations[right])?;
            let operation_index = group.multiplication_table()[left][right];
            let representative = &operations[operation_index];
            if !product.equivalent_mod_lattice(representative)? {
                return Err(SymmetryError::new(
                    "SeitzTableMismatch",
                    format!("explicit product ({left}, {right}) has the wrong spatial action"),
                ));
            }
            let lattice_translation = product.lattice_difference(representative)?;
            if !integer_vector(&lattice_translation)? {
                return Err(SymmetryError::new(
                    "NonintegerLatticeTranslation",
                    "explicit Seitz product difference is not an exact integer vector",
                ));
            }
            row.push(SeitzProductRecord {
                operation_index,
                lattice_translation,
            });
        }
        product_records.push(row);
    }
    Ok(CompiledSeitzGroup {
        operations,
        group,
        product_records,
    })
}

/// Returns the canonical representative of an exact rational fractional
/// coordinate in the half-open interval `[0, 1)`.
///
/// Fractional crystallographic translations must be rational.  Refusing a
/// non-rational field element here keeps site-orbit generation exact and
/// prevents an ordering-dependent numerical reduction.
pub fn fractional_part(value: &Element) -> SymmetryResult<Element> {
    if value
        .coefficients()
        .iter()
        .skip(1)
        .any(|entry| !entry.is_zero())
    {
        return Err(SymmetryError::new(
            "NonrationalLatticeTranslation",
            "lattice translation normalization requires an exact rational value",
        ));
    }
    let coefficient = &value.coefficients()[0];
    let numerator = BigInt::from_str(&coefficient.numerator_string()).map_err(|_| {
        SymmetryError::new(
            "InvalidLatticeTranslation",
            "rational numerator could not be normalized",
        )
    })?;
    let denominator = BigInt::from_str(&coefficient.denominator_string()).map_err(|_| {
        SymmetryError::new(
            "InvalidLatticeTranslation",
            "rational denominator could not be normalized",
        )
    })?;
    let mut remainder = numerator % &denominator;
    if remainder < BigInt::from(0) {
        remainder += &denominator;
    }
    Element::from_polynomial(
        Arc::clone(value.context()),
        &[Rational::new(remainder, denominator)?],
    )
    .map_err(Into::into)
}

fn normalize_spatial_mod_lattice(operation: &SpatialOperation) -> SymmetryResult<SpatialOperation> {
    SpatialOperation::new(
        operation.rotation.clone(),
        operation
            .translation
            .iter()
            .map(fractional_part)
            .collect::<SymmetryResult<Vec<_>>>()?,
    )
}

fn compare_element(left: &Element, right: &Element) -> Ordering {
    left.coefficients().cmp(right.coefficients())
}

fn compare_matrix(left: &ExactMatrix, right: &ExactMatrix) -> Ordering {
    left.entries()
        .iter()
        .zip(right.entries())
        .map(|(left, right)| compare_element(left, right))
        .find(|order| *order != Ordering::Equal)
        .unwrap_or(Ordering::Equal)
}

fn compare_spin_space(left: &SpinSpaceOperation, right: &SpinSpaceOperation) -> Ordering {
    compare_matrix(left.spatial.rotation(), right.spatial.rotation())
        .then_with(|| {
            left.spatial
                .translation()
                .iter()
                .zip(right.spatial.translation())
                .map(|(left, right)| compare_element(left, right))
                .find(|order| *order != Ordering::Equal)
                .unwrap_or(Ordering::Equal)
        })
        .then_with(|| compare_matrix(&left.spin_rotation, &right.spin_rotation))
        .then_with(|| left.antiunitary.cmp(&right.antiunitary))
}

/// Expands exact Seitz generators with the stable Mathematica layer ordering
/// and compiles the resulting ordered quotient group.
pub fn generate_seitz_group(generators: &[SeitzOperation]) -> SymmetryResult<CompiledSeitzGroup> {
    let context = generators
        .first()
        .ok_or_else(|| SymmetryError::new("EmptySeitzGenerators", "generator list is empty"))?
        .spatial()
        .rotation()
        .context();
    let dimension = generators[0].spatial().rotation().rows();
    let identity = SeitzOperation::new(
        "generated-identity",
        "GeneratedIdentity",
        SpatialOperation::new(
            ExactMatrix::identity(context, dimension)?,
            vec![Element::zero(context)?; dimension],
        )?,
        false,
    );
    let generated = generate_group_by(
        generators,
        &identity,
        |left, right| -> SymmetryResult<SeitzOperation> {
            let product = left.product(right)?;
            Ok(SeitzOperation::new(
                "",
                "",
                normalize_spatial_mod_lattice(product.spatial())?,
                product.antiunitary(),
            ))
        },
        |left, right| left.equivalent_mod_lattice(right).unwrap_or(false),
    )?;
    let operations = generated
        .into_iter()
        .enumerate()
        .map(|(index, operation)| {
            SeitzOperation::new(
                format!("generated-{index}"),
                format!("Generated{}", index + 1),
                operation.spatial,
                operation.antiunitary,
            )
        })
        .collect();
    compile_seitz_group(operations)
}
