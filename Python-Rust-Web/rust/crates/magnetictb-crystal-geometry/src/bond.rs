use std::cmp::Ordering;
use std::collections::BTreeMap;

use cyclotomic_nullspace::Rational;
use magnetictb_data::ExactExpression;

use crate::{GeometryError, GeometryResult};

/// Exact element of Q(sqrt(2), sqrt(3)), stored in the ordered basis
/// 1, sqrt(2), sqrt(3), sqrt(6).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct QuadraticReal {
    coefficients: [Rational; 4],
}

impl QuadraticReal {
    #[must_use]
    pub fn zero() -> Self {
        Self::from_rational(Rational::zero())
    }

    #[must_use]
    pub fn one() -> Self {
        Self::from_rational(Rational::one())
    }

    #[must_use]
    pub fn from_rational(value: Rational) -> Self {
        Self {
            coefficients: [value, Rational::zero(), Rational::zero(), Rational::zero()],
        }
    }

    #[must_use]
    pub fn sqrt_two() -> Self {
        Self {
            coefficients: [
                Rational::zero(),
                Rational::one(),
                Rational::zero(),
                Rational::zero(),
            ],
        }
    }

    #[must_use]
    pub fn sqrt_three() -> Self {
        Self {
            coefficients: [
                Rational::zero(),
                Rational::zero(),
                Rational::one(),
                Rational::zero(),
            ],
        }
    }

    #[must_use]
    pub fn coefficients(&self) -> &[Rational; 4] {
        &self.coefficients
    }

    #[must_use]
    pub fn as_rational(&self) -> Option<&Rational> {
        self.coefficients[1..]
            .iter()
            .all(Rational::is_zero)
            .then_some(&self.coefficients[0])
    }

    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.coefficients.iter().all(Rational::is_zero)
    }

    #[must_use]
    pub fn add(&self, right: &Self) -> Self {
        Self {
            coefficients: std::array::from_fn(|index| {
                self.coefficients[index].add(&right.coefficients[index])
            }),
        }
    }

    #[must_use]
    pub fn negate(&self) -> Self {
        Self {
            coefficients: std::array::from_fn(|index| self.coefficients[index].negate()),
        }
    }

    #[must_use]
    pub fn subtract(&self, right: &Self) -> Self {
        self.add(&right.negate())
    }

    #[must_use]
    pub fn scale(&self, scalar: &Rational) -> Self {
        Self {
            coefficients: std::array::from_fn(|index| self.coefficients[index].multiply(scalar)),
        }
    }

    #[must_use]
    pub fn multiply(&self, right: &Self) -> Self {
        let mut result = Self::zero();
        for left_basis in 0..4 {
            for right_basis in 0..4 {
                let common = left_basis & right_basis;
                let mut factor =
                    self.coefficients[left_basis].multiply(&right.coefficients[right_basis]);
                if common & 1 != 0 {
                    factor = factor.multiply(&Rational::from_i64(2));
                }
                if common & 2 != 0 {
                    factor = factor.multiply(&Rational::from_i64(3));
                }
                let basis = left_basis ^ right_basis;
                result.coefficients[basis] = result.coefficients[basis].add(&factor);
            }
        }
        result
    }

    /// Convert an exact quadratic real to the nearest finite binary64 value.
    ///
    /// This is an explicit numerical boundary used by geometry and by
    /// numerical evaluation of an already verified symbolic Hamiltonian.
    pub fn to_f64(&self) -> GeometryResult<f64> {
        let basis = [1.0, 2.0_f64.sqrt(), 3.0_f64.sqrt(), 6.0_f64.sqrt()];
        let mut result = 0.0;
        for (coefficient, basis_value) in self.coefficients.iter().zip(basis) {
            let value = coefficient.to_f64().ok_or_else(|| {
                GeometryError::new(
                    "NonFiniteBondGeometry",
                    "an exact bond coordinate is outside the finite f64 geometry range",
                )
            })?;
            result += value * basis_value;
        }
        if result.is_finite() {
            Ok(result)
        } else {
            Err(GeometryError::new(
                "NonFiniteBondGeometry",
                "an exact bond coordinate evaluated to a non-finite f64 value",
            ))
        }
    }

    fn sign(&self) -> Ordering {
        let first = QuadraticPair {
            rational: self.coefficients[0].clone(),
            sqrt_two: self.coefficients[1].clone(),
        };
        let second = QuadraticPair {
            rational: self.coefficients[2].clone(),
            sqrt_two: self.coefficients[3].clone(),
        };
        let first_sign = first.sign();
        let second_sign = second.sign();
        if first_sign == Ordering::Equal {
            return second_sign;
        }
        if second_sign == Ordering::Equal || first_sign == second_sign {
            return first_sign;
        }
        let difference = first.square().subtract(&second.square().scale_i64(3));
        if first_sign == Ordering::Greater {
            difference.sign()
        } else {
            difference.sign().reverse()
        }
    }

    pub fn from_expression(
        expression: &ExactExpression,
        bindings: &BTreeMap<String, Self>,
    ) -> GeometryResult<Self> {
        match expression {
            ExactExpression::Integer { value } => Ok(Self::from_rational(
                Rational::parse(value, "1").map_err(GeometryError::from)?,
            )),
            ExactExpression::Rational {
                numerator,
                denominator,
            } => Ok(Self::from_rational(
                Rational::parse(numerator, denominator).map_err(GeometryError::from)?,
            )),
            ExactExpression::Symbol { name } => bindings.get(name).cloned().ok_or_else(|| {
                GeometryError::new(
                    "MissingExactSymbol",
                    format!("no quadratic-real binding was supplied for {name}"),
                )
            }),
            ExactExpression::Call { head, arguments } if head == "System`Plus" => {
                let mut result = Self::zero();
                for argument in arguments {
                    result = result.add(&Self::from_expression(argument, bindings)?);
                }
                Ok(result)
            }
            ExactExpression::Call { head, arguments } if head == "System`Times" => {
                let mut result = Self::one();
                for argument in arguments {
                    result = result.multiply(&Self::from_expression(argument, bindings)?);
                }
                Ok(result)
            }
            ExactExpression::Call { head, arguments }
                if head == "System`Power" && arguments.len() == 2 =>
            {
                canonical_square_root_power(&arguments[0], &arguments[1])
            }
            ExactExpression::RootOfUnity { .. } | ExactExpression::Call { .. } => {
                Err(GeometryError::new(
                    "UnsupportedBondExactExpression",
                    "bond geometry supports rational values and canonical sqrt(2)/sqrt(3) expressions",
                ))
            }
        }
    }
}

impl Ord for QuadraticReal {
    fn cmp(&self, other: &Self) -> Ordering {
        self.subtract(other).sign()
    }
}

impl PartialOrd for QuadraticReal {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Clone)]
struct QuadraticPair {
    rational: Rational,
    sqrt_two: Rational,
}

impl QuadraticPair {
    fn sign(&self) -> Ordering {
        let rational_sign = self.rational.cmp(&Rational::zero());
        let radical_sign = self.sqrt_two.cmp(&Rational::zero());
        if rational_sign == Ordering::Equal {
            return radical_sign;
        }
        if radical_sign == Ordering::Equal || rational_sign == radical_sign {
            return rational_sign;
        }
        let difference = self.rational.multiply(&self.rational).subtract(
            &self
                .sqrt_two
                .multiply(&self.sqrt_two)
                .multiply(&Rational::from_i64(2)),
        );
        if rational_sign == Ordering::Greater {
            difference.cmp(&Rational::zero())
        } else {
            difference.cmp(&Rational::zero()).reverse()
        }
    }

    fn square(&self) -> Self {
        Self {
            rational: self.rational.multiply(&self.rational).add(
                &self
                    .sqrt_two
                    .multiply(&self.sqrt_two)
                    .multiply(&Rational::from_i64(2)),
            ),
            sqrt_two: self
                .rational
                .multiply(&self.sqrt_two)
                .multiply(&Rational::from_i64(2)),
        }
    }

    fn subtract(&self, right: &Self) -> Self {
        Self {
            rational: self.rational.subtract(&right.rational),
            sqrt_two: self.sqrt_two.subtract(&right.sqrt_two),
        }
    }

    fn scale_i64(&self, scalar: i64) -> Self {
        Self {
            rational: self.rational.multiply(&Rational::from_i64(scalar)),
            sqrt_two: self.sqrt_two.multiply(&Rational::from_i64(scalar)),
        }
    }
}

fn canonical_square_root_power(
    base: &ExactExpression,
    exponent: &ExactExpression,
) -> GeometryResult<QuadraticReal> {
    let ExactExpression::Integer { value: radicand } = base else {
        return Err(GeometryError::new(
            "UnsupportedBondExactExpression",
            "canonical radical base must be integer 2 or 3",
        ));
    };
    let (numerator, denominator) = match exponent {
        ExactExpression::Rational {
            numerator,
            denominator,
        } => (numerator.as_str(), denominator.as_str()),
        _ => {
            return Err(GeometryError::new(
                "UnsupportedBondExactExpression",
                "canonical radical exponent must be 1/2 or -1/2",
            ));
        }
    };
    if denominator != "2" || !matches!(numerator, "1" | "-1") {
        return Err(GeometryError::new(
            "UnsupportedBondExactExpression",
            "canonical radical exponent must be 1/2 or -1/2",
        ));
    }
    let root = match radicand.as_str() {
        "2" => QuadraticReal::sqrt_two(),
        "3" => QuadraticReal::sqrt_three(),
        _ => {
            return Err(GeometryError::new(
                "UnsupportedBondExactExpression",
                "canonical radical base must be integer 2 or 3",
            ));
        }
    };
    if numerator == "1" {
        Ok(root)
    } else {
        Ok(root.scale(&Rational::parse("1", radicand).map_err(GeometryError::from)?))
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PeriodicBondRecord {
    source_site: usize,
    target_site: usize,
    translation: Vec<i64>,
    displacement: Vec<QuadraticReal>,
    source_endpoint: Vec<QuadraticReal>,
    target_endpoint: Vec<QuadraticReal>,
    squared_distance: QuadraticReal,
}

impl PeriodicBondRecord {
    #[must_use]
    pub const fn source_site(&self) -> usize {
        self.source_site
    }

    #[must_use]
    pub const fn target_site(&self) -> usize {
        self.target_site
    }

    #[must_use]
    pub fn translation(&self) -> &[i64] {
        &self.translation
    }

    #[must_use]
    pub fn displacement(&self) -> &[QuadraticReal] {
        &self.displacement
    }

    #[must_use]
    pub fn source_endpoint(&self) -> &[QuadraticReal] {
        &self.source_endpoint
    }

    #[must_use]
    pub fn target_endpoint(&self) -> &[QuadraticReal] {
        &self.target_endpoint
    }

    #[must_use]
    pub const fn squared_distance(&self) -> &QuadraticReal {
        &self.squared_distance
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct PeriodicBondShell {
    squared_distance: QuadraticReal,
    numeric_distance: f64,
    bonds: Vec<PeriodicBondRecord>,
}

impl PeriodicBondShell {
    #[must_use]
    pub const fn squared_distance(&self) -> &QuadraticReal {
        &self.squared_distance
    }

    #[must_use]
    pub const fn numeric_distance(&self) -> f64 {
        self.numeric_distance
    }

    #[must_use]
    pub fn bonds(&self) -> &[PeriodicBondRecord] {
        &self.bonds
    }
}

/// Stable numerical options for adaptive periodic-bond search.
#[derive(Clone, Debug, PartialEq)]
pub struct BondSearchOptions {
    pub minimum_neighbors_per_site: usize,
    /// `None` is Mathematica's `Automatic`.
    pub initial_radius: Option<f64>,
    pub growth_factor: f64,
    pub maximum_iterations: usize,
    /// `None` is Mathematica's `Infinity`.
    pub maximum_radius: Option<f64>,
    pub distance_tolerance: f64,
}

impl Default for BondSearchOptions {
    fn default() -> Self {
        Self {
            minimum_neighbors_per_site: 0,
            initial_radius: None,
            growth_factor: 5.0 / 4.0,
            maximum_iterations: 24,
            maximum_radius: None,
            distance_tolerance: 1.0e-10,
        }
    }
}

/// Result and stopping evidence for the stable adaptive numerical search.
#[derive(Clone, Debug, PartialEq)]
pub struct PeriodicBondSearchResult {
    shells: Vec<PeriodicBondShell>,
    radius: f64,
    neighbor_counts_per_site: Vec<usize>,
    candidate_bond_count: usize,
    enumerated_translation_count: usize,
    requested_shells: usize,
    minimum_neighbors_per_site: usize,
    iterations: usize,
}

impl PeriodicBondSearchResult {
    #[must_use]
    pub fn shells(&self) -> &[PeriodicBondShell] {
        &self.shells
    }

    #[must_use]
    pub const fn radius(&self) -> f64 {
        self.radius
    }

    #[must_use]
    pub fn neighbor_counts_per_site(&self) -> &[usize] {
        &self.neighbor_counts_per_site
    }

    #[must_use]
    pub const fn candidate_bond_count(&self) -> usize {
        self.candidate_bond_count
    }

    #[must_use]
    pub const fn enumerated_translation_count(&self) -> usize {
        self.enumerated_translation_count
    }

    #[must_use]
    pub const fn requested_shells(&self) -> usize {
        self.requested_shells
    }

    #[must_use]
    pub const fn minimum_neighbors_per_site(&self) -> usize {
        self.minimum_neighbors_per_site
    }

    #[must_use]
    pub const fn iterations(&self) -> usize {
        self.iterations
    }
}

struct NumericBondRecord {
    squared_distance: f64,
    exact: PeriodicBondRecord,
}

struct RadiusSearchData {
    shells: Vec<PeriodicBondShell>,
    neighbor_counts_per_site: Vec<usize>,
    candidate_bond_count: usize,
    enumerated_translation_count: usize,
}

fn validate_periodic_bond_inputs(
    sites: &[Vec<QuadraticReal>],
    lattice: &[Vec<QuadraticReal>],
) -> GeometryResult<(usize, usize)> {
    let dimension = sites.first().map_or(0, Vec::len);
    let cartesian_dimension = lattice.first().map_or(0, Vec::len);
    if dimension == 0
        || cartesian_dimension < dimension
        || sites.iter().any(|site| site.len() != dimension)
        || lattice.len() != dimension
        || lattice
            .iter()
            .any(|vector| vector.len() != cartesian_dimension)
    {
        return Err(GeometryError::new(
            "InvalidPeriodicBondInput",
            "sites and lattice must define aligned nonempty real fractional geometry",
        ));
    }
    Ok((dimension, cartesian_dimension))
}

fn numeric_rows(rows: &[Vec<QuadraticReal>]) -> GeometryResult<Vec<Vec<f64>>> {
    rows.iter()
        .map(|row| row.iter().map(QuadraticReal::to_f64).collect())
        .collect()
}

fn gram_matrix(lattice: &[Vec<f64>]) -> Vec<Vec<f64>> {
    lattice
        .iter()
        .map(|left| {
            lattice
                .iter()
                .map(|right| left.iter().zip(right).map(|(a, b)| a * b).sum())
                .collect()
        })
        .collect()
}

/// Upper-triangular factor `U` satisfying `U^T U = matrix`.
fn cholesky_upper(matrix: &[Vec<f64>]) -> GeometryResult<Vec<Vec<f64>>> {
    let dimension = matrix.len();
    let mut upper = vec![vec![0.0; dimension]; dimension];
    for row in 0..dimension {
        for column in row..dimension {
            let prior = (0..row)
                .map(|index| upper[index][row] * upper[index][column])
                .sum::<f64>();
            let residual = matrix[row][column] - prior;
            if row == column {
                if !residual.is_finite() || residual <= 0.0 {
                    return Err(GeometryError::new(
                        "InvalidPeriodicBondInput",
                        "the numeric lattice rows must be linearly independent",
                    ));
                }
                upper[row][column] = residual.sqrt();
            } else {
                upper[row][column] = residual / upper[row][row];
            }
        }
    }
    Ok(upper)
}

const MAX_EXACT_F64_INTEGER: f64 = 9_007_199_254_740_992.0;
const MAX_EXACT_F64_INTEGER_I64: i64 = 9_007_199_254_740_992;

#[allow(clippy::cast_possible_truncation)]
fn bounded_integral_f64_to_i64(value: f64) -> i64 {
    debug_assert!(value.is_finite());
    debug_assert!(value.abs() <= MAX_EXACT_F64_INTEGER);
    value as i64
}

#[allow(clippy::cast_precision_loss)]
fn exact_i64_to_f64(value: i64) -> GeometryResult<f64> {
    if !(-MAX_EXACT_F64_INTEGER_I64..=MAX_EXACT_F64_INTEGER_I64).contains(&value) {
        return Err(GeometryError::new(
            "BondTranslationOverflow",
            "a bond translation is outside the exactly representable f64 integer range",
        ));
    }
    Ok(value as f64)
}

fn integer_translations_in_ball(
    shift: &[f64],
    triangular: &[Vec<f64>],
    radius: f64,
    tolerance: f64,
) -> GeometryResult<Vec<Vec<i64>>> {
    struct Enumerator<'a> {
        shift: &'a [f64],
        triangular: &'a [Vec<f64>],
        tolerance: f64,
        absolute_tolerance: f64,
        vector: Vec<i64>,
        shifted_vector: Vec<f64>,
        result: Vec<Vec<i64>>,
    }

    impl Enumerator<'_> {
        fn recurse(&mut self, level: usize, remaining: f64) -> GeometryResult<()> {
            if level == 0 {
                self.result.push(self.vector.clone());
                return Ok(());
            }
            let index = level - 1;
            let tail = ((index + 1)..self.shift.len())
                .map(|column| self.triangular[index][column] * self.shifted_vector[column])
                .sum::<f64>();
            let diagonal = self.triangular[index][index];
            let center = -self.shift[index] - tail / diagonal;
            let half_width = (remaining + self.absolute_tolerance).max(0.0).sqrt() / diagonal.abs();
            let lower_value = (center - half_width - 10.0 * self.tolerance).ceil();
            let upper_value = (center + half_width + 10.0 * self.tolerance).floor();
            if !lower_value.is_finite()
                || !upper_value.is_finite()
                || lower_value < -MAX_EXACT_F64_INTEGER
                || upper_value > MAX_EXACT_F64_INTEGER
            {
                return Err(GeometryError::new(
                    "BondTranslationOverflow",
                    "a numerical bond-search translation is outside the i64 range",
                ));
            }
            let lower = bounded_integral_f64_to_i64(lower_value);
            let upper = bounded_integral_f64_to_i64(upper_value);
            if lower > upper {
                return Ok(());
            }
            for integer in lower..=upper {
                self.vector[index] = integer;
                self.shifted_vector[index] = exact_i64_to_f64(integer)? + self.shift[index];
                let term = diagonal * self.shifted_vector[index] + tail;
                let next_remaining = remaining - term * term;
                if next_remaining >= -self.absolute_tolerance {
                    self.recurse(index, next_remaining.max(0.0))?;
                }
            }
            Ok(())
        }
    }

    let radius_squared = radius * radius;
    let mut enumerator = Enumerator {
        shift,
        triangular,
        tolerance,
        absolute_tolerance: tolerance * radius_squared.max(1.0),
        vector: vec![0; shift.len()],
        shifted_vector: vec![0.0; shift.len()],
        result: Vec::new(),
    };
    enumerator.recurse(shift.len(), radius_squared)?;
    Ok(enumerator.result)
}

fn numeric_squared_distance(
    displacement: &[f64],
    lattice: &[Vec<f64>],
    cartesian_dimension: usize,
) -> f64 {
    (0..cartesian_dimension)
        .map(|component| {
            displacement
                .iter()
                .zip(lattice)
                .map(|(coordinate, row)| coordinate * row[component])
                .sum::<f64>()
        })
        .map(|coordinate| coordinate * coordinate)
        .sum()
}

fn exact_bond_record(
    source_site: usize,
    target_site: usize,
    source: &[QuadraticReal],
    target: &[QuadraticReal],
    translation: &[i64],
    lattice: &[Vec<QuadraticReal>],
    cartesian_dimension: usize,
) -> PeriodicBondRecord {
    let target_endpoint = target
        .iter()
        .zip(translation)
        .map(|(coordinate, &integer)| {
            coordinate.add(&QuadraticReal::from_rational(Rational::from_i64(integer)))
        })
        .collect::<Vec<_>>();
    let displacement = target_endpoint
        .iter()
        .zip(source)
        .map(|(target, source)| target.subtract(source))
        .collect::<Vec<_>>();
    let mut cartesian = vec![QuadraticReal::zero(); cartesian_dimension];
    for fractional in 0..displacement.len() {
        for component in 0..cartesian_dimension {
            cartesian[component] = cartesian[component]
                .add(&displacement[fractional].multiply(&lattice[fractional][component]));
        }
    }
    let squared_distance = cartesian
        .iter()
        .fold(QuadraticReal::zero(), |total, value| {
            total.add(&value.multiply(value))
        });
    PeriodicBondRecord {
        source_site,
        target_site,
        translation: translation.to_vec(),
        displacement,
        source_endpoint: source.to_vec(),
        target_endpoint,
        squared_distance,
    }
}

#[allow(clippy::too_many_lines)]
fn bond_search_data_within_radius(
    sites: &[Vec<QuadraticReal>],
    lattice: &[Vec<QuadraticReal>],
    radius: f64,
    tolerance: f64,
) -> GeometryResult<RadiusSearchData> {
    let (dimension, cartesian_dimension) = validate_periodic_bond_inputs(sites, lattice)?;
    if !radius.is_finite() || radius < 0.0 || !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(GeometryError::new(
            "InvalidPeriodicBondInput",
            "radius must be finite and nonnegative and distance tolerance must be positive",
        ));
    }
    let numeric_sites = numeric_rows(sites)?;
    let numeric_lattice = numeric_rows(lattice)?;
    let triangular = cholesky_upper(&gram_matrix(&numeric_lattice))?;
    let radius_squared = radius * radius;
    let absolute_tolerance = tolerance * radius_squared.max(1.0);
    let mut records = Vec::new();
    let mut enumerated_translation_count = 0_usize;
    for (source_site, source) in sites.iter().enumerate() {
        for (target_site, target) in sites.iter().enumerate() {
            let shift = numeric_sites[target_site]
                .iter()
                .zip(&numeric_sites[source_site])
                .map(|(target, source)| target - source)
                .collect::<Vec<_>>();
            let translations =
                integer_translations_in_ball(&shift, &triangular, radius, tolerance)?;
            enumerated_translation_count = enumerated_translation_count
                .checked_add(translations.len())
                .ok_or_else(|| {
                    GeometryError::new(
                        "BondCandidateOverflow",
                        "the numerical bond-search candidate counter overflowed",
                    )
                })?;
            for translation in translations {
                let displacement = shift
                    .iter()
                    .zip(&translation)
                    .map(|(coordinate, &integer)| {
                        exact_i64_to_f64(integer).map(|integer| coordinate + integer)
                    })
                    .collect::<GeometryResult<Vec<_>>>()?;
                let squared_distance =
                    numeric_squared_distance(&displacement, &numeric_lattice, cartesian_dimension);
                if squared_distance <= radius_squared + absolute_tolerance {
                    records.push(NumericBondRecord {
                        squared_distance,
                        exact: exact_bond_record(
                            source_site,
                            target_site,
                            source,
                            target,
                            &translation,
                            lattice,
                            cartesian_dimension,
                        ),
                    });
                }
            }
        }
    }
    records.sort_by(|left, right| {
        left.squared_distance
            .total_cmp(&right.squared_distance)
            .then_with(|| left.exact.source_site.cmp(&right.exact.source_site))
            .then_with(|| left.exact.target_site.cmp(&right.exact.target_site))
            .then_with(|| left.exact.target_endpoint.cmp(&right.exact.target_endpoint))
    });

    let candidate_bond_count = records.len();
    let mut neighbor_counts_per_site = vec![0; sites.len()];
    for record in &records {
        if record.squared_distance > absolute_tolerance {
            neighbor_counts_per_site[record.exact.source_site] += 1;
        }
    }

    let mut raw_shells: Vec<(f64, Vec<PeriodicBondRecord>)> = Vec::new();
    for record in records {
        if let Some((reference, bonds)) = raw_shells.last_mut()
            && (record.squared_distance - *reference).abs()
                <= tolerance
                    * 1.0_f64
                        .max(record.squared_distance.abs())
                        .max(reference.abs())
        {
            bonds.push(record.exact);
        } else {
            raw_shells.push((record.squared_distance, vec![record.exact]));
        }
    }
    let shells = raw_shells
        .into_iter()
        .map(|(squared_distance, mut bonds)| {
            let exact_representative = bonds[0].squared_distance.clone();
            // Stable BondSearch first uses the fully sorted records to split
            // shells, then groups by source and sorts each source's endpoint
            // pairs by the target endpoint.
            bonds.sort_by(|left, right| {
                left.source_site
                    .cmp(&right.source_site)
                    .then_with(|| left.target_endpoint.cmp(&right.target_endpoint))
                    .then_with(|| left.target_site.cmp(&right.target_site))
            });
            PeriodicBondShell {
                squared_distance: exact_representative,
                numeric_distance: squared_distance.max(0.0).sqrt(),
                bonds,
            }
        })
        .collect::<Vec<_>>();
    debug_assert_eq!(dimension, triangular.len());
    Ok(RadiusSearchData {
        shells,
        neighbor_counts_per_site,
        candidate_bond_count,
        enumerated_translation_count,
    })
}

fn initial_bond_search_radius(
    sites: &[Vec<QuadraticReal>],
    lattice: &[Vec<QuadraticReal>],
    tolerance: f64,
) -> GeometryResult<f64> {
    let (dimension, cartesian_dimension) = validate_periodic_bond_inputs(sites, lattice)?;
    let numeric_sites = numeric_rows(sites)?;
    let numeric_lattice = numeric_rows(lattice)?;
    let mut distances = integer_box(dimension, 1)
        .into_iter()
        .filter(|translation| translation.iter().any(|&value| value != 0))
        .map(|translation| {
            let displacement = translation
                .into_iter()
                .map(exact_i64_to_f64)
                .collect::<GeometryResult<Vec<_>>>()?;
            Ok(
                numeric_squared_distance(&displacement, &numeric_lattice, cartesian_dimension)
                    .sqrt(),
            )
        })
        .collect::<GeometryResult<Vec<_>>>()?;
    for source in &numeric_sites {
        for target in &numeric_sites {
            distances.push(
                numeric_squared_distance(
                    &target
                        .iter()
                        .zip(source)
                        .map(|(target, source)| target - source)
                        .collect::<Vec<_>>(),
                    &numeric_lattice,
                    cartesian_dimension,
                )
                .sqrt(),
            );
        }
    }
    Ok(distances
        .into_iter()
        .filter(|distance| *distance > 10.0 * tolerance)
        .min_by(f64::total_cmp)
        .unwrap_or(1.0))
}

pub fn find_periodic_bond_shells(
    sites: &[Vec<QuadraticReal>],
    lattice: &[Vec<QuadraticReal>],
    requested_shells: usize,
    options: &BondSearchOptions,
) -> GeometryResult<PeriodicBondSearchResult> {
    validate_periodic_bond_inputs(sites, lattice)?;
    if requested_shells == 0
        || options
            .initial_radius
            .is_some_and(|radius| !radius.is_finite() || radius < 0.0)
        || !options.growth_factor.is_finite()
        || options.growth_factor <= 1.0
        || options.maximum_iterations == 0
        || options
            .maximum_radius
            .is_some_and(|radius| !radius.is_finite() || radius <= 0.0)
        || !options.distance_tolerance.is_finite()
        || options.distance_tolerance <= 0.0
    {
        return Err(GeometryError::new(
            "InvalidPeriodicBondInput",
            "requested shells and adaptive bond-search options are invalid",
        ));
    }
    let mut radius = if let Some(radius) = options.initial_radius {
        radius
    } else if requested_shells == 1 && options.minimum_neighbors_per_site == 0 {
        0.0
    } else {
        initial_bond_search_radius(sites, lattice, options.distance_tolerance)?
    };
    if let Some(maximum_radius) = options.maximum_radius {
        radius = radius.min(maximum_radius);
    }

    let mut last_data = None;
    let mut last_iteration = 0;
    for iteration in 1..=options.maximum_iterations {
        last_iteration = iteration;
        let data =
            bond_search_data_within_radius(sites, lattice, radius, options.distance_tolerance)?;
        let success = data.shells.len() >= requested_shells
            && data
                .neighbor_counts_per_site
                .iter()
                .copied()
                .min()
                .is_some_and(|count| count >= options.minimum_neighbors_per_site);
        if success {
            return Ok(PeriodicBondSearchResult {
                shells: data.shells,
                radius,
                neighbor_counts_per_site: data.neighbor_counts_per_site,
                candidate_bond_count: data.candidate_bond_count,
                enumerated_translation_count: data.enumerated_translation_count,
                requested_shells,
                minimum_neighbors_per_site: options.minimum_neighbors_per_site,
                iterations: iteration,
            });
        }
        last_data = Some(data);
        if options
            .maximum_radius
            .is_some_and(|maximum_radius| radius >= maximum_radius)
        {
            break;
        }
        radius = options.growth_factor * radius.max(10.0 * options.distance_tolerance.sqrt());
        if let Some(maximum_radius) = options.maximum_radius {
            radius = radius.min(maximum_radius);
        }
    }
    let Some(data) = last_data else {
        return Err(GeometryError::new(
            "SearchLimitReached",
            "adaptive search did not execute an iteration",
        ));
    };
    Err(GeometryError::new(
        "SearchLimitReached",
        format!(
            "adaptive search reached {} shells and {:?} neighbors per site at radius {radius} after {last_iteration} iterations",
            data.shells.len(),
            data.neighbor_counts_per_site
        ),
    ))
}

pub fn periodic_bond_shells_in_box(
    sites: &[Vec<QuadraticReal>],
    lattice: &[Vec<QuadraticReal>],
    translation_bound: i64,
) -> GeometryResult<Vec<PeriodicBondShell>> {
    let (dimension, cartesian_dimension) = validate_periodic_bond_inputs(sites, lattice)?;
    if translation_bound < 0 {
        return Err(GeometryError::new(
            "InvalidPeriodicBondInput",
            "sites, lattice, and explicit translation bound do not align",
        ));
    }
    let translations = integer_box(dimension, translation_bound);
    let mut records = Vec::new();
    for (source_site, source) in sites.iter().enumerate() {
        for (target_site, target) in sites.iter().enumerate() {
            for translation in &translations {
                records.push(exact_bond_record(
                    source_site,
                    target_site,
                    source,
                    target,
                    translation,
                    lattice,
                    cartesian_dimension,
                ));
            }
        }
    }
    records.sort_by(|left, right| {
        left.squared_distance
            .cmp(&right.squared_distance)
            .then_with(|| left.source_site.cmp(&right.source_site))
            .then_with(|| left.target_site.cmp(&right.target_site))
            .then_with(|| left.target_endpoint.cmp(&right.target_endpoint))
    });
    let mut shells: Vec<PeriodicBondShell> = Vec::new();
    for record in records {
        if let Some(shell) = shells.last_mut()
            && shell.squared_distance == record.squared_distance
        {
            shell.bonds.push(record);
        } else {
            shells.push(PeriodicBondShell {
                squared_distance: record.squared_distance.clone(),
                numeric_distance: record.squared_distance.to_f64()?.max(0.0).sqrt(),
                bonds: vec![record],
            });
        }
    }
    Ok(shells)
}

fn integer_box(dimension: usize, bound: i64) -> Vec<Vec<i64>> {
    fn append(result: &mut Vec<Vec<i64>>, current: &mut Vec<i64>, dimension: usize, bound: i64) {
        if current.len() == dimension {
            result.push(current.clone());
            return;
        }
        for value in -bound..=bound {
            current.push(value);
            append(result, current, dimension, bound);
            current.pop();
        }
    }
    let mut result = Vec::new();
    append(
        &mut result,
        &mut Vec::with_capacity(dimension),
        dimension,
        bound,
    );
    result
}

#[cfg(test)]
mod tests {
    use super::{
        BondSearchOptions, QuadraticReal, bond_search_data_within_radius, exact_i64_to_f64,
        find_periodic_bond_shells, periodic_bond_shells_in_box,
    };
    use cyclotomic_nullspace::Rational;

    fn rational(numerator: i64, denominator: i64) -> QuadraticReal {
        QuadraticReal::from_rational(
            Rational::parse(&numerator.to_string(), &denominator.to_string()).expect("rational"),
        )
    }

    #[test]
    fn exact_quadratic_order_distinguishes_radicals_without_floats() {
        assert!(QuadraticReal::sqrt_two() > rational(7, 5));
        assert!(QuadraticReal::sqrt_three() < rational(7, 4));
        assert_eq!(
            QuadraticReal::sqrt_two().multiply(&QuadraticReal::sqrt_three()),
            QuadraticReal::sqrt_two().multiply(&QuadraticReal::sqrt_three())
        );
    }

    #[test]
    fn bounded_one_dimensional_bonds_form_deterministic_shells() {
        let sites = vec![vec![rational(0, 1)], vec![rational(1, 2)]];
        let lattice = vec![vec![rational(1, 1)]];
        let shells = periodic_bond_shells_in_box(&sites, &lattice, 1).expect("bond shells");
        assert_eq!(shells[0].bonds().len(), 2);
        assert_eq!(shells[1].squared_distance(), &rational(1, 4));
        assert_eq!(shells[1].bonds().len(), 4);
        assert_eq!(shells[1].bonds()[0].translation(), [-1]);
        assert_eq!(shells[1].bonds()[3].translation(), [1]);
    }

    #[test]
    fn adaptive_zero_radius_and_nonzero_neighbor_stop_match_stable_semantics() {
        let sites = vec![vec![rational(0, 1)]];
        let lattice = vec![vec![rational(1, 1)]];
        let onsite = find_periodic_bond_shells(&sites, &lattice, 1, &BondSearchOptions::default())
            .expect("onsite shell");
        assert!(onsite.radius().abs() < f64::EPSILON);
        assert_eq!(onsite.iterations(), 1);
        assert_eq!(onsite.shells().len(), 1);
        assert_eq!(onsite.shells()[0].bonds().len(), 1);
        assert_eq!(onsite.neighbor_counts_per_site(), [0]);

        let neighbors = find_periodic_bond_shells(
            &sites,
            &lattice,
            1,
            &BondSearchOptions {
                minimum_neighbors_per_site: 1,
                ..BondSearchOptions::default()
            },
        )
        .expect("nonzero neighbors");
        assert!((neighbors.radius() - 1.0).abs() < f64::EPSILON);
        assert_eq!(neighbors.neighbor_counts_per_site(), [2]);
        assert_eq!(neighbors.candidate_bond_count(), 3);
    }

    #[test]
    fn fincke_pohst_nonorthogonal_ball_matches_small_box_and_avoids_box_candidates() {
        let sites = vec![vec![rational(0, 1), rational(0, 1)]];
        let lattice = vec![
            vec![rational(1, 1), rational(0, 1)],
            vec![
                rational(1, 2),
                QuadraticReal::sqrt_three().scale(&Rational::parse("1", "2").expect("half")),
            ],
        ];
        let numeric = bond_search_data_within_radius(&sites, &lattice, 1.01, 1.0e-10)
            .expect("nonorthogonal numeric ball");
        assert_eq!(numeric.candidate_bond_count, 7);
        assert_eq!(numeric.enumerated_translation_count, 7);
        assert_eq!(numeric.shells.len(), 2);
        assert_eq!(numeric.shells[0].bonds().len(), 1);
        assert_eq!(numeric.shells[1].bonds().len(), 6);

        let diagnostic = periodic_bond_shells_in_box(&sites, &lattice, 1).expect("box");
        let diagnostic_nearest = diagnostic
            .iter()
            .find(|shell| shell.squared_distance() == &rational(1, 1))
            .expect("nearest box shell");
        assert_eq!(
            numeric.shells[1]
                .bonds()
                .iter()
                .map(|bond| bond.translation().to_vec())
                .collect::<Vec<_>>(),
            diagnostic_nearest
                .bonds()
                .iter()
                .map(|bond| bond.translation().to_vec())
                .collect::<Vec<_>>()
        );
        assert!(numeric.candidate_bond_count < 3_usize.pow(2));
    }

    #[test]
    fn numerical_shell_tolerance_uses_first_sorted_distance_as_reference() {
        let epsilon = rational(1, 100_000_000);
        let sites = vec![
            vec![rational(0, 1)],
            vec![rational(1, 2)],
            vec![rational(1, 2).add(&epsilon)],
        ];
        let lattice = vec![vec![rational(1, 1)]];
        let within = bond_search_data_within_radius(&sites, &lattice, 0.6, 1.0e-6)
            .expect("within tolerance");
        let outside = bond_search_data_within_radius(&sites, &lattice, 0.6, 1.0e-10)
            .expect("outside tolerance");
        assert!(within.shells.len() < outside.shells.len());
        assert!(within.shells[0].numeric_distance().abs() < f64::EPSILON);
        assert!(
            outside
                .shells
                .windows(2)
                .all(|pair| pair[0].numeric_distance() <= pair[1].numeric_distance())
        );
    }

    #[test]
    fn adaptive_search_limit_is_explicit() {
        let error = find_periodic_bond_shells(
            &[vec![rational(0, 1)]],
            &[vec![rational(1, 1)]],
            2,
            &BondSearchOptions {
                initial_radius: Some(0.0),
                maximum_iterations: 1,
                maximum_radius: Some(0.1),
                ..BondSearchOptions::default()
            },
        )
        .expect_err("search must stop");
        assert_eq!(error.tag(), "SearchLimitReached");
    }

    #[test]
    fn exact_to_numeric_overflow_fails_instead_of_saturating() {
        assert_eq!(
            exact_i64_to_f64(9_007_199_254_740_993)
                .expect_err("2^53+1 must not round")
                .tag(),
            "BondTranslationOverflow"
        );
        let huge = QuadraticReal::from_rational(
            Rational::parse(&format!("1{}", "0".repeat(400)), "1").expect("huge integer"),
        );
        let error = find_periodic_bond_shells(
            &[vec![rational(0, 1)]],
            &[vec![huge]],
            1,
            &BondSearchOptions::default(),
        )
        .expect_err("nonfinite conversion must fail");
        assert_eq!(error.tag(), "NonFiniteBondGeometry");
    }
}
