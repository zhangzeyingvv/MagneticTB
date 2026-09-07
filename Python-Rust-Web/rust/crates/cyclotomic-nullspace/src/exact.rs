//! Direct Rust port of the Mathematica modules
//! `RationalArithmetic.wl`, `PolynomialArithmetic.wl`, and `CyclotomicField.wl`.

use crate::error::{ExactResult, fail};
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, ToPrimitive, Zero};
use std::collections::BTreeMap;
use std::str::FromStr;
use std::sync::Arc;

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct Rational(BigRational);

impl Rational {
    #[must_use]
    pub fn zero() -> Self {
        Self(BigRational::zero())
    }

    #[must_use]
    pub fn one() -> Self {
        Self(BigRational::one())
    }

    #[must_use]
    pub fn from_i64(value: i64) -> Self {
        Self(BigRational::from_integer(BigInt::from(value)))
    }

    pub fn new(numerator: BigInt, denominator: BigInt) -> ExactResult<Self> {
        if denominator.is_zero() {
            return fail("DivisionByZero", "rational denominator is zero");
        }
        Ok(Self(BigRational::new(numerator, denominator)))
    }

    pub fn parse(numerator: &str, denominator: &str) -> ExactResult<Self> {
        fn valid_integer(text: &str) -> bool {
            let digits = text.strip_prefix('-').unwrap_or(text);
            !digits.is_empty() && digits.bytes().all(|value| value.is_ascii_digit())
        }
        if !valid_integer(numerator) || !valid_integer(denominator) {
            return fail("MalformedSerialization", "invalid integer string");
        }
        let parsed_numerator = BigInt::from_str(numerator)
            .map_err(|_| crate::ExactError::new("MalformedSerialization", "invalid numerator"))?;
        let parsed_denominator = BigInt::from_str(denominator)
            .map_err(|_| crate::ExactError::new("MalformedSerialization", "invalid denominator"))?;
        Self::new(parsed_numerator, parsed_denominator)
    }

    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.0.is_zero()
    }

    #[must_use]
    pub fn is_integer(&self) -> bool {
        self.0.is_integer()
    }

    #[must_use]
    pub fn numerator_string(&self) -> String {
        self.0.numer().to_string()
    }

    #[must_use]
    pub fn denominator_string(&self) -> String {
        self.0.denom().to_string()
    }

    /// Convert an exact rational to the nearest finite binary64 value.
    ///
    /// This is intentionally an explicit boundary operation. Exact algebra
    /// never calls it; numerical geometry may use it for topology searches.
    #[must_use]
    pub fn to_f64(&self) -> Option<f64> {
        self.0.to_f64().filter(|value| value.is_finite())
    }

    #[must_use]
    pub(crate) fn numerator(&self) -> &BigInt {
        self.0.numer()
    }

    #[must_use]
    pub(crate) fn denominator(&self) -> &BigInt {
        self.0.denom()
    }

    #[must_use]
    pub fn add(&self, right: &Self) -> Self {
        Self(&self.0 + &right.0)
    }

    #[must_use]
    pub fn subtract(&self, right: &Self) -> Self {
        Self(&self.0 - &right.0)
    }

    #[must_use]
    pub fn multiply(&self, right: &Self) -> Self {
        Self(&self.0 * &right.0)
    }

    pub fn divide(&self, right: &Self) -> ExactResult<Self> {
        if right.is_zero() {
            return fail("DivisionByZero", "division by zero rational");
        }
        Ok(Self(&self.0 / &right.0))
    }

    #[must_use]
    pub fn negate(&self) -> Self {
        Self(-&self.0)
    }
}

impl From<i64> for Rational {
    fn from(value: i64) -> Self {
        Self::from_i64(value)
    }
}

pub type Polynomial = Vec<Rational>;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PolynomialDivision {
    pub quotient: Polynomial,
    pub remainder: Polynomial,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PolynomialExtendedGcd {
    pub gcd: Polynomial,
    pub left_coefficient: Polynomial,
    pub right_coefficient: Polynomial,
}

/// Port of `cqPolynomialNormalize`.
#[must_use]
pub fn polynomial_normalize(mut coefficients: Polynomial) -> Polynomial {
    if coefficients.is_empty() {
        return vec![Rational::zero()];
    }
    while coefficients.len() > 1 && coefficients.last().is_some_and(Rational::is_zero) {
        coefficients.pop();
    }
    coefficients
}

#[must_use]
pub fn polynomial_is_zero(polynomial: &[Rational]) -> bool {
    polynomial.iter().all(Rational::is_zero)
}

#[must_use]
pub fn polynomial_degree(polynomial: &[Rational]) -> Option<usize> {
    polynomial
        .iter()
        .rposition(|coefficient| !coefficient.is_zero())
}

pub fn polynomial_pad(polynomial: &[Rational], length: usize) -> ExactResult<Polynomial> {
    let mut normalized = polynomial_normalize(polynomial.to_vec());
    if normalized.len() > length {
        return fail(
            "MalformedSerialization",
            "polynomial is longer than requested padding",
        );
    }
    normalized.resize(length, Rational::zero());
    Ok(normalized)
}

#[must_use]
pub fn polynomial_add(left: &[Rational], right: &[Rational]) -> Polynomial {
    let size = left.len().max(right.len());
    let mut result = vec![Rational::zero(); size];
    for (index, value) in left.iter().enumerate() {
        result[index] = result[index].add(value);
    }
    for (index, value) in right.iter().enumerate() {
        result[index] = result[index].add(value);
    }
    polynomial_normalize(result)
}

#[must_use]
pub fn polynomial_negate(polynomial: &[Rational]) -> Polynomial {
    polynomial_normalize(polynomial.iter().map(Rational::negate).collect())
}

#[must_use]
pub fn polynomial_subtract(left: &[Rational], right: &[Rational]) -> Polynomial {
    polynomial_add(left, &polynomial_negate(right))
}

#[must_use]
pub fn polynomial_scale(scalar: &Rational, polynomial: &[Rational]) -> Polynomial {
    polynomial_normalize(
        polynomial
            .iter()
            .map(|coefficient| scalar.multiply(coefficient))
            .collect(),
    )
}

/// Direct convolution used by `cqPolynomialMultiply`.
#[must_use]
pub fn polynomial_multiply(left: &[Rational], right: &[Rational]) -> Polynomial {
    let left = polynomial_normalize(left.to_vec());
    let right = polynomial_normalize(right.to_vec());
    if polynomial_is_zero(&left) || polynomial_is_zero(&right) {
        return vec![Rational::zero()];
    }
    let mut result = vec![Rational::zero(); left.len() + right.len() - 1];
    for (left_index, left_value) in left.iter().enumerate() {
        for (right_index, right_value) in right.iter().enumerate() {
            let product = left_value.multiply(right_value);
            result[left_index + right_index] = result[left_index + right_index].add(&product);
        }
    }
    polynomial_normalize(result)
}

/// Port of `cqPolynomialDivRem`.
pub fn polynomial_div_rem(
    dividend: &[Rational],
    divisor: &[Rational],
) -> ExactResult<PolynomialDivision> {
    let mut remainder = polynomial_normalize(dividend.to_vec());
    let normalized_divisor = polynomial_normalize(divisor.to_vec());
    let Some(divisor_degree) = polynomial_degree(&normalized_divisor) else {
        return fail("DivisionByZero", "polynomial divisor is zero");
    };
    let mut quotient =
        vec![Rational::zero(); remainder.len().saturating_sub(divisor_degree).max(1)];
    let lead = &normalized_divisor[divisor_degree];
    while let Some(remainder_degree) = polynomial_degree(&remainder) {
        if remainder_degree < divisor_degree {
            break;
        }
        let shift = remainder_degree - divisor_degree;
        let factor = remainder[remainder_degree].divide(lead)?;
        quotient[shift] = quotient[shift].add(&factor);
        let mut term = vec![Rational::zero(); shift];
        term.extend(polynomial_scale(&factor, &normalized_divisor));
        remainder = polynomial_subtract(&remainder, &term);
    }
    Ok(PolynomialDivision {
        quotient: polynomial_normalize(quotient),
        remainder: polynomial_normalize(remainder),
    })
}

pub fn polynomial_exact_quotient(
    dividend: &[Rational],
    divisor: &[Rational],
) -> ExactResult<Polynomial> {
    let division = polynomial_div_rem(dividend, divisor)?;
    if !polynomial_is_zero(&division.remainder) {
        return fail(
            "SingularPolynomialElement",
            "polynomial division has nonzero remainder",
        );
    }
    Ok(division.quotient)
}

pub fn polynomial_mod(polynomial: &[Rational], modulus: &[Rational]) -> ExactResult<Polynomial> {
    Ok(polynomial_div_rem(polynomial, modulus)?.remainder)
}

pub fn polynomial_multiply_mod(
    left: &[Rational],
    right: &[Rational],
    modulus: &[Rational],
) -> ExactResult<Polynomial> {
    polynomial_mod(&polynomial_multiply(left, right), modulus)
}

pub fn polynomial_power_mod(
    base: &[Rational],
    exponent: usize,
    modulus: &[Rational],
) -> ExactResult<Polynomial> {
    let mut factor = polynomial_mod(base, modulus)?;
    let mut result = vec![Rational::one()];
    let mut power = exponent;
    while power > 0 {
        if power % 2 == 1 {
            result = polynomial_multiply_mod(&result, &factor, modulus)?;
        }
        power /= 2;
        if power > 0 {
            factor = polynomial_multiply_mod(&factor, &factor, modulus)?;
        }
    }
    Ok(polynomial_normalize(result))
}

/// Port of `cqPolynomialExtendedGCD`, including monic gcd normalization.
pub fn polynomial_extended_gcd(
    left: &[Rational],
    right: &[Rational],
) -> ExactResult<PolynomialExtendedGcd> {
    let mut old_remainder = polynomial_normalize(left.to_vec());
    let mut remainder = polynomial_normalize(right.to_vec());
    let mut old_left = vec![Rational::one()];
    let mut left_coefficient = vec![Rational::zero()];
    let mut old_right = vec![Rational::zero()];
    let mut right_coefficient = vec![Rational::one()];
    while !polynomial_is_zero(&remainder) {
        let division = polynomial_div_rem(&old_remainder, &remainder)?;
        let next_left = polynomial_subtract(
            &old_left,
            &polynomial_multiply(&division.quotient, &left_coefficient),
        );
        let next_right = polynomial_subtract(
            &old_right,
            &polynomial_multiply(&division.quotient, &right_coefficient),
        );
        old_remainder = remainder;
        remainder = division.remainder;
        old_left = left_coefficient;
        left_coefficient = next_left;
        old_right = right_coefficient;
        right_coefficient = next_right;
    }
    if polynomial_is_zero(&old_remainder) {
        return Ok(PolynomialExtendedGcd {
            gcd: vec![Rational::zero()],
            left_coefficient: vec![Rational::zero()],
            right_coefficient: vec![Rational::zero()],
        });
    }
    let scale = Rational::one().divide(old_remainder.last().expect("nonzero polynomial"))?;
    Ok(PolynomialExtendedGcd {
        gcd: polynomial_scale(&scale, &old_remainder),
        left_coefficient: polynomial_scale(&scale, &old_left),
        right_coefficient: polynomial_scale(&scale, &old_right),
    })
}

pub fn integer_factorization(value: usize) -> ExactResult<Vec<(usize, usize)>> {
    if value < 1 {
        return fail(
            "InvalidRootOfUnity",
            "integer factorization requires positive input",
        );
    }
    let mut remaining = value;
    let mut divisor = 2;
    let mut factors = Vec::new();
    while divisor <= remaining / divisor {
        let mut exponent = 0;
        while remaining % divisor == 0 {
            remaining /= divisor;
            exponent += 1;
        }
        if exponent > 0 {
            factors.push((divisor, exponent));
        }
        divisor = if divisor == 2 { 3 } else { divisor + 2 };
    }
    if remaining > 1 {
        factors.push((remaining, 1));
    }
    Ok(factors)
}

pub fn euler_phi(value: usize) -> ExactResult<usize> {
    let mut result = value;
    for (prime, _) in integer_factorization(value)? {
        result = (result / prime) * (prime - 1);
    }
    Ok(result)
}

pub fn integer_divisors(value: usize) -> ExactResult<Vec<usize>> {
    let mut divisors = vec![1];
    for (prime, exponent) in integer_factorization(value)? {
        let previous = divisors;
        divisors = Vec::new();
        let mut power = 1;
        for current_exponent in 0..=exponent {
            divisors.extend(previous.iter().map(|divisor| divisor * power));
            if current_exponent < exponent {
                power *= prime;
            }
        }
    }
    divisors.sort_unstable();
    Ok(divisors)
}

/// Port of `cqCyclotomicPolynomial` using x^N-1 factorization.
pub fn cyclotomic_polynomial(conductor: usize) -> ExactResult<Polynomial> {
    if conductor < 1 {
        return fail("InvalidRootOfUnity", "conductor must be positive");
    }
    let divisors = integer_divisors(conductor)?;
    let mut computed: BTreeMap<usize, Polynomial> = BTreeMap::new();
    for &current in &divisors {
        let mut polynomial = vec![Rational::zero(); current + 1];
        polynomial[0] = Rational::from(-1);
        polynomial[current] = Rational::one();
        for &proper in &divisors {
            if proper >= current {
                break;
            }
            if current % proper == 0 {
                polynomial = polynomial_exact_quotient(
                    &polynomial,
                    computed
                        .get(&proper)
                        .expect("proper divisor already computed"),
                )?;
            }
        }
        computed.insert(current, polynomial_normalize(polynomial));
    }
    Ok(computed
        .remove(&conductor)
        .expect("conductor is one of its divisors"))
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CyclotomicContext {
    pub conductor: usize,
    pub degree: usize,
    pub cyclotomic_polynomial: Polynomial,
}

/// Port of `cqCreateContext` (the Mathematica cache is an optimization, not semantics).
pub fn make_context(conductor: usize, max_degree: usize) -> ExactResult<Arc<CyclotomicContext>> {
    if conductor < 1 {
        return fail("InvalidRootOfUnity", "conductor must be positive");
    }
    if max_degree < 1 {
        return fail(
            "ConductorDegreeLimitExceeded",
            "maximum cyclotomic degree must be positive",
        );
    }
    let degree = euler_phi(conductor)?;
    if degree > max_degree {
        return fail(
            "ConductorDegreeLimitExceeded",
            "cyclotomic degree exceeds configured limit",
        );
    }
    let modulus = cyclotomic_polynomial(conductor)?;
    if modulus.len() != degree + 1 || modulus.last() != Some(&Rational::one()) {
        return fail("SingularPolynomialElement", "invalid cyclotomic polynomial");
    }
    Ok(Arc::new(CyclotomicContext {
        conductor,
        degree,
        cyclotomic_polynomial: modulus,
    }))
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Element {
    context: Arc<CyclotomicContext>,
    coefficients: Polynomial,
}

impl Element {
    /// Port of `cqElementFromPolynomial`.
    pub fn from_polynomial(
        context: Arc<CyclotomicContext>,
        polynomial: &[Rational],
    ) -> ExactResult<Self> {
        let reduced = polynomial_mod(polynomial, &context.cyclotomic_polynomial)?;
        let coefficients = polynomial_pad(&reduced, context.degree)?;
        Ok(Self {
            context,
            coefficients,
        })
    }

    pub fn zero(context: &Arc<CyclotomicContext>) -> ExactResult<Self> {
        Self::from_polynomial(Arc::clone(context), &[Rational::zero()])
    }

    pub fn one(context: &Arc<CyclotomicContext>) -> ExactResult<Self> {
        Self::from_polynomial(Arc::clone(context), &[Rational::one()])
    }

    /// Port of `cqRootOfUnityElement`.
    pub fn root_of_unity(
        context: &Arc<CyclotomicContext>,
        order: usize,
        power: i64,
    ) -> ExactResult<Self> {
        if order < 1 {
            return fail("InvalidRootOfUnity", "root order must be positive");
        }
        if context.conductor % order != 0 {
            return fail(
                "ConductorMismatch",
                "root order does not divide target conductor",
            );
        }
        let signed_order = i64::try_from(order)
            .map_err(|_| crate::ExactError::new("InvalidRootOfUnity", "root order exceeds i64"))?;
        let reduced_power = power.rem_euclid(signed_order);
        let exponent = usize::try_from(reduced_power).expect("remainder is nonnegative")
            * (context.conductor / order);
        let polynomial = polynomial_power_mod(
            &[Rational::zero(), Rational::one()],
            exponent,
            &context.cyclotomic_polynomial,
        )?;
        Self::from_polynomial(Arc::clone(context), &polynomial)
    }

    #[must_use]
    pub fn context(&self) -> &Arc<CyclotomicContext> {
        &self.context
    }

    #[must_use]
    pub fn coefficients(&self) -> &[Rational] {
        &self.coefficients
    }

    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.coefficients.iter().all(Rational::is_zero)
    }

    fn check_pair(&self, right: &Self) -> ExactResult<()> {
        if self.context != right.context {
            return fail(
                "ConductorMismatch",
                "cyclotomic elements have different contexts",
            );
        }
        Ok(())
    }

    pub fn add(&self, right: &Self) -> ExactResult<Self> {
        self.check_pair(right)?;
        Self::from_polynomial(
            Arc::clone(&self.context),
            &polynomial_add(&self.coefficients, &right.coefficients),
        )
    }

    pub fn negate(&self) -> ExactResult<Self> {
        Self::from_polynomial(
            Arc::clone(&self.context),
            &polynomial_negate(&self.coefficients),
        )
    }

    pub fn subtract(&self, right: &Self) -> ExactResult<Self> {
        self.add(&right.negate()?)
    }

    pub fn scale(&self, scalar: &Rational) -> ExactResult<Self> {
        Self::from_polynomial(
            Arc::clone(&self.context),
            &polynomial_scale(scalar, &self.coefficients),
        )
    }

    pub fn multiply(&self, right: &Self) -> ExactResult<Self> {
        self.check_pair(right)?;
        Self::from_polynomial(
            Arc::clone(&self.context),
            &polynomial_multiply(&self.coefficients, &right.coefficients),
        )
    }

    pub fn inverse(&self) -> ExactResult<Self> {
        if self.is_zero() {
            return fail("DivisionByZero", "cannot invert zero cyclotomic element");
        }
        let extended =
            polynomial_extended_gcd(&self.coefficients, &self.context.cyclotomic_polynomial)?;
        if extended.gcd != vec![Rational::one()] {
            return fail(
                "SingularPolynomialElement",
                "cyclotomic element is not invertible",
            );
        }
        Self::from_polynomial(Arc::clone(&self.context), &extended.left_coefficient)
    }

    pub fn divide(&self, right: &Self) -> ExactResult<Self> {
        self.multiply(&right.inverse()?)
    }

    pub fn power(&self, exponent: i64) -> ExactResult<Self> {
        let mut factor = if exponent < 0 {
            self.inverse()?
        } else {
            self.clone()
        };
        let mut power = exponent.unsigned_abs();
        let mut result = Self::one(&self.context)?;
        while power > 0 {
            if power % 2 == 1 {
                result = result.multiply(&factor)?;
            }
            power /= 2;
            if power > 0 {
                factor = factor.multiply(&factor)?;
            }
        }
        Ok(result)
    }

    /// Port of `cqElementConjugate`.
    pub fn conjugate(&self) -> ExactResult<Self> {
        let mut result = Self::zero(&self.context)?;
        for (index, coefficient) in self.coefficients.iter().enumerate() {
            if coefficient.is_zero() {
                continue;
            }
            let exponent =
                (self.context.conductor - index % self.context.conductor) % self.context.conductor;
            let basis_image = polynomial_power_mod(
                &[Rational::zero(), Rational::one()],
                exponent,
                &self.context.cyclotomic_polynomial,
            )?;
            let image = Self::from_polynomial(Arc::clone(&self.context), &basis_image)?;
            result = result.add(&image.scale(coefficient)?)?;
        }
        Ok(result)
    }
}

/// Port of `cqEmbedElement`.
pub fn embed_element(
    source_context: &Arc<CyclotomicContext>,
    target_context: &Arc<CyclotomicContext>,
    element: &Element,
) -> ExactResult<Element> {
    if element.context != *source_context
        || target_context.conductor % source_context.conductor != 0
    {
        return fail(
            "ConductorMismatch",
            "source conductor does not divide target conductor",
        );
    }
    let mut result = Element::zero(target_context)?;
    for (index, coefficient) in element.coefficients.iter().enumerate() {
        if coefficient.is_zero() {
            continue;
        }
        let power = i64::try_from(index).map_err(|_| {
            crate::ExactError::new("InvalidRootOfUnity", "embedding exponent exceeds i64")
        })?;
        let image = Element::root_of_unity(target_context, source_context.conductor, power)?;
        result = result.add(&image.scale(coefficient)?)?;
    }
    Ok(result)
}
