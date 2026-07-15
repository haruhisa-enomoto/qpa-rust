//! QPA's predefined Nakayama, canonical, and Kronecker algebras.
//!
//! This is a line-by-line translation of the first constructor block in
//! QPA 1.37 `lib/pathalgpredef.gi`.  In particular, nonempty presentations
//! create an ideal from the original relations, attach the configured GBNP
//! basis, and only then form the quotient.  The constructor filters set by QPA
//! are retained by [`PredefinedAlgebra`] instead of being inferred from the
//! resulting presentation.

use numerica::domains::{Field, integer::Integer};
use qpa_core::{
    ArrowSpec, NamedArrowSpec, PathAlgebra, PathAlgebraElement, PathAlgebraError, PathAlgebraIdeal,
    Quiver, QuiverAlgebra, QuiverError, QuotientError, QuotientPathAlgebra, VertexIndex,
};
use thiserror::Error;

use crate::{QpaGbnpError, gbnp_groebner_basis};

/// The diagnostic printed by QPA before it returns an invalid admissible
/// sequence unchanged.
pub const INVALID_ADMISSIBLE_SEQUENCE_DIAGNOSTIC: &str =
    "The admissible sequence is not a valid one.\n";

/// The constructor filter installed on a predefined QPA algebra.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum PredefinedAlgebraKind {
    /// `IsNakayamaAlgebra`.
    Nakayama,
    /// `IsCanonicalAlgebra`.
    Canonical,
    /// `IsKroneckerAlgebra`.
    Kronecker,
}

/// A path algebra or admissible quotient carrying QPA's constructor filter.
///
/// QPA does not recognize these three classes from arbitrary presentations in
/// this source block.  It sets a filter on the object returned by the
/// constructor.  This wrapper makes that provenance observable without adding
/// an independent recognition algorithm.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PredefinedAlgebra<F: Field> {
    kind: PredefinedAlgebraKind,
    algebra: QuiverAlgebra<F>,
}

impl<F: Field> PredefinedAlgebra<F> {
    fn new(kind: PredefinedAlgebraKind, algebra: QuiverAlgebra<F>) -> Self {
        Self { kind, algebra }
    }

    /// Return the constructor filter set by QPA.
    #[must_use]
    pub const fn kind(&self) -> PredefinedAlgebraKind {
        self.kind
    }

    /// Borrow the constructed path algebra or admissible quotient.
    #[must_use]
    pub const fn algebra(&self) -> &QuiverAlgebra<F> {
        &self.algebra
    }

    /// Consume the constructor tag and return the underlying quiver algebra.
    #[must_use]
    pub fn into_algebra(self) -> QuiverAlgebra<F> {
        self.algebra
    }

    /// Return whether the `IsNakayamaAlgebra` filter is installed.
    #[must_use]
    pub const fn has_nakayama_algebra_filter(&self) -> bool {
        matches!(self.kind, PredefinedAlgebraKind::Nakayama)
    }

    /// Return QPA's constructor-known `IsNakayamaAlgebra` value.
    ///
    /// `None` means the constructor did not install the filter; QPA has no
    /// recognition method here from which to infer `false`.
    #[must_use]
    pub const fn is_nakayama_algebra(&self) -> Option<bool> {
        if self.has_nakayama_algebra_filter() {
            Some(true)
        } else {
            None
        }
    }

    /// Return whether the `IsCanonicalAlgebra` filter is installed.
    #[must_use]
    pub const fn has_canonical_algebra_filter(&self) -> bool {
        matches!(self.kind, PredefinedAlgebraKind::Canonical)
    }

    /// Return QPA's constructor-known `IsCanonicalAlgebra` value.
    ///
    /// `None` means the constructor did not install the filter; QPA has no
    /// recognition method here from which to infer `false`.
    #[must_use]
    pub const fn is_canonical_algebra(&self) -> Option<bool> {
        if self.has_canonical_algebra_filter() {
            Some(true)
        } else {
            None
        }
    }

    /// Return whether the `IsKroneckerAlgebra` filter is installed.
    #[must_use]
    pub const fn has_kronecker_algebra_filter(&self) -> bool {
        matches!(self.kind, PredefinedAlgebraKind::Kronecker)
    }

    /// Return QPA's constructor-known `IsKroneckerAlgebra` value.
    ///
    /// `None` means the constructor did not install the filter; QPA has no
    /// recognition method here from which to infer `false`.
    #[must_use]
    pub const fn is_kronecker_algebra(&self) -> Option<bool> {
        if self.has_kronecker_algebra_filter() {
            Some(true)
        } else {
            None
        }
    }

    /// Return QPA's constructor-implied `IsFiniteTypeAlgebra` value.
    ///
    /// `None` means this source block installs no implication for the
    /// constructor; it does not assert mathematical falsity.
    #[must_use]
    pub const fn is_finite_type_algebra(&self) -> Option<bool> {
        if self.has_nakayama_algebra_filter() {
            Some(true)
        } else {
            None
        }
    }

    /// Return QPA's constructor-implied finite-global-dimension value.
    ///
    /// `None` means this source block installs no implication for the
    /// constructor; it does not assert mathematical falsity.
    #[must_use]
    pub const fn is_finite_global_dimension_algebra(&self) -> Option<bool> {
        if self.has_canonical_algebra_filter() || self.has_kronecker_algebra_filter() {
            Some(true)
        } else {
            None
        }
    }

    /// Return QPA's constructor-implied `IsHereditaryAlgebra` value.
    ///
    /// `None` means this source block installs no implication for the
    /// constructor; it does not assert mathematical falsity.
    #[must_use]
    pub const fn is_hereditary_algebra(&self) -> Option<bool> {
        if self.has_kronecker_algebra_filter() {
            Some(true)
        } else {
            None
        }
    }
}

/// QPA's two possible return types from `NakayamaAlgebra`.
///
/// QPA returns an invalid sequence unchanged after printing
/// [`INVALID_ADMISSIBLE_SEQUENCE_DIAGNOSTIC`].  That branch is not a GAP
/// error, so Rust retains it as a successful typed result.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum NakayamaAlgebraResult<F: Field> {
    /// A path algebra or admissible quotient tagged as Nakayama.
    Algebra(PredefinedAlgebra<F>),
    /// The original invalid sequence returned unchanged by QPA.
    InvalidSequence(Vec<Integer>),
}

impl<F: Field> NakayamaAlgebraResult<F> {
    /// Borrow the constructed algebra, if the sequence was admissible.
    #[must_use]
    pub const fn as_algebra(&self) -> Option<&PredefinedAlgebra<F>> {
        match self {
            Self::Algebra(algebra) => Some(algebra),
            Self::InvalidSequence(_) => None,
        }
    }

    /// Return the exact QPA diagnostic for the invalid-sequence branch.
    #[must_use]
    pub const fn diagnostic(&self) -> Option<&'static str> {
        match self {
            Self::Algebra(_) => None,
            Self::InvalidSequence(_) => Some(INVALID_ADMISSIBLE_SEQUENCE_DIAGNOSTIC),
        }
    }
}

/// Translate `NakayamaAlgebra(K, admiss_seq)` using QPA's default GBNP
/// backend.
///
/// The input uses Numerica's arbitrary-precision [`Integer`] so the intended
/// GAP integer domain is not artificially narrowed.  An invalid nonempty
/// sequence is returned in [`NakayamaAlgebraResult::InvalidSequence`], exactly
/// as in QPA.
///
/// # Errors
///
/// The empty list reproduces QPA's out-of-range list access as a typed error.
/// Valid entries that cannot be represented by the platform's allocation
/// index type are also rejected.  Quiver, path-algebra, GBNP, and quotient
/// failures are propagated.
pub fn nakayama_algebra<F: Field>(
    field: F,
    admissible_sequence: &[Integer],
) -> Result<NakayamaAlgebraResult<F>, PredefinedAlgebraError> {
    if admissible_sequence.is_empty() {
        return Err(PredefinedAlgebraError::EmptyAdmissibleSequence);
    }
    if !is_admissible_sequence(admissible_sequence) {
        return Ok(NakayamaAlgebraResult::InvalidSequence(
            admissible_sequence.to_vec(),
        ));
    }

    // For a singleton list QPA skips the loop containing its positivity test.
    // Consequently `[0]` and negative singletons reach the cyclic branch,
    // generate `One(KQ)` as their sole relation, and fail in the configured
    // Gröbner backend. Preserve that call path instead of underflowing a Rust
    // unsigned length conversion.
    if admissible_sequence.len() == 1 && admissible_sequence[0] < 1 {
        let only_vertex = vertex(1)?;
        let quiver = Quiver::from_indexed(1, vec![ArrowSpec::unnamed(only_vertex, only_vertex)])?;
        let path_algebra = PathAlgebra::new(field, quiver)?;
        let relation = path_algebra.one();
        let algebra = factor_predefined_presentation(path_algebra, &[relation])?;
        return Ok(NakayamaAlgebraResult::Algebra(PredefinedAlgebra::new(
            PredefinedAlgebraKind::Nakayama,
            algebra,
        )));
    }

    let sequence = admissible_sequence
        .iter()
        .enumerate()
        .map(|(position, value)| {
            usize::try_from(value.clone()).map_err(|_| {
                PredefinedAlgebraError::SizeNotRepresentable {
                    parameter: "admissible-sequence entry",
                    position: position + 1,
                    value: value.clone(),
                }
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let number_of_vertices = sequence.len();

    let (path_algebra, relations) = if sequence[number_of_vertices - 1] == 1 {
        linear_nakayama_presentation(field, &sequence)?
    } else {
        cyclic_nakayama_presentation(field, &sequence)?
    };
    let algebra = factor_predefined_presentation(path_algebra, &relations)?;
    Ok(NakayamaAlgebraResult::Algebra(PredefinedAlgebra::new(
        PredefinedAlgebraKind::Nakayama,
        algebra,
    )))
}

/// Translate QPA's sequence-first `NakayamaAlgebra(admiss_seq, K)` overload.
///
/// # Errors
///
/// Delegates to [`nakayama_algebra`].
pub fn nakayama_algebra_sequence_first<F: Field>(
    admissible_sequence: &[Integer],
    field: F,
) -> Result<NakayamaAlgebraResult<F>, PredefinedAlgebraError> {
    nakayama_algebra(field, admissible_sequence)
}

/// Translate `CanonicalAlgebra(K, weights, relcoeff)` using QPA's default
/// GBNP backend.
///
/// Coefficients have type `F::Element`, so QPA's membership check is enforced
/// by Rust's type system.  The number, order, and names of all vertices,
/// arrows, arm products, and relations follow the pinned source.
///
/// # Errors
///
/// Returns QPA-compatible typed errors for mismatched input lengths, weights
/// below two, and zero relation coefficients.  Structural and GBNP failures
/// are propagated.
pub fn canonical_algebra<F: Field>(
    field: F,
    weights: &[Integer],
    relation_coefficients: &[F::Element],
) -> Result<PredefinedAlgebra<F>, PredefinedAlgebraError> {
    if weights.len() != relation_coefficients.len() + 2 {
        return Err(PredefinedAlgebraError::CanonicalArmCoefficientMismatch {
            arms: weights.len(),
            coefficients: relation_coefficients.len(),
        });
    }
    if let Some((position, weight)) = weights.iter().enumerate().find(|(_, weight)| *weight < &2) {
        return Err(PredefinedAlgebraError::CanonicalWeightBelowTwo {
            position: position + 1,
            weight: weight.clone(),
        });
    }
    if let Some(position) = relation_coefficients
        .iter()
        .position(|coefficient| field.is_zero(coefficient))
    {
        return Err(PredefinedAlgebraError::ZeroCanonicalRelationCoefficient {
            position: position + 1,
        });
    }

    let weights = weights
        .iter()
        .enumerate()
        .map(|(position, value)| {
            usize::try_from(value.clone()).map_err(|_| {
                PredefinedAlgebraError::SizeNotRepresentable {
                    parameter: "canonical weight",
                    position: position + 1,
                    value: value.clone(),
                }
            })
        })
        .collect::<Result<Vec<_>, _>>()?;

    let mut vertices = vec![Box::<str>::from("v")];
    for (arm, &weight) in weights.iter().enumerate() {
        for position in 2..=weight {
            vertices.push(format!("v{}{position}", arm + 1).into_boxed_str());
        }
    }
    vertices.push(Box::<str>::from("w"));

    let mut arrows = Vec::new();
    for (arm, &weight) in weights.iter().enumerate() {
        let arm = arm + 1;
        arrows.push(NamedArrowSpec::named(
            "v",
            format!("v{arm}2"),
            format!("a{arm}"),
        ));
        for position in 2..weight {
            arrows.push(NamedArrowSpec::named(
                format!("v{arm}{position}"),
                format!("v{arm}{}", position + 1),
                format!("a{arm}{position}"),
            ));
        }
        arrows.push(NamedArrowSpec::named(
            format!("v{arm}{weight}"),
            "w",
            format!("b{arm}"),
        ));
    }

    let quiver = Quiver::from_named(vertices, arrows)?;
    let path_algebra = PathAlgebra::new(field, quiver)?;
    if weights.len() == 2 {
        return Ok(PredefinedAlgebra::new(
            PredefinedAlgebraKind::Canonical,
            QuiverAlgebra::Path(path_algebra),
        ));
    }

    let generators = path_algebra.generators();
    let arrow_generators = &generators[path_algebra.quiver().number_of_vertices()..];
    let mut arms = Vec::with_capacity(weights.len());
    let mut start = 0_usize;
    for &weight in &weights {
        let end = start
            .checked_add(weight)
            .ok_or(PredefinedAlgebraError::GeneratedIndexOverflow)?;
        arms.push(product_of_arrows(
            &path_algebra,
            arrow_generators,
            start..end,
        )?);
        start = end;
    }

    let mut relations = Vec::with_capacity(weights.len() - 2);
    for arm in 2..weights.len() {
        let relation = arms[arm]
            .sub(&arms[1])?
            .add(&arms[0].scale_left(&relation_coefficients[arm - 2])?)?;
        relations.push(relation);
    }
    let algebra = factor_predefined_presentation(path_algebra, &relations)?;
    Ok(PredefinedAlgebra::new(
        PredefinedAlgebraKind::Canonical,
        algebra,
    ))
}

/// Translate QPA's two-argument `CanonicalAlgebra(K, weights)` overload.
///
/// # Errors
///
/// QPA accepts this overload only for exactly two weights; construction then
/// delegates to [`canonical_algebra`] with an empty coefficient list.
pub fn canonical_algebra_two_weights<F: Field>(
    field: F,
    weights: &[Integer],
) -> Result<PredefinedAlgebra<F>, PredefinedAlgebraError> {
    if weights.len() != 2 {
        return Err(PredefinedAlgebraError::CanonicalTwoArgumentWeightCount {
            weights: weights.len(),
        });
    }
    canonical_algebra(field, weights, &[])
}

/// Translate `KroneckerAlgebra(K, n)`.
///
/// # Errors
///
/// QPA rejects fewer than two arrows.  Counts outside the platform allocation
/// index range and quiver/path-algebra construction failures are also
/// reported.
pub fn kronecker_algebra<F: Field>(
    field: F,
    number_of_arrows: &Integer,
) -> Result<PredefinedAlgebra<F>, PredefinedAlgebraError> {
    if number_of_arrows < &2 {
        return Err(PredefinedAlgebraError::KroneckerNeedsTwoArrows {
            number_of_arrows: number_of_arrows.clone(),
        });
    }
    let count = usize::try_from(number_of_arrows.clone()).map_err(|_| {
        PredefinedAlgebraError::SizeNotRepresentable {
            parameter: "Kronecker arrow count",
            position: 1,
            value: number_of_arrows.clone(),
        }
    })?;
    let source = vertex(1)?;
    let target = vertex(2)?;
    let arrows = (1..=count)
        .map(|position| ArrowSpec::named(source, target, format!("a{position}")))
        .collect();
    let quiver = Quiver::from_indexed(2, arrows)?;
    let path_algebra = PathAlgebra::new(field, quiver)?;
    Ok(PredefinedAlgebra::new(
        PredefinedAlgebraKind::Kronecker,
        QuiverAlgebra::Path(path_algebra),
    ))
}

fn is_admissible_sequence(sequence: &[Integer]) -> bool {
    let one = Integer::from(1);
    for pair in sequence.windows(2) {
        let previous_minus_one = &pair[0] - &one;
        if pair[1] < previous_minus_one || previous_minus_one < one {
            return false;
        }
    }
    sequence[0] >= &sequence[sequence.len() - 1] - &one
}

fn linear_nakayama_presentation<F: Field>(
    field: F,
    sequence: &[usize],
) -> Result<(PathAlgebra<F>, Vec<PathAlgebraElement<F>>), PredefinedAlgebraError> {
    let number_of_vertices = sequence.len();
    let mut arrows = Vec::with_capacity(number_of_vertices.saturating_sub(1));
    for position in 1..number_of_vertices {
        arrows.push(ArrowSpec::unnamed(vertex(position)?, vertex(position + 1)?));
    }
    let quiver = Quiver::from_indexed(number_of_vertices, arrows)?;
    let path_algebra = PathAlgebra::new(field, quiver)?;
    let generators = path_algebra.generators();
    let arrow_generators = &generators[number_of_vertices..];
    let mut relations = Vec::new();
    for position in 0..number_of_vertices.saturating_sub(1) {
        if sequence[position] < number_of_vertices - position
            && sequence[position] - 1 != sequence[position + 1]
        {
            let end = position
                .checked_add(sequence[position])
                .ok_or(PredefinedAlgebraError::GeneratedIndexOverflow)?;
            relations.push(product_of_arrows(
                &path_algebra,
                arrow_generators,
                position..end,
            )?);
        }
    }
    Ok((path_algebra, relations))
}

fn cyclic_nakayama_presentation<F: Field>(
    field: F,
    sequence: &[usize],
) -> Result<(PathAlgebra<F>, Vec<PathAlgebraElement<F>>), PredefinedAlgebraError> {
    let number_of_vertices = sequence.len();
    let mut arrows = Vec::with_capacity(number_of_vertices);
    for position in 1..number_of_vertices {
        arrows.push(ArrowSpec::unnamed(vertex(position)?, vertex(position + 1)?));
    }
    arrows.push(ArrowSpec::unnamed(vertex(number_of_vertices)?, vertex(1)?));
    let quiver = Quiver::from_indexed(number_of_vertices, arrows)?;
    let path_algebra = PathAlgebra::new(field, quiver)?;
    let generators = path_algebra.generators();
    let arrow_generators = &generators[number_of_vertices..];
    let mut relations = Vec::new();

    for position in 0..number_of_vertices {
        let next = (position + 1) % number_of_vertices;
        if sequence[position] - 1 == sequence[next] {
            continue;
        }
        let quotient = sequence[position] / number_of_vertices;
        let remainder = sequence[position] % number_of_vertices;
        let cycle = product_of_arrows(
            &path_algebra,
            arrow_generators,
            (position..number_of_vertices).chain(0..position),
        )?;
        let mut relation = path_algebra.one();
        for _ in 0..quotient {
            relation = relation.multiply(&cycle)?;
        }
        relation = append_cyclic_remainder(relation, arrow_generators, position, remainder)?;
        relations.push(relation);
    }
    Ok((path_algebra, relations))
}

fn append_cyclic_remainder<F: Field>(
    mut product: PathAlgebraElement<F>,
    arrows: &[PathAlgebraElement<F>],
    start: usize,
    length: usize,
) -> Result<PathAlgebraElement<F>, PredefinedAlgebraError> {
    for offset in 0..length {
        let index = (start + offset) % arrows.len();
        let factor = arrows
            .get(index)
            .ok_or(PredefinedAlgebraError::GeneratedArrowIndexOutOfRange { index })?;
        product = product.multiply(factor)?;
    }
    Ok(product)
}

fn product_of_arrows<F: Field>(
    path_algebra: &PathAlgebra<F>,
    arrows: &[PathAlgebraElement<F>],
    indices: impl IntoIterator<Item = usize>,
) -> Result<PathAlgebraElement<F>, PredefinedAlgebraError> {
    let mut product = path_algebra.one();
    for index in indices {
        let factor = arrows
            .get(index)
            .ok_or(PredefinedAlgebraError::GeneratedArrowIndexOutOfRange { index })?;
        product = product.multiply(factor)?;
    }
    Ok(product)
}

fn factor_predefined_presentation<F: Field>(
    path_algebra: PathAlgebra<F>,
    relations: &[PathAlgebraElement<F>],
) -> Result<QuiverAlgebra<F>, PredefinedAlgebraError> {
    if relations.is_empty() {
        return Ok(QuiverAlgebra::Path(path_algebra));
    }
    let complete_relations = gbnp_groebner_basis(relations, &path_algebra)?;
    let ideal = PathAlgebraIdeal::new(path_algebra, relations.to_vec())?;
    ideal.with_complete_groebner_basis(&complete_relations)?;
    let quotient = QuotientPathAlgebra::by_ideal(&ideal)?;
    let _ = quotient.is_admissible()?;
    Ok(QuiverAlgebra::Quotient(quotient))
}

fn vertex(index: usize) -> Result<VertexIndex, PredefinedAlgebraError> {
    VertexIndex::new(index).ok_or(PredefinedAlgebraError::GeneratedIndexOverflow)
}

/// Failure in QPA's predefined algebra constructors.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum PredefinedAlgebraError {
    /// QPA indexes both the first and last entry of an empty list.
    #[error("the admissible sequence must not be empty")]
    EmptyAdmissibleSequence,
    /// Canonical algebras have two more arms than relation coefficients.
    #[error(
        "number of canonical arms ({arms}) does not match relation coefficients ({coefficients})"
    )]
    CanonicalArmCoefficientMismatch {
        /// Number of supplied weights.
        arms: usize,
        /// Number of supplied relation coefficients.
        coefficients: usize,
    },
    /// Every canonical weight must be at least two.
    #[error("canonical weight {position} is below two: {weight}")]
    CanonicalWeightBelowTwo {
        /// One-based list position.
        position: usize,
        /// Rejected weight.
        weight: Integer,
    },
    /// Every canonical relation coefficient must be nonzero.
    #[error("canonical relation coefficient {position} is zero")]
    ZeroCanonicalRelationCoefficient {
        /// One-based list position.
        position: usize,
    },
    /// QPA's two-argument canonical overload accepts exactly two weights.
    #[error("the two-argument canonical constructor needs two weights, received {weights}")]
    CanonicalTwoArgumentWeightCount {
        /// Number of supplied weights.
        weights: usize,
    },
    /// A Kronecker quiver must have at least two parallel arrows.
    #[error("the Kronecker arrow count must be at least two: {number_of_arrows}")]
    KroneckerNeedsTwoArrows {
        /// Rejected arrow count.
        number_of_arrows: Integer,
    },
    /// An arbitrary-precision source count cannot index a Rust allocation.
    #[error("{parameter} at position {position} cannot be represented as a platform size: {value}")]
    SizeNotRepresentable {
        /// Kind of count being converted.
        parameter: &'static str,
        /// One-based position, or one for a scalar count.
        position: usize,
        /// Rejected arbitrary-precision value.
        value: Integer,
    },
    /// A generated endpoint or range overflowed the platform index type.
    #[error("a generated predefined-algebra index overflowed")]
    GeneratedIndexOverflow,
    /// A translated source range did not address a generated arrow.
    #[error("generated arrow index {index} is outside the predefined quiver")]
    GeneratedArrowIndexOutOfRange {
        /// Zero-based internal arrow index.
        index: usize,
    },
    /// Quiver construction failed.
    #[error(transparent)]
    Quiver(#[from] QuiverError),
    /// Path-algebra construction or arithmetic failed.
    #[error(transparent)]
    PathAlgebra(#[from] PathAlgebraError),
    /// QPA-to-GBNP conversion or `SGrobner` failed.
    #[error(transparent)]
    Gbnp(#[from] QpaGbnpError),
    /// Ideal or quotient construction failed.
    #[error(transparent)]
    Quotient(#[from] QuotientError),
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        Field, Ring,
        finite_field::Z2,
        integer::Integer,
        rational::{Q, RationalField},
    };
    use qpa_core::{AlgebraDimension, PathAlgebraElement, QuiverAlgebra};
    use serde_json::{Value, json};

    use super::{
        INVALID_ADMISSIBLE_SEQUENCE_DIAGNOSTIC, NakayamaAlgebraResult, PredefinedAlgebra,
        PredefinedAlgebraError, canonical_algebra, canonical_algebra_two_weights,
        kronecker_algebra, nakayama_algebra, nakayama_algebra_sequence_first,
    };

    fn integers(values: &[i64]) -> Vec<Integer> {
        values.iter().copied().map(Integer::from).collect()
    }

    fn small_integer(value: &Integer) -> i64 {
        value
            .to_string()
            .parse()
            .expect("the predefined-algebra fixture uses only small integers")
    }

    fn rational_coefficient(value: &<RationalField as numerica::domains::Set>::Element) -> Value {
        json!([
            small_integer(value.numerator_ref()),
            small_integer(value.denominator_ref())
        ])
    }

    fn element_json<F: Field>(
        element: &PathAlgebraElement<F>,
        coefficient_json: &impl Fn(&F::Element) -> Value,
    ) -> Value {
        Value::Array(
            element
                .terms()
                .iter()
                .map(|term| {
                    json!({
                        "path": term.path().external_representation(),
                        "coefficient": coefficient_json(term.coefficient())
                    })
                })
                .collect(),
        )
    }

    fn algebra_json<F: Field>(
        algebra: &PredefinedAlgebra<F>,
        coefficient_json: &impl Fn(&F::Element) -> Value,
    ) -> Value {
        let original = algebra.algebra().original_path_algebra();
        let (admissible_quotient, dimension) = match algebra.algebra() {
            QuiverAlgebra::Path(path_algebra) => (
                false,
                path_algebra
                    .canonical_basis()
                    .expect("oracle predefined path algebra is finite")
                    .len(),
            ),
            QuiverAlgebra::Quotient(quotient) => {
                let AlgebraDimension::Finite(dimension) = quotient.dimension().unwrap() else {
                    panic!("oracle predefined quotient is finite");
                };
                (
                    quotient.is_admissible().unwrap(),
                    usize::try_from(dimension).unwrap(),
                )
            }
        };
        json!({
            "path_algebra": algebra.algebra().is_path_algebra(),
            "admissible_quotient": admissible_quotient,
            "dimension": dimension,
            "vertices": original
                .quiver()
                .vertices()
                .map(|vertex| vertex.name().to_owned())
                .collect::<Vec<_>>(),
            "arrows": original
                .quiver()
                .arrows()
                .map(|arrow| json!({
                    "name": arrow.name(),
                    "source": arrow.source().index().get(),
                    "target": arrow.target().index().get()
                }))
                .collect::<Vec<_>>(),
            "relations": algebra
                .algebra()
                .relations()
                .iter()
                .map(|relation| element_json(relation, coefficient_json))
                .collect::<Vec<_>>()
        })
    }

    fn constructed<F: Field>(result: NakayamaAlgebraResult<F>) -> PredefinedAlgebra<F> {
        let NakayamaAlgebraResult::Algebra(algebra) = result else {
            panic!("fixture sequence must be admissible");
        };
        algebra
    }

    #[test]
    fn invalid_and_error_branches_follow_qpa_dispatch() {
        assert_eq!(
            nakayama_algebra(Q, &[]),
            Err(PredefinedAlgebraError::EmptyAdmissibleSequence)
        );

        let invalid_sequence = integers(&[1, 1]);
        let invalid = nakayama_algebra(Q, &invalid_sequence).unwrap();
        assert_eq!(
            invalid,
            NakayamaAlgebraResult::InvalidSequence(invalid_sequence.clone())
        );
        assert_eq!(
            invalid.diagnostic(),
            Some(INVALID_ADMISSIBLE_SEQUENCE_DIAGNOSTIC)
        );
        assert!(invalid.as_algebra().is_none());
        assert!(matches!(
            nakayama_algebra_sequence_first(&invalid_sequence, Q).unwrap(),
            NakayamaAlgebraResult::InvalidSequence(sequence) if sequence == invalid_sequence
        ));
        assert!(matches!(
            nakayama_algebra(Q, &integers(&[0])),
            Err(PredefinedAlgebraError::Gbnp(_))
        ));
        assert!(matches!(
            nakayama_algebra(Q, &integers(&[-2])),
            Err(PredefinedAlgebraError::Gbnp(_))
        ));

        let two = Q.nth(Integer::from(2));
        assert_eq!(
            canonical_algebra(Q, &integers(&[2, 2]), std::slice::from_ref(&two)),
            Err(PredefinedAlgebraError::CanonicalArmCoefficientMismatch {
                arms: 2,
                coefficients: 1
            })
        );
        assert_eq!(
            canonical_algebra(Q, &integers(&[2, 1]), &[]),
            Err(PredefinedAlgebraError::CanonicalWeightBelowTwo {
                position: 2,
                weight: Integer::from(1)
            })
        );
        assert_eq!(
            canonical_algebra(Q, &integers(&[2, 2, 2]), &[Q.zero()]),
            Err(PredefinedAlgebraError::ZeroCanonicalRelationCoefficient { position: 1 })
        );
        assert_eq!(
            canonical_algebra_two_weights(Q, &integers(&[2, 2, 2])),
            Err(PredefinedAlgebraError::CanonicalTwoArgumentWeightCount { weights: 3 })
        );
        assert_eq!(
            kronecker_algebra(Q, &Integer::from(1)),
            Err(PredefinedAlgebraError::KroneckerNeedsTwoArrows {
                number_of_arrows: Integer::from(1)
            })
        );

        let nakayama = constructed(nakayama_algebra(Q, &integers(&[1])).unwrap());
        assert_eq!(nakayama.is_nakayama_algebra(), Some(true));
        assert_eq!(nakayama.is_canonical_algebra(), None);
        assert_eq!(nakayama.is_kronecker_algebra(), None);
        assert_eq!(nakayama.is_finite_global_dimension_algebra(), None);
        assert_eq!(nakayama.is_hereditary_algebra(), None);
    }

    #[test]
    fn pinned_qpa_predefined_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_predefined.schema-v1.json"
        ))
        .unwrap();

        let nakayama_singleton = constructed(nakayama_algebra(Q, &integers(&[1])).unwrap());
        let nakayama_singleton_quotient =
            constructed(nakayama_algebra(Q, &integers(&[2])).unwrap());
        let nakayama_linear_path = constructed(nakayama_algebra(Q, &integers(&[3, 2, 1])).unwrap());
        let nakayama_linear_quotient =
            constructed(nakayama_algebra(Q, &integers(&[3, 3, 2, 1])).unwrap());
        let nakayama_cyclic = constructed(nakayama_algebra(Q, &integers(&[3, 3, 3])).unwrap());
        let nakayama_long_cyclic =
            constructed(nakayama_algebra(Q, &integers(&[7, 6, 7, 7])).unwrap());
        let canonical_two = canonical_algebra_two_weights(Q, &integers(&[2, 3])).unwrap();
        let canonical_four = canonical_algebra(
            Q,
            &integers(&[2, 2, 3, 4]),
            &[Q.nth(Integer::from(2)), Q.nth(Integer::from(-3))],
        )
        .unwrap();

        let rational = json!({
            "nakayama_singleton": algebra_json(&nakayama_singleton, &rational_coefficient),
            "nakayama_singleton_quotient": algebra_json(
                &nakayama_singleton_quotient,
                &rational_coefficient
            ),
            "nakayama_linear_path": algebra_json(&nakayama_linear_path, &rational_coefficient),
            "nakayama_linear_quotient": algebra_json(
                &nakayama_linear_quotient,
                &rational_coefficient
            ),
            "nakayama_cyclic": algebra_json(&nakayama_cyclic, &rational_coefficient),
            "nakayama_long_cyclic": algebra_json(&nakayama_long_cyclic, &rational_coefficient),
            "canonical_two": algebra_json(&canonical_two, &rational_coefficient),
            "canonical_four": algebra_json(&canonical_four, &rational_coefficient),
            "properties": {
                "nakayama": [
                    nakayama_cyclic.is_nakayama_algebra().unwrap(),
                    nakayama_cyclic.is_finite_type_algebra().unwrap()
                ],
                "canonical": [
                    canonical_four.is_canonical_algebra().unwrap(),
                    canonical_four
                        .is_finite_global_dimension_algebra()
                        .unwrap()
                ]
            }
        });

        let nakayama_cyclic_f2 = constructed(nakayama_algebra(Z2, &integers(&[3, 3, 3])).unwrap());
        let canonical_f2 =
            canonical_algebra(Z2, &integers(&[2, 3, 4]), std::slice::from_ref(&Z2.one())).unwrap();
        let kronecker_f2 = kronecker_algebra(Z2, &Integer::from(3)).unwrap();
        let z2_coefficient = |coefficient: &u8| Value::from(*coefficient);
        let f2 = json!({
            "nakayama_cyclic": algebra_json(&nakayama_cyclic_f2, &z2_coefficient),
            "canonical": algebra_json(&canonical_f2, &z2_coefficient),
            "kronecker": algebra_json(&kronecker_f2, &z2_coefficient),
            "properties": {
                "nakayama": [
                    nakayama_cyclic_f2.is_nakayama_algebra().unwrap(),
                    nakayama_cyclic_f2.is_finite_type_algebra().unwrap()
                ],
                "canonical": [
                    canonical_f2.is_canonical_algebra().unwrap(),
                    canonical_f2
                        .is_finite_global_dimension_algebra()
                        .unwrap()
                ],
                "kronecker": [
                    kronecker_f2.is_kronecker_algebra().unwrap(),
                    kronecker_f2
                        .is_finite_global_dimension_algebra()
                        .unwrap(),
                    kronecker_f2.is_hereditary_algebra().unwrap()
                ]
            }
        });

        assert_eq!(rational, expected["rational"]);
        assert_eq!(f2, expected["f2"]);
    }
}
