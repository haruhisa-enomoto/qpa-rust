//! QPA path-algebra constructors whose configured default backend is GBNP.
//!
//! `qpa-core` owns quotient representation and backend-neutral factor logic;
//! this module follows QPA's default call chain through `GBNPGroebnerBasis`.

use numerica::domains::Field;
use qpa_core::{
    AlgebraDimension, Arrow, CentreError, ElementInverseError, GroebnerBasis, PathAlgebra,
    PathAlgebraElement, PathAlgebraError, PathAlgebraIdeal, Quiver, QuiverAlgebra, QuiverError,
    QuotientError, QuotientPathAlgebra, QuotientPathAlgebraCentre, QuotientPathAlgebraElement,
    TriangularReducedError, VertexQuotientPresentation,
};
use thiserror::Error;

use crate::{
    InvLenLexError, QpaGbnpError, biggest_homogeneous_parts_inv_len_inv_lex, gbnp_groebner_basis,
    groebner_basis_inv_len_inv_lex,
};

/// Translate QPA's `/` method for a path algebra and relation list using the
/// default `GBNPGroebnerBasis` backend.
///
/// # Errors
///
/// Propagates the checked QPA-to-GBNP conversion, `SGrobner`, and quotient
/// construction failures.
pub fn quotient_by_relations<F: Field>(
    parent: PathAlgebra<F>,
    relators: &[PathAlgebraElement<F>],
) -> Result<QuotientPathAlgebra<F>, PathalgGbnpError> {
    let complete_relations = gbnp_groebner_basis(relators, &parent)?;
    Ok(QuotientPathAlgebra::by_complete_groebner_relations(
        parent,
        &complete_relations,
    )?)
}

/// Translate `QuiverAlgebraOfAmodAeA(A, elist)` for a path algebra through
/// QPA's default GBNP backend.
///
/// QPA returns the rebuilt path algebra itself because its relation list is
/// empty; the [`QuiverAlgebra`] result preserves that type distinction.
///
/// # Errors
///
/// Propagates full-subquiver reconstruction and path conversion failures.
pub fn quiver_algebra_of_path_algebra_mod_aea<F: Field>(
    algebra: &PathAlgebra<F>,
    vertex_positions: &[usize],
) -> Result<QuiverAlgebra<F>, PathalgGbnpError> {
    quotient_vertex_presentation(algebra.a_mod_aea_presentation(vertex_positions)?)
}

/// Translate `QuiverAlgebraOfAmodAeA(A, elist)` for a quotient through QPA's
/// default GBNP backend.
///
/// # Errors
///
/// Propagates `PathRemoval`, named-quiver conversion, GBNP, and quotient
/// construction failures.
pub fn quiver_algebra_of_quotient_mod_aea<F: Field>(
    algebra: &QuotientPathAlgebra<F>,
    vertex_positions: &[usize],
) -> Result<QuiverAlgebra<F>, PathalgGbnpError> {
    quotient_vertex_presentation(algebra.a_mod_aea_presentation(vertex_positions)?)
}

fn quotient_vertex_presentation<F: Field>(
    presentation: VertexQuotientPresentation<F>,
) -> Result<QuiverAlgebra<F>, PathalgGbnpError> {
    let (path_algebra, relations) = presentation.into_parts();
    if relations.is_empty() {
        Ok(QuiverAlgebra::Path(path_algebra))
    } else {
        Ok(QuiverAlgebra::Quotient(quotient_by_relations(
            path_algebra,
            &relations,
        )?))
    }
}

/// Compute and attach QPA's configured default Gröbner basis when the
/// quotient's defining ideal does not already carry one.
///
/// The attachment mutates the identity-preserving ideal cache. It deliberately
/// does not retrofit a normal-form function onto the quotient family, matching
/// QPA's lazy `Dimension` and `IsFiniteDimensional` methods.
///
/// # Errors
///
/// Propagates the checked QPA-to-GBNP conversion and basis construction.
pub fn ensure_quotient_groebner_basis<F: Field>(
    quotient: &QuotientPathAlgebra<F>,
) -> Result<GroebnerBasis<F>, PathalgGbnpError> {
    if let Some(basis) = quotient.ideal().groebner_basis() {
        return Ok(basis);
    }
    let relations = gbnp_groebner_basis(quotient.relations(), quotient.original_path_algebra())?;
    quotient.ideal().with_complete_groebner_basis(&relations)?;
    quotient
        .ideal()
        .groebner_basis()
        .ok_or(PathalgGbnpError::MissingAttachedBasis)
}

/// Translate QPA's lazy default `IsFiniteDimensional` method.
///
/// # Errors
///
/// Propagates default GBNP and nontip-dictionary failures.
pub fn quotient_is_finite_dimensional<F: Field>(
    quotient: &QuotientPathAlgebra<F>,
) -> Result<bool, PathalgGbnpError> {
    ensure_quotient_groebner_basis(quotient)?;
    Ok(quotient.is_finite_dimensional()?)
}

/// Translate QPA's lazy default `Dimension` method.
///
/// # Errors
///
/// Propagates default GBNP and nontip-dictionary failures.
pub fn quotient_dimension<F: Field>(
    quotient: &QuotientPathAlgebra<F>,
) -> Result<AlgebraDimension, PathalgGbnpError> {
    ensure_quotient_groebner_basis(quotient)?;
    Ok(quotient.dimension()?)
}

/// Translate QPA's `Inverse(elm)` through its lazy default GBNP
/// finite-dimensionality call and source polynomial.
///
/// `Ok(None)` is QPA's `fail` result for a nonunit.
///
/// # Errors
///
/// Propagates default GBNP, basis-coordinate, and quotient-arithmetic failures.
pub fn inverse_quotient_element<F: Field>(
    element: &QuotientPathAlgebraElement<F>,
) -> Result<Option<QuotientPathAlgebraElement<F>>, PathalgGbnpError> {
    ensure_quotient_groebner_basis(element.algebra())?;
    Ok(element.inverse()?)
}

/// Translate QPA's negative and zero powers for quotient elements.
///
/// `Ok(None)` is QPA's `fail` for a negative power of a nonunit.
///
/// # Errors
///
/// Propagates default GBNP, inversion, and quotient-arithmetic failures.
pub fn power_quotient_element<F: Field>(
    element: &QuotientPathAlgebraElement<F>,
    exponent: i64,
) -> Result<Option<QuotientPathAlgebraElement<F>>, PathalgGbnpError> {
    if exponent < 0 {
        ensure_quotient_groebner_basis(element.algebra())?;
    }
    Ok(element.pow_integer(exponent)?)
}

/// Translate QPA's quotient `Centre` method through its lazy default GBNP
/// dimensionality call.
///
/// # Errors
///
/// Propagates default GBNP, canonical-basis, commutator, and nullspace
/// failures.
pub fn quotient_centre<F: Field>(
    quotient: &QuotientPathAlgebra<F>,
) -> Result<QuotientPathAlgebraCentre<F>, PathalgGbnpError> {
    ensure_quotient_groebner_basis(quotient)?;
    Ok(quotient.centre()?)
}

/// Translate QPA's quotient `IsTriangularReduced` method through its lazy
/// default GBNP finite-dimensionality call.
///
/// # Errors
///
/// Propagates default GBNP, canonical-basis, and quotient-arithmetic failures.
pub fn quotient_is_triangular_reduced<F: Field>(
    quotient: &QuotientPathAlgebra<F>,
) -> Result<bool, PathalgGbnpError> {
    ensure_quotient_groebner_basis(quotient)?;
    Ok(quotient.is_triangular_reduced()?)
}

/// Translate QPA's `IsMonomialAlgebra` method for a path algebra.
///
/// QPA returns `true` without inspecting the quiver or computing a basis.
#[must_use]
pub const fn path_algebra_is_monomial<F: Field>(_: &PathAlgebra<F>) -> bool {
    true
}

/// Translate QPA's `IsMonomialAlgebra` method for a quotient through its
/// configured default GBNP backend.
///
/// QPA first accepts a trivial ideal or monomial defining generators without
/// computing a Gröbner basis. Only a nonmonomial generator list reaches the
/// lazy default backend.
///
/// # Errors
///
/// Propagates default GBNP and source-translated ideal failures.
pub fn quotient_is_monomial<F: Field>(
    quotient: &QuotientPathAlgebra<F>,
) -> Result<bool, PathalgGbnpError> {
    if quotient.ideal().is_trivial() || quotient.ideal().generators_are_monomial() {
        return Ok(true);
    }
    ensure_quotient_groebner_basis(quotient)?;
    Ok(quotient.ideal().is_monomial()?)
}

/// Translate QPA's lazy `IsAdmissibleQuotientOfPathAlgebra` property through
/// the configured default GBNP backend.
///
/// # Errors
///
/// Propagates default GBNP and the source-translated admissible-ideal test.
pub fn quotient_is_admissible<F: Field>(
    quotient: &QuotientPathAlgebra<F>,
) -> Result<bool, PathalgGbnpError> {
    ensure_quotient_groebner_basis(quotient)?;
    Ok(quotient.is_admissible()?)
}

/// Translate QPA's `MinimalGeneratingSetOfIdeal` through the
/// `GBNPGroebnerBasis` function installed by the default path-algebra
/// constructor.
///
/// # Errors
///
/// Propagates default GBNP, admissibility, quotient, basis, and exact
/// linear-algebra failures.
pub fn ideal_minimal_generating_set<F: Field>(
    ideal: &PathAlgebraIdeal<F>,
) -> Result<Vec<PathAlgebraElement<F>>, PathalgGbnpError> {
    if ideal.groebner_basis().is_none() {
        let relations = gbnp_groebner_basis(ideal.generators(), ideal.algebra())?;
        ideal.with_complete_groebner_basis(&relations)?;
    }
    ideal.minimal_generating_set_with(quotient_by_relations)
}

/// Translate QPA's `IsSpecialBiserialQuiver` degree test.
///
/// Every vertex must have at most two incoming and at most two outgoing
/// arrows. Parallel arrows and loops are counted exactly as QPA counts them.
#[must_use]
pub fn is_special_biserial_quiver(quiver: &Quiver) -> bool {
    quiver
        .vertices()
        .all(|vertex| vertex.in_degree() <= 2 && vertex.out_degree() <= 2)
}

/// Translate QPA's `IsSpecialBiserialAlgebra` method for a path algebra.
///
/// In addition to the quiver degree bound, every arrow has at most one
/// continuation on either side and the quiver must be acyclic.
#[must_use]
pub fn path_algebra_is_special_biserial<F: Field>(algebra: &PathAlgebra<F>) -> bool {
    let quiver = algebra.quiver();
    is_special_biserial_quiver(quiver)
        && quiver
            .arrows()
            .all(|arrow| arrow.target().out_degree() <= 1 && arrow.source().in_degree() <= 1)
        && quiver.is_acyclic()
}

/// Translate QPA's `IsSpecialBiserialAlgebra` method for a quotient through
/// its configured default GBNP backend.
///
/// # Errors
///
/// Propagates default admissibility, ideal membership, and path-arithmetic
/// failures.
pub fn quotient_is_special_biserial<F: Field>(
    quotient: &QuotientPathAlgebra<F>,
) -> Result<bool, PathalgGbnpError> {
    let quiver = quotient.original_path_algebra().quiver();
    if !is_special_biserial_quiver(quiver) || !quotient_is_admissible(quotient)? {
        return Ok(false);
    }
    for arrow in quiver.arrows() {
        let mut surviving_right = 0_usize;
        for right in quiver.outgoing_arrows(&arrow.target())? {
            if !arrow_product_is_in_ideal(quotient, &arrow, &right)? {
                surviving_right += 1;
            }
        }
        if surviving_right > 1 {
            return Ok(false);
        }

        let mut surviving_left = 0_usize;
        for left in quiver.incoming_arrows(&arrow.source())? {
            if !arrow_product_is_in_ideal(quotient, &left, &arrow)? {
                surviving_left += 1;
            }
        }
        if surviving_left > 1 {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Translate QPA's `IsStringAlgebra` method for a path algebra.
#[must_use]
pub fn path_algebra_is_string<F: Field>(algebra: &PathAlgebra<F>) -> bool {
    path_algebra_is_special_biserial(algebra)
}

/// Translate QPA's `IsStringAlgebra` method for a quotient.
///
/// # Errors
///
/// Propagates the default special-biserial and monomial-ideal computations.
pub fn quotient_is_string<F: Field>(
    quotient: &QuotientPathAlgebra<F>,
) -> Result<bool, PathalgGbnpError> {
    if !quotient_is_special_biserial(quotient)? {
        return Ok(false);
    }
    quotient_is_monomial(quotient)
}

/// Translate QPA's `IsGentleAlgebra` method for a path algebra.
#[must_use]
pub fn path_algebra_is_gentle<F: Field>(algebra: &PathAlgebra<F>) -> bool {
    path_algebra_is_special_biserial(algebra)
}

/// Translate QPA's `IsGentleAlgebra` method for a quotient through its default
/// GBNP quotient backend.
///
/// QPA requires a special-biserial admissible quotient whose minimal ideal
/// generators are monomial paths of length two, with at most one forbidden
/// continuation on either side of each arrow.
///
/// # Errors
///
/// Propagates default GBNP, admissibility, minimal-generator, membership, and
/// path-arithmetic failures.
pub fn quotient_is_gentle<F: Field>(
    quotient: &QuotientPathAlgebra<F>,
) -> Result<bool, PathalgGbnpError> {
    if !quotient_is_special_biserial(quotient)? || !quotient_is_admissible(quotient)? {
        return Ok(false);
    }
    let minimal_generators = ideal_minimal_generating_set(quotient.ideal())?;
    if minimal_generators.iter().any(|relation| {
        relation.terms().len() != 1 || relation.terms()[0].path().length() != Some(2)
    }) {
        return Ok(false);
    }

    let quiver = quotient.original_path_algebra().quiver();
    for arrow in quiver.arrows() {
        let mut forbidden_right = 0_usize;
        for right in quiver.outgoing_arrows(&arrow.target())? {
            if arrow_product_is_in_ideal(quotient, &arrow, &right)? {
                forbidden_right += 1;
            }
        }
        if forbidden_right > 1 {
            return Ok(false);
        }

        let mut forbidden_left = 0_usize;
        for left in quiver.incoming_arrows(&arrow.source())? {
            if arrow_product_is_in_ideal(quotient, &left, &arrow)? {
                forbidden_left += 1;
            }
        }
        if forbidden_left > 1 {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Translate QPA's true-method implication from `IsGentleAlgebra` to
/// `IsSpecialBiserialAlgebra`.
///
/// # Errors
///
/// Returns a typed dispatch error when the algebra is not gentle and
/// propagates the quotient predicate's default GBNP failures.
pub fn gentle_algebra_is_special_biserial<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<bool, PathalgGbnpError> {
    require_gentle_algebra(algebra)?;
    Ok(true)
}

/// Translate QPA's true-method implication from `IsGentleAlgebra` to
/// `IsGorensteinAlgebra`.
///
/// # Errors
///
/// Returns a typed dispatch error when the algebra is not gentle and
/// propagates the quotient predicate's default GBNP failures.
pub fn gentle_algebra_is_gorenstein<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<bool, PathalgGbnpError> {
    require_gentle_algebra(algebra)?;
    Ok(true)
}

fn require_gentle_algebra<F: Field>(algebra: &QuiverAlgebra<F>) -> Result<(), PathalgGbnpError> {
    let gentle = match algebra {
        QuiverAlgebra::Path(path_algebra) => path_algebra_is_gentle(path_algebra),
        QuiverAlgebra::Quotient(quotient) => quotient_is_gentle(quotient)?,
    };
    if gentle {
        Ok(())
    } else {
        Err(PathalgGbnpError::NonGentleAlgebraProperty)
    }
}

fn arrow_product_is_in_ideal<F: Field>(
    quotient: &QuotientPathAlgebra<F>,
    left: &Arrow,
    right: &Arrow,
) -> Result<bool, PathalgGbnpError> {
    let parent = quotient.original_path_algebra();
    let product = parent
        .element(left.path())?
        .multiply(&parent.element(right.path())?)?;
    Ok(quotient.ideal().contains(product)?)
}

fn require_admissible_quotient<F: Field>(
    quotient: &QuotientPathAlgebra<F>,
) -> Result<(), PathalgGbnpError> {
    if quotient_is_admissible(quotient)? {
        Ok(())
    } else {
        Err(PathalgGbnpError::NonAdmissibleQuotientProperty)
    }
}

/// Translate QPA's specialized `IsSemisimpleAlgebra` method selected by the
/// `IsAdmissibleQuotientOfPathAlgebra` filter.
///
/// Once that dispatch filter is established, QPA returns whether the original
/// quiver has no arrows. This entry point validates the filter before applying
/// the same test; QPA's separate radical-based generic overload is not used.
///
/// # Errors
///
/// Returns a typed dispatch error for a nonadmissible quotient and propagates
/// the default admissibility computation.
pub fn admissible_quotient_is_semisimple<F: Field>(
    quotient: &QuotientPathAlgebra<F>,
) -> Result<bool, PathalgGbnpError> {
    require_admissible_quotient(quotient)?;
    Ok(quotient.original_path_algebra().quiver().number_of_arrows() == 0)
}

/// Translate QPA's true-method implication from an admissible path-algebra
/// quotient to `IsElementaryAlgebra`.
///
/// # Errors
///
/// Returns a typed dispatch error for a nonadmissible quotient and propagates
/// the default admissibility computation.
pub fn admissible_quotient_is_elementary<F: Field>(
    quotient: &QuotientPathAlgebra<F>,
) -> Result<bool, PathalgGbnpError> {
    require_admissible_quotient(quotient)?;
    Ok(true)
}

/// Translate QPA's true-method implication `IsElementaryAlgebra` to
/// `IsBasicAlgebra` for an admissible path-algebra quotient.
///
/// # Errors
///
/// Returns a typed dispatch error for a nonadmissible quotient and propagates
/// the default admissibility computation.
pub fn admissible_quotient_is_basic<F: Field>(
    quotient: &QuotientPathAlgebra<F>,
) -> Result<bool, PathalgGbnpError> {
    admissible_quotient_is_elementary(quotient)
}

/// Translate `AssociatedMonomialAlgebra(A)` for a path algebra.
///
/// QPA returns the identical path-algebra object.
#[must_use]
pub fn associated_monomial_path_algebra<F: Field>(algebra: &PathAlgebra<F>) -> QuiverAlgebra<F> {
    QuiverAlgebra::Path(algebra.clone())
}

/// Translate `AssociatedMonomialAlgebra(A)` for a quotient through QPA's
/// configured default GBNP backend.
///
/// An already monomial quotient is returned identically. Otherwise QPA takes
/// `Tip` of every relation in the attached default Gröbner-basis object and
/// factors the original path algebra by those tips.
///
/// # Errors
///
/// Propagates the monomial test, default GBNP, and quotient construction.
pub fn associated_monomial_quotient<F: Field>(
    quotient: &QuotientPathAlgebra<F>,
) -> Result<QuiverAlgebra<F>, PathalgGbnpError> {
    if quotient_is_monomial(quotient)? {
        return Ok(QuiverAlgebra::Quotient(quotient.clone()));
    }
    let basis = ensure_quotient_groebner_basis(quotient)?;
    let relations = basis
        .relations()
        .iter()
        .map(PathAlgebraElement::leading_term)
        .collect::<Vec<_>>();
    Ok(QuiverAlgebra::Quotient(quotient_by_relations(
        quotient.original_path_algebra().clone(),
        &relations,
    )?))
}

/// Translate `AssociatedGradedAlgebra(A)` for a path algebra.
///
/// QPA checks finite dimensionality before its path-algebra identity branch.
/// Thus an acyclic parent is returned identically and a cyclic parent fails.
///
/// # Errors
///
/// Returns QPA's finite-dimensionality failure for a cyclic quiver.
pub fn associated_graded_path_algebra<F: Field>(
    algebra: &PathAlgebra<F>,
) -> Result<QuiverAlgebra<F>, PathalgGbnpError> {
    if !algebra.has_finite_canonical_basis() {
        return Err(PathalgGbnpError::InfiniteDimensionalAssociatedAlgebra);
    }
    Ok(QuiverAlgebra::Path(algebra.clone()))
}

/// Translate `AssociatedGradedAlgebra(A)` for a quotient through QPA's
/// configured default GBNP and inverse-length/inverse-lexicographic backends.
///
/// QPA requires a finite-dimensional admissible quotient, computes
/// `GroebnerBasisInvLenInvLex` from `RelationsOfAlgebra(A)`, retains the first
/// homogeneous block of every resulting relation, and factors the original
/// path algebra by that list. The source's cached-false admissibility bug is
/// retained: a bound property skips the rejection branch regardless of its
/// value and proceeds to the inverse-length backend.
///
/// # Errors
///
/// Propagates the finite/admissible preconditions, both source-translated
/// Gröbner computations, and quotient construction.
pub fn associated_graded_quotient<F: Field>(
    quotient: &QuotientPathAlgebra<F>,
) -> Result<QuiverAlgebra<F>, PathalgGbnpError> {
    ensure_quotient_groebner_basis(quotient)?;
    if !quotient.is_finite_dimensional()? {
        return Err(PathalgGbnpError::InfiniteDimensionalAssociatedAlgebra);
    }
    if quotient.cached_admissibility().is_none() && !quotient.is_admissible()? {
        return Err(PathalgGbnpError::NonAdmissibleAssociatedAlgebra);
    }
    let parent = quotient.original_path_algebra();
    let inverse_basis = groebner_basis_inv_len_inv_lex(quotient.relations(), parent)?;
    let relations = biggest_homogeneous_parts_inv_len_inv_lex(&inverse_basis)?;
    Ok(QuiverAlgebra::Quotient(quotient_by_relations(
        parent.clone(),
        &relations,
    )?))
}

/// Translate QPA's public `TruncatedPathAlgebra(K, Q, n)` through its default
/// GBNP backend.
///
/// # Errors
///
/// Propagates parent construction, `NthPowerOfArrowIdeal`, GBNP, and quotient
/// failures.
pub fn truncated_path_algebra<F: Field>(
    field: F,
    quiver: Quiver,
    exponent: i64,
) -> Result<QuotientPathAlgebra<F>, PathalgGbnpError> {
    let parent = PathAlgebra::new(field, quiver)?;
    let relators = parent.nth_power_of_arrow_ideal(exponent)?;
    quotient_by_relations(parent, &relators)
}

/// Translate QPA's `OppositePathAlgebra` for a quotient using the default
/// GBNP backend of the newly created opposite parent.
///
/// # Errors
///
/// Propagates opposite-relator conversion, GBNP, quotient, and reciprocal
/// cache failures.
pub fn opposite_path_algebra<F: Field>(
    quotient: &QuotientPathAlgebra<F>,
) -> Result<QuotientPathAlgebra<F>, PathalgGbnpError> {
    if let Some(opposite) = quotient.cached_opposite() {
        return Ok(opposite);
    }
    let opposite_parent = quotient.original_path_algebra().opposite()?;
    let opposite_relators = quotient.opposite_defining_relations()?;
    let complete_relations = gbnp_groebner_basis(&opposite_relators, &opposite_parent)?;
    let opposite = quotient
        .install_opposite_from_complete_relations(&opposite_relators, &complete_relations)?;
    if quotient.is_admissible()? {
        opposite.mark_admissible();
    }
    Ok(opposite)
}

/// Translate `OppositePathAlgebraElement` for a quotient element through the
/// default GBNP-backed opposite quotient.
///
/// # Errors
///
/// Propagates opposite quotient, termwise opposite, and normal-form failures.
pub fn opposite_quotient_element<F: Field>(
    element: &QuotientPathAlgebraElement<F>,
) -> Result<QuotientPathAlgebraElement<F>, PathalgGbnpError> {
    let opposite = opposite_path_algebra(element.algebra())?;
    Ok(opposite.element(element.representative().opposite()?)?)
}

/// Failure in QPA path-algebra constructors backed by the GBNP translation.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum PathalgGbnpError {
    /// An attached basis disappeared after the identity-preserving update.
    #[error("the quotient ideal did not retain its computed Groebner basis")]
    MissingAttachedBasis,
    /// `AssociatedGradedAlgebra` rejects an infinite-dimensional algebra.
    #[error("the entered algebra is not finite dimensional")]
    InfiniteDimensionalAssociatedAlgebra,
    /// `AssociatedGradedAlgebra` rejects a nonadmissible quotient.
    #[error("the entered algebra is not an admissible quotient of a path algebra")]
    NonAdmissibleAssociatedAlgebra,
    /// An admissible-quotient-only property received a quotient outside that
    /// GAP dispatch filter.
    #[error("the quotient is not an admissible quotient of a path algebra")]
    NonAdmissibleQuotientProperty,
    /// A gentle-algebra-only true-method implication received another algebra.
    #[error("the algebra is not a gentle algebra")]
    NonGentleAlgebraProperty,
    /// QPA's source-translated center construction failed.
    #[error(transparent)]
    Centre(#[from] CentreError),
    /// QPA's finite-dimensional inverse polynomial failed.
    #[error(transparent)]
    ElementInverse(#[from] ElementInverseError),
    /// QPA's triangular-reduced subset computation failed.
    #[error(transparent)]
    TriangularReduced(#[from] TriangularReducedError),
    /// QPA's inverse-length/inverse-lexicographic backend failed.
    #[error(transparent)]
    InvLenLex(#[from] InvLenLexError),
    /// Parent or opposite path-algebra construction failed.
    #[error(transparent)]
    PathAlgebra(#[from] PathAlgebraError),
    /// A source-translated quiver incidence lookup failed.
    #[error(transparent)]
    Quiver(#[from] QuiverError),
    /// QPA's checked GBNP adapter or `SGrobner` failed.
    #[error(transparent)]
    Gbnp(#[from] QpaGbnpError),
    /// Backend-neutral quotient construction failed.
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
    use qpa_core::{
        AlgebraDimension, ArrowIndex, ArrowSpec, PathAlgebra, PathAlgebraElement, PathAlgebraIdeal,
        Quiver, QuiverAlgebra, QuotientPathAlgebra, VertexIndex,
    };
    use serde_json::{Value, json};

    use super::{
        PathalgGbnpError, admissible_quotient_is_basic, admissible_quotient_is_elementary,
        admissible_quotient_is_semisimple, associated_graded_path_algebra,
        associated_graded_quotient, associated_monomial_path_algebra, associated_monomial_quotient,
        gentle_algebra_is_gorenstein, gentle_algebra_is_special_biserial,
        ideal_minimal_generating_set, inverse_quotient_element, is_special_biserial_quiver,
        opposite_path_algebra, opposite_quotient_element, path_algebra_is_gentle,
        path_algebra_is_monomial, path_algebra_is_special_biserial, path_algebra_is_string,
        power_quotient_element, quiver_algebra_of_path_algebra_mod_aea,
        quiver_algebra_of_quotient_mod_aea, quotient_by_relations, quotient_centre,
        quotient_dimension, quotient_is_admissible, quotient_is_finite_dimensional,
        quotient_is_gentle, quotient_is_monomial, quotient_is_special_biserial, quotient_is_string,
        quotient_is_triangular_reduced, truncated_path_algebra,
    };

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn arrow(index: usize) -> ArrowIndex {
        ArrowIndex::new(index).unwrap()
    }

    fn quiver() -> Quiver {
        Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(1), "a"),
                ArrowSpec::named(vertex(1), vertex(2), "b"),
                ArrowSpec::named(vertex(2), vertex(3), "c"),
                ArrowSpec::named(vertex(3), vertex(1), "d"),
                ArrowSpec::named(vertex(2), vertex(2), "e"),
            ],
        )
        .unwrap()
    }

    fn small_integer(value: &Integer) -> i64 {
        value
            .to_string()
            .parse()
            .expect("the pathalg fixture uses only small integers")
    }

    fn rational_element_json(element: &PathAlgebraElement<RationalField>) -> Value {
        Value::Array(
            element
                .terms()
                .iter()
                .map(|term| {
                    json!({
                        "path": term.path().external_representation(),
                        "coefficient": [
                            small_integer(term.coefficient().numerator_ref()),
                            small_integer(term.coefficient().denominator_ref())
                        ]
                    })
                })
                .collect(),
        )
    }

    fn z2_element_json(element: &PathAlgebraElement<Z2>) -> Value {
        Value::Array(
            element
                .terms()
                .iter()
                .map(|term| {
                    json!({
                        "path": term.path().external_representation(),
                        "coefficient": term.coefficient()
                    })
                })
                .collect(),
        )
    }

    fn quiver_json(quiver: &Quiver) -> Value {
        json!({
            "vertices": quiver.vertices().map(|vertex| vertex.name().to_owned()).collect::<Vec<_>>(),
            "arrows": quiver.arrows().map(|arrow| json!({
                "source": arrow.source().name(),
                "target": arrow.target().name(),
                "name": arrow.name()
            })).collect::<Vec<_>>()
        })
    }

    fn finite_dimension<F: numerica::domains::Field>(algebra: &QuiverAlgebra<F>) -> i64 {
        let AlgebraDimension::Finite(dimension) = algebra
            .as_quotient()
            .expect("the oracle presentation has nonempty relations")
            .dimension()
            .unwrap()
        else {
            panic!("the oracle vertex quotient must be finite");
        };
        small_integer(&dimension)
    }

    fn special_predicates<F: Field>(algebra: &QuotientPathAlgebra<F>) -> Value {
        json!([
            quotient_is_special_biserial(algebra).unwrap(),
            quotient_is_string(algebra).unwrap(),
            quotient_is_gentle(algebra).unwrap()
        ])
    }

    #[allow(clippy::many_single_char_names)]
    fn special_field_fixture<F: Field>(field: F) -> Value {
        let two_loop_quiver = Quiver::from_indexed(
            1,
            vec![
                ArrowSpec::named(vertex(1), vertex(1), "a"),
                ArrowSpec::named(vertex(1), vertex(1), "b"),
            ],
        )
        .unwrap();
        let parent = PathAlgebra::new(field.clone(), two_loop_quiver).unwrap();
        let generators = parent.generators();
        let a = &generators[1];
        let b = &generators[2];
        let ab = a.multiply(b).unwrap();
        let ba = b.multiply(a).unwrap();
        let a2 = a.pow(2).unwrap();
        let b2 = b.pow(2).unwrap();
        let b3 = b.pow(3).unwrap();
        let q1 =
            quotient_by_relations(parent.clone(), &[ab.clone(), a2.clone(), b2.clone()]).unwrap();
        let q2 = quotient_by_relations(
            parent.clone(),
            &[ab.sub(&ba).unwrap(), a2.clone(), b2.clone()],
        )
        .unwrap();
        let q3 = quotient_by_relations(parent, &[ab.add(&ba).unwrap(), a2.clone(), b2.clone(), ba])
            .unwrap();
        let q4 = quotient_by_relations(a.algebra().clone(), &[ab, a2, b3]).unwrap();

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let loop_parent = PathAlgebra::new(field, loop_quiver).unwrap();
        let x = loop_parent.generators()[1].clone();
        let square = quotient_by_relations(loop_parent.clone(), &[x.pow(2).unwrap()]).unwrap();
        let cubic = quotient_by_relations(loop_parent.clone(), &[x.pow(3).unwrap()]).unwrap();
        let linear = quotient_by_relations(loop_parent, &[x]).unwrap();
        let gentle = quotient_is_gentle(&square).unwrap();
        let square_algebra = QuiverAlgebra::Quotient(square.clone());
        let q1_quiver = q1.original_path_algebra().quiver();
        let q1_minimal_walks = ideal_minimal_generating_set(q1.ideal())
            .unwrap()
            .into_iter()
            .map(|relation| {
                relation.terms()[0]
                    .path()
                    .walk()
                    .unwrap()
                    .iter()
                    .map(|index| q1_quiver.arrow(*index).unwrap().name().to_owned())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        json!({
            "doc_q1": special_predicates(&q1),
            "doc_q1_minimal_walks": q1_minimal_walks,
            "doc_q2": special_predicates(&q2),
            "doc_q3": special_predicates(&q3),
            "doc_q4": special_predicates(&q4),
            "square_loop": special_predicates(&square),
            "cubic_loop": special_predicates(&cubic),
            "linear_loop": special_predicates(&linear),
            "gentle_implication": [
                false,
                gentle,
                gentle_algebra_is_gorenstein(&square_algebra).unwrap()
            ]
        })
    }

    #[test]
    fn pinned_qpa_special_biserial_string_and_gentle_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_special_biserial.schema-v1.json"
        ))
        .unwrap();
        let one_loop =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let three_out = Quiver::from_indexed(
            4,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(1), vertex(3), "b"),
                ArrowSpec::named(vertex(1), vertex(4), "c"),
            ],
        )
        .unwrap();
        let chain = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(2), vertex(3), "b"),
            ],
        )
        .unwrap();
        let chain_path = PathAlgebra::new(Q, chain).unwrap();
        let loop_path = PathAlgebra::new(Q, one_loop.clone()).unwrap();
        let chain_algebra = QuiverAlgebra::Path(chain_path.clone());

        assert_eq!(
            json!({
                "one_loop": is_special_biserial_quiver(&one_loop),
                "three_out": is_special_biserial_quiver(&three_out)
            }),
            expected["quivers"]
        );
        assert_eq!(
            json!({
                "chain": [
                    path_algebra_is_special_biserial(&chain_path),
                    path_algebra_is_string(&chain_path),
                    path_algebra_is_gentle(&chain_path)
                ],
                "one_loop": [
                    path_algebra_is_special_biserial(&loop_path),
                    path_algebra_is_string(&loop_path),
                    path_algebra_is_gentle(&loop_path)
                ]
            }),
            expected["paths"]
        );
        assert_eq!(special_field_fixture(Q), expected["rational"]);
        assert_eq!(special_field_fixture(Z2), expected["f2"]);
        assert!(gentle_algebra_is_special_biserial(&chain_algebra).unwrap());
        assert!(gentle_algebra_is_gorenstein(&chain_algebra).unwrap());
        assert_eq!(
            gentle_algebra_is_gorenstein(&QuiverAlgebra::Path(loop_path)),
            Err(PathalgGbnpError::NonGentleAlgebraProperty)
        );
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn pinned_qpa_admissible_property_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_admissible_properties.schema-v1.json"
        ))
        .unwrap();
        let chain_quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(2), vertex(3), "b"),
            ],
        )
        .unwrap();
        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();

        let chain_parent = PathAlgebra::new(Q, chain_quiver.clone()).unwrap();
        let chain_generators = chain_parent.generators();
        let chain_relation = chain_generators[3].multiply(&chain_generators[4]).unwrap();
        let chain = quotient_by_relations(chain_parent, &[chain_relation]).unwrap();
        let loop_algebra = truncated_path_algebra(Q, loop_quiver.clone(), 3).unwrap();
        let nonadmissible_parent = PathAlgebra::new(Q, loop_quiver.clone()).unwrap();
        let linear_relation = nonadmissible_parent.generators()[1].clone();
        let nonadmissible =
            quotient_by_relations(nonadmissible_parent, &[linear_relation]).unwrap();
        let nonadmissible_cached_before = nonadmissible.cached_admissibility().is_some();
        let nonadmissible_value = quotient_is_admissible(&nonadmissible).unwrap();
        let nonadmissible_cached_after = nonadmissible.cached_admissibility().is_some();

        assert_eq!(
            json!({
                "chain": {
                    "admissible": quotient_is_admissible(&chain).unwrap(),
                    "cached": chain.cached_admissibility().is_some(),
                    "semisimple": admissible_quotient_is_semisimple(&chain).unwrap(),
                    "elementary": admissible_quotient_is_elementary(&chain).unwrap(),
                    "basic": admissible_quotient_is_basic(&chain).unwrap()
                },
                "truncated_loop": {
                    "admissible": quotient_is_admissible(&loop_algebra).unwrap(),
                    "semisimple": admissible_quotient_is_semisimple(&loop_algebra).unwrap(),
                    "elementary": admissible_quotient_is_elementary(&loop_algebra).unwrap(),
                    "basic": admissible_quotient_is_basic(&loop_algebra).unwrap()
                },
                "nonadmissible": {
                    "value": nonadmissible_value,
                    "cached_before": nonadmissible_cached_before,
                    "cached_after": nonadmissible_cached_after
                }
            }),
            expected["rational"]
        );
        assert_eq!(
            admissible_quotient_is_semisimple(&nonadmissible),
            Err(PathalgGbnpError::NonAdmissibleQuotientProperty)
        );
        assert_eq!(
            associated_graded_quotient(&nonadmissible),
            Err(PathalgGbnpError::InvLenLex(
                crate::InvLenLexError::NonAdmissibleDefaultQuotient
            ))
        );

        let chain_parent2 = PathAlgebra::new(Z2, chain_quiver).unwrap();
        let chain_generators2 = chain_parent2.generators();
        let chain_relation2 = chain_generators2[3]
            .multiply(&chain_generators2[4])
            .unwrap();
        let chain2 = quotient_by_relations(chain_parent2, &[chain_relation2]).unwrap();
        let loop_algebra2 = truncated_path_algebra(Z2, loop_quiver, 3).unwrap();
        assert_eq!(
            json!({
                "chain": {
                    "admissible": quotient_is_admissible(&chain2).unwrap(),
                    "cached": chain2.cached_admissibility().is_some(),
                    "semisimple": admissible_quotient_is_semisimple(&chain2).unwrap(),
                    "elementary": admissible_quotient_is_elementary(&chain2).unwrap(),
                    "basic": admissible_quotient_is_basic(&chain2).unwrap()
                },
                "truncated_loop": {
                    "admissible": quotient_is_admissible(&loop_algebra2).unwrap(),
                    "semisimple": admissible_quotient_is_semisimple(&loop_algebra2).unwrap(),
                    "elementary": admissible_quotient_is_elementary(&loop_algebra2).unwrap(),
                    "basic": admissible_quotient_is_basic(&loop_algebra2).unwrap()
                }
            }),
            expected["f2"]
        );

        let discrete_parent =
            PathAlgebra::new(Q, Quiver::from_indexed(2, Vec::new()).unwrap()).unwrap();
        let discrete_ideal = PathAlgebraIdeal::new(discrete_parent, Vec::new()).unwrap();
        let discrete = QuotientPathAlgebra::by_ideal(&discrete_ideal).unwrap();
        assert!(admissible_quotient_is_semisimple(&discrete).unwrap());
    }

    #[test]
    #[allow(clippy::many_single_char_names, clippy::too_many_lines)]
    fn pinned_qpa_associated_algebra_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_associated.schema-v1.json"
        ))
        .unwrap();
        let finite_quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let loop_quiver = Quiver::from_indexed(
            1,
            vec![
                ArrowSpec::named(vertex(1), vertex(1), "x"),
                ArrowSpec::named(vertex(1), vertex(1), "y"),
            ],
        )
        .unwrap();

        let finite = PathAlgebra::new(Q, finite_quiver.clone()).unwrap();
        let finite_monomial = associated_monomial_path_algebra(&finite);
        let finite_graded = associated_graded_path_algebra(&finite).unwrap();
        let parent = PathAlgebra::new(Q, loop_quiver.clone()).unwrap();
        let generators = parent.generators();
        let mut relations = vec![
            generators[1]
                .pow(2)
                .unwrap()
                .add(&generators[2].pow(3).unwrap())
                .unwrap(),
        ];
        parent
            .add_nth_power_to_relations(&mut relations, 4)
            .unwrap();
        let algebra = quotient_by_relations(parent, &relations).unwrap();
        let monomial = associated_monomial_quotient(&algebra).unwrap();
        let graded = associated_graded_quotient(&algebra).unwrap();
        let truncated = truncated_path_algebra(Q, loop_quiver.clone(), 3).unwrap();
        let truncated_monomial = associated_monomial_quotient(&truncated).unwrap();
        let monomial_quotient = monomial.as_quotient().unwrap();

        assert_eq!(
            json!({
                "finite_path_is_monomial": path_algebra_is_monomial(&finite),
                "finite_path_identity": matches!(finite_monomial, QuiverAlgebra::Path(ref result) if result.is_identical(&finite)),
                "finite_path_graded_identity": matches!(finite_graded, QuiverAlgebra::Path(ref result) if result.is_identical(&finite)),
                "nonmonomial": quotient_is_monomial(&algebra).unwrap(),
                "monomial_identity": monomial_quotient.is_identical(&algebra),
                "monomial_result_is_monomial": quotient_is_monomial(monomial_quotient).unwrap(),
                "monomial_dimension": finite_dimension(&monomial),
                "monomial_relations": monomial.relations().iter().map(rational_element_json).collect::<Vec<_>>(),
                "graded_identity": graded.as_quotient().unwrap().is_identical(&algebra),
                "graded_dimension": finite_dimension(&graded),
                "graded_relations": graded.relations().iter().map(rational_element_json).collect::<Vec<_>>(),
                "already_monomial_identity": truncated_monomial.as_quotient().unwrap().is_identical(&truncated)
            }),
            expected["rational"]
        );

        let finite2 = PathAlgebra::new(Z2, finite_quiver).unwrap();
        let finite_monomial2 = associated_monomial_path_algebra(&finite2);
        let finite_graded2 = associated_graded_path_algebra(&finite2).unwrap();
        let parent2 = PathAlgebra::new(Z2, loop_quiver.clone()).unwrap();
        let generators2 = parent2.generators();
        let mut relations2 = vec![
            generators2[1]
                .pow(2)
                .unwrap()
                .add(&generators2[2].pow(3).unwrap())
                .unwrap(),
        ];
        parent2
            .add_nth_power_to_relations(&mut relations2, 4)
            .unwrap();
        let algebra2 = quotient_by_relations(parent2, &relations2).unwrap();
        let monomial2 = associated_monomial_quotient(&algebra2).unwrap();
        let graded2 = associated_graded_quotient(&algebra2).unwrap();
        let truncated2 = truncated_path_algebra(Z2, loop_quiver.clone(), 3).unwrap();
        let truncated_monomial2 = associated_monomial_quotient(&truncated2).unwrap();
        let monomial_quotient2 = monomial2.as_quotient().unwrap();

        assert_eq!(
            json!({
                "finite_path_is_monomial": path_algebra_is_monomial(&finite2),
                "finite_path_identity": matches!(finite_monomial2, QuiverAlgebra::Path(ref result) if result.is_identical(&finite2)),
                "finite_path_graded_identity": matches!(finite_graded2, QuiverAlgebra::Path(ref result) if result.is_identical(&finite2)),
                "nonmonomial": quotient_is_monomial(&algebra2).unwrap(),
                "monomial_identity": monomial_quotient2.is_identical(&algebra2),
                "monomial_result_is_monomial": quotient_is_monomial(monomial_quotient2).unwrap(),
                "monomial_dimension": finite_dimension(&monomial2),
                "monomial_relations": monomial2.relations().iter().map(z2_element_json).collect::<Vec<_>>(),
                "graded_identity": graded2.as_quotient().unwrap().is_identical(&algebra2),
                "graded_dimension": finite_dimension(&graded2),
                "graded_relations": graded2.relations().iter().map(z2_element_json).collect::<Vec<_>>(),
                "already_monomial_identity": truncated_monomial2.as_quotient().unwrap().is_identical(&truncated2)
            }),
            expected["f2"]
        );

        let cyclic_parent = PathAlgebra::new(Z2, loop_quiver.clone()).unwrap();
        assert_eq!(
            associated_graded_path_algebra(&cyclic_parent),
            Err(PathalgGbnpError::InfiniteDimensionalAssociatedAlgebra)
        );
        let linear_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let linear_parent = PathAlgebra::new(Z2, linear_quiver).unwrap();
        let x = linear_parent.generators()[1].clone();
        let nonadmissible = quotient_by_relations(linear_parent, &[x]).unwrap();
        assert_eq!(
            associated_graded_quotient(&nonadmissible),
            Err(PathalgGbnpError::NonAdmissibleAssociatedAlgebra)
        );
    }

    #[test]
    fn monomial_predicate_preserves_qpa_lazy_basis_branches() {
        let quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let parent = PathAlgebra::new(Q, quiver).unwrap();
        let x = parent.generators()[1].clone();

        let monomial_ideal =
            PathAlgebraIdeal::new(parent.clone(), vec![x.pow(2).unwrap()]).unwrap();
        let monomial = QuotientPathAlgebra::by_ideal(&monomial_ideal).unwrap();
        assert!(quotient_is_monomial(&monomial).unwrap());
        assert!(monomial_ideal.groebner_basis().is_none());

        let nonmonomial_ideal = PathAlgebraIdeal::new(
            parent,
            vec![x.pow(2).unwrap().add(&x.pow(3).unwrap()).unwrap()],
        )
        .unwrap();
        let nonmonomial = QuotientPathAlgebra::by_ideal(&nonmonomial_ideal).unwrap();
        assert!(!quotient_is_monomial(&nonmonomial).unwrap());
        assert!(nonmonomial_ideal.groebner_basis().is_some());
    }

    #[test]
    fn pinned_qpa_triangular_reduced_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_triangular_reduced.schema-v1.json"
        ))
        .unwrap();
        let one_vertex_quiver = Quiver::from_indexed(1, vec![]).unwrap();
        let arrow_quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "u")]).unwrap();
        let cycle_quiver = Quiver::from_indexed(
            2,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(2), vertex(1), "b"),
            ],
        )
        .unwrap();

        let one_vertex = PathAlgebra::new(Q, one_vertex_quiver.clone()).unwrap();
        let arrow = PathAlgebra::new(Q, arrow_quiver.clone()).unwrap();
        let truncated = truncated_path_algebra(Q, cycle_quiver.clone(), 2).unwrap();
        let cycle = PathAlgebra::new(Q, cycle_quiver.clone()).unwrap();
        let killed_reverse =
            quotient_by_relations(cycle.clone(), &[cycle.generators()[3].clone()]).unwrap();
        assert_eq!(
            json!({
                "one_vertex_path": one_vertex.is_triangular_reduced().unwrap(),
                "arrow_path": arrow.is_triangular_reduced().unwrap(),
                "two_cycle_truncated": quotient_is_triangular_reduced(&truncated).unwrap(),
                "killed_reverse_arrow": quotient_is_triangular_reduced(&killed_reverse).unwrap()
            }),
            expected["rational"]
        );

        let one_vertex2 = PathAlgebra::new(Z2, one_vertex_quiver).unwrap();
        let arrow2 = PathAlgebra::new(Z2, arrow_quiver).unwrap();
        let truncated2 = truncated_path_algebra(Z2, cycle_quiver.clone(), 2).unwrap();
        let cycle2 = PathAlgebra::new(Z2, cycle_quiver).unwrap();
        let killed_reverse2 =
            quotient_by_relations(cycle2.clone(), &[cycle2.generators()[3].clone()]).unwrap();
        assert_eq!(
            json!({
                "one_vertex_path": one_vertex2.is_triangular_reduced().unwrap(),
                "arrow_path": arrow2.is_triangular_reduced().unwrap(),
                "two_cycle_truncated": quotient_is_triangular_reduced(&truncated2).unwrap(),
                "killed_reverse_arrow": quotient_is_triangular_reduced(&killed_reverse2).unwrap()
            }),
            expected["f2"]
        );
    }

    #[test]
    #[allow(clippy::many_single_char_names, clippy::too_many_lines)]
    fn pinned_qpa_inverse_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_inverse.schema-v1.json"
        ))
        .unwrap();
        let quiver = Quiver::from_indexed(
            2,
            vec![
                ArrowSpec::named(vertex(1), vertex(1), "x"),
                ArrowSpec::named(vertex(2), vertex(2), "y"),
                ArrowSpec::named(vertex(1), vertex(2), "a"),
            ],
        )
        .unwrap();
        let algebra = truncated_path_algebra(Q, quiver.clone(), 2).unwrap();
        let generators = algebra.generators().unwrap();
        let q2 = Q.nth(Integer::from(2));
        let q3 = Q.nth(Integer::from(3));
        let unit = generators[0]
            .scale_left(&q2)
            .unwrap()
            .add(&generators[1].scale_left(&q3).unwrap())
            .unwrap()
            .add(&generators[2])
            .unwrap()
            .add(&generators[3])
            .unwrap()
            .add(&generators[4])
            .unwrap();
        let nonunit = generators[0].add(&generators[4]).unwrap();
        let inverse = inverse_quotient_element(&unit).unwrap().unwrap();
        let negative_two = power_quotient_element(&unit, -2).unwrap().unwrap();
        let zero_power = power_quotient_element(&nonunit, 0).unwrap().unwrap();

        assert_eq!(
            json!({
                "inverse": rational_element_json(inverse.representative()),
                "product": rational_element_json(unit.multiply(&inverse).unwrap().representative()),
                "negative_two": rational_element_json(negative_two.representative()),
                "nonunit_inverse_fail": inverse_quotient_element(&nonunit).unwrap().is_none(),
                "nonunit_negative_fail": power_quotient_element(&nonunit, -1).unwrap().is_none(),
                "zero_power": rational_element_json(zero_power.representative())
            }),
            expected["rational"]
        );

        let algebra2 = truncated_path_algebra(Z2, quiver, 2).unwrap();
        let generators2 = algebra2.generators().unwrap();
        let unit2 = generators2[0]
            .add(&generators2[1])
            .unwrap()
            .add(&generators2[2])
            .unwrap()
            .add(&generators2[3])
            .unwrap()
            .add(&generators2[4])
            .unwrap();
        let nonunit2 = generators2[2].clone();
        let inverse2 = inverse_quotient_element(&unit2).unwrap().unwrap();
        let negative_two2 = power_quotient_element(&unit2, -2).unwrap().unwrap();
        let zero_power2 = power_quotient_element(&nonunit2, 0).unwrap().unwrap();
        assert_eq!(
            json!({
                "inverse": z2_element_json(inverse2.representative()),
                "product": z2_element_json(unit2.multiply(&inverse2).unwrap().representative()),
                "negative_two": z2_element_json(negative_two2.representative()),
                "nonunit_inverse_fail": inverse_quotient_element(&nonunit2).unwrap().is_none(),
                "zero_power": z2_element_json(zero_power2.representative())
            }),
            expected["f2"]
        );
    }

    #[test]
    #[allow(clippy::many_single_char_names, clippy::too_many_lines)]
    fn pinned_qpa_vertex_quotient_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_vertex_quotient.schema-v1.json"
        ))
        .unwrap();
        let quiver = Quiver::from_indexed(
            4,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(2), vertex(3), "b"),
                ArrowSpec::named(vertex(1), vertex(3), "c"),
                ArrowSpec::named(vertex(3), vertex(4), "d"),
                ArrowSpec::named(vertex(2), vertex(4), "e"),
                ArrowSpec::named(vertex(1), vertex(4), "f"),
            ],
        )
        .unwrap();
        let algebra = PathAlgebra::new(Q, quiver.clone()).unwrap();
        let generators = algebra.generators();
        let a = &generators[4];
        let b = &generators[5];
        let c = &generators[6];
        let d = &generators[7];
        let e = &generators[8];
        let f = &generators[9];
        let relations = [
            a.multiply(b).unwrap().sub(c).unwrap(),
            b.multiply(d).unwrap().sub(e).unwrap(),
            c.multiply(d).unwrap().sub(f).unwrap(),
        ];
        let quotient = quotient_by_relations(algebra.clone(), &relations).unwrap();
        let remove_two = quiver_algebra_of_quotient_mod_aea(&quotient, &[2]).unwrap();
        let remove_three = quiver_algebra_of_quotient_mod_aea(&quotient, &[3]).unwrap();
        let none = quiver_algebra_of_quotient_mod_aea(&quotient, &[]).unwrap();
        let ignored = quiver_algebra_of_quotient_mod_aea(&quotient, &[0, 5, 99]).unwrap();
        let path_input = quiver_algebra_of_path_algebra_mod_aea(&algebra, &[2]).unwrap();

        assert_eq!(
            json!({
                "remove_two": {
                    "quiver": quiver_json(remove_two.original_path_algebra().quiver()),
                    "relations": remove_two.relations().iter().map(rational_element_json).collect::<Vec<_>>(),
                    "dimension": finite_dimension(&remove_two)
                },
                "remove_three": {
                    "quiver": quiver_json(remove_three.original_path_algebra().quiver()),
                    "relations": remove_three.relations().iter().map(rational_element_json).collect::<Vec<_>>(),
                    "dimension": finite_dimension(&remove_three)
                },
                "none": {
                    "quiver": quiver_json(none.original_path_algebra().quiver()),
                    "relations": none.relations().iter().map(rational_element_json).collect::<Vec<_>>(),
                    "dimension": finite_dimension(&none)
                },
                "ignored_same_quiver": quiver_json(ignored.original_path_algebra().quiver())
                    == quiver_json(none.original_path_algebra().quiver()),
                "path_algebra_input": {
                    "is_quotient": !path_input.is_path_algebra(),
                    "relations": path_input.relations().iter().map(rational_element_json).collect::<Vec<_>>()
                }
            }),
            expected["rational"]
        );

        let algebra2 = PathAlgebra::new(Z2, quiver).unwrap();
        let generators2 = algebra2.generators();
        let relations2 = [
            generators2[4]
                .multiply(&generators2[5])
                .unwrap()
                .add(&generators2[6])
                .unwrap(),
            generators2[5]
                .multiply(&generators2[7])
                .unwrap()
                .add(&generators2[8])
                .unwrap(),
            generators2[6]
                .multiply(&generators2[7])
                .unwrap()
                .add(&generators2[9])
                .unwrap(),
        ];
        let quotient2 = quotient_by_relations(algebra2, &relations2).unwrap();
        let remove_two2 = quiver_algebra_of_quotient_mod_aea(&quotient2, &[2]).unwrap();
        assert_eq!(
            json!({
                "remove_two": {
                    "quiver": quiver_json(remove_two2.original_path_algebra().quiver()),
                    "relations": remove_two2.relations().iter().map(z2_element_json).collect::<Vec<_>>(),
                    "dimension": finite_dimension(&remove_two2)
                }
            }),
            expected["f2"]
        );
    }

    #[test]
    fn qpa_default_truncation_and_opposite_use_gbnp() {
        let quiver = quiver();
        let truncated = truncated_path_algebra(Q, quiver.clone(), 3).unwrap();
        assert_eq!(
            truncated.dimension().unwrap(),
            AlgebraDimension::Finite(Integer::from(17))
        );
        let opposite = opposite_path_algebra(&truncated).unwrap();
        assert_eq!(
            opposite.dimension().unwrap(),
            AlgebraDimension::Finite(Integer::from(17))
        );
        assert!(
            opposite_path_algebra(&opposite)
                .unwrap()
                .is_identical(&truncated)
        );

        let a = truncated
            .path(quiver.arrow(arrow(1)).unwrap().path())
            .unwrap();
        let b = truncated
            .path(quiver.arrow(arrow(2)).unwrap().path())
            .unwrap();
        let c = truncated
            .path(quiver.arrow(arrow(3)).unwrap().path())
            .unwrap();
        let element = a
            .add(&a.pow(2).unwrap())
            .unwrap()
            .add(&b.multiply(&c).unwrap())
            .unwrap();
        assert_eq!(
            opposite_quotient_element(&element)
                .unwrap()
                .representative()
                .terms()
                .iter()
                .map(|term| term.path().external_representation())
                .collect::<Vec<_>>(),
            [vec![4], vec![4, 4], vec![6, 5]]
        );

        let binary = truncated_path_algebra(Z2, quiver, 3).unwrap();
        assert_eq!(
            binary.dimension().unwrap(),
            AlgebraDimension::Finite(Integer::from(17))
        );
    }

    #[test]
    fn lazy_dimension_attaches_default_gbnp_basis_without_changing_the_family() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_conversion.schema-v1.json"
        ))
        .unwrap();
        let quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let algebra = PathAlgebra::new(Q, quiver.clone()).unwrap();
        let x = algebra
            .element(quiver.arrow(arrow(1)).unwrap().path())
            .unwrap();
        let square = x.pow(2).unwrap();
        let cached_ideal = PathAlgebraIdeal::new(algebra.clone(), vec![square.clone()]).unwrap();
        let cached_quotient = QuotientPathAlgebra::by_ideal(&cached_ideal).unwrap();
        let cached_basis_before = cached_quotient
            .canonical_basis_object()
            .unwrap()
            .has_basis_vectors();
        quotient_dimension(&cached_quotient).unwrap();
        let cached_basis_after = cached_quotient
            .canonical_basis_object()
            .unwrap()
            .has_basis_vectors();

        let ideal = PathAlgebraIdeal::new(algebra, vec![square.clone()]).unwrap();
        let quotient = QuotientPathAlgebra::by_ideal(&ideal).unwrap();
        let before = quotient.element(x.add(&square).unwrap()).unwrap();
        let before_basis = ideal.groebner_basis().is_some();
        let before_normal_function = quotient.has_normal_form_function();
        let dimension = quotient_dimension(&quotient).unwrap();
        assert!(quotient_is_finite_dimensional(&quotient).unwrap());
        let after = quotient.element(x.add(&square).unwrap()).unwrap();
        let basis = quotient.canonical_basis().unwrap().unwrap();

        assert_eq!(
            json!({
                "before_basis": before_basis,
                "before_normal_function": before_normal_function,
                "before_element": rational_element_json(before.representative()),
                "before_element_is_normal": before.is_normal_form(),
                "cached_basis_before": cached_basis_before,
                "cached_basis_after": cached_basis_after,
                "dimension": match dimension {
                    AlgebraDimension::Finite(value) => small_integer(&value),
                    AlgebraDimension::Infinite => panic!("oracle quotient must be finite"),
                },
                "after_basis": ideal.groebner_basis().is_some(),
                "tip_reduced": ideal
                    .groebner_basis()
                    .unwrap()
                    .tip_reduced_property()
                    == Some(true),
                "after_normal_function": quotient.has_normal_form_function(),
                "after_element": rational_element_json(after.representative()),
                "after_element_is_normal": after.is_normal_form(),
                "basis": basis
                    .iter()
                    .map(|element| rational_element_json(element.representative()))
                    .collect::<Vec<_>>(),
                "basis_normal_flags": basis
                    .iter()
                    .map(qpa_core::QuotientPathAlgebraElement::is_normal_form)
                    .collect::<Vec<_>>()
            }),
            expected["lazy_default"]
        );

        let algebra2 = PathAlgebra::new(Z2, quiver).unwrap();
        let x2 = algebra2
            .element(algebra2.quiver().arrow(arrow(1)).unwrap().path())
            .unwrap();
        let ideal2 = PathAlgebraIdeal::new(algebra2, vec![x2.pow(2).unwrap()]).unwrap();
        let quotient2 = QuotientPathAlgebra::by_ideal(&ideal2).unwrap();
        assert_eq!(
            quotient_dimension(&quotient2).unwrap(),
            AlgebraDimension::Finite(Integer::from(2))
        );
        assert!(!quotient2.has_normal_form_function());
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn pinned_qpa_pathalg_centre_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_centre.schema-v1.json"
        ))
        .unwrap();
        let quiver = Quiver::from_indexed(
            2,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(2), vertex(1), "b"),
            ],
        )
        .unwrap();
        let parent = PathAlgebra::new(Q, quiver.clone()).unwrap();
        let parent_centre = parent.centre().unwrap();
        assert_eq!(parent_centre.dimension(), &AlgebraDimension::Infinite);
        assert!(
            parent_centre
                .contains(&parent_centre.generators()[1])
                .unwrap()
        );
        assert!(!parent_centre.contains(&parent.generators()[2]).unwrap());
        assert_eq!(
            json!({
                "generators": parent_centre
                    .generators()
                    .iter()
                    .map(rational_element_json)
                    .collect::<Vec<_>>()
            }),
            expected["parent_cycle"]
        );

        let acyclic_quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "u"),
                ArrowSpec::named(vertex(2), vertex(3), "v"),
            ],
        )
        .unwrap();
        let acyclic = PathAlgebra::new(Q, acyclic_quiver).unwrap();
        let acyclic_centre = acyclic.centre().unwrap();
        assert_eq!(
            acyclic_centre.dimension(),
            &AlgebraDimension::Finite(Integer::from(1))
        );
        assert_eq!(
            json!({
                "generators": acyclic_centre
                    .generators()
                    .iter()
                    .map(rational_element_json)
                    .collect::<Vec<_>>()
            }),
            expected["parent_acyclic"]
        );

        let discrete = PathAlgebra::new(Q, Quiver::from_indexed(2, Vec::new()).unwrap()).unwrap();
        let discrete_centre = discrete.centre().unwrap();
        assert_eq!(
            discrete_centre.dimension(),
            &AlgebraDimension::Finite(Integer::from(2))
        );
        assert_eq!(
            json!({
                "generators": discrete_centre
                    .generators()
                    .iter()
                    .map(rational_element_json)
                    .collect::<Vec<_>>()
            }),
            expected["parent_discrete"]
        );

        let quotient = truncated_path_algebra(Q, quiver.clone(), 4).unwrap();
        let centre = quotient_centre(&quotient).unwrap();
        let generators = quotient.generators().unwrap();
        let central = generators[2]
            .multiply(&generators[3])
            .unwrap()
            .add(&generators[3].multiply(&generators[2]).unwrap())
            .unwrap();
        let AlgebraDimension::Finite(ambient_dimension) = quotient.dimension().unwrap() else {
            panic!("oracle quotient must be finite");
        };
        assert_eq!(
            json!({
                "ambient_dimension": small_integer(&ambient_dimension),
                "basis": centre
                    .basis()
                    .iter()
                    .map(|element| rational_element_json(element.representative()))
                    .collect::<Vec<_>>(),
                "dimension": centre.dimension(),
                "central_member": centre.contains(&central).unwrap(),
                "arrow_member": centre.contains(&generators[2]).unwrap()
            }),
            expected["quotient"]
        );

        let quotient2 = truncated_path_algebra(Z2, quiver, 4).unwrap();
        let centre2 = quotient_centre(&quotient2).unwrap();
        assert_eq!(
            json!({
                "basis": centre2
                    .basis()
                    .iter()
                    .map(|element| z2_element_json(element.representative()))
                    .collect::<Vec<_>>(),
                "dimension": centre2.dimension()
            }),
            expected["f2"]
        );
    }
}
