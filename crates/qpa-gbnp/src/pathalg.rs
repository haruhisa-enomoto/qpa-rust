//! QPA path-algebra constructors whose configured default backend is GBNP.
//!
//! `qpa-core` owns quotient representation and backend-neutral factor logic;
//! this module follows QPA's default call chain through `GBNPGroebnerBasis`.

use numerica::domains::Field;
use qpa_core::{
    AlgebraDimension, GroebnerBasis, PathAlgebra, PathAlgebraElement, PathAlgebraError, Quiver,
    QuotientError, QuotientPathAlgebra, QuotientPathAlgebraElement,
};
use thiserror::Error;

use crate::{QpaGbnpError, gbnp_groebner_basis};

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
    Ok(quotient
        .install_opposite_from_complete_relations(&opposite_relators, &complete_relations)?)
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
    /// Parent or opposite path-algebra construction failed.
    #[error(transparent)]
    PathAlgebra(#[from] PathAlgebraError),
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
        finite_field::Z2,
        integer::Integer,
        rational::{Q, RationalField},
    };
    use qpa_core::{
        AlgebraDimension, ArrowIndex, ArrowSpec, PathAlgebra, PathAlgebraElement, PathAlgebraIdeal,
        Quiver, QuotientPathAlgebra, VertexIndex,
    };
    use serde_json::{Value, json};

    use super::{
        opposite_path_algebra, opposite_quotient_element, quotient_dimension,
        quotient_is_finite_dimensional, truncated_path_algebra,
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
}
