//! QPA path-algebra constructors whose configured default backend is GBNP.
//!
//! `qpa-core` owns quotient representation and backend-neutral factor logic;
//! this module follows QPA's default call chain through `GBNPGroebnerBasis`.

use numerica::domains::Field;
use qpa_core::{
    PathAlgebra, PathAlgebraElement, PathAlgebraError, Quiver, QuotientError, QuotientPathAlgebra,
    QuotientPathAlgebraElement,
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
    use numerica::domains::{finite_field::Z2, integer::Integer, rational::Q};
    use qpa_core::{AlgebraDimension, ArrowIndex, ArrowSpec, Quiver, VertexIndex};

    use super::{opposite_path_algebra, opposite_quotient_element, truncated_path_algebra};

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
}
