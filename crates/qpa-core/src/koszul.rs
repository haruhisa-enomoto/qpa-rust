//! QPA's quadratic-ideal and Koszul-dual helpers.
//!
//! This module translates `lib/pathalgkoszul.gi` from QPA 1.37.  In
//! particular, the quadratic perpendicular is obtained from GAP's
//! `NullspaceMat(TransposedMat(A))` construction, using the same ordered basis
//! of length-two paths and the same opposite-path name map.

use numerica::domains::Field;
use thiserror::Error;

use crate::{
    DenseMatrix, MatrixError, Path, PathAlgebra, PathAlgebraElement, PathAlgebraError, PathError,
    Quiver, QuiverError,
};

/// QPA's returned opposite path algebra and basis of the quadratic
/// perpendicular.
#[derive(Clone, Debug)]
pub struct QuadraticPerpendicular<F: Field> {
    opposite_algebra: PathAlgebra<F>,
    basis: Vec<PathAlgebraElement<F>>,
}

impl<F: Field> QuadraticPerpendicular<F> {
    /// Return the newly constructed path algebra on QPA's opposite quiver.
    #[must_use]
    pub const fn opposite_algebra(&self) -> &PathAlgebra<F> {
        &self.opposite_algebra
    }

    /// Return QPA's `basisOfIperp` list.
    #[must_use]
    pub fn basis(&self) -> &[PathAlgebraElement<F>] {
        &self.basis
    }

    /// Consume the result into QPA's two returned components.
    #[must_use]
    pub fn into_parts(self) -> (PathAlgebra<F>, Vec<PathAlgebraElement<F>>) {
        (self.opposite_algebra, self.basis)
    }
}

/// Translate QPA's `PathsOfLengthTwo` nested arrow loops and final `Sort`.
///
/// # Errors
///
/// Propagates a path-family or quiver-family invariant failure.
pub fn paths_of_length_two(quiver: &Quiver) -> Result<Vec<Path>, KoszulError> {
    let mut paths = Vec::new();
    for source in quiver.vertices() {
        for first in quiver.outgoing_arrows(&source)? {
            for second in quiver.outgoing_arrows(&first.target())? {
                paths.push(first.path().multiply(&second.path())?);
            }
        }
    }

    let mut comparison_error = None;
    paths.sort_by(|left, right| match quiver.compare_paths(left, right) {
        Ok(ordering) => ordering,
        Err(error) => {
            comparison_error = Some(error);
            std::cmp::Ordering::Equal
        }
    });
    comparison_error.map_or(Ok(paths), |error| Err(error.into()))
}

/// Translate QPA's `IsQuadraticIdeal` method for path-algebra elements.
///
/// The empty list and a list containing only zero elements are quadratic, as
/// in QPA. Rust's element type excludes QPA's separate quotient-element `fail`
/// branch at compile time.
///
/// # Errors
///
/// QPA dispatches this method only for a homogeneous list. Rust reports a list
/// containing elements from different path-algebra families explicitly.
pub fn is_quadratic_ideal<F: Field>(
    generators: &[PathAlgebraElement<F>],
) -> Result<bool, KoszulError> {
    ensure_homogeneous(generators)?;
    Ok(generators.iter().all(|generator| {
        generator
            .terms()
            .iter()
            .all(|term| term.path().length() == Some(2))
    }))
}

/// Translate QPA's `QuadraticPerpOfPathAlgebraIdeal`.
///
/// The input rows are coefficients in QPA's sorted length-two path basis. The
/// returned basis is the row nullspace of the transposed coefficient matrix,
/// translated term-by-term to the opposite quiver.
///
/// # Errors
///
/// Preserves QPA's empty-list defect, rejects a nonquadratic or heterogeneous
/// generator list, and propagates exact matrix/path construction failures.
pub fn quadratic_perpendicular_of_path_algebra_ideal<F: Field>(
    generators: &[PathAlgebraElement<F>],
) -> Result<QuadraticPerpendicular<F>, KoszulError> {
    if !is_quadratic_ideal(generators)? {
        return Err(KoszulError::NotQuadratic);
    }
    let algebra = generators
        .first()
        .map(PathAlgebraElement::algebra)
        .ok_or(KoszulError::EmptyGeneratorListBeforeDocumentedBranch)?;
    let degree_two_paths = paths_of_length_two(algebra.quiver())?;
    if degree_two_paths.is_empty() {
        return Err(KoszulError::CollapsedEmptyNullspaceMatrix);
    }

    let mut rows = Vec::with_capacity(generators.len());
    for generator in generators {
        let mut row = vec![algebra.field().zero(); degree_two_paths.len()];
        for term in generator.terms() {
            let position = degree_two_paths
                .iter()
                .position(|path| path == term.path())
                .ok_or(KoszulError::TermOutsideDegreeTwoBasis)?;
            row[position] = term.coefficient().clone();
        }
        rows.push(row);
    }

    let coefficients = DenseMatrix::from_rows(rows, algebra.field().clone())?;
    let perpendicular_rows = coefficients.transpose().nullspace_rows();
    let opposite_quiver = algebra.quiver().opposite()?;
    // QPA calls `PathAlgebra(K, Q_op)` here, not its cached
    // `OppositePathAlgebra(A)` operation, so each invocation creates a new
    // path-algebra object on the cached opposite quiver.
    let opposite_algebra = PathAlgebra::new(algebra.field().clone(), opposite_quiver)?;
    let opposite_paths = degree_two_paths
        .iter()
        .map(|path| algebra.quiver().opposite_path(path))
        .collect::<Result<Vec<_>, _>>()?;

    let mut basis = Vec::with_capacity(perpendicular_rows.nrows());
    for row in 0..perpendicular_rows.nrows() {
        basis.push(
            opposite_algebra.from_terms(
                opposite_paths
                    .iter()
                    .cloned()
                    .zip(perpendicular_rows.row(row).iter().cloned()),
            )?,
        );
    }

    Ok(QuadraticPerpendicular {
        opposite_algebra,
        basis,
    })
}

fn ensure_homogeneous<F: Field>(generators: &[PathAlgebraElement<F>]) -> Result<(), KoszulError> {
    let Some(first) = generators.first() else {
        return Ok(());
    };
    if generators
        .iter()
        .all(|generator| generator.algebra().is_identical(first.algebra()))
    {
        Ok(())
    } else {
        Err(KoszulError::HeterogeneousList)
    }
}

/// Failure in QPA's quadratic-ideal and Koszul-dual helpers.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum KoszulError {
    /// QPA's method filter requires one homogeneous element family.
    #[error("the generator list contains different path-algebra families")]
    HeterogeneousList,
    /// QPA returns `fail` after printing this diagnostic.
    #[error("BasisOfIperp(): ideal is not quadratic")]
    NotQuadratic,
    /// QPA indexes `gens[1]` before reaching its documented empty-list branch.
    #[error("QPA accesses the first generator before its empty-list branch")]
    EmptyGeneratorListBeforeDocumentedBranch,
    /// GAP loses the `n x 0` shape and cannot dispatch `NullspaceMat`.
    #[error("QPA cannot dispatch NullspaceMat after its empty matrix transpose")]
    CollapsedEmptyNullspaceMatrix,
    /// A supposedly quadratic monomial was absent from the parent quiver basis.
    #[error("a relation term is absent from the ordered length-two path basis")]
    TermOutsideDegreeTwoBasis,
    /// QPA quiver traversal failed.
    #[error(transparent)]
    Quiver(#[from] QuiverError),
    /// QPA path construction or comparison failed.
    #[error(transparent)]
    Path(#[from] PathError),
    /// QPA path-algebra construction failed.
    #[error(transparent)]
    PathAlgebra(#[from] PathAlgebraError),
    /// GAP-compatible nullspace construction failed.
    #[error(transparent)]
    Matrix(#[from] MatrixError),
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        Field,
        finite_field::{FiniteFieldCore, Z2, Zp},
        integer::Integer,
        rational::{Q, RationalField},
    };
    use serde_json::{Value, json};

    use crate::{ArrowIndex, ArrowSpec, PathAlgebra, PathAlgebraElement, Quiver, VertexIndex};

    use super::{
        KoszulError, is_quadratic_ideal, paths_of_length_two,
        quadratic_perpendicular_of_path_algebra_ideal,
    };

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn arrow(index: usize) -> ArrowIndex {
        ArrowIndex::new(index).unwrap()
    }

    fn loop_quiver() -> Quiver {
        Quiver::from_indexed(
            1,
            vec![
                ArrowSpec::named(vertex(1), vertex(1), "x"),
                ArrowSpec::named(vertex(1), vertex(1), "y"),
            ],
        )
        .unwrap()
    }

    fn small_integer(value: &Integer) -> i64 {
        value
            .to_string()
            .parse()
            .expect("the Koszul fixture uses only small integers")
    }

    fn element_json<F: Field>(
        element: &PathAlgebraElement<F>,
        coefficient_json: impl Fn(&F::Element) -> Value,
    ) -> Value {
        Value::Array(
            element
                .terms()
                .iter()
                .map(|term| {
                    json!({
                        "path": term.path().to_string(),
                        "coefficient": coefficient_json(term.coefficient())
                    })
                })
                .collect(),
        )
    }

    fn rational_json(element: &PathAlgebraElement<RationalField>) -> Value {
        element_json(element, |coefficient| {
            json!([
                small_integer(coefficient.numerator_ref()),
                small_integer(coefficient.denominator_ref())
            ])
        })
    }

    fn zp_json(element: &PathAlgebraElement<Zp>, field: &Zp) -> Value {
        element_json(element, |coefficient| {
            small_integer(&field.to_integer(coefficient)).into()
        })
    }

    fn z2_json(element: &PathAlgebraElement<Z2>) -> Value {
        element_json(element, |coefficient| {
            small_integer(&Z2.to_integer(coefficient)).into()
        })
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn pinned_qpa_pathalg_koszul_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_koszul.schema-v1.json"
        ))
        .unwrap();
        assert_eq!(
            expected["oracle"],
            json!({
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1"
            })
        );

        let quiver = loop_quiver();
        assert_eq!(
            json!(
                paths_of_length_two(&quiver)
                    .unwrap()
                    .iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>()
            ),
            expected["paths_loop"]
        );
        let multi = Quiver::from_indexed(
            2,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(2), vertex(1), "b"),
                ArrowSpec::named(vertex(1), vertex(1), "c"),
            ],
        )
        .unwrap();
        assert_eq!(
            json!(
                paths_of_length_two(&multi)
                    .unwrap()
                    .iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>()
            ),
            expected["paths_multi"]
        );

        let algebra = PathAlgebra::new(Q, quiver.clone()).unwrap();
        let x = algebra
            .element(quiver.arrow(arrow(1)).unwrap().path())
            .unwrap();
        let y = algebra
            .element(quiver.arrow(arrow(2)).unwrap().path())
            .unwrap();
        let rational_generators = vec![
            x.pow(2).unwrap().add(&y.pow(2).unwrap()).unwrap(),
            x.multiply(&y)
                .unwrap()
                .sub(
                    &y.multiply(&x)
                        .unwrap()
                        .scale_left(&Q.to_element(2.into(), 1.into(), true))
                        .unwrap(),
                )
                .unwrap(),
        ];
        assert_eq!(
            json!({
                "empty": is_quadratic_ideal::<RationalField>(&[]).unwrap(),
                "relations": is_quadratic_ideal(&rational_generators).unwrap(),
                "zero": is_quadratic_ideal(&[algebra.zero()]).unwrap(),
                "linear": is_quadratic_ideal(std::slice::from_ref(&x)).unwrap()
            }),
            expected["quadratic"]
        );
        let perpendicular =
            quadratic_perpendicular_of_path_algebra_ideal(&rational_generators).unwrap();
        assert_eq!(
            json!({
                "vertices": perpendicular.opposite_algebra().quiver().vertices().map(|vertex| vertex.name().to_owned()).collect::<Vec<_>>(),
                "arrows": perpendicular.opposite_algebra().quiver().arrows().map(|arrow| arrow.name().to_owned()).collect::<Vec<_>>()
            }),
            expected["opposite_quiver"]
        );
        assert_eq!(
            json!(
                perpendicular
                    .basis()
                    .iter()
                    .map(rational_json)
                    .collect::<Vec<_>>()
            ),
            expected["rational_perpendicular"]
        );
        let zero_perpendicular =
            quadratic_perpendicular_of_path_algebra_ideal(&[algebra.zero()]).unwrap();
        assert_eq!(
            json!(
                zero_perpendicular
                    .basis()
                    .iter()
                    .map(rational_json)
                    .collect::<Vec<_>>()
            ),
            expected["zero_perpendicular"]
        );
        let full_perpendicular = quadratic_perpendicular_of_path_algebra_ideal(&[
            x.pow(2).unwrap(),
            x.multiply(&y).unwrap(),
            y.multiply(&x).unwrap(),
            y.pow(2).unwrap(),
        ])
        .unwrap();
        assert_eq!(
            json!(
                full_perpendicular
                    .basis()
                    .iter()
                    .map(rational_json)
                    .collect::<Vec<_>>()
            ),
            expected["full_perpendicular"]
        );

        let algebra2 = PathAlgebra::new(Z2, quiver.clone()).unwrap();
        let x2 = algebra2
            .element(quiver.arrow(arrow(1)).unwrap().path())
            .unwrap();
        let y2 = algebra2
            .element(quiver.arrow(arrow(2)).unwrap().path())
            .unwrap();
        let perpendicular2 = quadratic_perpendicular_of_path_algebra_ideal(&[
            x2.pow(2).unwrap().add(&y2.pow(2).unwrap()).unwrap(),
            x2.multiply(&y2).unwrap(),
        ])
        .unwrap();
        assert_eq!(
            json!(
                perpendicular2
                    .basis()
                    .iter()
                    .map(z2_json)
                    .collect::<Vec<_>>()
            ),
            expected["f2_perpendicular"]
        );

        let field5 = Zp::new(5);
        let algebra5 = PathAlgebra::new(field5.clone(), quiver.clone()).unwrap();
        let x5 = algebra5
            .element(quiver.arrow(arrow(1)).unwrap().path())
            .unwrap();
        let y5 = algebra5
            .element(quiver.arrow(arrow(2)).unwrap().path())
            .unwrap();
        let relation5 = x5
            .multiply(&y5)
            .unwrap()
            .sub(
                &y5.multiply(&x5)
                    .unwrap()
                    .scale_left(&field5.to_element(2))
                    .unwrap(),
            )
            .unwrap();
        let perpendicular5 = quadratic_perpendicular_of_path_algebra_ideal(&[
            x5.pow(2).unwrap().add(&y5.pow(2).unwrap()).unwrap(),
            relation5,
        ])
        .unwrap();
        assert_eq!(
            json!(
                perpendicular5
                    .basis()
                    .iter()
                    .map(|element| zp_json(element, &field5))
                    .collect::<Vec<_>>()
            ),
            expected["f5_perpendicular"]
        );
    }

    #[test]
    fn source_failures_and_homogeneous_dispatch_are_explicit() {
        assert_eq!(
            quadratic_perpendicular_of_path_algebra_ideal::<RationalField>(&[]).unwrap_err(),
            KoszulError::EmptyGeneratorListBeforeDocumentedBranch
        );

        let quiver = loop_quiver();
        let algebra = PathAlgebra::new(Q, quiver.clone()).unwrap();
        let x = algebra
            .element(quiver.arrow(arrow(1)).unwrap().path())
            .unwrap();
        assert_eq!(
            quadratic_perpendicular_of_path_algebra_ideal(std::slice::from_ref(&x)).unwrap_err(),
            KoszulError::NotQuadratic
        );
        assert_eq!(
            KoszulError::NotQuadratic.to_string(),
            "BasisOfIperp(): ideal is not quadratic"
        );

        let foreign = PathAlgebra::new(Q, quiver).unwrap();
        assert_eq!(
            is_quadratic_ideal(&[algebra.zero(), foreign.zero()]).unwrap_err(),
            KoszulError::HeterogeneousList
        );

        let arrowless_quiver = Quiver::from_indexed(1, vec![]).unwrap();
        let arrowless = PathAlgebra::new(Q, arrowless_quiver).unwrap();
        assert_eq!(
            quadratic_perpendicular_of_path_algebra_ideal(&[arrowless.zero()]).unwrap_err(),
            KoszulError::CollapsedEmptyNullspaceMatrix
        );
    }
}
