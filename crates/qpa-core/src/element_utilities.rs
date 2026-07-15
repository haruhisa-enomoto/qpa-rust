//! QPA utilities acting on individual path- and quotient-algebra elements.
//!
//! This translates the `PathRemoval` and quiver-algebra support methods in
//! `lib/pathalg.gi`.  In particular, support is computed by multiplying every
//! vertex idempotent on the indicated side, rather than inferred from the
//! stored paths.

use numerica::domains::Field;

use crate::{
    PathAlgebra, PathAlgebraElement, PathAlgebraError, QuotientError, QuotientPathAlgebra,
    QuotientPathAlgebraElement,
};

impl<F: Field> PathAlgebraElement<F> {
    /// Translate QPA's `PathRemoval(x, elist)`.
    ///
    /// Every nontrivial monomial visiting a vertex whose one-based position is
    /// in `vertex_positions` is discarded. Positions outside the quiver are
    /// ignored, as they are by GAP's `Position(...) in elist` test.
    ///
    /// # Errors
    ///
    /// The pinned QPA method indexes the last arrow of every walk. It therefore
    /// fails on a vertex monomial, whose walk is empty; Rust reports the same
    /// source boundary as [`PathAlgebraError::VertexTermInPathRemoval`].
    pub fn path_removal(&self, vertex_positions: &[usize]) -> Result<Self, PathAlgebraError> {
        let mut result = self.algebra().zero();
        for term in self.terms() {
            let walk = term
                .path()
                .walk()
                .ok_or(PathAlgebraError::InvalidPathTermInPathRemoval)?;
            let Some(last_arrow_index) = walk.last() else {
                return Err(PathAlgebraError::VertexTermInPathRemoval);
            };
            let quiver = self.algebra().quiver();
            let mut visits_removed_vertex = false;
            for arrow_index in walk {
                let arrow = quiver
                    .arrow(*arrow_index)
                    .ok_or(PathAlgebraError::InvalidPathTermInPathRemoval)?;
                if vertex_positions.contains(&arrow.source().index().get()) {
                    visits_removed_vertex = true;
                    break;
                }
            }
            if !visits_removed_vertex {
                let last_arrow = quiver
                    .arrow(*last_arrow_index)
                    .ok_or(PathAlgebraError::InvalidPathTermInPathRemoval)?;
                visits_removed_vertex =
                    vertex_positions.contains(&last_arrow.target().index().get());
            }
            if !visits_removed_vertex {
                let monomial = self.algebra().element(term.path().clone())?;
                let rebuilt = self
                    .algebra()
                    .one()
                    .multiply(&monomial)?
                    .scale_left(term.coefficient())?;
                result = result.add(&rebuilt)?;
            }
        }
        Ok(result)
    }
}

impl<F: Field> PathAlgebra<F> {
    /// Translate QPA's `LeftSupportOfQuiverAlgebraElement(A, m)`.
    ///
    /// # Errors
    ///
    /// The element must belong to this identical path-algebra family.
    pub fn left_support_of_quiver_algebra_element(
        &self,
        element: &PathAlgebraElement<F>,
    ) -> Result<Vec<usize>, PathAlgebraError> {
        if !self.is_identical(element.algebra()) {
            return Err(PathAlgebraError::DifferentAlgebras);
        }
        let mut support = Vec::new();
        for vertex in self.quiver().vertices() {
            let idempotent = self.one().multiply(&self.element(vertex.path())?)?;
            if !idempotent.multiply(element)?.is_zero() {
                support.push(vertex.index().get());
            }
        }
        Ok(support)
    }

    /// Translate QPA's `RightSupportOfQuiverAlgebraElement(A, m)`.
    ///
    /// # Errors
    ///
    /// The element must belong to this identical path-algebra family.
    pub fn right_support_of_quiver_algebra_element(
        &self,
        element: &PathAlgebraElement<F>,
    ) -> Result<Vec<usize>, PathAlgebraError> {
        if !self.is_identical(element.algebra()) {
            return Err(PathAlgebraError::DifferentAlgebras);
        }
        let mut support = Vec::new();
        for vertex in self.quiver().vertices() {
            let idempotent = self.one().multiply(&self.element(vertex.path())?)?;
            if !element.multiply(&idempotent)?.is_zero() {
                support.push(vertex.index().get());
            }
        }
        Ok(support)
    }

    /// Return QPA's pair `[left support, right support]`.
    ///
    /// # Errors
    ///
    /// Propagates the two source-translated support computations.
    pub fn support_of_quiver_algebra_element(
        &self,
        element: &PathAlgebraElement<F>,
    ) -> Result<[Vec<usize>; 2], PathAlgebraError> {
        Ok([
            self.left_support_of_quiver_algebra_element(element)?,
            self.right_support_of_quiver_algebra_element(element)?,
        ])
    }
}

impl<F: Field> QuotientPathAlgebra<F> {
    /// Translate QPA's left-support method using actual quotient products.
    ///
    /// # Errors
    ///
    /// The element must belong to this identical quotient family; normal-form
    /// and multiplication failures are propagated.
    pub fn left_support_of_quiver_algebra_element(
        &self,
        element: &QuotientPathAlgebraElement<F>,
    ) -> Result<Vec<usize>, QuotientError> {
        if !self.is_identical(element.algebra()) {
            return Err(QuotientError::DifferentQuotients);
        }
        let mut support = Vec::new();
        for vertex in self.original_path_algebra().quiver().vertices() {
            let idempotent = self.one()?.multiply(&self.path(vertex.path())?)?;
            if !idempotent.multiply(element)?.is_zero() {
                support.push(vertex.index().get());
            }
        }
        Ok(support)
    }

    /// Translate QPA's right-support method using actual quotient products.
    ///
    /// # Errors
    ///
    /// The element must belong to this identical quotient family; normal-form
    /// and multiplication failures are propagated.
    pub fn right_support_of_quiver_algebra_element(
        &self,
        element: &QuotientPathAlgebraElement<F>,
    ) -> Result<Vec<usize>, QuotientError> {
        if !self.is_identical(element.algebra()) {
            return Err(QuotientError::DifferentQuotients);
        }
        let mut support = Vec::new();
        for vertex in self.original_path_algebra().quiver().vertices() {
            let idempotent = self.one()?.multiply(&self.path(vertex.path())?)?;
            if !element.multiply(&idempotent)?.is_zero() {
                support.push(vertex.index().get());
            }
        }
        Ok(support)
    }

    /// Return QPA's pair `[left support, right support]`.
    ///
    /// # Errors
    ///
    /// Propagates the two source-translated quotient support computations.
    pub fn support_of_quiver_algebra_element(
        &self,
        element: &QuotientPathAlgebraElement<F>,
    ) -> Result<[Vec<usize>; 2], QuotientError> {
        Ok([
            self.left_support_of_quiver_algebra_element(element)?,
            self.right_support_of_quiver_algebra_element(element)?,
        ])
    }
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        Ring,
        finite_field::Z2,
        integer::Integer,
        rational::{Q, RationalField},
    };
    use serde_json::{Value, json};

    use crate::{
        ArrowSpec, PathAlgebra, PathAlgebraElement, PathAlgebraError, Quiver, QuotientPathAlgebra,
        VertexIndex,
    };

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn fixture_quiver() -> Quiver {
        Quiver::from_indexed(
            4,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(2), vertex(3), "b"),
                ArrowSpec::named(vertex(1), vertex(3), "c"),
                ArrowSpec::named(vertex(3), vertex(4), "d"),
                ArrowSpec::named(vertex(2), vertex(4), "e"),
            ],
        )
        .unwrap()
    }

    fn small_integer(value: &Integer) -> i64 {
        value
            .to_string()
            .parse()
            .expect("the element-utilities fixture uses only small integers")
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

    #[test]
    #[allow(clippy::many_single_char_names, clippy::too_many_lines)]
    fn pinned_qpa_element_utilities_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_element_utilities.schema-v1.json"
        ))
        .unwrap();
        let quiver = fixture_quiver();
        let algebra = PathAlgebra::new(Q, quiver.clone()).unwrap();
        let generators = algebra.generators();
        let v1 = &generators[0];
        let a = &generators[4];
        let b = &generators[5];
        let c = &generators[6];
        let d = &generators[7];
        let ab = a.multiply(b).unwrap();
        let bd = b.multiply(d).unwrap();
        let q2 = Q.nth(Integer::from(2));
        let q3 = Q.nth(Integer::from(3));
        let q5 = Q.nth(Integer::from(5));
        let q7 = Q.nth(Integer::from(7));
        let x = ab
            .scale_left(&q2)
            .unwrap()
            .add(&c.scale_left(&q3).unwrap())
            .unwrap()
            .add(&bd.scale_left(&q5).unwrap())
            .unwrap()
            .add(&d.scale_left(&q7).unwrap())
            .unwrap();
        let m = v1
            .add(&a.scale_left(&q2).unwrap())
            .unwrap()
            .add(&ab.scale_left(&q3).unwrap())
            .unwrap()
            .add(&b.scale_left(&q5).unwrap())
            .unwrap()
            .add(&d.scale_left(&q7).unwrap())
            .unwrap();

        let relations = [ab.sub(c).unwrap(), bd.clone()];
        let quotient =
            QuotientPathAlgebra::by_relations_high_level(algebra.clone(), &relations).unwrap();
        let quotient_generators = quotient.generators().unwrap();
        let qm = quotient_generators[0]
            .add(&quotient_generators[4].scale_left(&q2).unwrap())
            .unwrap()
            .add(&quotient_generators[5].scale_left(&q3).unwrap())
            .unwrap()
            .add(&quotient_generators[7].scale_left(&q5).unwrap())
            .unwrap();
        let qab = quotient_generators[4]
            .multiply(&quotient_generators[5])
            .unwrap();
        let qmix = qab.add(&quotient_generators[6]).unwrap();
        let qcancelled = qab.sub(&quotient_generators[6]).unwrap();

        assert_eq!(
            json!({
                "path_removal": {
                    "none": rational_element_json(&x.path_removal(&[]).unwrap()),
                    "two": rational_element_json(&x.path_removal(&[2]).unwrap()),
                    "four": rational_element_json(&x.path_removal(&[4]).unwrap()),
                    "ignored": rational_element_json(&x.path_removal(&[0, 5, 99]).unwrap())
                },
                "support": {
                    "left": algebra.left_support_of_quiver_algebra_element(&m).unwrap(),
                    "right": algebra.right_support_of_quiver_algebra_element(&m).unwrap(),
                    "both": algebra.support_of_quiver_algebra_element(&m).unwrap(),
                    "zero": algebra.support_of_quiver_algebra_element(&algebra.zero()).unwrap()
                },
                "quotient": {
                    "support": quotient.support_of_quiver_algebra_element(&qm).unwrap(),
                    "mix_support": quotient.support_of_quiver_algebra_element(&qmix).unwrap(),
                    "cancelled_support": quotient.support_of_quiver_algebra_element(&qcancelled).unwrap()
                }
            }),
            expected["rational"]
        );

        let algebra2 = PathAlgebra::new(Z2, quiver).unwrap();
        let generators2 = algebra2.generators();
        let x2 = generators2[4]
            .add(&generators2[5].multiply(&generators2[7]).unwrap())
            .unwrap()
            .add(&generators2[7])
            .unwrap();
        assert_eq!(
            json!({
                "removed_two": z2_element_json(&x2.path_removal(&[2]).unwrap()),
                "removed_four": z2_element_json(&x2.path_removal(&[4]).unwrap()),
                "support": algebra2.support_of_quiver_algebra_element(&x2).unwrap()
            }),
            expected["f2"]
        );
    }

    #[test]
    fn path_removal_preserves_qpa_vertex_term_failure() {
        let algebra = PathAlgebra::new(Q, Quiver::from_indexed(1, vec![]).unwrap()).unwrap();
        assert_eq!(
            algebra.generators()[0].path_removal(&[1]),
            Err(PathAlgebraError::VertexTermInPathRemoval)
        );
    }

    #[test]
    fn support_rejects_foreign_algebra_families() {
        let first = PathAlgebra::new(Q, fixture_quiver()).unwrap();
        let second = PathAlgebra::new(Q, fixture_quiver()).unwrap();
        assert_eq!(
            first.left_support_of_quiver_algebra_element(&second.one()),
            Err(PathAlgebraError::DifferentAlgebras)
        );
    }
}
