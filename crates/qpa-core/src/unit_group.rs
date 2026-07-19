//! Finite unit groups of QPA quiver algebras.
//!
//! This is a direct translation of `Units` in `lib/pathalg.gi`: enumerate the
//! nonzero diagonal vertex combinations, add every element of the algebra
//! radical, apply GAP's stable `Unique`, construct the multiplication table,
//! and cache the ordered algebra elements as `ElementsOfUnits`.

use std::cmp::Ordering;

use numerica::domains::Field;
use thiserror::Error;

use crate::{
    CoefficientDomainError, QpaCoefficientField, QuiverAlgebra, QuiverAlgebraElement,
    QuiverAlgebraElementError, QuotientError, StructureConstantsAlgebra, StructureConstantsError,
};

/// QPA's finite table-group representation returned by `Units(A)`.
///
/// Elements retain QPA's `Unique` order. The raw multiplication table stores
/// one-based positions, exactly as the list passed to GAP's
/// `MagmaByMultiplicationTable`; the index helpers use ordinary zero-based
/// Rust indices.
#[derive(Clone, Debug)]
pub struct QuiverAlgebraUnitGroup<F: Field> {
    algebra: QuiverAlgebra<F>,
    elements: Vec<QuiverAlgebraElement<F>>,
    coordinates: Vec<Vec<F::Element>>,
    multiplication_table: Vec<Box<[usize]>>,
}

impl<F: Field> QuiverAlgebraUnitGroup<F> {
    /// Return the quiver algebra whose unit group this represents.
    #[must_use]
    pub const fn algebra(&self) -> &QuiverAlgebra<F> {
        &self.algebra
    }

    /// Return `ElementsOfUnits(A)` in QPA's cached order.
    #[must_use]
    pub fn elements(&self) -> &[QuiverAlgebraElement<F>] {
        &self.elements
    }

    /// Return unit coordinates in QPA's canonical algebra basis.
    #[must_use]
    pub fn coordinates(&self) -> &[Vec<F::Element>] {
        &self.coordinates
    }

    /// Return the order of the finite unit group.
    #[must_use]
    pub fn order(&self) -> usize {
        self.elements.len()
    }

    /// Return QPA's one-based multiplication table.
    #[must_use]
    pub fn multiplication_table(&self) -> &[Box<[usize]>] {
        &self.multiplication_table
    }

    /// Return the zero-based index of a product, when both indices exist.
    #[must_use]
    pub fn product_index(&self, left: usize, right: usize) -> Option<usize> {
        self.multiplication_table
            .get(left)?
            .get(right)?
            .checked_sub(1)
    }

    /// Return the zero-based identity index.
    ///
    /// QPA's field enumeration starts with zero and one, so the all-one
    /// diagonal plus the zero radical element is always the first unit.
    #[must_use]
    pub const fn identity_index(&self) -> usize {
        0
    }

    /// Return the zero-based inverse index of a group element.
    #[must_use]
    pub fn inverse_index(&self, index: usize) -> Option<usize> {
        (index < self.order()).then_some(())?;
        (0..self.order()).find(|&candidate| {
            self.product_index(index, candidate) == Some(self.identity_index())
                && self.product_index(candidate, index) == Some(self.identity_index())
        })
    }
}

impl<F: Field> QuiverAlgebra<F> {
    /// Return GAP's `HasElementsOfUnits(A)` attribute state.
    #[must_use]
    pub fn has_elements_of_units(&self) -> bool {
        match self {
            Self::Path(algebra) => algebra.cached_elements_of_units_coordinates().is_some(),
            Self::Quotient(algebra) => algebra.cached_elements_of_units_coordinates().is_some(),
        }
    }

    /// Return the cached `ElementsOfUnits(A)` list.
    ///
    /// # Errors
    ///
    /// Returns [`QuiverAlgebraUnitGroupError::UnboundElementsOfUnits`] before
    /// `Units(A)` has run, and propagates canonical-basis reconstruction
    /// failures.
    pub fn elements_of_units(
        &self,
    ) -> Result<Vec<QuiverAlgebraElement<F>>, QuiverAlgebraUnitGroupError> {
        let coordinates = self
            .cached_elements_of_units_coordinates()
            .ok_or(QuiverAlgebraUnitGroupError::UnboundElementsOfUnits)?;
        let basis = finite_canonical_basis(self)?;
        coordinates
            .iter()
            .map(|row| element_from_coordinates(self, &basis, row))
            .collect()
    }

    fn cached_elements_of_units_coordinates(&self) -> Option<&[Vec<F::Element>]> {
        match self {
            Self::Path(algebra) => algebra.cached_elements_of_units_coordinates(),
            Self::Quotient(algebra) => algebra.cached_elements_of_units_coordinates(),
        }
    }

    fn cache_elements_of_units_coordinates(
        &self,
        coordinates: Vec<Vec<F::Element>>,
    ) -> &[Vec<F::Element>] {
        match self {
            Self::Path(algebra) => algebra.cache_elements_of_units_coordinates(coordinates),
            Self::Quotient(algebra) => algebra.cache_elements_of_units_coordinates(coordinates),
        }
    }
}

impl<F: QpaCoefficientField> QuiverAlgebra<F> {
    /// Translate QPA's finite-dimensional, finite-field `Units(A)` method.
    ///
    /// Quotients must already carry the complete Gröbner basis used by QPA's
    /// canonical-basis machinery. The `qpa-gbnp` adapter installs that lazy
    /// default basis before delegating here.
    ///
    /// # Errors
    ///
    /// Preserves QPA's finite-dimensional check before its finite-field check,
    /// and propagates exact canonical-basis, radical, arithmetic, and table
    /// construction failures.
    #[allow(clippy::too_many_lines)]
    pub fn units(&self) -> Result<QuiverAlgebraUnitGroup<F>, QuiverAlgebraUnitGroupError> {
        let basis = finite_canonical_basis(self)?;
        let Some(field_elements) = self.field().qpa_elements()? else {
            return Err(QuiverAlgebraUnitGroupError::InfiniteCoefficientField);
        };
        let dimension = basis.len();

        // QPA's generic `RadicalOfAlgebra(A)` sees the canonical basis and its
        // multiplication table. Keep that dependency instead of substituting
        // a path-length description of the radical.
        let mut products = Vec::new();
        for left in &basis {
            for right in &basis {
                let product = left.multiply(right)?;
                products.extend(canonical_coordinates(self, &basis, &product)?);
            }
        }
        let structure = StructureConstantsAlgebra::new(self.field().clone(), dimension, products)?;

        let number_of_vertices = self.quiver().number_of_vertices();
        let generators = self.generators_with_one()?;
        let idempotent_coordinates = generators
            .iter()
            .take(number_of_vertices)
            .map(|element| canonical_coordinates(self, &basis, element))
            .collect::<Result<Vec<_>, QuiverAlgebraUnitGroupError>>()?;

        // `Filtered(Elements(k^n), x -> Product(x) <> Zero(k))`.
        let nonzero_field_elements = field_elements
            .iter()
            .filter(|element| !self.field().is_zero(element))
            .cloned()
            .collect::<Vec<_>>();
        let diagonal_scalars = cartesian_power(&nonzero_field_elements, number_of_vertices)?;
        let mut invertibles = Vec::with_capacity(diagonal_scalars.len());
        for scalars in diagonal_scalars {
            let mut coordinates = vec![self.field().zero(); dimension];
            for (scalar, idempotent) in scalars.iter().zip(&idempotent_coordinates) {
                for (coordinate, coefficient) in coordinates.iter_mut().zip(idempotent) {
                    self.field().add_mul_assign(coordinate, scalar, coefficient);
                }
            }
            invertibles.push(coordinates);
        }

        let radical = structure.radical()?;
        let radical_scalars = cartesian_power(&field_elements, radical.basis().nrows())?;
        let mut radical_elements = Vec::with_capacity(radical_scalars.len());
        for scalars in radical_scalars {
            let coordinates = radical.lift_coordinates(&scalars)?;
            let element = element_from_coordinates(self, &basis, &coordinates)?;
            insert_in_qpa_set_order(&mut radical_elements, (coordinates, element))?;
        }
        let radical_plus_initial = radical_elements
            .len()
            .checked_add(1)
            .ok_or(QuiverAlgebraUnitGroupError::EnumerationSizeOverflow)?;
        let generator_capacity = invertibles
            .len()
            .checked_mul(radical_plus_initial)
            .ok_or(QuiverAlgebraUnitGroupError::EnumerationSizeOverflow)?;
        let mut generator_coordinates = Vec::with_capacity(generator_capacity);
        generator_coordinates.extend(invertibles.iter().cloned());
        for (radical_element, _) in radical_elements {
            for invertible in &invertibles {
                generator_coordinates.push(
                    invertible
                        .iter()
                        .zip(&radical_element)
                        .map(|(left, right)| self.field().add(left, right))
                        .collect(),
                );
            }
        }
        let coordinates = stable_unique(generator_coordinates);

        // QPA's `MultiplicationTable(generators)` stores one-based list
        // positions before `MagmaByMultiplicationTable` and `AsGroup`.
        let mut multiplication_table = Vec::with_capacity(coordinates.len());
        for left in &coordinates {
            let mut row = Vec::with_capacity(coordinates.len());
            for right in &coordinates {
                let product = structure.multiply_coordinates(left, right)?;
                let position = coordinates
                    .iter()
                    .position(|candidate| candidate == &product)
                    .ok_or(QuiverAlgebraUnitGroupError::ProductOutsideUnitList)?;
                row.push(
                    position
                        .checked_add(1)
                        .ok_or(QuiverAlgebraUnitGroupError::EnumerationSizeOverflow)?,
                );
            }
            multiplication_table.push(row.into_boxed_slice());
        }

        let elements = coordinates
            .iter()
            .map(|row| element_from_coordinates(self, &basis, row))
            .collect::<Result<Vec<_>, _>>()?;
        self.cache_elements_of_units_coordinates(coordinates.clone());
        Ok(QuiverAlgebraUnitGroup {
            algebra: self.clone(),
            elements,
            coordinates,
            multiplication_table,
        })
    }
}

fn finite_canonical_basis<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Vec<QuiverAlgebraElement<F>>, QuiverAlgebraUnitGroupError> {
    match algebra {
        QuiverAlgebra::Path(path_algebra) => path_algebra
            .canonical_basis()
            .ok_or(QuiverAlgebraUnitGroupError::InfiniteDimensionalAlgebra)
            .map(|basis| basis.into_iter().map(QuiverAlgebraElement::Path).collect()),
        QuiverAlgebra::Quotient(quotient) => {
            if !quotient.is_finite_dimensional()? {
                return Err(QuiverAlgebraUnitGroupError::InfiniteDimensionalAlgebra);
            }
            quotient
                .canonical_basis()?
                .ok_or(QuiverAlgebraUnitGroupError::MissingFiniteCanonicalBasis)
                .map(|basis| {
                    basis
                        .into_iter()
                        .map(QuiverAlgebraElement::Quotient)
                        .collect()
                })
        }
    }
}

fn canonical_coordinates<F: Field>(
    algebra: &QuiverAlgebra<F>,
    basis: &[QuiverAlgebraElement<F>],
    element: &QuiverAlgebraElement<F>,
) -> Result<Vec<F::Element>, QuiverAlgebraUnitGroupError> {
    match (algebra, element) {
        (QuiverAlgebra::Path(path_algebra), QuiverAlgebraElement::Path(element)) => {
            if !path_algebra.is_identical(element.algebra()) {
                return Err(QuiverAlgebraUnitGroupError::DifferentAlgebraFamily);
            }
            basis
                .iter()
                .map(|basis_element| {
                    let QuiverAlgebraElement::Path(basis_element) = basis_element else {
                        return Err(QuiverAlgebraUnitGroupError::InvalidCanonicalBasis);
                    };
                    let [basis_term] = basis_element.terms() else {
                        return Err(QuiverAlgebraUnitGroupError::InvalidCanonicalBasis);
                    };
                    Ok(element
                        .terms()
                        .iter()
                        .find(|term| term.path() == basis_term.path())
                        .map_or_else(
                            || path_algebra.field().zero(),
                            |term| term.coefficient().clone(),
                        ))
                })
                .collect()
        }
        (QuiverAlgebra::Quotient(quotient), QuiverAlgebraElement::Quotient(element)) => {
            Ok(quotient.canonical_basis_coefficients(element)?)
        }
        _ => Err(QuiverAlgebraUnitGroupError::DifferentAlgebraVariant),
    }
}

fn element_from_coordinates<F: Field>(
    algebra: &QuiverAlgebra<F>,
    basis: &[QuiverAlgebraElement<F>],
    coordinates: &[F::Element],
) -> Result<QuiverAlgebraElement<F>, QuiverAlgebraUnitGroupError> {
    if coordinates.len() != basis.len() {
        return Err(QuiverAlgebraUnitGroupError::CoordinateLength {
            actual: coordinates.len(),
            expected: basis.len(),
        });
    }
    let mut result = algebra.zero();
    for (basis_element, coefficient) in basis.iter().zip(coordinates) {
        if !algebra.field().is_zero(coefficient) {
            result = result.add(&basis_element.scale_left(coefficient)?)?;
        }
    }
    Ok(result)
}

fn cartesian_power<T: Clone>(
    elements: &[T],
    exponent: usize,
) -> Result<Vec<Vec<T>>, QuiverAlgebraUnitGroupError> {
    let mut size = 1usize;
    for _ in 0..exponent {
        size = size
            .checked_mul(elements.len())
            .ok_or(QuiverAlgebraUnitGroupError::EnumerationSizeOverflow)?;
    }
    let mut rows = Vec::with_capacity(size);
    rows.push(Vec::new());
    for _ in 0..exponent {
        let mut expanded = Vec::with_capacity(
            rows.len()
                .checked_mul(elements.len())
                .ok_or(QuiverAlgebraUnitGroupError::EnumerationSizeOverflow)?,
        );
        for row in rows {
            for element in elements {
                let mut next = row.clone();
                next.push(element.clone());
                expanded.push(next);
            }
        }
        rows = expanded;
    }
    Ok(rows)
}

fn stable_unique<T: Eq>(rows: Vec<T>) -> Vec<T> {
    let mut result = Vec::new();
    for row in rows {
        if !result.contains(&row) {
            result.push(row);
        }
    }
    result
}

fn insert_in_qpa_set_order<F: Field>(
    elements: &mut Vec<(Vec<F::Element>, QuiverAlgebraElement<F>)>,
    candidate: (Vec<F::Element>, QuiverAlgebraElement<F>),
) -> Result<(), QuiverAlgebraUnitGroupError> {
    let mut position = 0;
    while position < elements.len()
        && qpa_element_compare(&elements[position].1, &candidate.1)? != Ordering::Greater
    {
        position += 1;
    }
    elements.insert(position, candidate);
    Ok(())
}

fn qpa_element_compare<F: Field>(
    left: &QuiverAlgebraElement<F>,
    right: &QuiverAlgebraElement<F>,
) -> Result<Ordering, QuiverAlgebraUnitGroupError> {
    match (left, right) {
        (QuiverAlgebraElement::Path(left), QuiverAlgebraElement::Path(right)) => left
            .compare(right)
            .map_err(QuiverAlgebraElementError::from)
            .map_err(QuiverAlgebraUnitGroupError::from),
        (QuiverAlgebraElement::Quotient(left), QuiverAlgebraElement::Quotient(right)) => {
            Ok(left.compare(right)?)
        }
        _ => Err(QuiverAlgebraUnitGroupError::DifferentAlgebraVariant),
    }
}

/// Failure in QPA's finite unit-group construction.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum QuiverAlgebraUnitGroupError {
    /// QPA rejects the algebra before inspecting its coefficient field.
    #[error("the entered algebra is not finite dimensional")]
    InfiniteDimensionalAlgebra,
    /// QPA's finite-dimensional quotient basis object had no vectors.
    #[error("the finite-dimensional quotient has no canonical basis vectors")]
    MissingFiniteCanonicalBasis,
    /// QPA rejects infinite coefficient fields.
    #[error("the algebra is not over a finite field")]
    InfiniteCoefficientField,
    /// `ElementsOfUnits` has no installed computation and `Units` has not set it.
    #[error("ElementsOfUnits is not bound for this algebra")]
    UnboundElementsOfUnits,
    /// An internal canonical-basis element was not a single path monomial.
    #[error("the canonical path-algebra basis contains a nonmonomial element")]
    InvalidCanonicalBasis,
    /// A coordinate conversion crossed distinct algebra variants.
    #[error("the element and canonical basis use different algebra variants")]
    DifferentAlgebraVariant,
    /// A coordinate conversion crossed distinct algebra objects.
    #[error("the element and canonical basis belong to different algebra families")]
    DifferentAlgebraFamily,
    /// A cached coordinate row no longer matches the immutable canonical basis.
    #[error("the unit has {actual} coordinates, expected {expected}")]
    CoordinateLength {
        /// Cached coordinate count.
        actual: usize,
        /// Canonical-basis dimension.
        expected: usize,
    },
    /// A finite Cartesian enumeration or one-based table position overflowed.
    #[error("the finite unit enumeration overflows addressable storage")]
    EnumerationSizeOverflow,
    /// Algebra multiplication unexpectedly left the enumerated unit list.
    #[error("a product is absent from the enumerated unit list")]
    ProductOutsideUnitList,
    /// Finite-field enumeration failed.
    #[error(transparent)]
    CoefficientDomain(#[from] CoefficientDomainError),
    /// Canonical-basis or quotient normal-form construction failed.
    #[error(transparent)]
    Quotient(#[from] QuotientError),
    /// Path- or quotient-algebra arithmetic failed.
    #[error(transparent)]
    Element(#[from] QuiverAlgebraElementError),
    /// Radical or structure-constant arithmetic failed.
    #[error(transparent)]
    StructureConstants(#[from] StructureConstantsError),
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        Field, Ring,
        finite_field::{FiniteFieldCore, Z2, Zp},
        rational::Q,
    };
    use serde_json::{Value, json};

    use crate::{ArrowSpec, PathAlgebra, Quiver, QuiverAlgebra, QuotientPathAlgebra, VertexIndex};

    use super::{QuiverAlgebraUnitGroup, QuiverAlgebraUnitGroupError};

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn table_json<F: Field>(group: &QuiverAlgebraUnitGroup<F>) -> Value {
        json!(
            group
                .multiplication_table()
                .iter()
                .map(|row| row.to_vec())
                .collect::<Vec<_>>()
        )
    }

    fn z2_coordinates_json(group: &QuiverAlgebraUnitGroup<Z2>) -> Value {
        json!(
            group
                .coordinates()
                .iter()
                .map(|row| {
                    row.iter()
                        .map(|coefficient| i64::from(!Z2.is_zero(coefficient)))
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        )
    }

    fn zp_coordinates_json(group: &QuiverAlgebraUnitGroup<Zp>) -> Value {
        let field = group.algebra().field();
        json!(
            group
                .coordinates()
                .iter()
                .map(|row| {
                    row.iter()
                        .map(|coefficient| {
                            field
                                .to_integer(coefficient)
                                .to_string()
                                .parse::<i64>()
                                .unwrap()
                        })
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        )
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn pinned_qpa_pathalg_units_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_units.schema-v1.json"
        ))
        .unwrap();
        assert_eq!(
            expected["oracle"],
            json!({
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            })
        );

        let binary_a2 = QuiverAlgebra::Path(
            PathAlgebra::new(
                Z2,
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        assert!(!binary_a2.has_elements_of_units());
        assert_eq!(
            binary_a2.elements_of_units().unwrap_err(),
            QuiverAlgebraUnitGroupError::UnboundElementsOfUnits
        );
        let binary_a2_group = binary_a2.units().unwrap();
        assert_eq!(binary_a2_group.identity_index(), 0);
        assert_eq!(binary_a2_group.inverse_index(1), Some(1));
        assert!(binary_a2.has_elements_of_units());
        assert_eq!(
            binary_a2_group.elements(),
            binary_a2.elements_of_units().unwrap()
        );
        assert_eq!(
            json!({
                "had_elements": false,
                "has_elements": binary_a2.has_elements_of_units(),
                "is_group": true,
                "size": binary_a2_group.order(),
                "elements": z2_coordinates_json(&binary_a2_group),
                "table": table_json(&binary_a2_group)
            }),
            expected["binary_a2"]
        );

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let field3 = Zp::new(3);
        let dual_parent = PathAlgebra::new(field3.clone(), loop_quiver.clone()).unwrap();
        let x = dual_parent.generators()[1].clone();
        let dual = QuiverAlgebra::Quotient(
            QuotientPathAlgebra::by_relations_high_level(dual_parent, &[x.pow(2).unwrap()])
                .unwrap(),
        );
        let dual_group = dual.units().unwrap();
        assert_eq!(
            json!({
                "is_group": true,
                "size": dual_group.order(),
                "elements": zp_coordinates_json(&dual_group),
                "table": table_json(&dual_group)
            }),
            expected["ternary_dual_numbers"]
        );

        let cubic_parent = PathAlgebra::new(Z2, loop_quiver).unwrap();
        let x = cubic_parent.generators()[1].clone();
        let cubic = QuiverAlgebra::Quotient(
            QuotientPathAlgebra::by_relations_high_level(cubic_parent, &[x.pow(3).unwrap()])
                .unwrap(),
        );
        let cubic_group = cubic.units().unwrap();
        assert_eq!(
            json!({
                "is_group": true,
                "size": cubic_group.order(),
                "elements": z2_coordinates_json(&cubic_group),
                "table": table_json(&cubic_group)
            }),
            expected["binary_cubic_numbers"]
        );

        let semisimple = QuiverAlgebra::Path(
            PathAlgebra::new(
                field3,
                Quiver::from_named(vec!["u".into(), "v".into()], vec![]).unwrap(),
            )
            .unwrap(),
        );
        let semisimple_group = semisimple.units().unwrap();
        assert_eq!(
            json!({
                "is_group": true,
                "size": semisimple_group.order(),
                "elements": zp_coordinates_json(&semisimple_group),
                "table": table_json(&semisimple_group)
            }),
            expected["ternary_semisimple"]
        );

        let quinary = QuiverAlgebra::Path(
            PathAlgebra::new(Zp::new(5), Quiver::from_indexed(1, vec![]).unwrap()).unwrap(),
        );
        let quinary_group = quinary.units().unwrap();
        assert_eq!(
            json!({
                "is_group": true,
                "size": quinary_group.order(),
                "elements": zp_coordinates_json(&quinary_group),
                "table": table_json(&quinary_group)
            }),
            expected["quinary_a1"]
        );

        let infinite = QuiverAlgebra::Path(
            PathAlgebra::new(
                Z2,
                Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "z")]).unwrap(),
            )
            .unwrap(),
        );
        let rational = QuiverAlgebra::Path(
            PathAlgebra::new(Q, Quiver::from_indexed(1, vec![]).unwrap()).unwrap(),
        );
        assert_eq!(
            infinite.units().unwrap_err(),
            QuiverAlgebraUnitGroupError::InfiniteDimensionalAlgebra
        );
        assert_eq!(
            rational.units().unwrap_err(),
            QuiverAlgebraUnitGroupError::InfiniteCoefficientField
        );
        assert_eq!(
            json!({
                "infinite_success": false,
                "rational_success": false
            }),
            expected["errors"]
        );
    }
}
