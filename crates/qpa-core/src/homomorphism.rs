//! Homomorphisms between QPA quiver algebras, stored by generator images.
//!
//! QPA tensor, opposite, and trivial-extension methods construct GAP algebra
//! homomorphisms and then explicitly retain `generators` and `genimages`.
//! This module gives those source-visible fields a typed Rust representation.

use numerica::domains::Field;
use thiserror::Error;

use crate::{QuiverAlgebra, QuiverAlgebraElement, QuiverAlgebraElementError};

/// A QPA algebra homomorphism specified on the ordinary quiver generators.
///
/// The stored generator order is vertices then arrows for path algebras, and
/// one then vertices then arrows for quotients, exactly QPA's
/// `GeneratorsOfAlgebra`. Images of arbitrary elements are evaluated by
/// replacing the quiver-generator indices in each path monomial, matching
/// GAP's mapping-by-images representation.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct QuiverAlgebraHomomorphism<F: Field> {
    source: QuiverAlgebra<F>,
    range: QuiverAlgebra<F>,
    generators: Vec<QuiverAlgebraElement<F>>,
    generator_images: Vec<QuiverAlgebraElement<F>>,
    generator_index_offset: usize,
}

impl<F: Field> QuiverAlgebraHomomorphism<F> {
    /// Construct the mapping on QPA's ordinary quiver-generator list.
    ///
    /// This is the typed analogue of the `AlgebraHomomorphismByImages` calls
    /// in QPA whose `gens` argument is `GeneratorsOfAlgebra(source)`.
    ///
    /// # Errors
    ///
    /// Rejects a list of the wrong length or an image outside the identical
    /// range algebra family, and propagates source-generator construction
    /// failures.
    pub fn by_generator_images(
        source: QuiverAlgebra<F>,
        range: QuiverAlgebra<F>,
        generator_images: Vec<QuiverAlgebraElement<F>>,
    ) -> Result<Self, QuiverAlgebraHomomorphismError> {
        let generators = source.generators()?;
        let generator_index_offset = usize::from(matches!(source, QuiverAlgebra::Quotient(_)));
        Self::by_explicit_generator_images(
            source,
            range,
            generators,
            generator_images,
            generator_index_offset,
        )
    }

    /// Construct a mapping on QPA's `GeneratorsOfAlgebraWithOne` list.
    ///
    /// For quotient path algebras this list is vertices followed by arrows,
    /// without the extra leading identity returned by ordinary
    /// `GeneratorsOfAlgebra`. Path algebras use the same list in both cases.
    ///
    /// # Errors
    ///
    /// Rejects a list of the wrong length or an image outside the identical
    /// range algebra family, and propagates source-generator construction
    /// failures.
    pub fn by_generators_with_one_images(
        source: QuiverAlgebra<F>,
        range: QuiverAlgebra<F>,
        generator_images: Vec<QuiverAlgebraElement<F>>,
    ) -> Result<Self, QuiverAlgebraHomomorphismError> {
        let generators = source.generators_with_one()?;
        Self::by_explicit_generator_images(source, range, generators, generator_images, 0)
    }

    fn by_explicit_generator_images(
        source: QuiverAlgebra<F>,
        range: QuiverAlgebra<F>,
        generators: Vec<QuiverAlgebraElement<F>>,
        generator_images: Vec<QuiverAlgebraElement<F>>,
        generator_index_offset: usize,
    ) -> Result<Self, QuiverAlgebraHomomorphismError> {
        if generator_images.len() != generators.len() {
            return Err(
                QuiverAlgebraHomomorphismError::GeneratorImageCountMismatch {
                    generators: generators.len(),
                    images: generator_images.len(),
                },
            );
        }
        for (position, image) in generator_images.iter().enumerate() {
            if !range.contains_element(image) {
                return Err(QuiverAlgebraHomomorphismError::ForeignRangeElement {
                    position: position + 1,
                });
            }
        }
        Ok(Self {
            source,
            range,
            generators,
            generator_images,
            generator_index_offset,
        })
    }

    /// Return QPA's `Source` algebra.
    #[must_use]
    pub const fn source(&self) -> &QuiverAlgebra<F> {
        &self.source
    }

    /// Return QPA's `Range` algebra.
    #[must_use]
    pub const fn range(&self) -> &QuiverAlgebra<F> {
        &self.range
    }

    /// Return the explicitly retained `generators` list.
    #[must_use]
    pub fn generators(&self) -> &[QuiverAlgebraElement<F>] {
        &self.generators
    }

    /// Return the explicitly retained `genimages` list.
    #[must_use]
    pub fn generator_images(&self) -> &[QuiverAlgebraElement<F>] {
        &self.generator_images
    }

    /// Return whether this map was entered on
    /// `GeneratorsOfAlgebraWithOne` rather than ordinary generators.
    #[must_use]
    pub const fn uses_generators_with_one(&self) -> bool {
        self.generator_index_offset == 0
    }

    /// Apply the homomorphism to one element of its identical source family.
    ///
    /// # Errors
    ///
    /// Rejects foreign source elements and propagates target arithmetic or
    /// quotient normal-form failures.
    pub fn image(
        &self,
        element: &QuiverAlgebraElement<F>,
    ) -> Result<QuiverAlgebraElement<F>, QuiverAlgebraHomomorphismError> {
        if !self.source.contains_element(element) {
            return Err(QuiverAlgebraHomomorphismError::ForeignSourceElement);
        }

        let terms = match element {
            QuiverAlgebraElement::Path(element) => element.terms(),
            QuiverAlgebraElement::Quotient(element) => element.representative().terms(),
        };
        let mut result = self.range.zero();
        for term in terms {
            let representation = term.path().external_representation();
            let Some((&first, rest)) = representation.split_first() else {
                return Err(QuiverAlgebraHomomorphismError::EmptyPathRepresentation);
            };
            let mut monomial_image = self.image_of_generator(first)?;
            for &generator in rest {
                monomial_image = monomial_image.multiply(&self.image_of_generator(generator)?)?;
            }
            result = result.add(&monomial_image.scale_left(term.coefficient())?)?;
        }
        Ok(result)
    }

    fn image_of_generator(
        &self,
        one_based_index: usize,
    ) -> Result<QuiverAlgebraElement<F>, QuiverAlgebraHomomorphismError> {
        self.generator_images
            .get(
                one_based_index
                    .wrapping_sub(1)
                    .saturating_add(self.generator_index_offset),
            )
            .cloned()
            .ok_or(QuiverAlgebraHomomorphismError::MissingGeneratorImage {
                index: one_based_index,
            })
    }
}

/// Failure while constructing or applying a QPA quiver-algebra homomorphism.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum QuiverAlgebraHomomorphismError {
    /// The ordinary generator and image lists must have equal length.
    #[error("received {images} images for {generators} algebra generators")]
    GeneratorImageCountMismatch {
        /// Number of source generators.
        generators: usize,
        /// Number of supplied images.
        images: usize,
    },
    /// Every stored generator image must belong to the range family.
    #[error("generator image {position} does not belong to the range algebra")]
    ForeignRangeElement {
        /// One-based position in the supplied image list.
        position: usize,
    },
    /// Mapping application is defined only on the identical source family.
    #[error("the element does not belong to the source algebra")]
    ForeignSourceElement,
    /// A nonzero source monomial must contain at least one quiver generator.
    #[error("a nonzero path has an empty external generator representation")]
    EmptyPathRepresentation,
    /// A path external representation referred outside the stored image list.
    #[error("the homomorphism has no image for generator index {index}")]
    MissingGeneratorImage {
        /// Missing one-based generator index.
        index: usize,
    },
    /// Source-generator construction or target arithmetic failed.
    #[error(transparent)]
    Element(#[from] QuiverAlgebraElementError),
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        finite_field::{FiniteFieldCore, Zp},
        rational::Q,
    };
    use serde_json::{Value, json};

    use crate::{
        GeneratorVariableError, NamedArrowSpec, PathAlgebra, Quiver, QuiverAlgebra,
        QuiverAlgebraElement, QuotientError, QuotientPathAlgebra,
    };

    use super::{QuiverAlgebraHomomorphism, QuiverAlgebraHomomorphismError};

    fn arrow_algebra() -> QuiverAlgebra<numerica::domains::rational::RationalField> {
        QuiverAlgebra::Path(
            PathAlgebra::new(
                Q,
                Quiver::from_named(
                    vec!["u".into(), "v".into()],
                    vec![NamedArrowSpec::named("u", "v", "a")],
                )
                .unwrap(),
            )
            .unwrap(),
        )
    }

    fn loop_algebra(field: Zp, first: &str, second: &str) -> QuiverAlgebra<Zp> {
        QuiverAlgebra::Path(
            PathAlgebra::new(
                field,
                Quiver::from_named(
                    vec!["v1".into()],
                    vec![
                        NamedArrowSpec::named("v1", "v1", first),
                        NamedArrowSpec::named("v1", "v1", second),
                    ],
                )
                .unwrap(),
            )
            .unwrap(),
        )
    }

    fn element_json(element: &QuiverAlgebraElement<Zp>, field: &Zp) -> Value {
        let representative = match element {
            QuiverAlgebraElement::Path(element) => element,
            QuiverAlgebraElement::Quotient(element) => element.representative(),
        };
        Value::Array(
            representative
                .terms()
                .iter()
                .map(|term| {
                    json!({
                        "path": term.path().to_string(),
                        "coefficient": field.to_integer(term.coefficient()).to_string().parse::<i64>().unwrap()
                    })
                })
                .collect(),
        )
    }

    #[test]
    fn identity_generator_images_evaluate_paths_and_sums() {
        let algebra = arrow_algebra();
        let generators = algebra.generators().unwrap();
        let map = QuiverAlgebraHomomorphism::by_generator_images(
            algebra.clone(),
            algebra.clone(),
            generators.clone(),
        )
        .unwrap();
        let element = generators[0].add(&generators[2]).unwrap();

        assert_eq!(map.generators(), generators);
        assert_eq!(map.generator_images(), generators);
        assert_eq!(map.image(&element).unwrap(), element);
        assert_eq!(map.image(&algebra.zero()).unwrap(), algebra.zero());
    }

    #[test]
    fn constructor_and_application_preserve_family_identity() {
        let source = arrow_algebra();
        let foreign = arrow_algebra();
        assert!(matches!(
            QuiverAlgebraHomomorphism::by_generator_images(
                source.clone(),
                source.clone(),
                foreign.generators().unwrap(),
            ),
            Err(QuiverAlgebraHomomorphismError::ForeignRangeElement { position: 1 })
        ));

        let map = QuiverAlgebraHomomorphism::by_generator_images(
            source.clone(),
            source.clone(),
            source.generators().unwrap(),
        )
        .unwrap();
        assert_eq!(
            map.image(&foreign.zero()).unwrap_err(),
            QuiverAlgebraHomomorphismError::ForeignSourceElement
        );
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn pinned_qpa_mapped_expression_and_ordered_by_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_mapping_ordering.schema-v1.json"
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

        let field = Zp::new(5);
        let source = loop_algebra(field.clone(), "a", "b");
        let range = loop_algebra(field.clone(), "x", "y");
        let source_generators = source.generators_with_one().unwrap();
        let range_generators = range.generators_with_one().unwrap();
        let range_images = vec![
            range_generators[0].clone(),
            range_generators[2].clone(),
            range_generators[1].clone(),
        ];
        let element = source_generators[0]
            .scale_left(&field.to_element(2))
            .unwrap()
            .add(
                &source_generators[1]
                    .multiply(&source_generators[1])
                    .unwrap()
                    .scale_left(&field.to_element(3))
                    .unwrap(),
            )
            .unwrap()
            .add(
                &source_generators[1]
                    .multiply(&source_generators[2])
                    .unwrap()
                    .scale_left(&field.to_element(4))
                    .unwrap(),
            )
            .unwrap()
            .add(&source_generators[2])
            .unwrap();
        let mapping = QuiverAlgebraHomomorphism::by_generators_with_one_images(
            source.clone(),
            range,
            range_images,
        )
        .unwrap();
        let mapped = mapping.image(&element).unwrap();
        assert_eq!(
            json!({
                "mapped": element_json(&mapped, &field),
                "mapped_zero": element_json(&mapping.image(&source.zero()).unwrap(), &field),
                "unretained_success": false,
                "retained_success": true,
                "retained_image": element_json(&mapped, &field)
            }),
            expected["path"]
        );

        let assignments = source.generator_variable_assignments().unwrap();
        let invalid = QuiverAlgebra::Path(
            PathAlgebra::new(
                field.clone(),
                Quiver::from_named(
                    vec!["v1".into()],
                    vec![NamedArrowSpec::named("v1", "v1", "not-valid")],
                )
                .unwrap(),
            )
            .unwrap(),
        );
        let read_only = QuiverAlgebra::Path(
            PathAlgebra::new(
                field.clone(),
                Quiver::from_named(
                    vec!["v1".into()],
                    vec![NamedArrowSpec::named("v1", "v1", "Z")],
                )
                .unwrap(),
            )
            .unwrap(),
        );
        let invalid_assignment = invalid.generator_variable_assignments();
        assert_eq!(
            invalid_assignment.as_ref().unwrap_err(),
            &GeneratorVariableError::InvalidIdentifier {
                name: "not-valid".into(),
            }
        );
        let read_only_assignment =
            read_only.generator_variable_assignments_with_read_only(|name| name == "Z");
        assert_eq!(
            read_only_assignment.as_ref().unwrap_err(),
            &GeneratorVariableError::ReadOnlyIdentifier { name: "Z".into() }
        );
        assert_eq!(
            json!({
                "success": true,
                "invalid_success": invalid_assignment.is_ok(),
                "read_only_success": read_only_assignment.is_ok(),
                "names": assignments
                    .iter()
                    .map(crate::GeneratorVariableAssignment::name)
                    .collect::<Vec<_>>(),
                "values": assignments
                    .iter()
                    .map(|assignment| element_json(assignment.element(), &field))
                    .collect::<Vec<_>>()
            }),
            expected["assignment"]
        );

        let QuiverAlgebra::Path(source_parent) = source else {
            unreachable!();
        };
        let source_parent_generators = source_parent.generators();
        let source_quotient = QuotientPathAlgebra::by_relations_high_level(
            source_parent,
            &[source_parent_generators[1]
                .multiply(&source_parent_generators[1])
                .unwrap()],
        )
        .unwrap();
        let range_parent = loop_algebra(field.clone(), "x", "y");
        let QuiverAlgebra::Path(range_parent) = range_parent else {
            unreachable!();
        };
        let range_parent_generators = range_parent.generators();
        let range_quotient = QuotientPathAlgebra::by_relations_high_level(
            range_parent,
            &[range_parent_generators[2]
                .multiply(&range_parent_generators[2])
                .unwrap()],
        )
        .unwrap();
        let source = QuiverAlgebra::Quotient(source_quotient.clone());
        let range = QuiverAlgebra::Quotient(range_quotient);
        let source_generators = source.generators_with_one().unwrap();
        let range_generators = range.generators_with_one().unwrap();
        let range_images = vec![
            range_generators[0].clone(),
            range_generators[2].clone(),
            range_generators[1].clone(),
        ];
        let element = source_generators[0]
            .scale_left(&field.to_element(2))
            .unwrap()
            .add(
                &source_generators[1]
                    .multiply(&source_generators[2])
                    .unwrap()
                    .scale_left(&field.to_element(3))
                    .unwrap(),
            )
            .unwrap()
            .add(
                &source_generators[2]
                    .multiply(&source_generators[2])
                    .unwrap()
                    .scale_left(&field.to_element(4))
                    .unwrap(),
            )
            .unwrap();
        let mapping =
            QuiverAlgebraHomomorphism::by_generators_with_one_images(source, range, range_images)
                .unwrap();
        let mapped = mapping.image(&element).unwrap();
        assert_eq!(
            json!({
                "mapped": element_json(&mapped, &field),
                "retained_success": true,
                "retained_image": element_json(&mapped, &field),
                "ordered_success": source_quotient
                    .ordered_by(source_quotient.ordering())
                    .is_ok()
            }),
            expected["quotient"]
        );
        assert_eq!(
            source_quotient
                .ordered_by(source_quotient.ordering())
                .unwrap_err(),
            QuotientError::UnboundRelatorsForOrderedBy
        );
    }
}
