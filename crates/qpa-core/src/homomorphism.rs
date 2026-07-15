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
        let quotient_unit_offset = usize::from(matches!(self.source, QuiverAlgebra::Quotient(_)));
        self.generator_images
            .get(
                one_based_index
                    .wrapping_sub(1)
                    .saturating_add(quotient_unit_offset),
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
    use numerica::domains::rational::Q;

    use crate::{NamedArrowSpec, PathAlgebra, Quiver, QuiverAlgebra};

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
}
