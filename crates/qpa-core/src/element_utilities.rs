//! QPA utilities acting on individual path- and quotient-algebra elements.
//!
//! This translates the `PathRemoval` and quiver-algebra support methods in
//! `lib/pathalg.gi`.  In particular, support is computed by multiplying every
//! vertex idempotent on the indicated side, rather than inferred from the
//! stored paths.

use numerica::domains::Field;
use thiserror::Error;

use crate::{
    NamedArrowSpec, Path, PathAlgebra, PathAlgebraElement, PathAlgebraError, PathAlgebraModule,
    Quiver, QuotientError, QuotientPathAlgebra, QuotientPathAlgebraElement, TensorProductAlgebra,
};

/// A QPA quiver algebra that may be a path algebra or a proper quotient.
///
/// `QuiverAlgebraOfAmodAeA` returns the rebuilt path algebra itself when no
/// converted relations survive, while nonempty relations produce a quotient.
/// This enum retains that observable type distinction.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum QuiverAlgebra<F: Field> {
    /// QPA's empty-relation result.
    Path(PathAlgebra<F>),
    /// A quotient by the computed relation ideal.
    Quotient(QuotientPathAlgebra<F>),
}

impl<F: Field> QuiverAlgebra<F> {
    /// Return whether QPA returned an actual path algebra.
    #[must_use]
    pub const fn is_path_algebra(&self) -> bool {
        matches!(self, Self::Path(_))
    }

    /// Return the path algebra containing the sparse defining relations.
    #[must_use]
    pub fn original_path_algebra(&self) -> &PathAlgebra<F> {
        match self {
            Self::Path(algebra) => algebra,
            Self::Quotient(algebra) => algebra.original_path_algebra(),
        }
    }

    /// Return the common coefficient field of this path algebra or quotient.
    #[must_use]
    pub fn field(&self) -> &F {
        self.original_path_algebra().field()
    }

    /// Return QPA's `QuiverOfPathAlgebra` value.
    #[must_use]
    pub fn quiver(&self) -> &Quiver {
        self.original_path_algebra().quiver()
    }

    /// Return QPA's `RelationsOfAlgebra` value.
    #[must_use]
    pub fn relations(&self) -> &[PathAlgebraElement<F>] {
        match self {
            Self::Path(_) => &[],
            Self::Quotient(algebra) => algebra.relations(),
        }
    }

    /// Borrow the quotient result, if nonempty relations survived.
    #[must_use]
    pub const fn as_quotient(&self) -> Option<&QuotientPathAlgebra<F>> {
        match self {
            Self::Path(_) => None,
            Self::Quotient(algebra) => Some(algebra),
        }
    }

    /// Return the additive zero in the observable algebra variant.
    #[must_use]
    pub fn zero(&self) -> QuiverAlgebraElement<F> {
        match self {
            Self::Path(algebra) => QuiverAlgebraElement::Path(algebra.zero()),
            Self::Quotient(algebra) => QuiverAlgebraElement::Quotient(algebra.zero()),
        }
    }

    /// Return QPA's identity element in the observable algebra variant.
    ///
    /// # Errors
    ///
    /// Propagates quotient normal-form failures.
    pub fn one(&self) -> Result<QuiverAlgebraElement<F>, QuiverAlgebraElementError> {
        match self {
            Self::Path(algebra) => Ok(QuiverAlgebraElement::Path(algebra.one())),
            Self::Quotient(algebra) => Ok(QuiverAlgebraElement::Quotient(algebra.one()?)),
        }
    }

    /// Return QPA's ordinary `GeneratorsOfAlgebra` list.
    ///
    /// Path algebras return vertices followed by arrows. QPA quotients are
    /// algebras-with-one built from that same list, so GAP's ordinary algebra
    /// generators prepend the identity: one, vertices, then arrows.
    ///
    /// # Errors
    ///
    /// Propagates quotient normal-form failures while embedding generators.
    pub fn generators(&self) -> Result<Vec<QuiverAlgebraElement<F>>, QuiverAlgebraElementError> {
        match self {
            Self::Path(algebra) => Ok(algebra
                .generators()
                .into_iter()
                .map(QuiverAlgebraElement::Path)
                .collect()),
            Self::Quotient(algebra) => {
                let mut generators = Vec::with_capacity(
                    algebra
                        .original_path_algebra()
                        .quiver()
                        .number_of_vertices()
                        + algebra.original_path_algebra().quiver().number_of_arrows()
                        + 1,
                );
                generators.push(QuiverAlgebraElement::Quotient(algebra.one()?));
                generators.extend(
                    algebra
                        .generators()?
                        .into_iter()
                        .map(QuiverAlgebraElement::Quotient),
                );
                Ok(generators)
            }
        }
    }

    /// Return QPA's `GeneratorsOfAlgebraWithOne` list.
    ///
    /// For both variants this is the vertex-then-arrow list; unlike
    /// [`Self::generators`], the quotient case does not prepend the identity.
    ///
    /// # Errors
    ///
    /// Propagates quotient normal-form failures while embedding generators.
    pub fn generators_with_one(
        &self,
    ) -> Result<Vec<QuiverAlgebraElement<F>>, QuiverAlgebraElementError> {
        match self {
            Self::Path(algebra) => Ok(algebra
                .generators()
                .into_iter()
                .map(QuiverAlgebraElement::Path)
                .collect()),
            Self::Quotient(algebra) => Ok(algebra
                .generators()?
                .into_iter()
                .map(QuiverAlgebraElement::Quotient)
                .collect()),
        }
    }

    /// Return whether an element belongs to this identical algebra family.
    #[must_use]
    pub fn contains_element(&self, element: &QuiverAlgebraElement<F>) -> bool {
        match (self, element) {
            (Self::Path(algebra), QuiverAlgebraElement::Path(element)) => {
                algebra.is_identical(element.algebra())
            }
            (Self::Quotient(algebra), QuiverAlgebraElement::Quotient(element)) => {
                algebra.is_identical(element.algebra())
            }
            _ => false,
        }
    }

    /// Return QPA's `IsEnvelopingAlgebra` property.
    ///
    /// The generic source method is false; only the quotient object explicitly
    /// marked by `EnvelopingAlgebra` returns true.
    #[must_use]
    pub fn is_enveloping_algebra(&self) -> bool {
        match self {
            Self::Path(_) => false,
            Self::Quotient(algebra) => algebra.is_enveloping_algebra(),
        }
    }

    /// Return an already-live cached `EnvelopingAlgebra` attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_enveloping_algebra(&self) -> Option<TensorProductAlgebra<F>> {
        match self {
            Self::Path(algebra) => algebra.cached_enveloping_algebra(),
            Self::Quotient(algebra) => algebra.cached_enveloping_algebra(),
        }
    }

    /// Install and return QPA's cached `EnvelopingAlgebra` attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_enveloping_algebra(
        &self,
        enveloping: &TensorProductAlgebra<F>,
    ) -> TensorProductAlgebra<F> {
        match self {
            Self::Path(algebra) => algebra.cache_enveloping_algebra(enveloping),
            Self::Quotient(algebra) => algebra.cache_enveloping_algebra(enveloping),
        }
    }

    /// Return an already-live cached
    /// `AlgebraAsModuleOverEnvelopingAlgebra` attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_regular_enveloping_module(&self) -> Option<PathAlgebraModule<F>> {
        match self {
            Self::Path(algebra) => algebra.cached_regular_enveloping_module(),
            Self::Quotient(algebra) => algebra.cached_regular_enveloping_module(),
        }
    }

    /// Return QPA's `HasAlgebraAsModuleOverEnvelopingAlgebra` state.
    #[doc(hidden)]
    #[must_use]
    pub fn has_regular_enveloping_module_attribute(&self) -> bool {
        match self {
            Self::Path(algebra) => algebra.has_regular_enveloping_module_attribute(),
            Self::Quotient(algebra) => algebra.has_regular_enveloping_module_attribute(),
        }
    }

    /// Return whether the bound attribute value is GAP's `fail`.
    #[doc(hidden)]
    #[must_use]
    pub fn regular_enveloping_module_attribute_is_fail(&self) -> bool {
        match self {
            Self::Path(algebra) => algebra.regular_enveloping_module_attribute_is_fail(),
            Self::Quotient(algebra) => algebra.regular_enveloping_module_attribute_is_fail(),
        }
    }

    /// Bind the regular enveloping-module attribute to GAP's `fail` value.
    #[doc(hidden)]
    pub fn mark_regular_enveloping_module_attribute_fail(&self) {
        match self {
            Self::Path(algebra) => algebra.mark_regular_enveloping_module_attribute_fail(),
            Self::Quotient(algebra) => algebra.mark_regular_enveloping_module_attribute_fail(),
        }
    }

    /// Install and return QPA's cached regular enveloping-module attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_regular_enveloping_module(
        &self,
        module: &PathAlgebraModule<F>,
    ) -> PathAlgebraModule<F> {
        match self {
            Self::Path(algebra) => algebra.cache_regular_enveloping_module(module),
            Self::Quotient(algebra) => algebra.cache_regular_enveloping_module(module),
        }
    }

    /// Look up a vertex or arrow using QPA's record-component name semantics
    /// and return its image in this algebra.
    ///
    /// Repeated names resolve to the last arrow, then the last vertex, exactly
    /// as for the source quiver record.
    ///
    /// # Errors
    ///
    /// Propagates path embedding and quotient normal-form failures.
    pub fn generator_element_by_name(
        &self,
        name: &str,
    ) -> Result<Option<QuiverAlgebraElement<F>>, QuiverAlgebraElementError> {
        let Some(generator) = self.quiver().generator_by_name(name) else {
            return Ok(None);
        };
        Ok(Some(match self {
            Self::Path(algebra) => QuiverAlgebraElement::Path(algebra.element(generator.path())?),
            Self::Quotient(algebra) => {
                QuiverAlgebraElement::Quotient(algebra.path(generator.path())?)
            }
        }))
    }
}

impl<F: Field> From<PathAlgebra<F>> for QuiverAlgebra<F> {
    fn from(algebra: PathAlgebra<F>) -> Self {
        Self::Path(algebra)
    }
}

impl<F: Field> From<QuotientPathAlgebra<F>> for QuiverAlgebra<F> {
    fn from(algebra: QuotientPathAlgebra<F>) -> Self {
        Self::Quotient(algebra)
    }
}

/// An element of QPA's `IsQuiverAlgebra` union, retaining whether its algebra
/// is a path algebra or a proper quotient.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum QuiverAlgebraElement<F: Field> {
    /// An element of a path algebra.
    Path(PathAlgebraElement<F>),
    /// An element of a quotient of a path algebra.
    Quotient(QuotientPathAlgebraElement<F>),
}

impl<F: Field> QuiverAlgebraElement<F> {
    /// Return the exact coefficient field.
    #[must_use]
    pub fn field(&self) -> &F {
        match self {
            Self::Path(element) => element.algebra().field(),
            Self::Quotient(element) => element.algebra().original_path_algebra().field(),
        }
    }

    /// Return whether this is the additive zero.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        match self {
            Self::Path(element) => element.is_zero(),
            Self::Quotient(element) => element.is_zero(),
        }
    }

    /// Return QPA's leading term in the same algebra variant.
    #[must_use]
    pub fn leading_term(&self) -> Self {
        match self {
            Self::Path(element) => Self::Path(element.leading_term()),
            Self::Quotient(element) => Self::Quotient(element.leading_term()),
        }
    }

    /// Return QPA's leading coefficient, including field zero for zero.
    #[must_use]
    pub fn leading_coefficient(&self) -> F::Element {
        match self {
            Self::Path(element) => element.leading_coefficient(),
            Self::Quotient(element) => element.leading_coefficient(),
        }
    }

    /// Return QPA's leading monomial, including the zero path for zero.
    #[must_use]
    pub fn leading_monomial(&self) -> Path {
        match self {
            Self::Path(element) => element.leading_monomial(),
            Self::Quotient(element) => element.leading_monomial(),
        }
    }

    /// Add two elements in one identical quiver-algebra family.
    ///
    /// # Errors
    ///
    /// Reports mixed variants, different algebra families, and quotient
    /// normal-form failures.
    pub fn add(&self, right: &Self) -> Result<Self, QuiverAlgebraElementError> {
        match (self, right) {
            (Self::Path(left), Self::Path(right)) => Ok(Self::Path(left.add(right)?)),
            (Self::Quotient(left), Self::Quotient(right)) => Ok(Self::Quotient(left.add(right)?)),
            _ => Err(QuiverAlgebraElementError::DifferentVariants),
        }
    }

    /// Subtract two elements in one identical quiver-algebra family.
    ///
    /// # Errors
    ///
    /// Reports mixed variants, different algebra families, and quotient
    /// normal-form failures.
    pub fn sub(&self, right: &Self) -> Result<Self, QuiverAlgebraElementError> {
        match (self, right) {
            (Self::Path(left), Self::Path(right)) => Ok(Self::Path(left.sub(right)?)),
            (Self::Quotient(left), Self::Quotient(right)) => Ok(Self::Quotient(left.sub(right)?)),
            _ => Err(QuiverAlgebraElementError::DifferentVariants),
        }
    }

    /// Multiply two elements in one identical quiver-algebra family.
    ///
    /// # Errors
    ///
    /// Reports mixed variants, different algebra families, and quotient
    /// normal-form failures.
    pub fn multiply(&self, right: &Self) -> Result<Self, QuiverAlgebraElementError> {
        match (self, right) {
            (Self::Path(left), Self::Path(right)) => Ok(Self::Path(left.multiply(right)?)),
            (Self::Quotient(left), Self::Quotient(right)) => {
                Ok(Self::Quotient(left.multiply(right)?))
            }
            _ => Err(QuiverAlgebraElementError::DifferentVariants),
        }
    }

    /// Multiply by a field scalar on the left.
    ///
    /// # Errors
    ///
    /// Propagates path- or quotient-algebra arithmetic failures.
    pub fn scale_left(&self, scalar: &F::Element) -> Result<Self, QuiverAlgebraElementError> {
        match self {
            Self::Path(element) => Ok(Self::Path(element.scale_left(scalar)?)),
            Self::Quotient(element) => Ok(Self::Quotient(element.scale_left(scalar)?)),
        }
    }
}

/// Failure in arithmetic over QPA's path-algebra/quotient union.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum QuiverAlgebraElementError {
    /// GAP has no common element family for a path-algebra element and a
    /// quotient element.
    #[error("the elements belong to different quiver-algebra variants")]
    DifferentVariants,
    /// GAP's generic `Product` has no identity for an empty element list in
    /// this method family.
    #[error("an empty quiver-algebra element product is undefined")]
    EmptyProduct,
    /// GAP's generic `Sum` has no additive identity supplied for an empty
    /// element list in this method family.
    #[error("an empty quiver-algebra element sum is undefined")]
    EmptySum,
    /// Path-algebra arithmetic failed.
    #[error(transparent)]
    PathAlgebra(#[from] PathAlgebraError),
    /// Quotient arithmetic or normal-form reduction failed.
    #[error(transparent)]
    Quotient(#[from] QuotientError),
}

/// The rebuilt named quiver and converted relation list immediately before
/// QPA applies its configured Gröbner backend in `QuiverAlgebraOfAmodAeA`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct VertexQuotientPresentation<F: Field> {
    path_algebra: PathAlgebra<F>,
    relations: Vec<PathAlgebraElement<F>>,
}

impl<F: Field> VertexQuotientPresentation<F> {
    /// Return the rebuilt path algebra.
    #[must_use]
    pub const fn path_algebra(&self) -> &PathAlgebra<F> {
        &self.path_algebra
    }

    /// Return relations after `PathRemoval` and arrow-name conversion.
    #[must_use]
    pub fn relations(&self) -> &[PathAlgebraElement<F>] {
        &self.relations
    }

    /// Consume the presentation into its parent and relation list.
    #[must_use]
    pub fn into_parts(self) -> (PathAlgebra<F>, Vec<PathAlgebraElement<F>>) {
        (self.path_algebra, self.relations)
    }
}

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
    /// Rebuild QPA's presentation for `A/AeA` before Gröbner computation.
    ///
    /// A path algebra has no defining relations, so this always returns an
    /// empty converted relation list on the full subquiver complementary to
    /// `vertex_positions`.
    ///
    /// # Errors
    ///
    /// Propagates named-quiver and path-algebra construction failures.
    pub fn a_mod_aea_presentation(
        &self,
        vertex_positions: &[usize],
    ) -> Result<VertexQuotientPresentation<F>, PathAlgebraError> {
        vertex_quotient_presentation(self, &[], vertex_positions)
    }

    /// Translate `QuiverAlgebraOfAmodAeA` through QPA's high-level backend.
    ///
    /// The empty-relation result remains a path algebra, matching QPA's
    /// dispatch for `KQ / []`.
    ///
    /// # Errors
    ///
    /// Propagates presentation and high-level Gröbner failures.
    pub fn quiver_algebra_of_a_mod_aea_high_level(
        &self,
        vertex_positions: &[usize],
    ) -> Result<QuiverAlgebra<F>, QuotientError> {
        quotient_vertex_presentation_high_level(self.a_mod_aea_presentation(vertex_positions)?)
    }

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
    /// Rebuild QPA's presentation for `A/AeA` before Gröbner computation.
    ///
    /// This uses `RelationsOfAlgebra(A)`, applies `PathRemoval` to each stored
    /// relation, removes zeros, and converts every surviving walk by matching
    /// arrow names in the rebuilt full subquiver.
    ///
    /// # Errors
    ///
    /// Propagates `PathRemoval`, quiver reconstruction, and sparse arithmetic
    /// failures.
    pub fn a_mod_aea_presentation(
        &self,
        vertex_positions: &[usize],
    ) -> Result<VertexQuotientPresentation<F>, PathAlgebraError> {
        vertex_quotient_presentation(
            self.original_path_algebra(),
            self.relations(),
            vertex_positions,
        )
    }

    /// Translate `QuiverAlgebraOfAmodAeA` through QPA's high-level backend.
    ///
    /// # Errors
    ///
    /// Propagates presentation and high-level Gröbner failures.
    pub fn quiver_algebra_of_a_mod_aea_high_level(
        &self,
        vertex_positions: &[usize],
    ) -> Result<QuiverAlgebra<F>, QuotientError> {
        quotient_vertex_presentation_high_level(self.a_mod_aea_presentation(vertex_positions)?)
    }

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

fn quotient_vertex_presentation_high_level<F: Field>(
    presentation: VertexQuotientPresentation<F>,
) -> Result<QuiverAlgebra<F>, QuotientError> {
    let (path_algebra, relations) = presentation.into_parts();
    if relations.is_empty() {
        Ok(QuiverAlgebra::Path(path_algebra))
    } else {
        Ok(QuiverAlgebra::Quotient(
            QuotientPathAlgebra::by_relations_high_level(path_algebra, &relations)?,
        ))
    }
}

fn vertex_quotient_presentation<F: Field>(
    parent: &PathAlgebra<F>,
    relations: &[PathAlgebraElement<F>],
    removed_vertex_positions: &[usize],
) -> Result<VertexQuotientPresentation<F>, PathAlgebraError> {
    let retained_vertices = parent
        .quiver()
        .vertices()
        .filter(|vertex| !removed_vertex_positions.contains(&vertex.index().get()))
        .collect::<Vec<_>>();
    let vertex_names = retained_vertices
        .iter()
        .map(|vertex| Box::<str>::from(vertex.name()))
        .collect::<Vec<_>>();
    let retained_arrows = parent
        .quiver()
        .arrows()
        .filter(|arrow| {
            retained_vertices.contains(&arrow.source())
                && retained_vertices.contains(&arrow.target())
        })
        .collect::<Vec<_>>();
    let arrow_specs = retained_arrows
        .iter()
        .map(|arrow| {
            NamedArrowSpec::named(arrow.source().name(), arrow.target().name(), arrow.name())
        })
        .collect();
    let quiver = Quiver::from_named(vertex_names, arrow_specs)?;
    let path_algebra = PathAlgebra::new(parent.field().clone(), quiver)?;

    let mut converted_relations = Vec::new();
    for relation in relations {
        let relation = relation.path_removal(removed_vertex_positions)?;
        if relation.is_zero() {
            continue;
        }
        let mut converted = path_algebra.zero();
        for term in relation.terms() {
            let walk = term
                .path()
                .walk()
                .ok_or(PathAlgebraError::InvalidPathTermInPathRemoval)?;
            let mut mapped_arrows = Vec::with_capacity(walk.len());
            for arrow_index in walk {
                let original_arrow = parent
                    .quiver()
                    .arrow(*arrow_index)
                    .ok_or(PathAlgebraError::InvalidPathTermInPathRemoval)?;
                let mapped_arrow = path_algebra
                    .quiver()
                    .arrows()
                    .find(|arrow| arrow.name() == original_arrow.name())
                    .ok_or_else(|| PathAlgebraError::MissingArrowInVertexQuotient {
                        name: original_arrow.name().into(),
                    })?;
                mapped_arrows.push(mapped_arrow);
            }
            let mut mapped_arrows = mapped_arrows.into_iter();
            let Some(first) = mapped_arrows.next() else {
                return Err(PathAlgebraError::VertexTermInPathRemoval);
            };
            let mut path = first.path();
            for arrow in mapped_arrows {
                path = path.multiply(&arrow.path())?;
            }
            let monomial = path_algebra.one().multiply(&path_algebra.element(path)?)?;
            converted = converted.add(&monomial.scale_left(term.coefficient())?)?;
        }
        converted_relations.push(converted);
    }
    Ok(VertexQuotientPresentation {
        path_algebra,
        relations: converted_relations,
    })
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
