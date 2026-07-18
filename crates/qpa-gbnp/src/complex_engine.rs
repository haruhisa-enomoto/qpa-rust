//! The category, complex, and chain-map machinery from QPA's `complex.gd` and
//! `complex.gi`.
//!
//! QPA represents a complex by an [`InfiniteList`] of differentials.  This
//! module keeps that representation: lazy position functions, automatically
//! detected repeating tails, and their callbacks are delegated directly to
//! the source-translated infinite-list layer in `qpa-core`.

// This module deliberately exposes the callback-rich, fallible API and owned
// lazy-list specifications used by QPA.  Flattening those signatures merely to
// satisfy style lints would obscure the correspondence with `complex.gd/.gi`.
#![allow(
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::needless_pass_by_value,
    clippy::similar_names,
    clippy::single_match_else,
    clippy::too_many_lines,
    clippy::type_complexity,
    clippy::unnecessary_literal_bound
)]

use std::{
    fmt::{self, Write as _},
    sync::{Arc, OnceLock, Weak},
};

use numerica::domains::{Field, integer::Integer};
use qpa_core::{
    ExtendedInteger, HalfInfiniteListSpecification, InfiniteList, InfiniteListDirection,
    InfiniteListError, PathAlgebraModule, PathAlgebraModuleError, PathAlgebraModuleHomomorphism,
    QpaCoefficientField, QuiverAlgebra, cokernel_projection, direct_sum_inclusions,
    direct_sum_of_qpa_modules, direct_sum_projections, image_inclusion,
    image_of_module_homomorphism, image_projection, image_projection_inclusion, kernel_inclusion,
    kernel_of_module_homomorphism, lifting_inclusion_morphisms, multiply_lists_of_maps,
    zero_module, zero_module_mapping,
};
use thiserror::Error;

use crate::module::{
    EnvelopingModuleError, injective_envelope, lifting_morphism_from_projective, projective_cover,
};

/// A category carrying precisely the operations stored by QPA's `Cat`
/// record, together with the ordinary source/range and additive-inverse
/// operations used by the installed complex methods.
pub trait QpaCategory: Clone + Eq + fmt::Debug + Send + Sync + 'static {
    /// Objects of the category.
    type Object: Clone + Eq + fmt::Debug + Send + Sync + 'static;
    /// Morphisms of the category.
    type Morphism: Clone + Eq + fmt::Debug + Send + Sync + 'static;

    /// QPA's printable category name.
    fn name(&self) -> &str;
    /// QPA's `zeroObj` property.
    fn zero_object(&self) -> Self::Object;
    /// QPA's `isZeroObj` property.
    fn is_zero_object(&self, object: &Self::Object) -> bool;
    /// QPA's `zeroMap` property.
    fn zero_map(
        &self,
        source: &Self::Object,
        range: &Self::Object,
    ) -> Result<Self::Morphism, ComplexError>;
    /// QPA's `isZeroMapping` property.
    fn is_zero_mapping(&self, morphism: &Self::Morphism) -> bool;
    /// QPA's `composeMaps(g,f)`, meaning first `f` and then `g`.
    fn compose_maps(
        &self,
        g: &Self::Morphism,
        f: &Self::Morphism,
    ) -> Result<Self::Morphism, ComplexError>;
    /// Source of a morphism.
    fn source(&self, morphism: &Self::Morphism) -> Self::Object;
    /// Range of a morphism.
    fn range(&self, morphism: &Self::Morphism) -> Self::Object;
    /// Additive inverse of a morphism.
    fn negate(&self, morphism: &Self::Morphism) -> Result<Self::Morphism, ComplexError>;
    /// QPA's `ker` property.
    fn kernel(&self, morphism: &Self::Morphism) -> Result<Self::Object, ComplexError>;
    /// QPA's `im` property.
    fn image(&self, morphism: &Self::Morphism) -> Result<Self::Object, ComplexError>;
    /// QPA's `isExact(g,f)` property.
    fn is_exact(&self, g: &Self::Morphism, f: &Self::Morphism) -> Result<bool, ComplexError>;
    /// QPA's `objStr` property.
    fn object_string(&self, object: &Self::Object) -> String;
}

/// Additional module-category operations used by QPA's homology and
/// truncation methods.  QPA installs those methods on every complex but calls
/// module-specific global operations internally; the separate trait preserves
/// that actual dispatch boundary.
pub trait QpaAbelianCategory: QpaCategory {
    /// `KernelInclusion`.
    fn kernel_inclusion(&self, morphism: &Self::Morphism) -> Result<Self::Morphism, ComplexError>;
    /// `ImageProjection`.
    fn image_projection(&self, morphism: &Self::Morphism) -> Result<Self::Morphism, ComplexError>;
    /// `ImageInclusion`.
    fn image_inclusion(&self, morphism: &Self::Morphism) -> Result<Self::Morphism, ComplexError>;
    /// `CoKernelProjection`.
    fn cokernel_projection(
        &self,
        morphism: &Self::Morphism,
    ) -> Result<Self::Morphism, ComplexError>;
    /// `LiftingInclusionMorphisms`, where `None` is GAP's `fail`.
    fn lifting_inclusion_morphisms(
        &self,
        first: &Self::Morphism,
        second: &Self::Morphism,
    ) -> Result<Option<Self::Morphism>, ComplexError>;
}

/// Callback collection accepted by [`CallbackCategory::new`], the typed Rust
/// counterpart of QPA's `properties` record passed to `Cat`.
pub struct CallbackCategoryProperties<O, M> {
    /// Printable `name` field.
    pub name: Arc<str>,
    /// `zeroObj` field.
    pub zero_object: O,
    /// `isZeroObj` field.
    pub is_zero_object: Arc<dyn Fn(&O) -> bool + Send + Sync>,
    /// `zeroMap` field.
    pub zero_map: Arc<dyn Fn(&O, &O) -> Result<M, ComplexError> + Send + Sync>,
    /// `isZeroMapping` field.
    pub is_zero_mapping: Arc<dyn Fn(&M) -> bool + Send + Sync>,
    /// `composeMaps` field.
    pub compose_maps: Arc<dyn Fn(&M, &M) -> Result<M, ComplexError> + Send + Sync>,
    /// Source accessor used by the installed methods.
    pub source: Arc<dyn Fn(&M) -> O + Send + Sync>,
    /// Range accessor used by the installed methods.
    pub range: Arc<dyn Fn(&M) -> O + Send + Sync>,
    /// Additive inverse used by signed shifts.
    pub negate: Arc<dyn Fn(&M) -> Result<M, ComplexError> + Send + Sync>,
    /// `ker` field.
    pub kernel: Arc<dyn Fn(&M) -> Result<O, ComplexError> + Send + Sync>,
    /// `im` field.
    pub image: Arc<dyn Fn(&M) -> Result<O, ComplexError> + Send + Sync>,
    /// `isExact` field.
    pub is_exact: Arc<dyn Fn(&M, &M) -> Result<bool, ComplexError> + Send + Sync>,
    /// `objStr` field.
    pub object_string: Arc<dyn Fn(&O) -> String + Send + Sync>,
}

impl<O: Clone, M> Clone for CallbackCategoryProperties<O, M> {
    fn clone(&self) -> Self {
        Self {
            name: Arc::clone(&self.name),
            zero_object: self.zero_object.clone(),
            is_zero_object: Arc::clone(&self.is_zero_object),
            zero_map: Arc::clone(&self.zero_map),
            is_zero_mapping: Arc::clone(&self.is_zero_mapping),
            compose_maps: Arc::clone(&self.compose_maps),
            source: Arc::clone(&self.source),
            range: Arc::clone(&self.range),
            negate: Arc::clone(&self.negate),
            kernel: Arc::clone(&self.kernel),
            image: Arc::clone(&self.image),
            is_exact: Arc::clone(&self.is_exact),
            object_string: Arc::clone(&self.object_string),
        }
    }
}

/// A runtime category constructed from an identity and callback record, as in
/// QPA's `Cat(identity, properties)`.
pub struct CallbackCategory<I, O, M> {
    identity: I,
    properties: CallbackCategoryProperties<O, M>,
}

impl<I: Clone, O: Clone, M> Clone for CallbackCategory<I, O, M> {
    fn clone(&self) -> Self {
        Self {
            identity: self.identity.clone(),
            properties: self.properties.clone(),
        }
    }
}

impl<I: fmt::Debug, O, M> fmt::Debug for CallbackCategory<I, O, M> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("CallbackCategory")
            .field("identity", &self.identity)
            .field("name", &self.properties.name)
            .finish_non_exhaustive()
    }
}

impl<I: PartialEq, O, M> PartialEq for CallbackCategory<I, O, M> {
    fn eq(&self, other: &Self) -> bool {
        self.identity == other.identity
    }
}

impl<I: Eq, O, M> Eq for CallbackCategory<I, O, M> {}

impl<I, O, M> CallbackCategory<I, O, M> {
    /// Translate QPA's `Cat(identity, properties)` constructor.
    #[must_use]
    pub const fn new(identity: I, properties: CallbackCategoryProperties<O, M>) -> Self {
        Self {
            identity,
            properties,
        }
    }

    /// Return QPA's category identity payload.
    #[must_use]
    pub const fn identity(&self) -> &I {
        &self.identity
    }
}

impl<I, O, M> QpaCategory for CallbackCategory<I, O, M>
where
    I: Clone + Eq + fmt::Debug + Send + Sync + 'static,
    O: Clone + Eq + fmt::Debug + Send + Sync + 'static,
    M: Clone + Eq + fmt::Debug + Send + Sync + 'static,
{
    type Object = O;
    type Morphism = M;

    fn name(&self) -> &str {
        &self.properties.name
    }

    fn zero_object(&self) -> O {
        self.properties.zero_object.clone()
    }

    fn is_zero_object(&self, object: &O) -> bool {
        (self.properties.is_zero_object)(object)
    }

    fn zero_map(&self, source: &O, range: &O) -> Result<M, ComplexError> {
        (self.properties.zero_map)(source, range)
    }

    fn is_zero_mapping(&self, morphism: &M) -> bool {
        (self.properties.is_zero_mapping)(morphism)
    }

    fn compose_maps(&self, g: &M, f: &M) -> Result<M, ComplexError> {
        (self.properties.compose_maps)(g, f)
    }

    fn source(&self, morphism: &M) -> O {
        (self.properties.source)(morphism)
    }

    fn range(&self, morphism: &M) -> O {
        (self.properties.range)(morphism)
    }

    fn negate(&self, morphism: &M) -> Result<M, ComplexError> {
        (self.properties.negate)(morphism)
    }

    fn kernel(&self, morphism: &M) -> Result<O, ComplexError> {
        (self.properties.kernel)(morphism)
    }

    fn image(&self, morphism: &M) -> Result<O, ComplexError> {
        (self.properties.image)(morphism)
    }

    fn is_exact(&self, g: &M, f: &M) -> Result<bool, ComplexError> {
        (self.properties.is_exact)(g, f)
    }

    fn object_string(&self, object: &O) -> String {
        (self.properties.object_string)(object)
    }
}

/// QPA's category of right modules over one path algebra or quotient.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RightModuleCategory<F: Field> {
    algebra: QuiverAlgebra<F>,
    zero: PathAlgebraModule<F>,
}

/// Translate `CatOfRightAlgebraModules(A)`.
///
/// # Errors
///
/// Propagates construction of QPA's zero module.
pub fn cat_of_right_algebra_modules<F>(
    algebra: &QuiverAlgebra<F>,
) -> Result<RightModuleCategory<F>, ComplexError>
where
    F: QpaCoefficientField + Send + Sync,
    F::Element: Send + Sync,
{
    Ok(RightModuleCategory {
        algebra: algebra.clone(),
        zero: zero_module(algebra)?,
    })
}

impl<F> RightModuleCategory<F>
where
    F: QpaCoefficientField + Send + Sync,
    F::Element: Send + Sync,
{
    /// Return the algebra retained in QPA's category identity.
    #[must_use]
    pub const fn algebra(&self) -> &QuiverAlgebra<F> {
        &self.algebra
    }
}

impl<F> QpaCategory for RightModuleCategory<F>
where
    F: QpaCoefficientField + Send + Sync + 'static,
    F::Element: Send + Sync,
{
    type Object = PathAlgebraModule<F>;
    type Morphism = PathAlgebraModuleHomomorphism<F>;

    fn name(&self) -> &str {
        "right modules over algebra"
    }

    fn zero_object(&self) -> Self::Object {
        self.zero.clone()
    }

    fn is_zero_object(&self, object: &Self::Object) -> bool {
        object.dimension() == 0
    }

    fn zero_map(
        &self,
        source: &Self::Object,
        range: &Self::Object,
    ) -> Result<Self::Morphism, ComplexError> {
        Ok(zero_module_mapping(source.clone(), range.clone())?)
    }

    fn is_zero_mapping(&self, morphism: &Self::Morphism) -> bool {
        morphism.is_zero()
    }

    fn compose_maps(
        &self,
        g: &Self::Morphism,
        f: &Self::Morphism,
    ) -> Result<Self::Morphism, ComplexError> {
        Ok(f.compose(g)?)
    }

    fn source(&self, morphism: &Self::Morphism) -> Self::Object {
        morphism.source().clone()
    }

    fn range(&self, morphism: &Self::Morphism) -> Self::Object {
        morphism.range().clone()
    }

    fn negate(&self, morphism: &Self::Morphism) -> Result<Self::Morphism, ComplexError> {
        Ok(morphism.negate()?)
    }

    fn kernel(&self, morphism: &Self::Morphism) -> Result<Self::Object, ComplexError> {
        Ok(kernel_of_module_homomorphism(morphism)?)
    }

    fn image(&self, morphism: &Self::Morphism) -> Result<Self::Object, ComplexError> {
        Ok(image_of_module_homomorphism(morphism)?)
    }

    fn is_exact(&self, g: &Self::Morphism, f: &Self::Morphism) -> Result<bool, ComplexError> {
        Ok(kernel_of_module_homomorphism(g)?.dimension()
            == image_of_module_homomorphism(f)?.dimension())
    }

    fn object_string(&self, object: &Self::Object) -> String {
        format!(
            "({})",
            object
                .dimension_vector()
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(",")
        )
    }
}

impl<F> QpaAbelianCategory for RightModuleCategory<F>
where
    F: QpaCoefficientField + Send + Sync + 'static,
    F::Element: Send + Sync,
{
    fn kernel_inclusion(&self, morphism: &Self::Morphism) -> Result<Self::Morphism, ComplexError> {
        Ok(kernel_inclusion(morphism)?)
    }

    fn image_projection(&self, morphism: &Self::Morphism) -> Result<Self::Morphism, ComplexError> {
        Ok(image_projection(morphism)?)
    }

    fn image_inclusion(&self, morphism: &Self::Morphism) -> Result<Self::Morphism, ComplexError> {
        Ok(image_inclusion(morphism)?)
    }

    fn cokernel_projection(
        &self,
        morphism: &Self::Morphism,
    ) -> Result<Self::Morphism, ComplexError> {
        Ok(cokernel_projection(morphism)?)
    }

    fn lifting_inclusion_morphisms(
        &self,
        first: &Self::Morphism,
        second: &Self::Morphism,
    ) -> Result<Option<Self::Morphism>, ComplexError> {
        Ok(lifting_inclusion_morphisms(first, second)?)
    }
}

/// QPA's five accepted descriptions of one side of a complex.
pub enum ComplexHalf<C: QpaCategory> {
    /// QPA's string `"zero"`.
    Zero,
    /// A repeating differential list.
    Repeat(Vec<C::Morphism>),
    /// A predecessor-generated half-list (`"next"` or `"next/repeat"`).
    Next {
        /// QPA's next-differential function.
        function:
            Arc<dyn Fn(&C::Morphism) -> Result<C::Morphism, ComplexError> + Send + Sync + 'static>,
        /// Value supplied before the first generated differential.
        initial_value: C::Morphism,
        /// Whether QPA's collision-based repeat conversion is enabled.
        repeatify: bool,
    },
    /// A degree-generated half-list (`"pos"`).
    Position {
        /// Function receiving the complex and requested degree.
        function: Arc<
            dyn Fn(&QpaComplex<C>, &Integer) -> Result<C::Morphism, ComplexError>
                + Send
                + Sync
                + 'static,
        >,
        /// Whether consecutively generated values are retained.
        store_values: bool,
    },
}

impl<C: QpaCategory> Clone for ComplexHalf<C> {
    fn clone(&self) -> Self {
        match self {
            Self::Zero => Self::Zero,
            Self::Repeat(values) => Self::Repeat(values.clone()),
            Self::Next {
                function,
                initial_value,
                repeatify,
            } => Self::Next {
                function: Arc::clone(function),
                initial_value: initial_value.clone(),
                repeatify: *repeatify,
            },
            Self::Position {
                function,
                store_values,
            } => Self::Position {
                function: Arc::clone(function),
                store_values: *store_values,
            },
        }
    }
}

struct QpaComplexData<C: QpaCategory> {
    category: C,
    differentials: OnceLock<InfiniteList<C::Morphism, ComplexError>>,
    zero_complex: bool,
    projective_complex: OnceLock<bool>,
    injective_complex: OnceLock<bool>,
    // Cuts, shifts, and splices can share lazy half-list generators with
    // source complexes. Keep those sources alive while this complex can ask
    // a shared generator for another differential.
    _dependencies: Vec<Arc<QpaComplexData<C>>>,
}

/// QPA's `IsQPAComplex` object.
pub struct QpaComplex<C: QpaCategory>(Arc<QpaComplexData<C>>);

impl<C: QpaCategory> Clone for QpaComplex<C> {
    fn clone(&self) -> Self {
        Self(Arc::clone(&self.0))
    }
}

impl<C: QpaCategory> fmt::Debug for QpaComplex<C> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("QpaComplex")
            .field("category", &self.0.category)
            .field("zero_complex", &self.0.zero_complex)
            .finish_non_exhaustive()
    }
}

impl<C: QpaCategory> PartialEq for QpaComplex<C> {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }
}

impl<C: QpaCategory> Eq for QpaComplex<C> {}

impl<C: QpaCategory> QpaComplex<C> {
    /// Return QPA's `CatOfComplex` attribute.
    #[must_use]
    pub fn category(&self) -> &C {
        &self.0.category
    }

    /// Return whether the object carries QPA's `IsZeroComplex` filter.
    #[must_use]
    pub fn is_zero_complex(&self) -> bool {
        self.0.zero_complex
    }

    /// Return QPA's cached `IsProjectiveComplex` property, if installed.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_is_projective_complex(&self) -> Option<bool> {
        self.0.projective_complex.get().copied()
    }

    /// Install QPA's `IsProjectiveComplex` property value.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_is_projective_complex(&self, value: bool) -> bool {
        *self.0.projective_complex.get_or_init(|| value)
    }

    /// Return QPA's cached `IsInjectiveComplex` property, if installed.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_is_injective_complex(&self) -> Option<bool> {
        self.0.injective_complex.get().copied()
    }

    /// Install QPA's `IsInjectiveComplex` property value.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_is_injective_complex(&self, value: bool) -> bool {
        *self.0.injective_complex.get_or_init(|| value)
    }

    /// Return QPA's `DifferentialsOfComplex` attribute.
    #[must_use]
    pub fn differentials(&self) -> InfiniteList<C::Morphism, ComplexError> {
        self.0
            .differentials
            .get()
            .expect("QPA complex differential list is initialized before exposure")
            .clone()
    }

    /// Translate `DifferentialOfComplex(C,i)` and `C^i`.
    ///
    /// # Errors
    ///
    /// Propagates lazy infinite-list generation and complex validation errors.
    pub fn differential(&self, degree: &Integer) -> Result<C::Morphism, ComplexError> {
        self.differentials().value(degree)
    }

    /// Translate `ObjectOfComplex(C,i)`.
    ///
    /// # Errors
    ///
    /// Propagates [`Self::differential`].
    pub fn object(&self, degree: &Integer) -> Result<C::Object, ComplexError> {
        Ok(self.category().source(&self.differential(degree)?))
    }

    /// Translate `CyclesOfComplex(C,i)`.
    pub fn cycles(&self, degree: &Integer) -> Result<C::Object, ComplexError> {
        self.category().kernel(&self.differential(degree)?)
    }

    /// Translate `BoundariesOfComplex(C,i)`.
    pub fn boundaries(&self, degree: &Integer) -> Result<C::Object, ComplexError> {
        self.category()
            .image(&self.differential(&(degree + Integer::from(1)))?)
    }

    /// Translate `HighestKnownDegree`.
    pub fn highest_known_degree(&self) -> Result<Option<ExtendedInteger>, ComplexError> {
        self.differentials().highest_known_position()
    }

    /// Translate `LowestKnownDegree`.
    pub fn lowest_known_degree(&self) -> Result<Option<ExtendedInteger>, ComplexError> {
        self.differentials().lowest_known_position()
    }

    /// Translate `IsPositiveRepeating`.
    pub fn is_positive_repeating(&self) -> Result<bool, ComplexError> {
        self.differentials()
            .positive_part()?
            .ok_or(ComplexError::MissingComplexPositivePart)?
            .is_repeating()
    }

    /// Translate `IsNegativeRepeating`.
    pub fn is_negative_repeating(&self) -> Result<bool, ComplexError> {
        self.differentials()
            .negative_part()?
            .ok_or(ComplexError::MissingComplexNegativePart)?
            .is_repeating()
    }

    /// Translate `PositiveRepeatDegrees`, returning `None` for GAP's `fail`.
    pub fn positive_repeat_degrees(&self) -> Result<Option<Vec<Integer>>, ComplexError> {
        let positive = self
            .differentials()
            .positive_part()?
            .ok_or(ComplexError::MissingComplexPositivePart)?;
        let Some(repeating) = positive.repeating_list()? else {
            return Ok(None);
        };
        Ok(Some(integer_interval(
            positive.start_position(),
            &(positive.start_position() + Integer::from(repeating.len()) - Integer::from(1)),
        )))
    }

    /// Translate `NegativeRepeatDegrees`, returning `None` for GAP's `fail`.
    pub fn negative_repeat_degrees(&self) -> Result<Option<Vec<Integer>>, ComplexError> {
        let negative = self
            .differentials()
            .negative_part()?
            .ok_or(ComplexError::MissingComplexNegativePart)?;
        let Some(repeating) = negative.repeating_list()? else {
            return Ok(None);
        };
        let last = negative.start_position() - Integer::from(repeating.len()) + Integer::from(1);
        Ok(Some(integer_interval(&last, negative.start_position())))
    }

    /// Translate `UpperBound`, returning `None` for GAP's `fail`.
    pub fn upper_bound(&self) -> Result<Option<ExtendedInteger>, ComplexError> {
        if self.is_zero_complex() {
            return Ok(Some(ExtendedInteger::NegativeInfinity));
        }
        let differentials = self.differentials();
        let positive = differentials
            .positive_part()?
            .ok_or(ComplexError::MissingComplexPositivePart)?;
        let Some(repeating) = positive.repeating_list()? else {
            return Ok(None);
        };
        if repeating.len() != 1
            || !self
                .category()
                .is_zero_object(&self.object(positive.start_position())?)
        {
            return Ok(Some(ExtendedInteger::PositiveInfinity));
        }
        let mut degree = differentials.middle_end()? - Integer::from(1);
        while self.category().is_zero_object(&self.object(&degree)?) {
            degree -= Integer::from(1);
            if degree < *differentials.middle_start() {
                return Ok(None);
            }
        }
        Ok(Some(ExtendedInteger::Finite(degree)))
    }

    /// Translate `LowerBound`, returning `None` for GAP's `fail`.
    pub fn lower_bound(&self) -> Result<Option<ExtendedInteger>, ComplexError> {
        if self.is_zero_complex() {
            return Ok(Some(ExtendedInteger::PositiveInfinity));
        }
        let differentials = self.differentials();
        let negative = differentials
            .negative_part()?
            .ok_or(ComplexError::MissingComplexNegativePart)?;
        let Some(repeating) = negative.repeating_list()? else {
            return Ok(None);
        };
        if repeating.len() != 1
            || !self
                .category()
                .is_zero_object(&self.object(negative.start_position())?)
        {
            return Ok(Some(ExtendedInteger::NegativeInfinity));
        }
        let mut degree = differentials.middle_start().clone();
        let middle_end = differentials.middle_end()?;
        while self.category().is_zero_object(&self.object(&degree)?) {
            degree += Integer::from(1);
            if degree > middle_end {
                return Ok(None);
            }
        }
        Ok(Some(ExtendedInteger::Finite(degree)))
    }

    /// Translate `IsFiniteComplex`, returning `None` for GAP's `fail`.
    pub fn is_finite_complex(&self) -> Result<Option<bool>, ComplexError> {
        let upper = self.upper_bound()?;
        let lower = self.lower_bound()?;
        if self.is_zero_complex()
            || matches!(upper, Some(ExtendedInteger::Finite(_)))
                && matches!(lower, Some(ExtendedInteger::Finite(_)))
        {
            Ok(Some(true))
        } else if matches!(
            upper,
            Some(ExtendedInteger::NegativeInfinity | ExtendedInteger::PositiveInfinity)
        ) || matches!(
            lower,
            Some(ExtendedInteger::NegativeInfinity | ExtendedInteger::PositiveInfinity)
        ) {
            Ok(Some(false))
        } else {
            Ok(None)
        }
    }

    /// Translate `LengthOfComplex`, returning `None` for GAP's `fail`.
    pub fn length(&self) -> Result<Option<ExtendedInteger>, ComplexError> {
        if self.is_zero_complex() {
            return Ok(Some(ExtendedInteger::Finite(Integer::from(0))));
        }
        match self.is_finite_complex()? {
            Some(true) => {
                let Some(ExtendedInteger::Finite(upper)) = self.upper_bound()? else {
                    return Err(ComplexError::InconsistentFiniteBounds);
                };
                let Some(ExtendedInteger::Finite(lower)) = self.lower_bound()? else {
                    return Err(ComplexError::InconsistentFiniteBounds);
                };
                Ok(Some(ExtendedInteger::Finite(
                    upper - lower + Integer::from(1),
                )))
            }
            Some(false) => Ok(Some(ExtendedInteger::PositiveInfinity)),
            None => Ok(None),
        }
    }

    /// Translate `IsExactInDegree`.
    pub fn is_exact_in_degree(&self, degree: &Integer) -> Result<bool, ComplexError> {
        self.category().is_exact(
            &self.differential(degree)?,
            &self.differential(&(degree + Integer::from(1)))?,
        )
    }

    /// Translate `ForEveryDegree`, returning `None` for GAP's `fail`.
    pub fn for_every_degree(
        &self,
        function: impl Fn(&C::Morphism, &C::Morphism) -> Result<bool, ComplexError>,
    ) -> Result<Option<bool>, ComplexError> {
        let differentials = self.differentials();
        let positive = differentials
            .positive_part()?
            .ok_or(ComplexError::MissingComplexPositivePart)?;
        let negative = differentials
            .negative_part()?
            .ok_or(ComplexError::MissingComplexNegativePart)?;
        let Some(positive_repeat) = positive.repeating_list()? else {
            return Ok(None);
        };
        let Some(negative_repeat) = negative.repeating_list()? else {
            return Ok(None);
        };
        let first = negative.start_position() - Integer::from(negative_repeat.len());
        let last = positive.start_position() + Integer::from(positive_repeat.len());
        let mut degree = first;
        while degree <= last {
            if !function(
                &self.differential(&degree)?,
                &self.differential(&(&degree + Integer::from(1)))?,
            )? {
                return Ok(Some(false));
            }
            degree += Integer::from(1);
        }
        Ok(Some(true))
    }

    /// Translate the `IsExactSequence` property.
    ///
    /// QPA's installed method attempts to return `fail` for an unknown
    /// nonrepeating complex. GAP rejects a non-boolean property result, so that
    /// observable branch is a typed error here.
    pub fn is_exact_sequence(&self) -> Result<bool, ComplexError> {
        self.for_every_degree(|g, f| self.category().is_exact(g, f))?
            .ok_or(ComplexError::PropertyReturnedNonBoolean {
                property: "IsExactSequence",
            })
    }

    /// Translate the `IsShortExactSequence` property, including GAP's error
    /// when its source method attempts to return `fail`.
    pub fn is_short_exact_sequence(&self) -> Result<bool, ComplexError> {
        let length = self
            .length()?
            .ok_or(ComplexError::PropertyReturnedNonBoolean {
                property: "IsShortExactSequence",
            })?;
        if length == ExtendedInteger::Finite(Integer::from(3)) {
            self.is_exact_sequence()
        } else {
            Ok(false)
        }
    }

    /// Translate signed `Shift(C,n)`.
    pub fn shift(&self, amount: &Integer) -> Result<Self, ComplexError> {
        let mut differentials = self.differentials().shift(amount)?;
        if (amount % Integer::from(2)) != 0 {
            let category = self.category().clone();
            differentials = differentials.map(move |map| category.negate(map))?;
        }
        Ok(Self::by_differential_list_with_dependencies(
            self.category().clone(),
            differentials,
            &[self],
        ))
    }

    /// Translate `ShiftUnsigned(C,n)`.
    pub fn shift_unsigned(&self, amount: &Integer) -> Result<Self, ComplexError> {
        Ok(Self::by_differential_list_with_dependencies(
            self.category().clone(),
            self.differentials().shift(amount)?,
            &[self],
        ))
    }

    /// Translate `ComplexByDifferentialList` without revalidating the list.
    #[must_use]
    pub fn by_differential_list(
        category: C,
        differentials: InfiniteList<C::Morphism, ComplexError>,
    ) -> Self {
        Self::by_differential_list_with_dependencies(category, differentials, &[])
    }

    fn by_differential_list_with_dependencies(
        category: C,
        differentials: InfiniteList<C::Morphism, ComplexError>,
        dependencies: &[&Self],
    ) -> Self {
        let data = Arc::new(QpaComplexData {
            category,
            differentials: OnceLock::new(),
            zero_complex: false,
            projective_complex: OnceLock::new(),
            injective_complex: OnceLock::new(),
            _dependencies: dependencies
                .iter()
                .map(|complex| Arc::clone(&complex.0))
                .collect(),
        });
        data.differentials
            .set(differentials)
            .expect("new QPA complex has an empty differential slot");
        Self(data)
    }

    pub(crate) fn with_differential_list(
        &self,
        category: C,
        differentials: InfiniteList<C::Morphism, ComplexError>,
    ) -> Self {
        Self::by_differential_list_with_dependencies(category, differentials, &[self])
    }

    /// Return QPA's printed complex representation.
    pub fn qpa_string(&self) -> Result<String, ComplexError> {
        if self.is_zero_complex() {
            return Ok("0 -> 0".to_owned());
        }
        let differentials = self.differentials();
        let upper = self.upper_bound()?;
        let lower = self.lower_bound()?;
        let top = if matches!(upper, Some(ExtendedInteger::Finite(_))) {
            Some(differentials.middle_end()? - Integer::from(1))
        } else if self.is_positive_repeating()? {
            Some(differentials.middle_end()?)
        } else {
            finite_known(self.highest_known_degree()?)
        };
        let bottom = if self.is_negative_repeating()? {
            Some(differentials.middle_start().clone())
        } else {
            finite_known(self.lowest_known_degree()?)
        };
        let mut output = if matches!(upper, Some(ExtendedInteger::Finite(_))) {
            "0 -> ".to_owned()
        } else {
            "--- -> ".to_owned()
        };
        if self.is_positive_repeating()?
            && upper == Some(ExtendedInteger::PositiveInfinity)
            && let Some(mut degrees) = self.positive_repeat_degrees()?
        {
            output.push_str("[ ");
            degrees.reverse();
            for degree in degrees {
                write!(
                    output,
                    "{}:{} -> ",
                    degree,
                    self.category().object_string(&self.object(&degree)?)
                )
                .expect("writing to a String cannot fail");
            }
            output.push_str("] ");
        }
        if let (Some(top), Some(bottom)) = (top, bottom) {
            let mut degree = top;
            while degree >= bottom {
                write!(
                    output,
                    "{}:{} -> ",
                    degree,
                    self.category().object_string(&self.object(&degree)?)
                )
                .expect("writing to a String cannot fail");
                degree -= Integer::from(1);
            }
        }
        if self.is_negative_repeating()?
            && lower == Some(ExtendedInteger::NegativeInfinity)
            && let Some(mut degrees) = self.negative_repeat_degrees()?
        {
            output.push_str("[ ");
            degrees.reverse();
            for degree in degrees {
                write!(
                    output,
                    "{}:{} -> ",
                    degree,
                    self.category().object_string(&self.object(&degree)?)
                )
                .expect("writing to a String cannot fail");
            }
            output.push_str("] ");
        }
        output.push_str(if matches!(lower, Some(ExtendedInteger::Finite(_))) {
            "0"
        } else {
            "---"
        });
        Ok(output)
    }
}

impl<C: QpaAbelianCategory> QpaComplex<C> {
    /// Translate `HomologyOfComplex(C,i)` through QPA's inclusion lifting and
    /// cokernel construction.
    pub fn homology(&self, degree: &Integer) -> Result<C::Object, ComplexError> {
        let kernel = self
            .category()
            .kernel_inclusion(&self.differential(degree)?)?;
        let image = self
            .category()
            .image_inclusion(&self.differential(&(degree + Integer::from(1)))?)?;
        let lifting = self
            .category()
            .lifting_inclusion_morphisms(&kernel, &image)?
            .ok_or(ComplexError::MissingInclusionLifting)?;
        let cokernel = self.category().cokernel_projection(&lifting)?;
        Ok(self.category().range(&cokernel))
    }

    /// Translate `GoodTruncationBelow(C,i)`.
    pub fn good_truncation_below(&self, degree: &Integer) -> Result<Self, ComplexError> {
        let category = self.category().clone();
        let differentials = self.differentials();
        let truncation = differentials.positive_part_from(&(degree + Integer::from(2)))?;
        let kernel = category.kernel_inclusion(&self.differential(degree)?)?;
        let image_projection =
            category.image_projection(&self.differential(&(degree + Integer::from(1)))?)?;
        let image_inclusion =
            category.image_inclusion(&self.differential(&(degree + Integer::from(1)))?)?;
        let lifting = category
            .lifting_inclusion_morphisms(&kernel, &image_inclusion)?
            .ok_or(ComplexError::MissingInclusionLifting)?;
        let new_part = InfiniteList::finite(
            degree.clone(),
            vec![
                category.zero_map(&category.source(&kernel), &category.zero_object())?,
                category.compose_maps(&lifting, &image_projection)?,
            ],
        );
        let zero_part = zero_complex(category.clone())?
            .differentials()
            .negative_part_from(&(degree - Integer::from(1)))?;
        Ok(Self::by_differential_list_with_dependencies(
            category,
            InfiniteList::concatenate(&[truncation, new_part, zero_part])?,
            &[self],
        ))
    }

    /// Translate QPA's explicitly documented "not working" implementation of
    /// `GoodTruncationAbove(C,i)` without repairing its source behavior.
    pub fn good_truncation_above(&self, degree: &Integer) -> Result<Self, ComplexError> {
        let category = self.category().clone();
        let differentials = self.differentials();
        let truncation = differentials.negative_part_from(&(degree - Integer::from(1)))?;
        let factor_inclusion = category.image_inclusion(&self.differential(degree)?)?;
        let factor = category.range(&factor_inclusion);
        let new_part = InfiniteList::finite(
            degree.clone(),
            vec![
                factor_inclusion,
                category.zero_map(&category.zero_object(), &factor)?,
            ],
        );
        let zero_part = zero_complex(category.clone())?
            .differentials()
            .positive_part_from(&(degree + Integer::from(2)))?;
        Ok(Self::by_differential_list_with_dependencies(
            category,
            InfiniteList::concatenate(&[zero_part, new_part, truncation])?,
            &[self],
        ))
    }

    /// Translate `GoodTruncation(C,upper,lower)`.
    pub fn good_truncation(&self, upper: &Integer, lower: &Integer) -> Result<Self, ComplexError> {
        if lower > upper {
            return Err(ComplexError::TruncationBoundsReversed);
        }
        self.good_truncation_below(lower)?
            .good_truncation_above(upper)
    }

    /// Translate `BrutalTruncationBelow(C,i)`.
    pub fn brutal_truncation_below(&self, degree: &Integer) -> Result<Self, ComplexError> {
        let category = self.category().clone();
        let truncation = self
            .differentials()
            .positive_part_from(&(degree + Integer::from(1)))?;
        let new_part = InfiniteList::finite(
            degree.clone(),
            vec![category.zero_map(&self.object(degree)?, &category.zero_object())?],
        );
        let zero_part = zero_complex(category.clone())?
            .differentials()
            .negative_part_from(&(degree - Integer::from(1)))?;
        Ok(Self::by_differential_list_with_dependencies(
            category,
            InfiniteList::concatenate(&[truncation, new_part, zero_part])?,
            &[self],
        ))
    }

    /// Translate `BrutalTruncationAbove(C,i)`.
    pub fn brutal_truncation_above(&self, degree: &Integer) -> Result<Self, ComplexError> {
        let category = self.category().clone();
        let truncation = self.differentials().negative_part_from(degree)?;
        let new_part = InfiniteList::finite(
            degree + Integer::from(1),
            vec![category.zero_map(&category.zero_object(), &self.object(degree)?)?],
        );
        let zero_part = zero_complex(category.clone())?
            .differentials()
            .positive_part_from(&(degree + Integer::from(2)))?;
        Ok(Self::by_differential_list_with_dependencies(
            category,
            InfiniteList::concatenate(&[zero_part, new_part, truncation])?,
            &[self],
        ))
    }

    /// Translate `BrutalTruncation(C,upper,lower)`.
    pub fn brutal_truncation(
        &self,
        upper: &Integer,
        lower: &Integer,
    ) -> Result<Self, ComplexError> {
        if lower > upper {
            return Err(ComplexError::TruncationBoundsReversed);
        }
        self.brutal_truncation_below(lower)?
            .brutal_truncation_above(upper)
    }

    /// Translate `SyzygyTruncation(C,i)`.
    pub fn syzygy_truncation(&self, degree: &Integer) -> Result<Self, ComplexError> {
        let category = self.category().clone();
        let truncation = self.differentials().negative_part_from(degree)?;
        let kernel_inclusion = category.kernel_inclusion(&self.differential(degree)?)?;
        let kernel = category.source(&kernel_inclusion);
        let new_part = InfiniteList::finite(
            degree + Integer::from(1),
            vec![
                kernel_inclusion,
                category.zero_map(&category.zero_object(), &kernel)?,
            ],
        );
        let zero_part = zero_complex(category.clone())?
            .differentials()
            .positive_part_from(&(degree + Integer::from(3)))?;
        Ok(Self::by_differential_list_with_dependencies(
            category,
            InfiniteList::concatenate(&[zero_part, new_part, truncation])?,
            &[self],
        ))
    }

    /// Translate `CosyzygyTruncation(C,i)`.
    pub fn cosyzygy_truncation(&self, degree: &Integer) -> Result<Self, ComplexError> {
        let category = self.category().clone();
        let truncation = self.differentials().positive_part_from(degree)?;
        let cokernel_projection = category.cokernel_projection(&self.differential(degree)?)?;
        let cokernel = category.range(&cokernel_projection);
        let new_part = InfiniteList::finite(
            degree - Integer::from(2),
            vec![
                category.zero_map(&cokernel, &category.zero_object())?,
                cokernel_projection,
            ],
        );
        let zero_part = zero_complex(category.clone())?
            .differentials()
            .negative_part_from(&(degree - Integer::from(3)))?;
        Ok(Self::by_differential_list_with_dependencies(
            category,
            InfiniteList::concatenate(&[truncation, new_part, zero_part])?,
            &[self],
        ))
    }

    /// Translate `SyzygyCosyzygyTruncation(C,upper,lower)` literally, including
    /// QPA's finite middle differential interval.
    pub fn syzygy_cosyzygy_truncation(
        &self,
        upper: &Integer,
        lower: &Integer,
    ) -> Result<Self, ComplexError> {
        if lower > upper {
            return Err(ComplexError::TruncationBoundsReversed);
        }
        let category = self.category().clone();
        let middle = self.differentials().finite_part_as_list(lower, upper)?;
        let truncation = InfiniteList::finite(lower.clone(), middle);

        let kernel_inclusion = category.kernel_inclusion(&self.differential(upper)?)?;
        let kernel = category.source(&kernel_inclusion);
        let upper_new = InfiniteList::finite(
            upper + Integer::from(1),
            vec![
                kernel_inclusion,
                category.zero_map(&category.zero_object(), &kernel)?,
            ],
        );
        let upper_zero = zero_complex(category.clone())?
            .differentials()
            .positive_part_from(&(upper + Integer::from(3)))?;

        let cokernel_projection = category.cokernel_projection(&self.differential(lower)?)?;
        let cokernel = category.range(&cokernel_projection);
        let lower_new = InfiniteList::finite(
            lower - Integer::from(2),
            vec![
                category.zero_map(&cokernel, &category.zero_object())?,
                cokernel_projection,
            ],
        );
        let lower_zero = zero_complex(category.clone())?
            .differentials()
            .negative_part_from(&(lower - Integer::from(3)))?;

        Ok(Self::by_differential_list_with_dependencies(
            category,
            InfiniteList::concatenate(&[upper_zero, upper_new, truncation, lower_new, lower_zero])?,
            &[self],
        ))
    }

    /// Translate `CutComplexAbove`.  Like QPA, this scans until it encounters
    /// the first zero object if the stored upper bound is not finite.
    pub fn cut_complex_above(&self) -> Result<Self, ComplexError> {
        if matches!(self.upper_bound()?, Some(ExtendedInteger::Finite(_))) {
            return Ok(self.clone());
        }
        let Some(ExtendedInteger::Finite(mut degree)) = self.lower_bound()? else {
            return Err(ComplexError::CutAboveNeedsFiniteLowerBound);
        };
        loop {
            if self.category().is_zero_object(&self.object(&degree)?) {
                return self.brutal_truncation_above(&(&degree - Integer::from(1)));
            }
            degree += Integer::from(1);
        }
    }

    /// Translate `CutComplexBelow`.  Like QPA, this scans until it encounters
    /// the first zero object if the stored lower bound is not finite.
    pub fn cut_complex_below(&self) -> Result<Self, ComplexError> {
        if matches!(self.lower_bound()?, Some(ExtendedInteger::Finite(_))) {
            return Ok(self.clone());
        }
        let Some(ExtendedInteger::Finite(mut degree)) = self.upper_bound()? else {
            return Err(ComplexError::CutBelowNeedsFiniteUpperBound);
        };
        loop {
            if self.category().is_zero_object(&self.object(&degree)?) {
                return self.brutal_truncation_below(&(&degree + Integer::from(1)));
            }
            degree -= Integer::from(1);
        }
    }
}

fn finite_known(value: Option<ExtendedInteger>) -> Option<Integer> {
    match value {
        Some(ExtendedInteger::Finite(value)) => Some(value),
        _ => None,
    }
}

fn integer_interval(first: &Integer, last: &Integer) -> Vec<Integer> {
    let mut result = Vec::new();
    let mut value = first.clone();
    while value <= *last {
        result.push(value.clone());
        value += Integer::from(1);
    }
    result
}

fn check_differential_pair<C: QpaCategory>(
    category: &C,
    top_degree: Integer,
    higher: &C::Morphism,
    lower: &C::Morphism,
    higher_name: &str,
    lower_name: &str,
) -> Result<(), ComplexError> {
    if category.range(higher) != category.source(lower) {
        return Err(ComplexError::DifferentialObjectMismatch {
            top_degree,
            higher: higher_name.to_owned(),
            lower: lower_name.to_owned(),
        });
    }
    if !category.is_zero_mapping(&category.compose_maps(lower, higher)?) {
        return Err(ComplexError::NonzeroDifferentialComposition {
            top_degree,
            higher: higher_name.to_owned(),
            lower: lower_name.to_owned(),
        });
    }
    Ok(())
}

fn check_differential_list<C: QpaCategory>(
    category: &C,
    list: &[C::Morphism],
    start_degree: &Integer,
    name: &str,
) -> Result<(), ComplexError> {
    for index in 1..list.len() {
        check_differential_pair(
            category,
            start_degree + Integer::from(index),
            &list[index],
            &list[index - 1],
            &format!("element {} in {name} list", index + 1),
            &format!("element {index} in {name} list"),
        )?;
    }
    Ok(())
}

fn check_repeating_differential_list<C: QpaCategory>(
    category: &C,
    list: &[C::Morphism],
    start_degree: &Integer,
    name: &str,
) -> Result<(), ComplexError> {
    check_differential_list(category, list, start_degree, name)?;
    let Some(first) = list.first() else {
        return Err(ComplexError::EmptyComplexRepeat {
            side: name.to_owned(),
        });
    };
    let last = list.last().expect("nonempty repeat has a final element");
    check_differential_pair(
        category,
        start_degree + Integer::from(list.len()),
        first,
        last,
        &format!("element 1 in {name} list"),
        &format!("element {} in {name} list", list.len()),
    )
}

fn normalize_half<C: QpaCategory>(
    half: ComplexHalf<C>,
    complex: &QpaComplex<C>,
    zero_map: C::Morphism,
) -> HalfInfiniteListSpecification<C::Morphism, ComplexError> {
    match half {
        ComplexHalf::Zero => HalfInfiniteListSpecification::Repeat(vec![zero_map]),
        ComplexHalf::Repeat(values) => HalfInfiniteListSpecification::Repeat(values),
        ComplexHalf::Next {
            function,
            initial_value,
            repeatify,
        } => HalfInfiniteListSpecification::Next {
            function,
            initial_value,
            repeatify,
        },
        ComplexHalf::Position {
            function,
            store_values,
        } => {
            let complex = complex.clone();
            HalfInfiniteListSpecification::Position {
                function: Arc::new(move |degree| function(&complex, degree)),
                store_values,
            }
        }
    }
}

/// Translate QPA's `Complex(cat, base, middle, positive, negative)`.
#[allow(clippy::too_many_lines)]
pub fn complex<C: QpaCategory>(
    category: C,
    base_position: Integer,
    middle: Vec<C::Morphism>,
    positive: ComplexHalf<C>,
    negative: ComplexHalf<C>,
) -> Result<QpaComplex<C>, ComplexError> {
    check_differential_list(&category, &middle, &base_position, "middle")?;

    if let ComplexHalf::Repeat(values) = &positive {
        check_repeating_differential_list(
            &category,
            values,
            &(&base_position + Integer::from(middle.len())),
            "positive",
        )?;
        if let (Some(first), Some(last_middle)) = (values.first(), middle.last()) {
            check_differential_pair(
                &category,
                &base_position + Integer::from(middle.len()),
                first,
                last_middle,
                "element 1 in positive list",
                &format!("element {} in middle list", middle.len()),
            )?;
        }
    }
    if let ComplexHalf::Repeat(values) = &negative {
        let mut reversed = values.clone();
        reversed.reverse();
        check_repeating_differential_list(
            &category,
            &reversed,
            &(&base_position + Integer::from(middle.len())),
            "negative",
        )?;
        if let (Some(first_middle), Some(last_negative)) = (middle.first(), reversed.last()) {
            check_differential_pair(
                &category,
                base_position.clone(),
                first_middle,
                last_negative,
                "element 1 in middle list",
                &format!("element {} in negative list", reversed.len()),
            )?;
        }
        if middle.is_empty()
            && let ComplexHalf::Repeat(positive_values) = &positive
            && let (Some(first_positive), Some(last_negative)) =
                (positive_values.first(), reversed.last())
        {
            check_differential_pair(
                &category,
                base_position.clone(),
                first_positive,
                last_negative,
                "element 1 in positive list",
                &format!("element {} in negative list", reversed.len()),
            )?;
        }
    }

    let positive_is_zero = matches!(positive, ComplexHalf::Zero);
    let negative_is_zero = matches!(negative, ComplexHalf::Zero);
    let last_middle_object = if positive_is_zero {
        Some(
            category.source(
                middle
                    .last()
                    .ok_or(ComplexError::EmptyMiddleWithZeroHalf { side: "positive" })?,
            ),
        )
    } else {
        None
    };
    let first_middle_object = if negative_is_zero {
        Some(
            category.range(
                middle
                    .first()
                    .ok_or(ComplexError::EmptyMiddleWithZeroHalf { side: "negative" })?,
            ),
        )
    } else {
        None
    };
    let mut normalized_middle = middle;
    let mut normalized_base = base_position;
    if positive_is_zero {
        let last_object = last_middle_object.expect("positive zero captured its boundary object");
        if !category.is_zero_object(&last_object) {
            normalized_middle.push(category.zero_map(&category.zero_object(), &last_object)?);
        }
        while normalized_middle
            .last()
            .is_some_and(|map| category.is_zero_object(&category.range(map)))
        {
            normalized_middle.pop();
        }
    }
    if negative_is_zero {
        let first_object = first_middle_object.expect("negative zero captured its boundary object");
        if !category.is_zero_object(&first_object) {
            normalized_middle.insert(
                0,
                category.zero_map(&first_object, &category.zero_object())?,
            );
            normalized_base -= Integer::from(1);
        }
        while normalized_middle
            .first()
            .is_some_and(|map| category.is_zero_object(&category.source(map)))
        {
            normalized_middle.remove(0);
            normalized_base += Integer::from(1);
        }
    }

    if positive_is_zero
        && normalized_middle.is_empty()
        && let ComplexHalf::Next {
            function,
            initial_value,
            repeatify: true,
        } = &negative
    {
        let second = function(initial_value)?;
        if category.is_zero_mapping(initial_value) && category.is_zero_mapping(&second) {
            return zero_complex(category);
        }
    }
    if negative_is_zero
        && normalized_middle.is_empty()
        && let ComplexHalf::Next {
            function,
            initial_value,
            repeatify: true,
        } = &positive
    {
        let second = function(initial_value)?;
        if category.is_zero_mapping(initial_value) && category.is_zero_mapping(&second) {
            return zero_complex(category);
        }
    }

    let data = Arc::new(QpaComplexData {
        category: category.clone(),
        differentials: OnceLock::new(),
        zero_complex: false,
        projective_complex: OnceLock::new(),
        injective_complex: OnceLock::new(),
        _dependencies: Vec::new(),
    });
    let result = QpaComplex(data);
    let zero_map = category.zero_map(&category.zero_object(), &category.zero_object())?;
    let positive = normalize_half(positive, &result, zero_map.clone());
    let negative = normalize_half(negative, &result, zero_map);
    let weak: Weak<QpaComplexData<C>> = Arc::downgrade(&result.0);
    let callback = Arc::new(move |degree: &Integer, direction, _list_type| {
        let complex = QpaComplex(
            weak.upgrade()
                .ok_or(ComplexError::ComplexDroppedDuringLazyEvaluation)?,
        );
        let (lower_degree, higher_degree) = match direction {
            InfiniteListDirection::Positive => (degree - Integer::from(1), degree.clone()),
            InfiniteListDirection::Negative => (degree.clone(), degree + Integer::from(1)),
        };
        let lower = complex.differential(&lower_degree)?;
        let higher = complex.differential(&higher_degree)?;
        check_differential_pair(
            complex.category(),
            higher_degree,
            &higher,
            &lower,
            "new lazy differential",
            "adjacent known differential",
        )
    });
    let differentials = InfiniteList::new(
        normalized_base,
        normalized_middle,
        positive,
        negative,
        Some(callback),
    );
    result
        .0
        .differentials
        .set(differentials)
        .expect("new QPA complex has an empty differential slot");
    Ok(result)
}

/// Translate `FiniteComplex`.
pub fn finite_complex<C: QpaCategory>(
    category: C,
    base_position: Integer,
    differentials: Vec<C::Morphism>,
) -> Result<QpaComplex<C>, ComplexError> {
    complex(
        category,
        base_position,
        differentials,
        ComplexHalf::Zero,
        ComplexHalf::Zero,
    )
}

/// Translate `ZeroComplex`.
pub fn zero_complex<C: QpaCategory>(category: C) -> Result<QpaComplex<C>, ComplexError> {
    let zero = category.zero_object();
    let zero_map = category.zero_map(&zero, &zero)?;
    let data = Arc::new(QpaComplexData {
        category,
        differentials: OnceLock::new(),
        zero_complex: true,
        projective_complex: OnceLock::new(),
        injective_complex: OnceLock::new(),
        _dependencies: Vec::new(),
    });
    data.differentials
        .set(InfiniteList::constant(zero_map))
        .expect("new QPA zero complex has an empty differential slot");
    Ok(QpaComplex(data))
}

/// Translate `StalkComplex`.
pub fn stalk_complex<C: QpaCategory>(
    category: C,
    object: &C::Object,
    degree: Integer,
) -> Result<QpaComplex<C>, ComplexError> {
    let zero = category.zero_object();
    finite_complex(
        category.clone(),
        degree,
        vec![
            category.zero_map(object, &zero)?,
            category.zero_map(&zero, object)?,
        ],
    )
}

/// Translate `ShortExactSequence` and its exactness check.
pub fn short_exact_sequence<C: QpaCategory>(
    category: C,
    first: C::Morphism,
    second: C::Morphism,
) -> Result<QpaComplex<C>, ComplexError> {
    let result = finite_complex(category, Integer::from(0), vec![second, first])?;
    if !result.is_short_exact_sequence()? {
        return Err(ComplexError::NotShortExact);
    }
    Ok(result)
}

/// Translate `YonedaProduct`.
pub fn yoneda_product<C: QpaCategory>(
    first: &QpaComplex<C>,
    second: &QpaComplex<C>,
) -> Result<QpaComplex<C>, ComplexError> {
    if first.is_zero_complex() {
        return Ok(second.clone());
    }
    let Some(ExtendedInteger::Finite(lower_first)) = first.lower_bound()? else {
        return Err(ComplexError::YonedaFirstNotBoundedBelow);
    };
    let Some(ExtendedInteger::Finite(upper_second)) = second.upper_bound()? else {
        return Err(ComplexError::YonedaSecondNotBoundedAbove);
    };
    if first.object(&lower_first)? != second.object(&upper_second)? {
        return Err(ComplexError::YonedaIncompatibleObjects);
    }
    let shifted = first
        .differentials()
        .shift(&(&lower_first - &upper_second + Integer::from(1)))?;
    let connection = first.category().compose_maps(
        &second.differential(&upper_second)?,
        &first.differential(&(&lower_first + Integer::from(1)))?,
    )?;
    let lists = vec![
        shifted.positive_part_from(&(&upper_second + Integer::from(1)))?,
        InfiniteList::finite(Integer::from(0), vec![connection]),
        second
            .differentials()
            .negative_part_from(&(&upper_second - Integer::from(1)))?,
    ];
    Ok(QpaComplex::by_differential_list_with_dependencies(
        first.category().clone(),
        InfiniteList::concatenate(&lists)?,
        &[first, second],
    ))
}

/// QPA's descriptions of one side of a chain map.
pub enum ChainMapHalf<C: QpaCategory> {
    /// QPA's string `"zero"`.
    Zero,
    /// A repeating morphism list.
    Repeat(Vec<C::Morphism>),
    /// A predecessor-generated side.
    Next {
        /// Next-morphism function.
        function:
            Arc<dyn Fn(&C::Morphism) -> Result<C::Morphism, ComplexError> + Send + Sync + 'static>,
        /// Initial predecessor value.
        initial_value: C::Morphism,
        /// Whether collision-based repeat conversion is enabled.
        repeatify: bool,
    },
    /// A degree-generated side.
    Position {
        /// Function receiving the chain map and requested degree.
        function: Arc<
            dyn Fn(&QpaChainMap<C>, &Integer) -> Result<C::Morphism, ComplexError>
                + Send
                + Sync
                + 'static,
        >,
        /// Whether generated values are retained.
        store_values: bool,
    },
}

impl<C: QpaCategory> Clone for ChainMapHalf<C> {
    fn clone(&self) -> Self {
        match self {
            Self::Zero => Self::Zero,
            Self::Repeat(values) => Self::Repeat(values.clone()),
            Self::Next {
                function,
                initial_value,
                repeatify,
            } => Self::Next {
                function: Arc::clone(function),
                initial_value: initial_value.clone(),
                repeatify: *repeatify,
            },
            Self::Position {
                function,
                store_values,
            } => Self::Position {
                function: Arc::clone(function),
                store_values: *store_values,
            },
        }
    }
}

struct QpaChainMapData<C: QpaCategory> {
    source: QpaComplex<C>,
    range: QpaComplex<C>,
    morphisms: OnceLock<InfiniteList<C::Morphism, ComplexError>>,
}

/// QPA's `IsChainMap` object.
pub struct QpaChainMap<C: QpaCategory>(Arc<QpaChainMapData<C>>);

impl<C: QpaCategory> Clone for QpaChainMap<C> {
    fn clone(&self) -> Self {
        Self(Arc::clone(&self.0))
    }
}

impl<C: QpaCategory> fmt::Debug for QpaChainMap<C> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("QpaChainMap")
            .field("source", &self.0.source)
            .field("range", &self.0.range)
            .finish_non_exhaustive()
    }
}

impl<C: QpaCategory> PartialEq for QpaChainMap<C> {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }
}

impl<C: QpaCategory> Eq for QpaChainMap<C> {}

impl<C: QpaCategory> QpaChainMap<C> {
    /// QPA's `Source` attribute.
    #[must_use]
    pub fn source(&self) -> &QpaComplex<C> {
        &self.0.source
    }

    /// QPA's `Range` attribute.
    #[must_use]
    pub fn range(&self) -> &QpaComplex<C> {
        &self.0.range
    }

    /// QPA's `MorphismsOfChainMap` attribute.
    #[must_use]
    pub fn morphisms(&self) -> InfiniteList<C::Morphism, ComplexError> {
        self.0
            .morphisms
            .get()
            .expect("QPA chain map morphism list is initialized before exposure")
            .clone()
    }

    /// Translate `MorphismOfChainMap(map,i)` and `map^i`.
    pub fn morphism(&self, degree: &Integer) -> Result<C::Morphism, ComplexError> {
        self.morphisms().value(degree)
    }

    /// Translate chain-map `HighestKnownDegree`.
    pub fn highest_known_degree(&self) -> Result<Option<ExtendedInteger>, ComplexError> {
        self.morphisms().highest_known_position()
    }

    /// Translate chain-map `LowestKnownDegree`.
    pub fn lowest_known_degree(&self) -> Result<Option<ExtendedInteger>, ComplexError> {
        self.morphisms().lowest_known_position()
    }

    /// QPA's printed representation.
    #[must_use]
    pub fn qpa_string(&self) -> &'static str {
        "<chain map>"
    }
}

fn normalize_chain_map_half<C: QpaCategory>(
    half: ChainMapHalf<C>,
    map: &QpaChainMap<C>,
    source: &QpaComplex<C>,
    range: &QpaComplex<C>,
) -> HalfInfiniteListSpecification<C::Morphism, ComplexError> {
    match half {
        ChainMapHalf::Zero => {
            let category = source.category().clone();
            let source = source.clone();
            let range = range.clone();
            HalfInfiniteListSpecification::Position {
                function: Arc::new(move |degree| {
                    category.zero_map(&source.object(degree)?, &range.object(degree)?)
                }),
                store_values: false,
            }
        }
        ChainMapHalf::Repeat(values) => HalfInfiniteListSpecification::Repeat(values),
        ChainMapHalf::Next {
            function,
            initial_value,
            repeatify,
        } => HalfInfiniteListSpecification::Next {
            function,
            initial_value,
            repeatify,
        },
        ChainMapHalf::Position {
            function,
            store_values,
        } => {
            let map = map.clone();
            HalfInfiniteListSpecification::Position {
                function: Arc::new(move |degree| function(&map, degree)),
                store_values,
            }
        }
    }
}

fn check_chain_map_domain<C: QpaCategory>(
    map: &QpaChainMap<C>,
    degree: &Integer,
) -> Result<(), ComplexError> {
    let category = map.source().category();
    let morphism = map.morphism(degree)?;
    if category.source(&morphism) != map.source().object(degree)? {
        return Err(ComplexError::ChainMapWrongSource {
            degree: degree.clone(),
        });
    }
    if category.range(&morphism) != map.range().object(degree)? {
        return Err(ComplexError::ChainMapWrongRange {
            degree: degree.clone(),
        });
    }
    Ok(())
}

fn check_chain_map_square<C: QpaCategory>(
    map: &QpaChainMap<C>,
    top_degree: &Integer,
) -> Result<(), ComplexError> {
    let category = map.source().category();
    let range_route = category.compose_maps(
        &map.range().differential(top_degree)?,
        &map.morphism(top_degree)?,
    )?;
    let source_route = category.compose_maps(
        &map.morphism(&(top_degree - Integer::from(1)))?,
        &map.source().differential(top_degree)?,
    )?;
    if range_route != source_route {
        return Err(ComplexError::ChainMapSquareDoesNotCommute {
            top_degree: top_degree.clone(),
        });
    }
    Ok(())
}

/// Translate QPA's six-argument `ChainMap` constructor.
#[allow(clippy::too_many_lines)]
pub fn chain_map<C: QpaCategory>(
    source: QpaComplex<C>,
    range: QpaComplex<C>,
    base_position: Integer,
    middle: Vec<C::Morphism>,
    positive: ChainMapHalf<C>,
    negative: ChainMapHalf<C>,
) -> Result<QpaChainMap<C>, ComplexError> {
    if source.category() != range.category() {
        return Err(ComplexError::ChainMapDifferentCategories);
    }
    let negative_repeat_length = match &negative {
        ChainMapHalf::Repeat(values) => Some(values.len()),
        _ => None,
    };
    let positive_repeat_length = match &positive {
        ChainMapHalf::Repeat(values) => Some(values.len()),
        _ => None,
    };
    let middle_length = middle.len();
    let data = Arc::new(QpaChainMapData {
        source: source.clone(),
        range: range.clone(),
        morphisms: OnceLock::new(),
    });
    let result = QpaChainMap(data);
    let positive = normalize_chain_map_half(positive, &result, &source, &range);
    let negative = normalize_chain_map_half(negative, &result, &source, &range);
    let weak: Weak<QpaChainMapData<C>> = Arc::downgrade(&result.0);
    let callback = Arc::new(move |degree: &Integer, direction, _list_type| {
        let map = QpaChainMap(
            weak.upgrade()
                .ok_or(ComplexError::ChainMapDroppedDuringLazyEvaluation)?,
        );
        check_chain_map_domain(&map, degree)?;
        let top_degree = match direction {
            InfiniteListDirection::Positive => degree.clone(),
            InfiniteListDirection::Negative => degree + Integer::from(1),
        };
        check_chain_map_square(&map, &top_degree)
    });
    let morphisms = InfiniteList::new(
        base_position.clone(),
        middle,
        positive,
        negative,
        Some(callback),
    );
    result
        .0
        .morphisms
        .set(morphisms)
        .expect("new QPA chain map has an empty morphism slot");

    let first = negative_repeat_length.map_or_else(
        || base_position.clone(),
        |length| &base_position - Integer::from(1 + length),
    );
    let last = positive_repeat_length.map_or_else(
        || &base_position + Integer::from(middle_length) - Integer::from(1),
        |length| &base_position + Integer::from(middle_length + length),
    );
    let mut degree = first.clone();
    while degree <= last {
        check_chain_map_domain(&result, &degree)?;
        degree += Integer::from(1);
    }
    let mut degree = first + Integer::from(1);
    while degree <= last {
        check_chain_map_square(&result, &degree)?;
        degree += Integer::from(1);
    }
    Ok(result)
}

/// Translate `FiniteChainMap`.
pub fn finite_chain_map<C: QpaCategory>(
    source: QpaComplex<C>,
    range: QpaComplex<C>,
    base_degree: Integer,
    morphisms: Vec<C::Morphism>,
) -> Result<QpaChainMap<C>, ComplexError> {
    chain_map(
        source,
        range,
        base_degree,
        morphisms,
        ChainMapHalf::Zero,
        ChainMapHalf::Zero,
    )
}

/// Translate `ZeroChainMap`.
pub fn zero_chain_map<C: QpaCategory>(
    source: QpaComplex<C>,
    range: QpaComplex<C>,
) -> Result<QpaChainMap<C>, ComplexError> {
    chain_map(
        source,
        range,
        Integer::from(0),
        Vec::new(),
        ChainMapHalf::Zero,
        ChainMapHalf::Zero,
    )
}

/// One lazy side supplied to [`complex_and_chain_maps`].
pub enum ComplexAndChainMapsHalf<C: QpaCategory> {
    /// QPA's string `"zero"`.
    Zero,
    /// Repeating tuples `[differential, incoming..., outgoing...]`.
    Repeat(Vec<Vec<C::Morphism>>),
    /// Predecessor-generated tuples.
    Next {
        /// Next-tuple function.
        function: Arc<
            dyn Fn(&Vec<C::Morphism>) -> Result<Vec<C::Morphism>, ComplexError>
                + Send
                + Sync
                + 'static,
        >,
        /// Initial predecessor tuple.
        initial_value: Vec<C::Morphism>,
        /// Whether collision-based repeat conversion is enabled.
        repeatify: bool,
    },
    /// Degree-generated tuples, with the new complex and every new map in the
    /// same argument order used by QPA's `CallFuncList` call.
    Position {
        /// Tuple generator.
        function: Arc<
            dyn Fn(
                    &QpaComplex<C>,
                    &[QpaChainMap<C>],
                    &[QpaChainMap<C>],
                    &Integer,
                ) -> Result<Vec<C::Morphism>, ComplexError>
                + Send
                + Sync
                + 'static,
        >,
        /// Whether generated tuples are retained.
        store_values: bool,
    },
}

struct CombinedContext<C: QpaCategory> {
    complex: OnceLock<QpaComplex<C>>,
    incoming: OnceLock<Vec<QpaChainMap<C>>>,
    outgoing: OnceLock<Vec<QpaChainMap<C>>>,
}

/// Return value of QPA's list-valued `ComplexAndChainMaps` constructor.
pub struct ComplexAndChainMaps<C: QpaCategory> {
    /// Newly constructed complex.
    pub complex: QpaComplex<C>,
    /// Maps from the entered source complexes to `complex`.
    pub incoming: Vec<QpaChainMap<C>>,
    /// Maps from `complex` to the entered range complexes.
    pub outgoing: Vec<QpaChainMap<C>>,
}

fn normalize_combined_half<C: QpaCategory>(
    half: ComplexAndChainMapsHalf<C>,
    context: &Arc<CombinedContext<C>>,
) -> (
    bool,
    HalfInfiniteListSpecification<Vec<C::Morphism>, ComplexError>,
) {
    match half {
        ComplexAndChainMapsHalf::Zero => (
            true,
            HalfInfiniteListSpecification::Repeat(vec![Vec::new()]),
        ),
        ComplexAndChainMapsHalf::Repeat(values) => {
            (false, HalfInfiniteListSpecification::Repeat(values))
        }
        ComplexAndChainMapsHalf::Next {
            function,
            initial_value,
            repeatify,
        } => (
            false,
            HalfInfiniteListSpecification::Next {
                function,
                initial_value,
                repeatify,
            },
        ),
        ComplexAndChainMapsHalf::Position {
            function,
            store_values,
        } => {
            let context = Arc::clone(context);
            (
                false,
                HalfInfiniteListSpecification::Position {
                    function: Arc::new(move |degree| {
                        function(
                            context
                                .complex
                                .get()
                                .ok_or(ComplexError::CombinedComplexNotInitialized)?,
                            context
                                .incoming
                                .get()
                                .ok_or(ComplexError::CombinedMapsNotInitialized)?,
                            context
                                .outgoing
                                .get()
                                .ok_or(ComplexError::CombinedMapsNotInitialized)?,
                            degree,
                        )
                    }),
                    store_values,
                },
            )
        }
    }
}

fn combined_component<C: QpaCategory>(
    list: &InfiniteList<Vec<C::Morphism>, ComplexError>,
    degree: &Integer,
    index: usize,
) -> Result<C::Morphism, ComplexError> {
    list.value(degree)?
        .get(index)
        .cloned()
        .ok_or(ComplexError::CombinedTupleTooShort {
            degree: degree.clone(),
            index: index + 1,
        })
}

/// Translate QPA's `ComplexAndChainMaps` circular lazy construction.
#[allow(clippy::too_many_lines)]
pub fn complex_and_chain_maps<C: QpaCategory>(
    source_complexes: &[QpaComplex<C>],
    range_complexes: &[QpaComplex<C>],
    base_position: Integer,
    middle: Vec<Vec<C::Morphism>>,
    positive: ComplexAndChainMapsHalf<C>,
    negative: ComplexAndChainMapsHalf<C>,
) -> Result<ComplexAndChainMaps<C>, ComplexError> {
    let category = source_complexes
        .first()
        .or_else(|| range_complexes.first())
        .ok_or(ComplexError::CombinedComplexListEmpty)?
        .category()
        .clone();
    if source_complexes
        .iter()
        .chain(range_complexes)
        .any(|complex| complex.category() != &category)
    {
        return Err(ComplexError::ChainMapDifferentCategories);
    }
    let expected_width = 1 + source_complexes.len() + range_complexes.len();
    if let Some((position, _tuple)) = middle
        .iter()
        .enumerate()
        .find(|(_, tuple)| tuple.len() < expected_width)
    {
        return Err(ComplexError::CombinedTupleTooShort {
            degree: &base_position + Integer::from(position),
            index: expected_width,
        });
    }

    let context = Arc::new(CombinedContext {
        complex: OnceLock::new(),
        incoming: OnceLock::new(),
        outgoing: OnceLock::new(),
    });
    let (positive_zero, positive) = normalize_combined_half(positive, &context);
    let (negative_zero, negative) = normalize_combined_half(negative, &context);
    let combined = InfiniteList::new(
        base_position.clone(),
        middle.clone(),
        positive,
        negative,
        None,
    );

    let middle_differentials = middle.iter().map(|tuple| tuple[0].clone()).collect();
    let positive_differentials = if positive_zero {
        ComplexHalf::Zero
    } else {
        let combined = combined.clone();
        ComplexHalf::Position {
            function: Arc::new(move |_complex, degree| {
                combined_component::<C>(&combined, degree, 0)
            }),
            store_values: true,
        }
    };
    let negative_differentials = if negative_zero {
        ComplexHalf::Zero
    } else {
        let combined = combined.clone();
        ComplexHalf::Position {
            function: Arc::new(move |_complex, degree| {
                combined_component::<C>(&combined, degree, 0)
            }),
            store_values: true,
        }
    };
    let new_complex = complex(
        category,
        base_position.clone(),
        middle_differentials,
        positive_differentials,
        negative_differentials,
    )?;
    context
        .complex
        .set(new_complex.clone())
        .map_err(|_| ComplexError::CombinedComplexAlreadyInitialized)?;

    let mut incoming = Vec::with_capacity(source_complexes.len());
    for (index, source) in source_complexes.iter().enumerate() {
        let component = 1 + index;
        let middle_maps = middle
            .iter()
            .map(|tuple| tuple[component].clone())
            .collect();
        let positive_maps = if positive_zero {
            ChainMapHalf::Zero
        } else {
            let combined = combined.clone();
            ChainMapHalf::Position {
                function: Arc::new(move |_map, degree| {
                    combined_component::<C>(&combined, degree, component)
                }),
                store_values: true,
            }
        };
        let negative_maps = if negative_zero {
            ChainMapHalf::Zero
        } else {
            let combined = combined.clone();
            ChainMapHalf::Position {
                function: Arc::new(move |_map, degree| {
                    combined_component::<C>(&combined, degree, component)
                }),
                store_values: true,
            }
        };
        incoming.push(chain_map(
            source.clone(),
            new_complex.clone(),
            base_position.clone(),
            middle_maps,
            positive_maps,
            negative_maps,
        )?);
    }

    let mut outgoing = Vec::with_capacity(range_complexes.len());
    for (index, range) in range_complexes.iter().enumerate() {
        let component = 1 + source_complexes.len() + index;
        let middle_maps = middle
            .iter()
            .map(|tuple| tuple[component].clone())
            .collect();
        let positive_maps = if positive_zero {
            ChainMapHalf::Zero
        } else {
            let combined = combined.clone();
            ChainMapHalf::Position {
                function: Arc::new(move |_map, degree| {
                    combined_component::<C>(&combined, degree, component)
                }),
                store_values: true,
            }
        };
        let negative_maps = if negative_zero {
            ChainMapHalf::Zero
        } else {
            let combined = combined.clone();
            ChainMapHalf::Position {
                function: Arc::new(move |_map, degree| {
                    combined_component::<C>(&combined, degree, component)
                }),
                store_values: true,
            }
        };
        outgoing.push(chain_map(
            new_complex.clone(),
            range.clone(),
            base_position.clone(),
            middle_maps,
            positive_maps,
            negative_maps,
        )?);
    }
    context
        .incoming
        .set(incoming.clone())
        .map_err(|_| ComplexError::CombinedMapsAlreadyInitialized)?;
    context
        .outgoing
        .set(outgoing.clone())
        .map_err(|_| ComplexError::CombinedMapsAlreadyInitialized)?;
    Ok(ComplexAndChainMaps {
        complex: new_complex,
        incoming,
        outgoing,
    })
}

/// A QPA complex in the category of right modules over a Numerica-backed
/// path algebra or quotient.
pub type ModuleComplex<F> = QpaComplex<RightModuleCategory<F>>;

/// A QPA chain map between right-module complexes.
pub type ModuleChainMap<F> = QpaChainMap<RightModuleCategory<F>>;

/// Construct QPA's lazy projective-resolution complex, including the entered
/// module in degree `-1`.
pub fn projective_resolution_complex<F>(
    module: &PathAlgebraModule<F>,
) -> Result<ModuleComplex<F>, ComplexError>
where
    F: QpaCoefficientField + Send + Sync + 'static,
    F::Element: Send + Sync,
{
    let category = cat_of_right_algebra_modules(module.right_acting_algebra())?;
    let cover = projective_cover(module)?;
    let initial = cover.clone();
    complex(
        category,
        Integer::from(0),
        vec![cover],
        ComplexHalf::Next {
            function: Arc::new(|differential| {
                let inclusion = kernel_inclusion(differential)?;
                let cover = projective_cover(inclusion.source())?;
                Ok(cover.compose(&inclusion)?)
            }),
            initial_value: initial,
            repeatify: true,
        },
        ComplexHalf::Zero,
    )
}

/// Construct QPA's lazy injective-resolution complex.
pub fn injective_resolution_complex<F>(
    module: &PathAlgebraModule<F>,
) -> Result<ModuleComplex<F>, ComplexError>
where
    F: QpaCoefficientField + Send + Sync + 'static,
    F::Element: Send + Sync,
{
    let category = cat_of_right_algebra_modules(module.right_acting_algebra())?;
    let envelope = injective_envelope(module)?;
    let initial = envelope.clone();
    complex(
        category,
        Integer::from(0),
        vec![envelope],
        ComplexHalf::Zero,
        ComplexHalf::Next {
            function: Arc::new(|differential| {
                let projection = cokernel_projection(differential)?;
                let envelope = injective_envelope(projection.range())?;
                Ok(projection.compose(&envelope)?)
            }),
            initial_value: initial,
            repeatify: true,
        },
    )
}

fn finite_module_bound(
    bound: Option<ExtendedInteger>,
    error: ComplexError,
) -> Result<Integer, ComplexError> {
    match bound {
        Some(ExtendedInteger::Finite(value)) => Ok(value),
        _ => Err(error),
    }
}

fn projective_lifting<F: Field>(
    surjection: &PathAlgebraModuleHomomorphism<F>,
    map: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, ComplexError> {
    lifting_morphism_from_projective(surjection, map)?
        .ok_or(ComplexError::ProjectiveLiftingReturnedFail)
}

/// Translate QPA's `ComparisonLifting(f,PC,EC)`.
#[allow(clippy::too_many_lines)]
pub fn comparison_lifting<F>(
    morphism: &PathAlgebraModuleHomomorphism<F>,
    projective_complex: &ModuleComplex<F>,
    exact_complex: &ModuleComplex<F>,
) -> Result<ModuleChainMap<F>, ComplexError>
where
    F: QpaCoefficientField + Send + Sync + 'static,
    F::Element: Send + Sync,
{
    if projective_complex.lower_bound()? != exact_complex.lower_bound()? {
        return Err(ComplexError::ComparisonDifferentLowerBounds);
    }
    let lower = finite_module_bound(
        projective_complex.lower_bound()?,
        ComplexError::ComparisonLowerBoundNotFinite,
    )?;
    if projective_complex.object(&lower)? != *morphism.source()
        || exact_complex.object(&lower)? != *morphism.range()
    {
        return Err(ComplexError::ComparisonMapWrongObjects);
    }

    let surjection = projective_complex
        .differential(&(&lower + Integer::from(1)))?
        .compose(morphism)?;
    let first_lifting = projective_lifting(
        &exact_complex.differential(&(&lower + Integer::from(1)))?,
        &surjection,
    )?;
    let mut middle = vec![morphism.clone(), first_lifting];

    let projective = projective_complex.clone();
    let exact = exact_complex.clone();
    let next_lifting = Arc::new(
        move |previous: &PathAlgebraModuleHomomorphism<F>, degree: &Integer| {
            if projective.object(degree)?.dimension() == 0 {
                return Ok(zero_module_mapping(
                    projective.object(degree)?,
                    exact.object(degree)?,
                )?);
            }
            let composite = projective.differential(degree)?.compose(previous)?;
            let (projection, inclusion) = image_projection_inclusion(&composite)?;
            let kernel = kernel_inclusion(&exact.differential(&(degree - Integer::from(1)))?)?;
            let into_kernel = lifting_inclusion_morphisms(&kernel, &inclusion)?
                .ok_or(ComplexError::MissingInclusionLifting)?;
            let map = projection.compose(&into_kernel)?;

            let (exact_projection, exact_inclusion) =
                image_projection_inclusion(&exact.differential(degree)?)?;
            let exact_into_kernel = lifting_inclusion_morphisms(&kernel, &exact_inclusion)?
                .ok_or(ComplexError::MissingInclusionLifting)?;
            let exact_surjection = exact_projection.compose(&exact_into_kernel)?;
            projective_lifting(&exact_surjection, &map)
        },
    );

    match (
        projective_complex.upper_bound()?,
        exact_complex.upper_bound()?,
    ) {
        (
            Some(ExtendedInteger::Finite(projective_upper)),
            Some(ExtendedInteger::Finite(exact_upper)),
        ) => {
            let mut degree = &lower + Integer::from(2);
            while degree <= projective_upper {
                let previous = middle
                    .last()
                    .ok_or(ComplexError::EmptyComparisonMiddle)?
                    .clone();
                middle.push(next_lifting(&previous, &degree)?);
                degree += Integer::from(1);
            }
            while degree <= exact_upper {
                middle.push(zero_module_mapping(
                    projective_complex.object(&degree)?,
                    exact_complex.object(&degree)?,
                )?);
                degree += Integer::from(1);
            }
            finite_chain_map(
                projective_complex.clone(),
                exact_complex.clone(),
                lower,
                middle,
            )
        }
        _ => {
            let next = Arc::clone(&next_lifting);
            chain_map(
                projective_complex.clone(),
                exact_complex.clone(),
                lower,
                middle,
                ChainMapHalf::Position {
                    function: Arc::new(move |map, degree| {
                        next(&map.morphism(&(degree - Integer::from(1)))?, degree)
                    }),
                    store_values: true,
                },
                ChainMapHalf::Zero,
            )
        }
    }
}

/// Translate `ComparisonLiftingToProjectiveResolution`.
pub fn comparison_lifting_to_projective_resolution<F>(
    morphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<ModuleChainMap<F>, ComplexError>
where
    F: QpaCoefficientField + Send + Sync + 'static,
    F::Element: Send + Sync,
{
    comparison_lifting(
        morphism,
        &projective_resolution_complex(morphism.source())?,
        &projective_resolution_complex(morphism.range())?,
    )
}

fn module_direct_sum<F: Field>(
    modules: &[PathAlgebraModule<F>],
) -> Result<PathAlgebraModule<F>, ComplexError> {
    direct_sum_of_qpa_modules(modules)?.ok_or(ComplexError::MappingConeDirectSumFailed)
}

fn cone_tuple<F>(
    map: &ModuleChainMap<F>,
    degree: &Integer,
    upper_object: &PathAlgebraModule<F>,
    lower_object: &PathAlgebraModule<F>,
) -> Result<Vec<PathAlgebraModuleHomomorphism<F>>, ComplexError>
where
    F: QpaCoefficientField + Send + Sync + 'static,
    F::Element: Send + Sync,
{
    let source = map.source();
    let range = map.range();
    let projections = direct_sum_projections(upper_object)?;
    let inclusions = direct_sum_inclusions(lower_object)?;
    let matrix = vec![
        vec![
            source
                .differential(&(degree - Integer::from(1)))?
                .negate()?,
            zero_module_mapping(
                range.object(degree)?,
                source.object(&(degree - Integer::from(2)))?,
            )?,
        ],
        vec![
            map.morphism(&(degree - Integer::from(1)))?.negate()?,
            range.differential(degree)?,
        ],
    ];
    let differential = multiply_lists_of_maps(&projections, &matrix, &inclusions)?;
    Ok(vec![
        differential,
        direct_sum_inclusions(upper_object)?[1].clone(),
        direct_sum_projections(upper_object)?[0].clone(),
    ])
}

/// Translate `MappingCone(f)`, returning QPA's `[cone, inclusion, projection]`
/// data with its circular lazy construction.
pub fn mapping_cone<F>(
    map: &ModuleChainMap<F>,
) -> Result<ComplexAndChainMaps<RightModuleCategory<F>>, ComplexError>
where
    F: QpaCoefficientField + Send + Sync + 'static,
    F::Element: Send + Sync,
{
    let source = map.source();
    let range = map.range();
    let base = [
        source.lower_bound()?,
        range.lower_bound()?,
        source.upper_bound()?,
        range.upper_bound()?,
    ]
    .into_iter()
    .find_map(|bound| match bound {
        Some(ExtendedInteger::Finite(value)) => Some(value),
        _ => None,
    })
    .unwrap_or_else(|| Integer::from(0));

    let upper_object = module_direct_sum(&[
        source.object(&(&base - Integer::from(1)))?,
        range.object(&base)?,
    ])?;
    let lower_object = module_direct_sum(&[
        source.object(&(&base - Integer::from(2)))?,
        range.object(&(&base - Integer::from(1)))?,
    ])?;
    let middle = vec![cone_tuple(map, &base, &upper_object, &lower_object)?];

    let positive_map = map.clone();
    let positive = ComplexAndChainMapsHalf::Position {
        function: Arc::new(
            move |complex: &ModuleComplex<F>,
                  _incoming: &[ModuleChainMap<F>],
                  _outgoing: &[ModuleChainMap<F>],
                  degree: &Integer| {
                let next_object = module_direct_sum(&[
                    positive_map.source().object(&(degree - Integer::from(1)))?,
                    positive_map.range().object(degree)?,
                ])?;
                let previous_object = complex
                    .category()
                    .source(&complex.differential(&(degree - Integer::from(1)))?);
                cone_tuple(&positive_map, degree, &next_object, &previous_object)
            },
        ),
        store_values: true,
    };

    let negative_map = map.clone();
    let negative = ComplexAndChainMapsHalf::Position {
        function: Arc::new(
            move |complex: &ModuleComplex<F>,
                  _incoming: &[ModuleChainMap<F>],
                  _outgoing: &[ModuleChainMap<F>],
                  degree: &Integer| {
                let next_object = module_direct_sum(&[
                    negative_map.source().object(&(degree - Integer::from(2)))?,
                    negative_map.range().object(&(degree - Integer::from(1)))?,
                ])?;
                let previous_object = complex
                    .category()
                    .range(&complex.differential(&(degree + Integer::from(1)))?);
                cone_tuple(&negative_map, degree, &previous_object, &next_object)
            },
        ),
        store_values: true,
    };
    complex_and_chain_maps(
        std::slice::from_ref(range),
        &[source.shift(&Integer::from(-1))?],
        base,
        middle,
        positive,
        negative,
    )
}

/// Typed failures exposed by the source-translated complex layer.
#[derive(Debug, Error)]
pub enum ComplexError {
    /// Infinite-list construction or lazy access failed.
    #[error(transparent)]
    InfiniteList(#[from] InfiniteListError),
    /// A foundational module or module-map operation failed.
    #[error(transparent)]
    Module(#[from] PathAlgebraModuleError),
    /// A higher module operation failed.
    #[error(transparent)]
    ModuleOperation(#[from] EnvelopingModuleError),
    /// A runtime `Cat` callback reported a failure.
    #[error("category callback failed: {0}")]
    CategoryCallback(String),
    /// Consecutive differentials have incompatible intermediate objects.
    #[error("range/source mismatch between {higher} and {lower} at degree {top_degree}")]
    DifferentialObjectMismatch {
        /// Higher differential degree.
        top_degree: Integer,
        /// Higher differential description.
        higher: String,
        /// Lower differential description.
        lower: String,
    },
    /// Consecutive differentials have nonzero composition.
    #[error("nonzero composition between {higher} and {lower} at degree {top_degree}")]
    NonzeroDifferentialComposition {
        /// Higher differential degree.
        top_degree: Integer,
        /// Higher differential description.
        higher: String,
        /// Lower differential description.
        lower: String,
    },
    /// QPA indexes a nonempty middle when a side is `"zero"`.
    #[error("{side} zero half requires a nonempty middle differential list")]
    EmptyMiddleWithZeroHalf {
        /// Side whose QPA code indexed the middle.
        side: &'static str,
    },
    /// QPA's repeating-complex check requires a nonempty period.
    #[error("the {side} repeating differential list is empty")]
    EmptyComplexRepeat {
        /// Repeating side.
        side: String,
    },
    /// Lazy evaluation outlived its owning complex.
    #[error("the owning complex was dropped during lazy differential evaluation")]
    ComplexDroppedDuringLazyEvaluation,
    /// A complex made without a positive half reached an operation requiring it.
    #[error("complex has no positive infinite-list part")]
    MissingComplexPositivePart,
    /// A complex made without a negative half reached an operation requiring it.
    #[error("complex has no negative infinite-list part")]
    MissingComplexNegativePart,
    /// Bounds contradicted a successful finite-complex result.
    #[error("finite complex has nonfinite or unknown bounds")]
    InconsistentFiniteBounds,
    /// A GAP property method attempted to return `fail` instead of a boolean.
    #[error("method for property {property} did not return true or false")]
    PropertyReturnedNonBoolean {
        /// Property whose installed source method returned `fail`.
        property: &'static str,
    },
    /// The entered pair did not form a short exact sequence.
    #[error("not exact")]
    NotShortExact,
    /// First Yoneda factor was not bounded below.
    #[error("first argument in Yoneda product must be bounded below")]
    YonedaFirstNotBoundedBelow,
    /// Second Yoneda factor was not bounded above.
    #[error("second argument in Yoneda product must be bounded above")]
    YonedaSecondNotBoundedAbove,
    /// Yoneda boundary objects did not agree.
    #[error("non-compatible complexes for Yoneda product")]
    YonedaIncompatibleObjects,
    /// QPA's inclusion-lifting call returned `fail` where a map was consumed.
    #[error("LiftingInclusionMorphisms returned fail")]
    MissingInclusionLifting,
    /// A two-sided truncation received its bounds in the forbidden order.
    #[error("First input integer must be greater than or equal to the second")]
    TruncationBoundsReversed,
    /// `CutComplexAbove` could not begin its source loop.
    #[error("CutComplexAbove requires a finite lower bound")]
    CutAboveNeedsFiniteLowerBound,
    /// `CutComplexBelow` could not begin its source loop.
    #[error("CutComplexBelow requires a finite upper bound")]
    CutBelowNeedsFiniteUpperBound,
    /// Source and range complexes have different category identities.
    #[error("source and range of chain map must be complexes over the same cat")]
    ChainMapDifferentCategories,
    /// One chain-map component has the wrong source object.
    #[error("incorrect source of chain map morphism in degree {degree}")]
    ChainMapWrongSource {
        /// Offending degree.
        degree: Integer,
    },
    /// One chain-map component has the wrong range object.
    #[error("incorrect range of chain map morphism in degree {degree}")]
    ChainMapWrongRange {
        /// Offending degree.
        degree: Integer,
    },
    /// A chain-map square failed QPA's structural equality test.
    #[error(
        "chain map morphisms at degree {top_degree} and its predecessor do not commute with complex differentials"
    )]
    ChainMapSquareDoesNotCommute {
        /// Top degree of the square.
        top_degree: Integer,
    },
    /// Lazy evaluation outlived its owning chain map.
    #[error("the owning chain map was dropped during lazy morphism evaluation")]
    ChainMapDroppedDuringLazyEvaluation,
    /// `ComplexAndChainMaps` received no complex from which to recover a cat.
    #[error("ComplexAndChainMaps requires at least one source or range complex")]
    CombinedComplexListEmpty,
    /// A combined tuple lacked a required differential or component map.
    #[error("combined tuple at degree {degree} has no element {index}")]
    CombinedTupleTooShort {
        /// Tuple degree.
        degree: Integer,
        /// Missing one-based component index.
        index: usize,
    },
    /// Circular combined evaluation occurred before the complex was assigned.
    #[error("ComplexAndChainMaps complex is not initialized")]
    CombinedComplexNotInitialized,
    /// Circular combined evaluation occurred before the chain maps were assigned.
    #[error("ComplexAndChainMaps chain maps are not initialized")]
    CombinedMapsNotInitialized,
    /// Internal combined context was assigned twice.
    #[error("ComplexAndChainMaps complex was already initialized")]
    CombinedComplexAlreadyInitialized,
    /// Internal combined map context was assigned twice.
    #[error("ComplexAndChainMaps chain maps were already initialized")]
    CombinedMapsAlreadyInitialized,
    /// QPA's projective lifting returned `fail` where a map was consumed.
    #[error("LiftingMorphismFromProjective returned fail")]
    ProjectiveLiftingReturnedFail,
    /// Comparison complexes did not have the same lower bound.
    #[error("The complexes must have the same lower bound")]
    ComparisonDifferentLowerBounds,
    /// Comparison lifting needs the common lower bound as an integer.
    #[error("the comparison-complex lower bound is not finite")]
    ComparisonLowerBoundNotFinite,
    /// The entered module map did not match the lower-degree objects.
    #[error("The map has wrong source and/or range")]
    ComparisonMapWrongObjects,
    /// Internal comparison-lifting state unexpectedly lost its seed maps.
    #[error("comparison lifting has an empty middle map list")]
    EmptyComparisonMiddle,
    /// QPA's nonempty mapping-cone direct sum unexpectedly returned `fail`.
    #[error("DirectSumOfQPAModules returned fail while constructing a mapping cone")]
    MappingConeDirectSumFailed,
}

#[cfg(test)]
mod tests {
    use std::sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    };

    use numerica::domains::{
        Field,
        finite_field::{FiniteFieldCore, Z2},
        integer::Integer,
        rational::{Q, Rational},
    };
    use qpa_core::{
        ArrowSpec, DenseMatrix, PathAlgebra, QpaCoefficientField, Quiver, QuiverAlgebra,
        identity_module_mapping, right_module_hom_over_algebra, right_module_over_path_algebra,
    };
    use serde_json::{Value, json};

    use super::*;

    fn map_json<F: Field>(
        map: &PathAlgebraModuleHomomorphism<F>,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        json!([
            map.source().dimension_vector(),
            map.range().dimension_vector(),
            map.vertex_matrices()
                .iter()
                .map(|matrix| (0..matrix.nrows())
                    .map(|row| matrix.row(row).iter().map(coefficient).collect::<Vec<_>>())
                    .collect::<Vec<_>>())
                .collect::<Vec<_>>()
        ])
    }

    fn complex_sample<F>(
        complex: &ModuleComplex<F>,
        degrees: &[i64],
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value
    where
        F: QpaCoefficientField + Send + Sync + 'static,
        F::Element: Send + Sync,
    {
        json!({
            "objects": degrees.iter().map(|degree| complex
                .object(&Integer::from(*degree)).unwrap().dimension_vector().to_vec())
                .collect::<Vec<_>>(),
            "differentials": degrees.iter().map(|degree| map_json(
                &complex.differential(&Integer::from(*degree)).unwrap(), coefficient))
                .collect::<Vec<_>>(),
        })
    }

    fn chain_map_sample<F>(
        map: &ModuleChainMap<F>,
        degrees: &[i64],
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value
    where
        F: QpaCoefficientField + Send + Sync + 'static,
        F::Element: Send + Sync,
    {
        json!({
            "morphisms": degrees.iter().map(|degree| map_json(
                &map.morphism(&Integer::from(*degree)).unwrap(), coefficient))
                .collect::<Vec<_>>(),
            "highest_known": known_string(map.highest_known_degree().unwrap()),
            "lowest_known": known_string(map.lowest_known_degree().unwrap()),
        })
    }

    fn extended_string(value: &ExtendedInteger) -> String {
        match value {
            ExtendedInteger::NegativeInfinity => "-infinity".to_owned(),
            ExtendedInteger::Finite(value) => value.to_string(),
            ExtendedInteger::PositiveInfinity => "infinity".to_owned(),
        }
    }

    fn known_string(value: Option<ExtendedInteger>) -> String {
        value.map_or_else(|| "none".to_owned(), |value| extended_string(&value))
    }

    fn bound_string(value: Option<ExtendedInteger>) -> String {
        value.map_or_else(|| "fail".to_owned(), |value| extended_string(&value))
    }

    fn repeat_json(value: Option<Vec<Integer>>) -> Value {
        value.map_or_else(
            || json!("fail"),
            |degrees| {
                json!(
                    degrees
                        .iter()
                        .map(|degree| degree.to_string().parse::<i64>().unwrap())
                        .collect::<Vec<_>>()
                )
            },
        )
    }

    fn properties_json<F>(complex: &ModuleComplex<F>) -> Value
    where
        F: QpaCoefficientField + Send + Sync + 'static,
        F::Element: Send + Sync,
    {
        let finiteness = complex.is_finite_complex().unwrap();
        let exact = if finiteness.is_none() {
            json!("property_error")
        } else {
            json!(complex.is_exact_sequence().unwrap())
        };
        let short_exact = if finiteness.is_none() {
            json!("property_error")
        } else {
            json!(complex.is_short_exact_sequence().unwrap())
        };
        json!({
            "upper": bound_string(complex.upper_bound().unwrap()),
            "lower": bound_string(complex.lower_bound().unwrap()),
            "finite": finiteness.map_or_else(|| json!("fail"), |value| json!(value)),
            "length": bound_string(complex.length().unwrap()),
            "exact": exact,
            "short_exact": short_exact,
            "positive_repeating": complex.is_positive_repeating().unwrap(),
            "negative_repeating": complex.is_negative_repeating().unwrap(),
            "positive_repeat_degrees": repeat_json(
                complex.positive_repeat_degrees().unwrap()),
            "negative_repeat_degrees": repeat_json(
                complex.negative_repeat_degrees().unwrap()),
            "highest_known": known_string(complex.highest_known_degree().unwrap()),
            "lowest_known": known_string(complex.lowest_known_degree().unwrap()),
        })
    }

    fn module_map<F: Field>(
        source: &PathAlgebraModule<F>,
        range: &PathAlgebraModule<F>,
        rows: usize,
        columns: usize,
        entries: Vec<F::Element>,
    ) -> PathAlgebraModuleHomomorphism<F> {
        right_module_hom_over_algebra(
            source.clone(),
            range.clone(),
            vec![
                DenseMatrix::from_flat(
                    rows,
                    columns,
                    entries,
                    source.right_acting_algebra().field().clone(),
                )
                .unwrap(),
            ],
        )
        .unwrap()
    }

    fn complex_case<F>(
        field: F,
        coefficient: impl Fn(&F::Element) -> Value + Copy + Send + Sync + 'static,
    ) -> Value
    where
        F: QpaCoefficientField + Send + Sync + 'static,
        F::Element: Send + Sync,
    {
        let algebra = QuiverAlgebra::Path(
            PathAlgebra::new(
                field.clone(),
                Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap(),
            )
            .unwrap(),
        );
        let category = cat_of_right_algebra_modules(&algebra).unwrap();
        let one = right_module_over_path_algebra(algebra.clone(), vec![1], &[])
            .unwrap()
            .unwrap();
        let two = right_module_over_path_algebra(algebra.clone(), vec![2], &[])
            .unwrap()
            .unwrap();
        let f = module_map(&one, &two, 1, 2, vec![field.one(), field.zero()]);
        let g = module_map(&two, &one, 2, 1, vec![field.zero(), field.one()]);
        let zero_endo = zero_module_mapping(one.clone(), one.clone()).unwrap();
        let ses = short_exact_sequence(category.clone(), f, g).unwrap();
        let stalk = stalk_complex(category.clone(), &one, Integer::from(2)).unwrap();
        let periodic = complex(
            category.clone(),
            Integer::from(0),
            vec![zero_endo.clone()],
            ComplexHalf::Repeat(vec![zero_endo.clone()]),
            ComplexHalf::Repeat(vec![zero_endo.clone()]),
        )
        .unwrap();

        let calls = Arc::new(AtomicUsize::new(0));
        let positive_calls = Arc::clone(&calls);
        let positive_zero = zero_endo.clone();
        let negative_calls = Arc::clone(&calls);
        let negative_zero = zero_endo.clone();
        let positioned = complex(
            category.clone(),
            Integer::from(0),
            vec![zero_endo],
            ComplexHalf::Position {
                function: Arc::new(move |_complex, _degree| {
                    positive_calls.fetch_add(1, Ordering::SeqCst);
                    Ok(positive_zero.clone())
                }),
                store_values: false,
            },
            ComplexHalf::Position {
                function: Arc::new(move |_complex, _degree| {
                    negative_calls.fetch_add(1, Ordering::SeqCst);
                    Ok(negative_zero.clone())
                }),
                store_values: false,
            },
        )
        .unwrap();
        positioned.differential(&Integer::from(5)).unwrap();
        positioned.differential(&Integer::from(-5)).unwrap();

        let shift_one = ses.shift(&Integer::from(1)).unwrap();
        let shift_minus_one = ses.shift(&Integer::from(-1)).unwrap();
        let brutal_below = ses.brutal_truncation_below(&Integer::from(0)).unwrap();
        let brutal_above = ses.brutal_truncation_above(&Integer::from(0)).unwrap();
        let brutal = ses
            .brutal_truncation(&Integer::from(0), &Integer::from(0))
            .unwrap();
        let good_below = ses.good_truncation_below(&Integer::from(0)).unwrap();
        let good_above = ses.good_truncation_above(&Integer::from(0)).unwrap();
        let good = ses
            .good_truncation(&Integer::from(0), &Integer::from(0))
            .unwrap();
        let syzygy = ses.syzygy_truncation(&Integer::from(0)).unwrap();
        let cosyzygy = ses.cosyzygy_truncation(&Integer::from(0)).unwrap();
        let both = ses
            .syzygy_cosyzygy_truncation(&Integer::from(0), &Integer::from(0))
            .unwrap();
        let yoneda = yoneda_product(&ses, &ses).unwrap();
        let identity_chain = finite_chain_map(
            ses.clone(),
            ses.clone(),
            Integer::from(-1),
            (-1..=1)
                .map(|degree| {
                    identity_module_mapping(ses.object(&Integer::from(degree)).unwrap()).unwrap()
                })
                .collect(),
        )
        .unwrap();
        let zero_chain = zero_chain_map(ses.clone(), ses.clone()).unwrap();
        let cone = mapping_cone(&identity_chain).unwrap();
        let comparison = comparison_lifting_to_projective_resolution(
            &identity_module_mapping(one.clone()).unwrap(),
        )
        .unwrap();
        let cut_above = projective_resolution_complex(&one)
            .unwrap()
            .cut_complex_above()
            .unwrap();
        let cut_below = injective_resolution_complex(&one)
            .unwrap()
            .cut_complex_below()
            .unwrap();
        let degrees = [-2, -1, 0, 1, 2, 3];
        let cone_degrees = [-3, -2, -1, 0, 1, 2, 3];

        json!({
            "category_equal": category == cat_of_right_algebra_modules(&algebra).unwrap(),
            "zero": complex_sample(&zero_complex(category.clone()).unwrap(), &[-1,0,1], coefficient),
            "ses": complex_sample(&ses, &degrees, coefficient),
            "ses_homology": (-2..=2).map(|degree| ses
                .homology(&Integer::from(degree)).unwrap().dimension_vector().to_vec())
                .collect::<Vec<_>>(),
            "stalk": complex_sample(&stalk, &[0,1,2,3,4], coefficient),
            "periodic": complex_sample(&periodic, &[-2,-1,0,1,2], coefficient),
            "position_calls": calls.load(Ordering::SeqCst),
            "position_known": [
                known_string(positioned.lowest_known_degree().unwrap()),
                known_string(positioned.highest_known_degree().unwrap())
            ],
            "shift_one": complex_sample(&shift_one, &degrees, coefficient),
            "shift_minus_one": complex_sample(&shift_minus_one, &degrees, coefficient),
            "brutal_below": complex_sample(&brutal_below, &degrees, coefficient),
            "brutal_above": complex_sample(&brutal_above, &degrees, coefficient),
            "brutal": complex_sample(&brutal, &degrees, coefficient),
            "good_below": complex_sample(&good_below, &degrees, coefficient),
            "good_above": complex_sample(&good_above, &degrees, coefficient),
            "good": complex_sample(&good, &degrees, coefficient),
            "syzygy": complex_sample(&syzygy, &degrees, coefficient),
            "cosyzygy": complex_sample(&cosyzygy, &degrees, coefficient),
            "syzygy_cosyzygy": complex_sample(&both, &degrees, coefficient),
            "yoneda": complex_sample(&yoneda, &cone_degrees, coefficient),
            "identity_chain": chain_map_sample(&identity_chain, &degrees, coefficient),
            "zero_chain": chain_map_sample(&zero_chain, &degrees, coefficient),
            "cone": complex_sample(&cone.complex, &cone_degrees, coefficient),
            "cone_inclusion": chain_map_sample(&cone.incoming[0], &cone_degrees, coefficient),
            "cone_projection": chain_map_sample(&cone.outgoing[0], &cone_degrees, coefficient),
            "comparison": chain_map_sample(&comparison, &[-2,-1,0,1,2], coefficient),
            "cut_above": complex_sample(&cut_above, &[-2,-1,0,1,2], coefficient),
            "cut_below": complex_sample(&cut_below, &[-2,-1,0,1,2], coefficient),
            "properties": {
                "zero": properties_json(&zero_complex(category).unwrap()),
                "ses": properties_json(&ses),
                "stalk": properties_json(&stalk),
                "periodic": properties_json(&periodic),
                "good": properties_json(&good),
                "yoneda": properties_json(&yoneda),
            }
        })
    }

    #[test]
    fn pinned_qpa_complex_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/complex.schema-v1.json"
        ))
        .unwrap();
        let actual = json!({
            "schema_version": 1,
            "kind": "complex",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            },
            "binary": complex_case(Z2, |entry| json!(
                Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
            )),
            "rational": complex_case(Q, |entry: &Rational| json!([
                entry.numerator_ref().to_string().parse::<i64>().unwrap(),
                entry.denominator_ref().to_string().parse::<i64>().unwrap(),
            ])),
        });
        assert_eq!(actual, expected);
    }

    #[test]
    fn qpa_complex_printing_surface_matches() {
        let algebra = QuiverAlgebra::Path(
            PathAlgebra::new(Q, Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap()).unwrap(),
        );
        let category = cat_of_right_algebra_modules(&algebra).unwrap();
        let one = right_module_over_path_algebra(algebra.clone(), vec![1], &[])
            .unwrap()
            .unwrap();
        let two = right_module_over_path_algebra(algebra, vec![2], &[])
            .unwrap()
            .unwrap();
        let f = module_map(&one, &two, 1, 2, vec![1.into(), 0.into()]);
        let g = module_map(&two, &one, 2, 1, vec![0.into(), 1.into()]);
        let ses = short_exact_sequence(category.clone(), f, g).unwrap();
        let zero_endo = zero_module_mapping(one.clone(), one).unwrap();
        let periodic = complex(
            category.clone(),
            Integer::from(0),
            vec![zero_endo.clone()],
            ComplexHalf::Repeat(vec![zero_endo.clone()]),
            ComplexHalf::Repeat(vec![zero_endo]),
        )
        .unwrap();
        let chain = finite_chain_map(
            ses.clone(),
            ses.clone(),
            Integer::from(-1),
            (-1..=1)
                .map(|degree| {
                    identity_module_mapping(ses.object(&Integer::from(degree)).unwrap()).unwrap()
                })
                .collect(),
        )
        .unwrap();
        let actual = format!(
            "<cat: {}>\n{}\n{}\n{}\n{}\n",
            category.name(),
            zero_complex(category.clone())
                .unwrap()
                .qpa_string()
                .unwrap(),
            ses.qpa_string().unwrap(),
            periodic.qpa_string().unwrap(),
            chain.qpa_string()
        );
        assert_eq!(
            actual,
            include_str!("../../../tests/fixtures/complex_print.txt")
        );
    }

    #[test]
    fn source_failures_remain_typed() {
        let algebra = QuiverAlgebra::Path(
            PathAlgebra::new(Q, Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap()).unwrap(),
        );
        let category = cat_of_right_algebra_modules(&algebra).unwrap();
        assert!(matches!(
            finite_complex(category.clone(), Integer::from(0), Vec::new()),
            Err(ComplexError::EmptyMiddleWithZeroHalf { side: "positive" })
        ));
        let one = right_module_over_path_algebra(algebra, vec![1], &[])
            .unwrap()
            .unwrap();
        let identity = identity_module_mapping(one.clone()).unwrap();
        assert!(matches!(
            finite_complex(
                category.clone(),
                Integer::from(0),
                vec![identity.clone(), identity]
            ),
            Err(ComplexError::NonzeroDifferentialComposition { .. })
        ));
        let zero = zero_module_mapping(one.clone(), one).unwrap();
        let invalid_short = short_exact_sequence(category, zero.clone(), zero);
        assert!(matches!(invalid_short, Err(ComplexError::NotShortExact)));
        let fixture: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/complex_errors.schema-v1.json"
        ))
        .unwrap();
        assert!(
            fixture["success"]
                .as_object()
                .unwrap()
                .values()
                .all(|value| value == &json!(false))
        );
    }
}
