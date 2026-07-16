//! Matrix representations of QPA right modules over quiver algebras.
//!
//! The checked constructor in this module translates the three-argument
//! `RightModuleOverPathAlgebra` method from `lib/module.gi`. In particular,
//! arrow matrices act on row vectors, omitted arrows receive GAP's physical
//! `1`-dimensional stand-ins at zero-dimensional vertices, and every defining
//! relation of a quotient is evaluated before a module is returned.

use std::{
    fmt,
    sync::{Arc, OnceLock, RwLock, Weak},
};

use numerica::domains::Field;
use thiserror::Error;

use crate::{
    DenseMatrix, Generator, MatrixError, Path, PathEndpoint, QuiverAlgebra, QuiverAlgebraElement,
    QuiverAlgebraElementError, QuotientError, TensorProductAlgebra,
};

/// One labeled arrow matrix in QPA's three-argument module constructor.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct NamedArrowMatrix<F: Field> {
    label: Box<str>,
    matrix: DenseMatrix<F>,
}

impl<F: Field> NamedArrowMatrix<F> {
    /// Pair a quiver arrow's string label with its proposed right action.
    #[must_use]
    pub fn new(label: impl Into<Box<str>>, matrix: DenseMatrix<F>) -> Self {
        Self {
            label: label.into(),
            matrix,
        }
    }

    /// Return the entered arrow label.
    #[must_use]
    pub fn label(&self) -> &str {
        &self.label
    }

    /// Return the entered matrix.
    #[must_use]
    pub const fn matrix(&self) -> &DenseMatrix<F> {
        &self.matrix
    }
}

/// One entry returned by QPA's `MatricesOfPathAlgebraModule` operation.
///
/// GAP returns the integer zero, rather than a physical matrix, for every
/// arrow of the zero module. The enum preserves that observable distinction.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum PathAlgebraModuleMatrix<F: Field> {
    /// GAP's scalar `0` entry for an arrow of the zero module.
    ZeroScalar,
    /// An arrow-action matrix for a nonzero module.
    Matrix(DenseMatrix<F>),
}

#[derive(Debug)]
struct PathAlgebraModuleData<F: Field> {
    acting_algebra: QuiverAlgebra<F>,
    dimension_vector: Vec<usize>,
    dimension: usize,
    arrow_matrices: Vec<DenseMatrix<F>>,
    retained_tensor: Option<TensorProductAlgebra<F>>,
    dual_weak: RwLock<Option<WeakPathAlgebraModule<F>>>,
    dual_strong: OnceLock<PathAlgebraModule<F>>,
    endomorphism_algebra: OnceLock<ModuleEndomorphismAlgebra<F>>,
    top_module: OnceLock<PathAlgebraModule<F>>,
    top_projection_weak: RwLock<Option<WeakPathAlgebraModuleHomomorphism<F>>>,
    direct_sum: OnceLock<DirectSumMetadata<F>>,
    minimal_generators: OnceLock<Vec<Vec<Vec<F::Element>>>>,
    is_projective: OnceLock<bool>,
    is_injective: OnceLock<bool>,
    is_indecomposable: OnceLock<bool>,
    is_simple_qpa: OnceLock<bool>,
    is_semisimple: OnceLock<bool>,
    projective_cover_weak: RwLock<Option<WeakPathAlgebraModuleHomomorphism<F>>>,
    injective_envelope_weak: RwLock<Option<WeakPathAlgebraModuleHomomorphism<F>>>,
    first_syzygy: OnceLock<PathAlgebraModule<F>>,
}

impl<F: Field> PathAlgebraModuleData<F> {
    fn new(
        acting_algebra: QuiverAlgebra<F>,
        dimension_vector: Vec<usize>,
        dimension: usize,
        arrow_matrices: Vec<DenseMatrix<F>>,
        retained_tensor: Option<TensorProductAlgebra<F>>,
    ) -> Self {
        Self {
            acting_algebra,
            dimension_vector,
            dimension,
            arrow_matrices,
            retained_tensor,
            dual_weak: RwLock::new(None),
            dual_strong: OnceLock::new(),
            endomorphism_algebra: OnceLock::new(),
            top_module: OnceLock::new(),
            top_projection_weak: RwLock::new(None),
            direct_sum: OnceLock::new(),
            minimal_generators: OnceLock::new(),
            is_projective: OnceLock::new(),
            is_injective: OnceLock::new(),
            is_indecomposable: OnceLock::new(),
            is_simple_qpa: OnceLock::new(),
            is_semisimple: OnceLock::new(),
            projective_cover_weak: RwLock::new(None),
            injective_envelope_weak: RwLock::new(None),
            first_syzygy: OnceLock::new(),
        }
    }
}

#[derive(Debug)]
struct DirectSumMetadata<F: Field> {
    summands: Vec<PathAlgebraModule<F>>,
    projections: RwLock<Vec<WeakPathAlgebraModuleHomomorphism<F>>>,
    inclusions: RwLock<Vec<WeakPathAlgebraModuleHomomorphism<F>>>,
}

/// A QPA right path-algebra matrix module.
///
/// Clones preserve GAP object identity. Structural equality follows QPA's
/// method: acting algebra, dimension vector, and reported arrow matrices must
/// agree.
#[derive(Clone)]
pub struct PathAlgebraModule<F: Field>(Arc<PathAlgebraModuleData<F>>);

impl<F: Field> PathAlgebraModule<F> {
    /// Return QPA's `RightActingAlgebra` value.
    #[must_use]
    pub fn right_acting_algebra(&self) -> &QuiverAlgebra<F> {
        &self.0.acting_algebra
    }

    /// Return QPA's cached `DimensionVector` value.
    #[must_use]
    pub fn dimension_vector(&self) -> &[usize] {
        &self.0.dimension_vector
    }

    /// Return the total vector-space dimension.
    #[must_use]
    pub fn dimension(&self) -> usize {
        self.0.dimension
    }

    /// Return QPA's additive zero module element.
    #[must_use]
    pub fn zero_element(&self) -> PathAlgebraModuleElement<F> {
        PathAlgebraModuleElement::zero(self.clone())
    }

    /// Return `BasisVectors(Basis(M))` in QPA's vertex-major order.
    #[must_use]
    pub fn basis_elements(&self) -> Vec<PathAlgebraModuleElement<F>> {
        let field = self.right_acting_algebra().field();
        let mut basis = Vec::with_capacity(self.dimension());
        for (vertex, &dimension) in self.dimension_vector().iter().enumerate() {
            for coordinate in 0..dimension {
                let mut element = self.zero_element();
                element.coordinates[vertex][coordinate] = field.one();
                basis.push(element);
            }
        }
        basis
    }

    /// Construct an element from one logical coordinate row at each vertex.
    ///
    /// Zero-dimensional vertices are entered as empty rows; internally they
    /// receive QPA's one-coordinate zero stand-in.
    ///
    /// # Errors
    ///
    /// Rejects a wrong vertex count or a row whose length differs from the
    /// corresponding dimension-vector entry.
    pub fn element_from_vertex_coordinates(
        &self,
        coordinates: Vec<Vec<F::Element>>,
    ) -> Result<PathAlgebraModuleElement<F>, PathAlgebraModuleError> {
        if coordinates.len() != self.dimension_vector().len() {
            return Err(PathAlgebraModuleError::ElementVertexCount {
                actual: coordinates.len(),
                expected: self.dimension_vector().len(),
            });
        }
        let mut physical = Vec::with_capacity(coordinates.len());
        for (vertex, (row, &dimension)) in coordinates
            .into_iter()
            .zip(self.dimension_vector())
            .enumerate()
        {
            if row.len() != dimension {
                return Err(PathAlgebraModuleError::ElementCoordinateCount {
                    vertex: vertex + 1,
                    actual: row.len(),
                    expected: dimension,
                });
            }
            physical.push(if dimension == 0 {
                vec![self.right_acting_algebra().field().zero()]
            } else {
                row
            });
        }
        Ok(PathAlgebraModuleElement {
            module: self.clone(),
            coordinates: physical,
        })
    }

    /// Return arrow actions in quiver arrow order.
    ///
    /// These are the physical matrices retained by QPA's element family. At a
    /// zero-dimensional endpoint, an omitted arrow uses a one-row or
    /// one-column stand-in so relation products remain representable in GAP.
    #[must_use]
    pub fn arrow_matrices(&self) -> &[DenseMatrix<F>] {
        &self.0.arrow_matrices
    }

    /// Translate `MatricesOfPathAlgebraModule` exactly, including its scalar
    /// zero entries for a total-zero module.
    #[must_use]
    pub fn matrices_of_path_algebra_module(&self) -> Vec<PathAlgebraModuleMatrix<F>> {
        if self.dimension() == 0 {
            vec![PathAlgebraModuleMatrix::ZeroScalar; self.0.arrow_matrices.len()]
        } else {
            self.0
                .arrow_matrices
                .iter()
                .cloned()
                .map(PathAlgebraModuleMatrix::Matrix)
                .collect()
        }
    }

    /// Return whether two handles denote the identical GAP module object.
    #[must_use]
    pub fn is_identical(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }

    /// Return QPA's installed `IsDirectSumOfModules` filter.
    #[must_use]
    pub fn is_direct_sum_of_modules(&self) -> bool {
        self.0.direct_sum.get().is_some()
    }

    /// Return the ordered summands retained by `DirectSumOfQPAModules`.
    #[must_use]
    pub fn direct_sum_summands(&self) -> Option<&[Self]> {
        self.0
            .direct_sum
            .get()
            .map(|metadata| metadata.summands.as_slice())
    }

    /// Return QPA's cached `IsProjectiveModule` value, if installed.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_is_projective_module(&self) -> Option<bool> {
        self.0.is_projective.get().copied()
    }

    /// Install QPA's `IsProjectiveModule` property value.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_is_projective_module(&self, value: bool) -> bool {
        *self.0.is_projective.get_or_init(|| value)
    }

    /// Return QPA's cached `IsInjectiveModule` value, if installed.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_is_injective_module(&self) -> Option<bool> {
        self.0.is_injective.get().copied()
    }

    /// Install QPA's `IsInjectiveModule` property value.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_is_injective_module(&self, value: bool) -> bool {
        *self.0.is_injective.get_or_init(|| value)
    }

    /// Return QPA's cached `IsIndecomposableModule` value, if installed.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_is_indecomposable_module(&self) -> Option<bool> {
        self.0.is_indecomposable.get().copied()
    }

    /// Install QPA's `IsIndecomposableModule` property value.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_is_indecomposable_module(&self, value: bool) -> bool {
        *self.0.is_indecomposable.get_or_init(|| value)
    }

    /// Return QPA's cached `IsSimpleQPAModule` value, if installed.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_is_simple_qpa_module(&self) -> Option<bool> {
        self.0.is_simple_qpa.get().copied()
    }

    /// Install QPA's `IsSimpleQPAModule` property and its true-method
    /// implications to indecomposability and semisimplicity.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_is_simple_qpa_module(&self, value: bool) -> bool {
        let installed = *self.0.is_simple_qpa.get_or_init(|| value);
        if installed {
            let _ = self.cache_is_indecomposable_module(true);
            let _ = self.cache_is_semisimple_module(true);
        }
        installed
    }

    /// Return QPA's cached `IsSemisimpleModule` value, if installed.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_is_semisimple_module(&self) -> Option<bool> {
        self.0.is_semisimple.get().copied()
    }

    /// Install QPA's `IsSemisimpleModule` property value.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_is_semisimple_module(&self, value: bool) -> bool {
        *self.0.is_semisimple.get_or_init(|| value)
    }

    /// Return QPA's installed minimal-generator attribute as module elements.
    #[doc(hidden)]
    pub fn cached_minimal_generating_set(
        &self,
    ) -> Result<Option<Vec<PathAlgebraModuleElement<F>>>, PathAlgebraModuleError> {
        self.0
            .minimal_generators
            .get()
            .map(|generators| {
                generators
                    .iter()
                    .map(|coordinates| self.element_from_vertex_coordinates(coordinates.clone()))
                    .collect()
            })
            .transpose()
    }

    /// Install QPA's minimal-generator attribute without retaining an `Arc`
    /// cycle through its elements.
    #[doc(hidden)]
    pub fn cache_minimal_generating_set(
        &self,
        generators: &[PathAlgebraModuleElement<F>],
    ) -> Result<Vec<PathAlgebraModuleElement<F>>, PathAlgebraModuleError> {
        let mut coordinates = Vec::with_capacity(generators.len());
        for generator in generators {
            if !generator.module().is_identical(self) {
                return Err(PathAlgebraModuleError::MinimalGeneratorInWrongModule);
            }
            coordinates.push(
                generator
                    .external_representation()
                    .iter()
                    .zip(self.dimension_vector())
                    .map(|(row, &dimension)| row[..dimension].to_vec())
                    .collect(),
            );
        }
        let _ = self.0.minimal_generators.set(coordinates);
        self.cached_minimal_generating_set()?
            .ok_or(PathAlgebraModuleError::MinimalGeneratingSetCacheFailed)
    }

    /// Return an already-live cached `ProjectiveCover` attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_projective_cover(&self) -> Option<PathAlgebraModuleHomomorphism<F>> {
        self.0
            .projective_cover_weak
            .read()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .as_ref()
            .and_then(WeakPathAlgebraModuleHomomorphism::upgrade)
    }

    /// Install QPA's `ProjectiveCover` attribute without creating an `Arc`
    /// cycle from the range module through the map back to itself.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_projective_cover(
        &self,
        cover: &PathAlgebraModuleHomomorphism<F>,
    ) -> PathAlgebraModuleHomomorphism<F> {
        let mut slot = self
            .0
            .projective_cover_weak
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if let Some(cached) = slot
            .as_ref()
            .and_then(WeakPathAlgebraModuleHomomorphism::upgrade)
        {
            return cached;
        }
        *slot = Some(cover.downgrade());
        cover.clone()
    }

    /// Return an already-live cached `InjectiveEnvelope` attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_injective_envelope(&self) -> Option<PathAlgebraModuleHomomorphism<F>> {
        self.0
            .injective_envelope_weak
            .read()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .as_ref()
            .and_then(WeakPathAlgebraModuleHomomorphism::upgrade)
    }

    /// Install QPA's `InjectiveEnvelope` attribute without retaining the map's
    /// source module through an `Arc` cycle.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_injective_envelope(
        &self,
        envelope: &PathAlgebraModuleHomomorphism<F>,
    ) -> PathAlgebraModuleHomomorphism<F> {
        let mut slot = self
            .0
            .injective_envelope_weak
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if let Some(cached) = slot
            .as_ref()
            .and_then(WeakPathAlgebraModuleHomomorphism::upgrade)
        {
            return cached;
        }
        *slot = Some(envelope.downgrade());
        envelope.clone()
    }

    /// Return QPA's cached `1stSyzygy` attribute, if installed.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_first_syzygy(&self) -> Option<Self> {
        self.0.first_syzygy.get().cloned()
    }

    /// Install and return QPA's `1stSyzygy` attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_first_syzygy(&self, syzygy: &Self) -> Self {
        self.0.first_syzygy.get_or_init(|| syzygy.clone()).clone()
    }

    /// Retain tensor-product metadata when an attribute-created module acts
    /// over an enveloping algebra. This keeps that cached GAP algebra object
    /// alive for the lifetime of the module.
    #[doc(hidden)]
    #[must_use]
    pub fn retained_tensor_product(&self) -> Option<&TensorProductAlgebra<F>> {
        self.0.retained_tensor.as_ref()
    }

    /// Return an already-live cached `DualOfModule` attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_dual_module(&self) -> Option<Self> {
        self.0.dual_strong.get().cloned().or_else(|| {
            self.0
                .dual_weak
                .read()
                .unwrap_or_else(std::sync::PoisonError::into_inner)
                .as_ref()
                .and_then(WeakPathAlgebraModule::upgrade)
        })
    }

    /// Install only this module's `DualOfModule` cache direction.
    ///
    /// QPA's zero-module branch returns before setting the reciprocal
    /// attribute on the result.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_dual_module(&self, dual: &Self) -> Self {
        let mut slot = self
            .0
            .dual_weak
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if let Some(cached) = slot.as_ref().and_then(WeakPathAlgebraModule::upgrade) {
            return cached;
        }
        *slot = Some(dual.downgrade());
        dual.clone()
    }

    /// Install QPA's reciprocal `DualOfModule` attribute pair.
    ///
    /// # Errors
    ///
    /// Reports a concurrent reciprocal-cache initialization race.
    #[doc(hidden)]
    pub fn install_dual_module(&self, dual: &Self) -> Result<Self, PathAlgebraModuleError> {
        if let Some(cached) = self.cached_dual_module() {
            return Ok(cached);
        }
        if dual.0.dual_strong.set(self.clone()).is_err() {
            return Err(PathAlgebraModuleError::DualCacheRace);
        }
        Ok(self.cache_dual_module(dual))
    }

    pub(crate) fn downgrade(&self) -> WeakPathAlgebraModule<F> {
        WeakPathAlgebraModule(Arc::downgrade(&self.0))
    }

    fn cached_endomorphism_algebra(&self) -> Option<ModuleEndomorphismAlgebra<F>> {
        self.0.endomorphism_algebra.get().cloned()
    }

    fn cache_endomorphism_algebra(
        &self,
        algebra: &ModuleEndomorphismAlgebra<F>,
    ) -> ModuleEndomorphismAlgebra<F> {
        self.0
            .endomorphism_algebra
            .get_or_init(|| algebra.clone())
            .clone()
    }

    fn cached_top_module(&self) -> Option<Self> {
        self.0.top_module.get().cloned()
    }

    fn cache_top_module(&self, top: &Self) -> Self {
        self.0.top_module.get_or_init(|| top.clone()).clone()
    }

    fn cached_top_projection(&self) -> Option<PathAlgebraModuleHomomorphism<F>> {
        self.0
            .top_projection_weak
            .read()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .as_ref()
            .and_then(WeakPathAlgebraModuleHomomorphism::upgrade)
    }

    fn cache_top_projection(
        &self,
        projection: &PathAlgebraModuleHomomorphism<F>,
    ) -> PathAlgebraModuleHomomorphism<F> {
        let mut slot = self
            .0
            .top_projection_weak
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if let Some(cached) = slot
            .as_ref()
            .and_then(WeakPathAlgebraModuleHomomorphism::upgrade)
        {
            return cached;
        }
        *slot = Some(projection.downgrade());
        projection.clone()
    }
}

/// An element of a QPA right path-algebra matrix module.
///
/// Coordinates are stored vertexwise.  A zero-dimensional vertex retains the
/// single zero coordinate created by QPA's `ZeroModElement`; logical basis and
/// coefficient lists omit those physical stand-ins.
#[derive(Clone, Debug)]
pub struct PathAlgebraModuleElement<F: Field> {
    module: PathAlgebraModule<F>,
    coordinates: Vec<Vec<F::Element>>,
}

impl<F: Field> PathAlgebraModuleElement<F> {
    fn zero(module: PathAlgebraModule<F>) -> Self {
        let zero = module.right_acting_algebra().field().zero();
        let coordinates = module
            .dimension_vector()
            .iter()
            .map(|&dimension| vec![zero.clone(); physical_dimension(dimension)])
            .collect();
        Self {
            module,
            coordinates,
        }
    }

    fn from_logical_flat(
        module: PathAlgebraModule<F>,
        coordinates: &[F::Element],
    ) -> Result<Self, PathAlgebraModuleError> {
        if coordinates.len() != module.dimension() {
            return Err(PathAlgebraModuleError::ElementFlatCoordinateCount {
                actual: coordinates.len(),
                expected: module.dimension(),
            });
        }
        let mut offset = 0;
        let rows = module
            .dimension_vector()
            .iter()
            .map(|&dimension| {
                if dimension == 0 {
                    vec![module.right_acting_algebra().field().zero()]
                } else {
                    let row = coordinates[offset..offset + dimension].to_vec();
                    offset += dimension;
                    row
                }
            })
            .collect();
        Ok(Self {
            module,
            coordinates: rows,
        })
    }

    /// Return the identical module family containing this element.
    #[must_use]
    pub const fn module(&self) -> &PathAlgebraModule<F> {
        &self.module
    }

    /// Return QPA's vertexwise `ExtRepOfObj` coordinates, including physical
    /// zero-space stand-ins.
    #[must_use]
    pub fn external_representation(&self) -> &[Vec<F::Element>] {
        &self.coordinates
    }

    /// Flatten the physical external representation exactly as GAP `Flat`.
    #[must_use]
    pub fn flat_external_representation(&self) -> Vec<F::Element> {
        self.coordinates.iter().flatten().cloned().collect()
    }

    /// Return coefficients in QPA's standard module basis.
    #[must_use]
    pub fn basis_coefficients(&self) -> Vec<F::Element> {
        self.coordinates
            .iter()
            .zip(self.module.dimension_vector())
            .flat_map(|(row, &dimension)| row[..dimension].iter().cloned())
            .collect()
    }

    /// Translate `SupportModuleElement` as one-based vertex positions.
    #[must_use]
    pub fn support_vertices(&self) -> Vec<usize> {
        let field = self.module.right_acting_algebra().field();
        self.coordinates
            .iter()
            .zip(self.module.dimension_vector())
            .enumerate()
            .filter_map(|(vertex, (row, &dimension))| {
                row[..dimension]
                    .iter()
                    .any(|entry| !field.is_zero(entry))
                    .then_some(vertex + 1)
            })
            .collect()
    }

    /// Return whether this is the additive zero.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.support_vertices().is_empty()
    }

    /// Add elements from the identical QPA element family.
    ///
    /// # Errors
    ///
    /// Rejects elements belonging to different module objects.
    pub fn add(&self, right: &Self) -> Result<Self, PathAlgebraModuleError> {
        self.ensure_same_module(right)?;
        let field = self.module.right_acting_algebra().field();
        let coordinates = self
            .coordinates
            .iter()
            .zip(&right.coordinates)
            .map(|(left, right)| {
                left.iter()
                    .zip(right)
                    .map(|(left, right)| field.add(left, right))
                    .collect()
            })
            .collect();
        Ok(Self {
            module: self.module.clone(),
            coordinates,
        })
    }

    /// Multiply this element by a field scalar.
    #[must_use]
    pub fn scale(&self, scalar: &F::Element) -> Self {
        let field = self.module.right_acting_algebra().field();
        Self {
            module: self.module.clone(),
            coordinates: self
                .coordinates
                .iter()
                .map(|row| row.iter().map(|entry| field.mul(entry, scalar)).collect())
                .collect(),
        }
    }

    /// Apply QPA's right action by a path-algebra or quotient element.
    ///
    /// # Errors
    ///
    /// Rejects an element outside the identical acting-algebra family and
    /// propagates exact matrix failures.
    pub fn act(&self, element: &QuiverAlgebraElement<F>) -> Result<Self, PathAlgebraModuleError> {
        let algebra = self.module.right_acting_algebra();
        if !algebra.contains_element(element) {
            return Err(PathAlgebraModuleError::ForeignActingElement);
        }
        let terms = match element {
            QuiverAlgebraElement::Path(element) => element.terms(),
            QuiverAlgebraElement::Quotient(element) => element.representative().terms(),
        };
        let mut result = self.module.zero_element();
        for term in terms {
            let PathEndpoint::Vertex(source) = term.path().source() else {
                return Err(PathAlgebraModuleError::ZeroPathAction);
            };
            let PathEndpoint::Vertex(target) = term.path().target() else {
                return Err(PathAlgebraModuleError::ZeroPathAction);
            };
            let source_index = source.index().zero_based();
            let target_index = target.index().zero_based();
            let action = path_action_matrix(
                term.path(),
                self.module.dimension_vector(),
                self.module.arrow_matrices(),
                algebra.field().clone(),
            )?;
            let row = DenseMatrix::from_flat(
                1,
                self.coordinates[source_index].len(),
                self.coordinates[source_index].clone(),
                algebra.field().clone(),
            )?;
            let image = row.multiply(&action)?.scale_left(term.coefficient());
            for (entry, value) in result.coordinates[target_index]
                .iter_mut()
                .zip(image.row(0))
            {
                algebra.field().add_assign(entry, value);
            }
        }
        Ok(result)
    }

    fn ensure_same_module(&self, right: &Self) -> Result<(), PathAlgebraModuleError> {
        if self.module.is_identical(&right.module) {
            Ok(())
        } else {
            Err(PathAlgebraModuleError::DifferentElementFamilies)
        }
    }
}

impl<F: Field> PartialEq for PathAlgebraModuleElement<F> {
    fn eq(&self, other: &Self) -> bool {
        self.module.is_identical(&other.module) && self.coordinates == other.coordinates
    }
}

impl<F: Field> Eq for PathAlgebraModuleElement<F> {}

impl<F: Field> fmt::Debug for PathAlgebraModule<F> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_tuple("PathAlgebraModule")
            .field(&self.dimension_vector())
            .finish()
    }
}

impl<F: Field> fmt::Display for PathAlgebraModule<F> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "<{:?}>", self.dimension_vector())
    }
}

impl<F: Field> PartialEq for PathAlgebraModule<F> {
    fn eq(&self, other: &Self) -> bool {
        self.right_acting_algebra() == other.right_acting_algebra()
            && self.dimension_vector() == other.dimension_vector()
            && self.matrices_of_path_algebra_module() == other.matrices_of_path_algebra_module()
    }
}

impl<F: Field> Eq for PathAlgebraModule<F> {}

#[derive(Clone, Debug)]
pub(crate) struct WeakPathAlgebraModule<F: Field>(Weak<PathAlgebraModuleData<F>>);

impl<F: Field> WeakPathAlgebraModule<F> {
    pub(crate) fn upgrade(&self) -> Option<PathAlgebraModule<F>> {
        self.0.upgrade().map(PathAlgebraModule)
    }
}

/// A cycle-free cache for a GAP algebra attribute whose value is a list of
/// modules over that same algebra.
#[derive(Debug)]
pub(crate) struct WeakPathAlgebraModuleListCache<F: Field> {
    modules: RwLock<Vec<WeakPathAlgebraModule<F>>>,
    bound: OnceLock<()>,
}

impl<F: Field> WeakPathAlgebraModuleListCache<F> {
    pub(crate) const fn new() -> Self {
        Self {
            modules: RwLock::new(Vec::new()),
            bound: OnceLock::new(),
        }
    }

    pub(crate) fn cached(&self) -> Option<Vec<PathAlgebraModule<F>>> {
        self.bound.get()?;
        self.modules
            .read()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .iter()
            .map(WeakPathAlgebraModule::upgrade)
            .collect()
    }

    pub(crate) fn merge(&self, generated: &[PathAlgebraModule<F>]) -> Vec<PathAlgebraModule<F>> {
        let mut slot = self
            .modules
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if self.bound.get().is_none() || slot.len() != generated.len() {
            *slot = generated.iter().map(PathAlgebraModule::downgrade).collect();
            let _ = self.bound.set(());
            return generated.to_vec();
        }
        generated
            .iter()
            .enumerate()
            .map(|(index, module)| {
                if let Some(cached) = slot[index].upgrade() {
                    cached
                } else {
                    slot[index] = module.downgrade();
                    module.clone()
                }
            })
            .collect()
    }
}

#[derive(Debug)]
struct PathAlgebraModuleHomomorphismData<F: Field> {
    source: PathAlgebraModule<F>,
    range: PathAlgebraModule<F>,
    vertex_matrices: Vec<DenseMatrix<F>>,
    kernel_inclusion: OnceLock<PathAlgebraModuleHomomorphism<F>>,
    image_projection_inclusion: OnceLock<(
        PathAlgebraModuleHomomorphism<F>,
        PathAlgebraModuleHomomorphism<F>,
    )>,
    cokernel_projection: OnceLock<PathAlgebraModuleHomomorphism<F>>,
    dual_homomorphism: OnceLock<PathAlgebraModuleHomomorphism<F>>,
    is_injective: OnceLock<bool>,
    is_surjective: OnceLock<bool>,
    is_split_monomorphism: OnceLock<bool>,
    is_split_epimorphism: OnceLock<bool>,
    is_right_minimal: OnceLock<bool>,
    is_left_minimal: OnceLock<bool>,
    left_inverse_weak: RwLock<Option<WeakPathAlgebraModuleHomomorphism<F>>>,
    left_inverse_fail: OnceLock<()>,
    right_inverse_weak: RwLock<Option<WeakPathAlgebraModuleHomomorphism<F>>>,
    right_inverse_fail: OnceLock<()>,
    right_minimal_version: OnceLock<MinimalVersionAttribute<F>>,
    left_minimal_version: OnceLock<MinimalVersionAttribute<F>>,
    kernel_of_what: RwLock<Option<WeakPathAlgebraModuleHomomorphism<F>>>,
    image_of_what: RwLock<Option<WeakPathAlgebraModuleHomomorphism<F>>>,
    cokernel_of_what: RwLock<Option<WeakPathAlgebraModuleHomomorphism<F>>>,
}

/// A QPA homomorphism between right path-algebra matrix modules.
///
/// The matrix at vertex `i` acts on row vectors from the source vertex space
/// to the range vertex space. Physical one-row/one-column zero-space
/// stand-ins are retained exactly as in `RightModuleHomOverAlgebra`. Clones
/// retain mapping identity so QPA's cached kernel, image, and cokernel
/// attributes can return the identical map object.
#[derive(Clone, Debug)]
pub struct PathAlgebraModuleHomomorphism<F: Field>(Arc<PathAlgebraModuleHomomorphismData<F>>);

/// An ordered pair of QPA module homomorphisms returned by categorical constructions.
pub type PathAlgebraModuleHomomorphismPair<F> = (
    PathAlgebraModuleHomomorphism<F>,
    PathAlgebraModuleHomomorphism<F>,
);

/// An ordered triple of QPA module homomorphisms returned by submodule constructions.
pub type PathAlgebraModuleHomomorphismTriple<F> = (
    PathAlgebraModuleHomomorphism<F>,
    PathAlgebraModuleHomomorphism<F>,
    PathAlgebraModuleHomomorphism<F>,
);

/// The heterogeneous result of QPA's one-step minimal-version operations.
#[derive(Clone, Debug)]
pub enum MoreMinimalVersionResult<F: Field> {
    /// QPA returned the identical entered map because no summand was removed.
    Unchanged(PathAlgebraModuleHomomorphism<F>),
    /// QPA returned the retained map followed by the removed zero component.
    Reduced {
        /// The restriction/corestriction that remains after one step.
        minimal_part: PathAlgebraModuleHomomorphism<F>,
        /// The zero component supported on the summand removed in this step.
        removed_part: PathAlgebraModuleHomomorphism<F>,
    },
}

/// QPA's minimal map together with the summands removed in successive steps.
pub type PathAlgebraModuleMinimalVersion<F> =
    (PathAlgebraModuleHomomorphism<F>, Vec<PathAlgebraModule<F>>);

#[derive(Debug)]
struct MinimalVersionAttribute<F: Field> {
    // `None` denotes the identical owning map and avoids an `Arc` self-cycle.
    other_homomorphism: Option<PathAlgebraModuleHomomorphism<F>>,
    removed_modules: Vec<PathAlgebraModule<F>>,
}

#[derive(Clone, Debug)]
struct WeakPathAlgebraModuleHomomorphism<F: Field>(Weak<PathAlgebraModuleHomomorphismData<F>>);

impl<F: Field> WeakPathAlgebraModuleHomomorphism<F> {
    fn upgrade(&self) -> Option<PathAlgebraModuleHomomorphism<F>> {
        self.0.upgrade().map(PathAlgebraModuleHomomorphism)
    }
}

fn cache_weak_homomorphism<F: Field>(
    slot: &RwLock<Option<WeakPathAlgebraModuleHomomorphism<F>>>,
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> PathAlgebraModuleHomomorphism<F> {
    let mut slot = slot
        .write()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    if let Some(cached) = slot
        .as_ref()
        .and_then(WeakPathAlgebraModuleHomomorphism::upgrade)
    {
        return cached;
    }
    *slot = Some(homomorphism.downgrade());
    homomorphism.clone()
}

/// QPA's `HomOverAlgebraWithBasisFunction` result.
///
/// The stored maps retain QPA's exact nullspace/full-matrix basis order. The
/// coordinate operation returns coefficients in that original map basis,
/// including QPA's one-entry zero vector for an empty Hom basis.
#[derive(Clone, Debug)]
pub struct PathAlgebraModuleHomomorphismBasis<F: Field> {
    basis: Vec<PathAlgebraModuleHomomorphism<F>>,
    flattened_basis: Option<DenseMatrix<F>>,
    field: F,
}

impl<F: Field> PathAlgebraModuleHomomorphismBasis<F> {
    /// Return the ordered basis produced by `HomOverAlgebra`.
    #[must_use]
    pub fn basis(&self) -> &[PathAlgebraModuleHomomorphism<F>] {
        &self.basis
    }

    /// Express a map in QPA's returned Hom basis.
    ///
    /// `Ok(None)` means the flattened map is not in the basis span. For an
    /// empty Hom basis, QPA's installed closure returns `[0]`.
    ///
    /// # Errors
    ///
    /// Propagates an incompatible flattened matrix length.
    pub fn coefficients(
        &self,
        homomorphism: &PathAlgebraModuleHomomorphism<F>,
    ) -> Result<Option<Vec<F::Element>>, MatrixError> {
        let Some(basis) = &self.flattened_basis else {
            return Ok(Some(vec![self.field.zero()]));
        };
        basis.solution_row(&flatten_module_homomorphism(homomorphism))
    }
}

#[derive(Debug)]
struct ModuleEndomorphismAlgebraData<F: Field> {
    field: F,
    basis_matrices: Vec<DenseMatrix<F>>,
    flattened_basis: DenseMatrix<F>,
    one: DenseMatrix<F>,
}

/// QPA's matrix-algebra realization of `EndOverAlgebra(M)`.
///
/// Elements use block-diagonal matrices on the vertex-major underlying vector
/// space of `M`. The basis is GAP's semi-echelon basis of the subalgebra with
/// one generated by the block matrices of `HomOverAlgebra(M, M)`.
#[derive(Clone, Debug)]
pub struct ModuleEndomorphismAlgebra<F: Field>(Arc<ModuleEndomorphismAlgebraData<F>>);

impl<F: Field> ModuleEndomorphismAlgebra<F> {
    /// Return whether two handles denote the identical cached QPA algebra.
    #[must_use]
    pub fn is_identical(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }

    /// Return the vector-space dimension.
    #[must_use]
    pub fn dimension(&self) -> usize {
        self.0.basis_matrices.len()
    }

    /// Return the exact coefficient field.
    #[must_use]
    pub fn field(&self) -> &F {
        &self.0.field
    }

    /// Return GAP's ordered semi-echelon matrix basis.
    #[must_use]
    pub fn basis_matrices(&self) -> &[DenseMatrix<F>] {
        &self.0.basis_matrices
    }

    /// Return the block identity matrix.
    #[must_use]
    pub fn one(&self) -> &DenseMatrix<F> {
        &self.0.one
    }

    /// Return the zero endomorphism matrix.
    #[must_use]
    pub fn zero(&self) -> DenseMatrix<F> {
        DenseMatrix::zero(self.0.one.nrows(), self.0.one.ncols(), self.0.field.clone())
    }

    /// Return coordinates in GAP's semi-echelon endomorphism-algebra basis.
    ///
    /// `Ok(None)` means the matrix is outside the represented subalgebra.
    ///
    /// # Errors
    ///
    /// Propagates a wrong flattened matrix length.
    pub fn coefficients(
        &self,
        matrix: &DenseMatrix<F>,
    ) -> Result<Option<Vec<F::Element>>, MatrixError> {
        if matrix.field() != &self.0.field {
            return Err(MatrixError::DifferentFields);
        }
        if (matrix.nrows(), matrix.ncols()) != (self.0.one.nrows(), self.0.one.ncols()) {
            return Err(MatrixError::DifferentShapes {
                left: (matrix.nrows(), matrix.ncols()),
                right: (self.0.one.nrows(), self.0.one.ncols()),
            });
        }
        self.0.flattened_basis.solution_row(&flatten_matrix(matrix))
    }

    /// Reconstruct a represented matrix from its ordered basis coefficients.
    ///
    /// # Errors
    ///
    /// Requires one coefficient per basis matrix and propagates exact matrix
    /// multiplication failures.
    pub fn matrix_from_coefficients(
        &self,
        coefficients: &[F::Element],
    ) -> Result<DenseMatrix<F>, MatrixError> {
        if coefficients.len() != self.dimension() {
            return Err(MatrixError::TargetLength {
                actual: coefficients.len(),
                expected: self.dimension(),
            });
        }
        let coefficient_row = DenseMatrix::from_flat(
            1,
            coefficients.len(),
            coefficients.to_vec(),
            self.0.field.clone(),
        )?;
        let flattened = coefficient_row.multiply(&self.0.flattened_basis)?;
        DenseMatrix::from_flat(
            self.0.one.nrows(),
            self.0.one.ncols(),
            flattened.row(0).to_vec(),
            self.0.field.clone(),
        )
    }

    /// Multiply two represented endomorphism matrices in QPA order.
    ///
    /// # Errors
    ///
    /// Propagates matrix field or shape mismatches.
    pub fn multiply(
        &self,
        left: &DenseMatrix<F>,
        right: &DenseMatrix<F>,
    ) -> Result<DenseMatrix<F>, MatrixError> {
        left.multiply(right)
    }
}

#[derive(Debug)]
struct ModuleEndomorphismQuotientAlgebraData<F: Field> {
    parent: ModuleEndomorphismAlgebra<F>,
    ideal_basis_matrices: Vec<DenseMatrix<F>>,
    projection: DenseMatrix<F>,
    lifts: DenseMatrix<F>,
}

/// QPA's quotient algebra `End(M) / P(M,M)` represented in natural-quotient
/// coordinates.
///
/// Its ordered basis is the basis returned for GAP's quotient algebra. A
/// zero-dimensional quotient deliberately has no multiplicative-one element,
/// matching GAP's result when the ideal is the whole endomorphism algebra.
#[derive(Clone, Debug)]
pub struct ModuleEndomorphismQuotientAlgebra<F: Field>(
    Arc<ModuleEndomorphismQuotientAlgebraData<F>>,
);

/// An element of [`ModuleEndomorphismQuotientAlgebra`] in its GAP basis.
#[derive(Clone, Debug)]
pub struct ModuleEndomorphismQuotientElement<F: Field> {
    algebra: ModuleEndomorphismQuotientAlgebra<F>,
    coordinates: Vec<F::Element>,
}

/// QPA's `NaturalHomomorphismByIdeal(End(M), P(M,M))` return value.
#[derive(Clone, Debug)]
pub struct ModuleEndomorphismQuotientHomomorphism<F: Field> {
    source: ModuleEndomorphismAlgebra<F>,
    range: ModuleEndomorphismQuotientAlgebra<F>,
}

impl<F: Field> ModuleEndomorphismQuotientAlgebra<F> {
    /// Return whether two handles denote the identical quotient algebra.
    #[must_use]
    pub fn is_identical(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }

    /// Return the vector-space dimension of the quotient.
    #[must_use]
    pub fn dimension(&self) -> usize {
        self.0.lifts.nrows()
    }

    /// Return GAP's ordered basis of the kernel ideal as parent matrices.
    #[must_use]
    pub fn ideal_basis_matrices(&self) -> &[DenseMatrix<F>] {
        &self.0.ideal_basis_matrices
    }

    /// Return the quotient basis as coordinate elements.
    #[must_use]
    pub fn basis_elements(&self) -> Vec<ModuleEndomorphismQuotientElement<F>> {
        (0..self.dimension())
            .map(|coordinate| {
                let mut coordinates = vec![self.0.parent.field().zero(); self.dimension()];
                coordinates[coordinate] = self.0.parent.field().one();
                ModuleEndomorphismQuotientElement {
                    algebra: self.clone(),
                    coordinates,
                }
            })
            .collect()
    }

    /// Return the quotient zero element.
    #[must_use]
    pub fn zero(&self) -> ModuleEndomorphismQuotientElement<F> {
        ModuleEndomorphismQuotientElement {
            algebra: self.clone(),
            coordinates: vec![self.0.parent.field().zero(); self.dimension()],
        }
    }

    /// Return the quotient one, or `None` for GAP's zero-dimensional quotient
    /// algebra where `One` has no applicable method.
    ///
    /// # Errors
    ///
    /// Propagates the exact natural projection of the parent identity.
    pub fn one(
        &self,
    ) -> Result<Option<ModuleEndomorphismQuotientElement<F>>, PathAlgebraModuleError> {
        if self.dimension() == 0 {
            return Ok(None);
        }
        Ok(Some(self.project_parent_matrix(self.0.parent.one())?))
    }

    /// Multiply quotient elements using lifted representatives and the parent
    /// endomorphism algebra.
    ///
    /// # Errors
    ///
    /// Both elements must belong to this identical quotient family; exact
    /// lift, multiplication, coordinate, and projection failures propagate.
    pub fn multiply(
        &self,
        left: &ModuleEndomorphismQuotientElement<F>,
        right: &ModuleEndomorphismQuotientElement<F>,
    ) -> Result<ModuleEndomorphismQuotientElement<F>, PathAlgebraModuleError> {
        if !self.is_identical(left.algebra()) || !self.is_identical(right.algebra()) {
            return Err(PathAlgebraModuleError::ForeignEndomorphismQuotientElement);
        }
        let left = self.lift_element(left)?;
        let right = self.lift_element(right)?;
        let product = self.0.parent.multiply(&left, &right)?;
        self.project_parent_matrix(&product)
    }

    fn project_parent_matrix(
        &self,
        matrix: &DenseMatrix<F>,
    ) -> Result<ModuleEndomorphismQuotientElement<F>, PathAlgebraModuleError> {
        let coordinates = self
            .0
            .parent
            .coefficients(matrix)?
            .ok_or(PathAlgebraModuleError::MatrixOutsideEndomorphismAlgebra)?;
        let row = DenseMatrix::from_flat(
            1,
            coordinates.len(),
            coordinates,
            self.0.parent.field().clone(),
        )?;
        let projected = row.multiply(&self.0.projection)?;
        Ok(ModuleEndomorphismQuotientElement {
            algebra: self.clone(),
            coordinates: projected.row(0).to_vec(),
        })
    }

    fn lift_element(
        &self,
        element: &ModuleEndomorphismQuotientElement<F>,
    ) -> Result<DenseMatrix<F>, PathAlgebraModuleError> {
        if element.coordinates.len() != self.dimension() {
            return Err(
                PathAlgebraModuleError::EndomorphismQuotientElementCoordinateCount {
                    actual: element.coordinates.len(),
                    expected: self.dimension(),
                },
            );
        }
        let row = DenseMatrix::from_flat(
            1,
            element.coordinates.len(),
            element.coordinates.clone(),
            self.0.parent.field().clone(),
        )?;
        let parent_coordinates = row.multiply(&self.0.lifts)?;
        Ok(self
            .0
            .parent
            .matrix_from_coefficients(parent_coordinates.row(0))?)
    }
}

impl<F: Field> ModuleEndomorphismQuotientElement<F> {
    /// Return the identical quotient algebra family.
    #[must_use]
    pub fn algebra(&self) -> &ModuleEndomorphismQuotientAlgebra<F> {
        &self.algebra
    }

    /// Return coefficients in GAP's quotient basis.
    #[must_use]
    pub fn coordinates(&self) -> &[F::Element] {
        &self.coordinates
    }

    /// Return whether all quotient coordinates vanish.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.coordinates
            .iter()
            .all(|entry| self.algebra.0.parent.field().is_zero(entry))
    }
}

impl<F: Field> ModuleEndomorphismQuotientHomomorphism<F> {
    /// Return QPA's source `EndOverAlgebra(M)`.
    #[must_use]
    pub fn source(&self) -> &ModuleEndomorphismAlgebra<F> {
        &self.source
    }

    /// Return QPA's stable endomorphism quotient algebra.
    #[must_use]
    pub fn range(&self) -> &ModuleEndomorphismQuotientAlgebra<F> {
        &self.range
    }

    /// Return the natural image of a represented parent matrix.
    ///
    /// # Errors
    ///
    /// The matrix must belong to the represented source algebra; exact
    /// coordinate and quotient projection failures propagate.
    pub fn image(
        &self,
        matrix: &DenseMatrix<F>,
    ) -> Result<ModuleEndomorphismQuotientElement<F>, PathAlgebraModuleError> {
        self.range.project_parent_matrix(matrix)
    }

    /// Return GAP's ordered basis for the kernel ideal.
    #[must_use]
    pub fn kernel_basis_matrices(&self) -> &[DenseMatrix<F>] {
        self.range.ideal_basis_matrices()
    }
}

impl<F: Field> PartialEq for PathAlgebraModuleHomomorphism<F> {
    fn eq(&self, other: &Self) -> bool {
        self.source() == other.source()
            && self.range() == other.range()
            && self.vertex_matrices() == other.vertex_matrices()
    }
}

impl<F: Field> Eq for PathAlgebraModuleHomomorphism<F> {}

impl<F: Field> PathAlgebraModuleHomomorphism<F> {
    /// Return QPA's `Source` module.
    #[must_use]
    pub fn source(&self) -> &PathAlgebraModule<F> {
        &self.0.source
    }

    /// Return QPA's `Range` module.
    #[must_use]
    pub fn range(&self) -> &PathAlgebraModule<F> {
        &self.0.range
    }

    /// Translate `MatricesOfPathAlgebraMatModuleHomomorphism`.
    #[must_use]
    pub fn vertex_matrices(&self) -> &[DenseMatrix<F>] {
        &self.0.vertex_matrices
    }

    /// Flatten QPA's internal `f!.maps` list in vertex, row, column order.
    ///
    /// Physical zero-space stand-ins are retained, as required by source
    /// methods that take vector-space spans of raw module-map data.
    #[must_use]
    pub fn flat_vertex_matrix_entries(&self) -> Vec<F::Element> {
        self.vertex_matrices()
            .iter()
            .flat_map(|matrix| matrix.as_flat().iter().cloned())
            .collect()
    }

    /// Return whether two handles denote the identical QPA mapping object.
    #[must_use]
    pub fn is_identical(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }

    /// Return QPA's cached `DualOfModuleHomomorphism` attribute, if set.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_dual_homomorphism(&self) -> Option<Self> {
        self.0.dual_homomorphism.get().cloned()
    }

    /// Install QPA's one-way `DualOfModuleHomomorphism` attribute cache.
    ///
    /// The source method does not install the reciprocal attribute on the
    /// generated dual map.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_dual_homomorphism(&self, dual: &Self) -> Self {
        self.0
            .dual_homomorphism
            .get_or_init(|| dual.clone())
            .clone()
    }

    /// Apply the mapping to a module element.
    ///
    /// # Errors
    ///
    /// Rejects an element outside the identical source family and propagates
    /// exact matrix failures.
    pub fn image(
        &self,
        element: &PathAlgebraModuleElement<F>,
    ) -> Result<PathAlgebraModuleElement<F>, PathAlgebraModuleError> {
        if !self.source().is_identical(element.module()) {
            return Err(PathAlgebraModuleError::ForeignMapSourceElement);
        }
        let mut image = self.range().zero_element();
        for (vertex, matrix) in self.vertex_matrices().iter().enumerate() {
            let row = DenseMatrix::from_flat(
                1,
                element.coordinates[vertex].len(),
                element.coordinates[vertex].clone(),
                matrix.field().clone(),
            )?;
            image.coordinates[vertex] = row.multiply(matrix)?.row(0).to_vec();
        }
        Ok(image)
    }

    /// Translate QPA's `PreImagesRepresentativeNC` for a module map.
    ///
    /// `Ok(None)` is GAP's `fail` when no preimage exists.
    ///
    /// # Errors
    ///
    /// Rejects an element outside the identical range family and propagates
    /// exact row-solution failures.
    pub fn preimage_representative(
        &self,
        element: &PathAlgebraModuleElement<F>,
    ) -> Result<Option<PathAlgebraModuleElement<F>>, PathAlgebraModuleError> {
        if !self.range().is_identical(element.module()) {
            return Err(PathAlgebraModuleError::ForeignMapRangeElement);
        }
        let mut preimage = self.source().zero_element();
        for (vertex, matrix) in self.vertex_matrices().iter().enumerate() {
            let Some(row) = matrix.solution_row(&element.coordinates[vertex])? else {
                return Ok(None);
            };
            preimage.coordinates[vertex] = row;
        }
        Ok(Some(preimage))
    }

    /// Translate QPA's `ImagesSet` list branch, omitting zero images while
    /// retaining the entered order and duplicates.
    ///
    /// # Errors
    ///
    /// Propagates [`Self::image`], including its identical-source check.
    pub fn images_set(
        &self,
        elements: &[PathAlgebraModuleElement<F>],
    ) -> Result<Vec<PathAlgebraModuleElement<F>>, PathAlgebraModuleError> {
        let mut images = Vec::new();
        for element in elements {
            let image = self.image(element)?;
            if !image.is_zero() {
                images.push(image);
            }
        }
        Ok(images)
    }

    /// Return QPA's zero mapping with the same source and range.
    ///
    /// # Errors
    ///
    /// Propagates the checked zero-map constructor.
    pub fn zero(&self) -> Result<Self, PathAlgebraModuleError> {
        zero_module_mapping(self.source().clone(), self.range().clone())
    }

    /// Add two maps in the same structural homomorphism set.
    ///
    /// # Errors
    ///
    /// Rejects different sources or ranges and propagates exact matrix or
    /// checked module-map failures.
    pub fn add(&self, right: &Self) -> Result<Self, PathAlgebraModuleError> {
        if self.source() != right.source() || self.range() != right.range() {
            return Err(PathAlgebraModuleError::DifferentHomomorphismSets);
        }
        let matrices = self
            .vertex_matrices()
            .iter()
            .zip(right.vertex_matrices())
            .map(|(left, right)| left.add(right))
            .collect::<Result<Vec<_>, _>>()?;
        right_module_hom_over_algebra(self.source().clone(), self.range().clone(), matrices)
    }

    /// Compose maps in QPA's row-action order: first `self`, then `right`.
    ///
    /// # Errors
    ///
    /// Rejects a structural range/source mismatch and propagates exact matrix
    /// or checked module-map failures.
    pub fn compose(&self, right: &Self) -> Result<Self, PathAlgebraModuleError> {
        if self.range() != right.source() {
            return Err(PathAlgebraModuleError::NoncomposableHomomorphisms);
        }
        let matrices = self
            .vertex_matrices()
            .iter()
            .zip(right.vertex_matrices())
            .map(|(left, right)| left.multiply(right))
            .collect::<Result<Vec<_>, _>>()?;
        right_module_hom_over_algebra(self.source().clone(), right.range().clone(), matrices)
    }

    /// Multiply every vertex matrix by a coefficient-field scalar.
    ///
    /// Numerica's element type makes the coefficient domain part of the Rust
    /// type for `QQ` and fixed fields; runtime prime-field descriptors are
    /// retained by each matrix and the checked constructor.
    ///
    /// # Errors
    ///
    /// Propagates the checked module-map constructor.
    pub fn scale(&self, scalar: &F::Element) -> Result<Self, PathAlgebraModuleError> {
        right_module_hom_over_algebra(
            self.source().clone(),
            self.range().clone(),
            self.vertex_matrices()
                .iter()
                .map(|matrix| matrix.scale_left(scalar))
                .collect(),
        )
    }

    /// Return QPA's additive inverse map.
    ///
    /// # Errors
    ///
    /// Propagates the checked module-map constructor.
    pub fn negate(&self) -> Result<Self, PathAlgebraModuleError> {
        let field = self.source().right_acting_algebra().field();
        self.scale(&field.neg(&field.one()))
    }

    /// Return whether every stored physical vertex matrix is zero.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.vertex_matrices().iter().all(DenseMatrix::is_zero)
    }

    /// Translate QPA's componentwise `IsInjective` method.
    #[must_use]
    pub fn is_injective(&self) -> bool {
        *self.0.is_injective.get_or_init(|| {
            self.source()
                .dimension_vector()
                .iter()
                .enumerate()
                .all(|(vertex, &dimension)| {
                    dimension == 0 || logical_vertex_map(self, vertex).rank() == dimension
                })
        })
    }

    /// Translate QPA's `IsSurjective` image-space dimension test.
    #[must_use]
    pub fn is_surjective(&self) -> bool {
        *self.0.is_surjective.get_or_init(|| {
            self.range()
                .dimension_vector()
                .iter()
                .enumerate()
                .all(|(vertex, &dimension)| logical_vertex_map(self, vertex).rank() == dimension)
        })
    }

    /// Translate QPA's short-circuit injective-and-surjective isomorphism
    /// predicate.
    #[must_use]
    pub fn is_isomorphism(&self) -> bool {
        self.is_injective() && self.is_surjective()
    }

    /// Return QPA's cached `IsSplitMonomorphism` value, if installed.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_is_split_monomorphism(&self) -> Option<bool> {
        self.0.is_split_monomorphism.get().copied()
    }

    /// Install QPA's `IsSplitMonomorphism` property value.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_is_split_monomorphism(&self, value: bool) -> bool {
        *self.0.is_split_monomorphism.get_or_init(|| value)
    }

    /// Return QPA's cached `IsSplitEpimorphism` value, if installed.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_is_split_epimorphism(&self) -> Option<bool> {
        self.0.is_split_epimorphism.get().copied()
    }

    /// Install QPA's `IsSplitEpimorphism` property value.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_is_split_epimorphism(&self, value: bool) -> bool {
        *self.0.is_split_epimorphism.get_or_init(|| value)
    }

    /// Return QPA's cached `IsRightMinimal` value, if installed.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_is_right_minimal(&self) -> Option<bool> {
        self.0.is_right_minimal.get().copied()
    }

    /// Install QPA's `IsRightMinimal` property value.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_is_right_minimal(&self, value: bool) -> bool {
        *self.0.is_right_minimal.get_or_init(|| value)
    }

    /// Return QPA's cached `IsLeftMinimal` value, if installed.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_is_left_minimal(&self) -> Option<bool> {
        self.0.is_left_minimal.get().copied()
    }

    /// Install QPA's `IsLeftMinimal` property value.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_is_left_minimal(&self, value: bool) -> bool {
        *self.0.is_left_minimal.get_or_init(|| value)
    }

    /// Return an already-live cached `LeftInverseOfHomomorphism` attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_left_inverse(&self) -> Option<Self> {
        self.0
            .left_inverse_weak
            .read()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .as_ref()
            .and_then(WeakPathAlgebraModuleHomomorphism::upgrade)
    }

    /// Return whether QPA's left-inverse attribute is bound to `false`.
    #[doc(hidden)]
    #[must_use]
    pub fn left_inverse_attribute_is_fail(&self) -> bool {
        self.0.left_inverse_fail.get().is_some()
    }

    /// Cache QPA's failed left-inverse attribute result.
    #[doc(hidden)]
    pub fn mark_left_inverse_attribute_fail(&self) {
        let _ = self.0.left_inverse_fail.set(());
    }

    /// Install a live `LeftInverseOfHomomorphism` attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_left_inverse(&self, inverse: &Self) -> Self {
        cache_weak_homomorphism(&self.0.left_inverse_weak, inverse)
    }

    /// Return an already-live cached `RightInverseOfHomomorphism` attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_right_inverse(&self) -> Option<Self> {
        self.0
            .right_inverse_weak
            .read()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .as_ref()
            .and_then(WeakPathAlgebraModuleHomomorphism::upgrade)
    }

    /// Return whether QPA's right-inverse attribute is bound to `false`.
    #[doc(hidden)]
    #[must_use]
    pub fn right_inverse_attribute_is_fail(&self) -> bool {
        self.0.right_inverse_fail.get().is_some()
    }

    /// Cache QPA's failed right-inverse attribute result.
    #[doc(hidden)]
    pub fn mark_right_inverse_attribute_fail(&self) {
        let _ = self.0.right_inverse_fail.set(());
    }

    /// Install a live `RightInverseOfHomomorphism` attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_right_inverse(&self, inverse: &Self) -> Self {
        cache_weak_homomorphism(&self.0.right_inverse_weak, inverse)
    }

    /// Return the map whose cached kernel inclusion this map is, if it is
    /// still live.
    #[must_use]
    pub fn kernel_of_what(&self) -> Option<Self> {
        self.0
            .kernel_of_what
            .read()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .as_ref()
            .and_then(WeakPathAlgebraModuleHomomorphism::upgrade)
    }

    /// Return the map whose cached image factor this map belongs to, if it is
    /// still live.
    #[must_use]
    pub fn image_of_what(&self) -> Option<Self> {
        self.0
            .image_of_what
            .read()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .as_ref()
            .and_then(WeakPathAlgebraModuleHomomorphism::upgrade)
    }

    /// Return the map whose cached cokernel projection this map is, if it is
    /// still live.
    #[must_use]
    pub fn cokernel_of_what(&self) -> Option<Self> {
        self.0
            .cokernel_of_what
            .read()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .as_ref()
            .and_then(WeakPathAlgebraModuleHomomorphism::upgrade)
    }

    fn set_is_injective(&self, value: bool) {
        let _ = self.0.is_injective.set(value);
    }

    fn set_is_surjective(&self, value: bool) {
        let _ = self.0.is_surjective.set(value);
    }

    fn downgrade(&self) -> WeakPathAlgebraModuleHomomorphism<F> {
        WeakPathAlgebraModuleHomomorphism(Arc::downgrade(&self.0))
    }

    fn cached_kernel_inclusion(&self) -> Option<Self> {
        self.0.kernel_inclusion.get().cloned()
    }

    fn cache_kernel_inclusion(&self, inclusion: &Self) -> Self {
        let _ = self.0.kernel_inclusion.set(inclusion.clone());
        self.0
            .kernel_inclusion
            .get()
            .expect("kernel inclusion was just initialized")
            .clone()
    }

    fn cached_image_projection_inclusion(&self) -> Option<(Self, Self)> {
        self.0.image_projection_inclusion.get().cloned()
    }

    fn cache_image_projection_inclusion(
        &self,
        projection: &Self,
        inclusion: &Self,
    ) -> (Self, Self) {
        let _ = self
            .0
            .image_projection_inclusion
            .set((projection.clone(), inclusion.clone()));
        self.0
            .image_projection_inclusion
            .get()
            .expect("image factorization was just initialized")
            .clone()
    }

    fn cached_cokernel_projection(&self) -> Option<Self> {
        self.0.cokernel_projection.get().cloned()
    }

    fn cache_cokernel_projection(&self, projection: &Self) -> Self {
        let _ = self.0.cokernel_projection.set(projection.clone());
        self.0
            .cokernel_projection
            .get()
            .expect("cokernel projection was just initialized")
            .clone()
    }

    fn cached_right_minimal_version(&self) -> Option<PathAlgebraModuleMinimalVersion<F>> {
        self.0.right_minimal_version.get().map(|attribute| {
            (
                attribute
                    .other_homomorphism
                    .clone()
                    .unwrap_or_else(|| self.clone()),
                attribute.removed_modules.clone(),
            )
        })
    }

    fn cache_right_minimal_version(
        &self,
        homomorphism: &Self,
        removed_modules: &[PathAlgebraModule<F>],
    ) -> PathAlgebraModuleMinimalVersion<F> {
        let attribute = self
            .0
            .right_minimal_version
            .get_or_init(|| MinimalVersionAttribute {
                other_homomorphism: (!self.is_identical(homomorphism))
                    .then(|| homomorphism.clone()),
                removed_modules: removed_modules.to_vec(),
            });
        (
            attribute
                .other_homomorphism
                .clone()
                .unwrap_or_else(|| self.clone()),
            attribute.removed_modules.clone(),
        )
    }

    fn cached_left_minimal_version(&self) -> Option<PathAlgebraModuleMinimalVersion<F>> {
        self.0.left_minimal_version.get().map(|attribute| {
            (
                attribute
                    .other_homomorphism
                    .clone()
                    .unwrap_or_else(|| self.clone()),
                attribute.removed_modules.clone(),
            )
        })
    }

    fn cache_left_minimal_version(
        &self,
        homomorphism: &Self,
        removed_modules: &[PathAlgebraModule<F>],
    ) -> PathAlgebraModuleMinimalVersion<F> {
        let attribute = self
            .0
            .left_minimal_version
            .get_or_init(|| MinimalVersionAttribute {
                other_homomorphism: (!self.is_identical(homomorphism))
                    .then(|| homomorphism.clone()),
                removed_modules: removed_modules.to_vec(),
            });
        (
            attribute
                .other_homomorphism
                .clone()
                .unwrap_or_else(|| self.clone()),
            attribute.removed_modules.clone(),
        )
    }

    fn set_kernel_of_what(&self, homomorphism: &Self) {
        *self
            .0
            .kernel_of_what
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner) = Some(homomorphism.downgrade());
    }

    fn set_image_of_what(&self, homomorphism: &Self) {
        *self
            .0
            .image_of_what
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner) = Some(homomorphism.downgrade());
    }

    fn set_cokernel_of_what(&self, homomorphism: &Self) {
        *self
            .0
            .cokernel_of_what
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner) = Some(homomorphism.downgrade());
    }
}

/// Translate QPA's checked `RightModuleHomOverAlgebra` constructor.
///
/// # Errors
///
/// Requires structurally equal acting algebras, one correctly shaped matrix
/// over the common field per vertex, and every arrow square to commute.
pub fn right_module_hom_over_algebra<F: Field>(
    source: PathAlgebraModule<F>,
    range: PathAlgebraModule<F>,
    vertex_matrices: Vec<DenseMatrix<F>>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    if source.right_acting_algebra() != range.right_acting_algebra() {
        return Err(PathAlgebraModuleError::DifferentHomomorphismAlgebras);
    }
    if vertex_matrices.len() != source.dimension_vector().len() {
        return Err(PathAlgebraModuleError::HomomorphismMatrixCount {
            actual: vertex_matrices.len(),
            expected: source.dimension_vector().len(),
        });
    }
    for (vertex, matrix) in vertex_matrices.iter().enumerate() {
        if matrix.field() != source.right_acting_algebra().field() {
            return Err(PathAlgebraModuleError::HomomorphismMatrixField { vertex: vertex + 1 });
        }
        let expected = (
            physical_dimension(source.dimension_vector()[vertex]),
            physical_dimension(range.dimension_vector()[vertex]),
        );
        let actual = (matrix.nrows(), matrix.ncols());
        if actual != expected
            || (source.dimension_vector()[vertex] == 0 || range.dimension_vector()[vertex] == 0)
                && !matrix.is_zero()
        {
            return Err(PathAlgebraModuleError::HomomorphismMatrixShape {
                vertex: vertex + 1,
                actual,
                expected,
            });
        }
    }

    for arrow in source.right_acting_algebra().quiver().arrows() {
        let source_vertex = arrow.source().index().zero_based();
        let target_vertex = arrow.target().index().zero_based();
        let arrow_index = arrow.index().zero_based();
        let left =
            source.arrow_matrices()[arrow_index].multiply(&vertex_matrices[target_vertex])?;
        let right =
            vertex_matrices[source_vertex].multiply(&range.arrow_matrices()[arrow_index])?;
        if left != right {
            return Err(PathAlgebraModuleError::NoncommutingHomomorphismSquare {
                arrow: arrow.name().into(),
            });
        }
    }

    Ok(PathAlgebraModuleHomomorphism(Arc::new(
        PathAlgebraModuleHomomorphismData {
            source,
            range,
            vertex_matrices,
            kernel_inclusion: OnceLock::new(),
            image_projection_inclusion: OnceLock::new(),
            cokernel_projection: OnceLock::new(),
            dual_homomorphism: OnceLock::new(),
            is_injective: OnceLock::new(),
            is_surjective: OnceLock::new(),
            is_split_monomorphism: OnceLock::new(),
            is_split_epimorphism: OnceLock::new(),
            is_right_minimal: OnceLock::new(),
            is_left_minimal: OnceLock::new(),
            left_inverse_weak: RwLock::new(None),
            left_inverse_fail: OnceLock::new(),
            right_inverse_weak: RwLock::new(None),
            right_inverse_fail: OnceLock::new(),
            right_minimal_version: OnceLock::new(),
            left_minimal_version: OnceLock::new(),
            kernel_of_what: RwLock::new(None),
            image_of_what: RwLock::new(None),
            cokernel_of_what: RwLock::new(None),
        },
    )))
}

/// Translate QPA's cached `ZeroModule(A)` attribute.
///
/// The module has one zero dimension for every vertex and no entered arrow
/// matrices. A weak algebra-side cache preserves GAP object identity while a
/// returned zero module remains live without leaking an algebra/module cycle.
///
/// # Errors
///
/// Reports the theoretically impossible rejection of the all-zero
/// representation by the checked module constructor.
pub fn zero_module<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<PathAlgebraModule<F>, PathAlgebraModuleError> {
    if let Some(cached) = algebra.cached_zero_module() {
        return Ok(cached);
    }
    let generated = right_module_over_path_algebra(
        algebra.clone(),
        vec![0; algebra.quiver().number_of_vertices()],
        &[],
    )?
    .ok_or(PathAlgebraModuleError::GeneratedZeroModuleViolatesRelations)?;
    Ok(algebra.cache_zero_module(&generated))
}

/// Translate QPA's zero mapping between two modules.
///
/// # Errors
///
/// Rejects modules over different algebras.
pub fn zero_module_mapping<F: Field>(
    source: PathAlgebraModule<F>,
    range: PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    if source.right_acting_algebra() != range.right_acting_algebra() {
        return Err(PathAlgebraModuleError::DifferentHomomorphismAlgebras);
    }
    let field = source.right_acting_algebra().field().clone();
    let matrices = source
        .dimension_vector()
        .iter()
        .zip(range.dimension_vector())
        .map(|(&source_dimension, &range_dimension)| {
            DenseMatrix::zero(
                physical_dimension(source_dimension),
                physical_dimension(range_dimension),
                field.clone(),
            )
        })
        .collect();
    right_module_hom_over_algebra(source, range, matrices)
}

/// Translate QPA's identity mapping on a module.
///
/// At a zero-dimensional vertex QPA stores a `1 x 1` zero matrix rather than
/// an identity matrix.
///
/// # Errors
///
/// Propagates the checked module-map constructor.
pub fn identity_module_mapping<F: Field>(
    module: PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    let field = module.right_acting_algebra().field().clone();
    let matrices = module
        .dimension_vector()
        .iter()
        .map(|&dimension| {
            if dimension == 0 {
                DenseMatrix::zero(1, 1, field.clone())
            } else {
                DenseMatrix::identity(dimension, field.clone())
            }
        })
        .collect();
    right_module_hom_over_algebra(module.clone(), module, matrices)
}

/// Translate QPA's `InverseOfIsomorphism(f)` operation.
///
/// QPA retains physical zero matrices at zero-dimensional vertices and takes
/// the ordinary exact inverse of every other vertex matrix.
///
/// # Errors
///
/// Rejects a map that is not an isomorphism and propagates exact matrix or
/// checked module-map construction failures.
pub fn inverse_of_isomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    if !homomorphism.is_isomorphism() {
        return Err(PathAlgebraModuleError::EnteredMapNotIsomorphism);
    }
    let inverse_matrices = homomorphism
        .vertex_matrices()
        .iter()
        .map(|matrix| {
            if matrix.is_zero() {
                Ok(matrix.clone())
            } else {
                matrix
                    .inverse()?
                    .ok_or(PathAlgebraModuleError::InternalIsomorphismMatrixSingular)
            }
        })
        .collect::<Result<Vec<_>, PathAlgebraModuleError>>()?;
    right_module_hom_over_algebra(
        homomorphism.range().clone(),
        homomorphism.source().clone(),
        inverse_matrices,
    )
}

/// Translate QPA's cached `RightInverseOfHomomorphism(f)` attribute.
///
/// `Ok(None)` is QPA's returned `false`: either `f` is not injective, the
/// reverse Hom basis is empty, or GAP's row solver finds no retraction.
/// Successful calls install the same reciprocal inverse attributes and split
/// properties as the source method.
///
/// # Errors
///
/// Propagates Hom-space construction, composition, exact row solving, and
/// checked linear-combination failures.
pub fn right_inverse_of_homomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<Option<PathAlgebraModuleHomomorphism<F>>, PathAlgebraModuleError> {
    if let Some(cached) = homomorphism.cached_right_inverse() {
        return Ok(Some(cached));
    }
    if homomorphism.right_inverse_attribute_is_fail() {
        return Ok(None);
    }
    if !homomorphism.is_injective() {
        homomorphism.mark_right_inverse_attribute_fail();
        return Ok(None);
    }
    let candidates = hom_over_algebra(homomorphism.range(), homomorphism.source())?
        .ok_or(PathAlgebraModuleError::InternalSplittingHomRejected)?;
    if candidates.is_empty() {
        homomorphism.mark_right_inverse_attribute_fail();
        return Ok(None);
    }
    let composites = candidates
        .iter()
        .map(|candidate| homomorphism.compose(candidate))
        .collect::<Result<Vec<_>, _>>()?;
    let identity = identity_module_mapping(homomorphism.source().clone())?;
    let Some(coefficients) = solve_map_linear_combination(&composites, &identity)? else {
        let _ = homomorphism.cache_is_split_monomorphism(false);
        homomorphism.mark_right_inverse_attribute_fail();
        return Ok(None);
    };
    let inverse = linear_combination_of_homomorphisms(&candidates, &coefficients)?;
    let _ = inverse.cache_is_split_epimorphism(true);
    let _ = homomorphism.cache_is_split_monomorphism(true);
    let _ = inverse.cache_left_inverse(homomorphism);
    Ok(Some(homomorphism.cache_right_inverse(&inverse)))
}

/// Translate QPA's cached `LeftInverseOfHomomorphism(f)` attribute.
///
/// `Ok(None)` is QPA's returned `false`: either `f` is not surjective, the
/// reverse Hom basis is empty, or GAP's row solver finds no section.
/// Successful calls install the source method's reciprocal attributes and
/// split properties.
///
/// # Errors
///
/// Propagates Hom-space construction, composition, exact row solving, and
/// checked linear-combination failures.
pub fn left_inverse_of_homomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<Option<PathAlgebraModuleHomomorphism<F>>, PathAlgebraModuleError> {
    if let Some(cached) = homomorphism.cached_left_inverse() {
        return Ok(Some(cached));
    }
    if homomorphism.left_inverse_attribute_is_fail() {
        return Ok(None);
    }
    if !homomorphism.is_surjective() {
        homomorphism.mark_left_inverse_attribute_fail();
        return Ok(None);
    }
    let candidates = hom_over_algebra(homomorphism.range(), homomorphism.source())?
        .ok_or(PathAlgebraModuleError::InternalSplittingHomRejected)?;
    if candidates.is_empty() {
        homomorphism.mark_left_inverse_attribute_fail();
        return Ok(None);
    }
    let composites = candidates
        .iter()
        .map(|candidate| candidate.compose(homomorphism))
        .collect::<Result<Vec<_>, _>>()?;
    let identity = identity_module_mapping(homomorphism.range().clone())?;
    let Some(coefficients) = solve_map_linear_combination(&composites, &identity)? else {
        let _ = homomorphism.cache_is_split_epimorphism(false);
        homomorphism.mark_left_inverse_attribute_fail();
        return Ok(None);
    };
    let inverse = linear_combination_of_homomorphisms(&candidates, &coefficients)?;
    let _ = inverse.cache_is_split_monomorphism(true);
    let _ = homomorphism.cache_is_split_epimorphism(true);
    let _ = inverse.cache_right_inverse(homomorphism);
    Ok(Some(homomorphism.cache_left_inverse(&inverse)))
}

/// Translate QPA's `IsSplitMonomorphism(f)` property method.
///
/// # Errors
///
/// Propagates [`right_inverse_of_homomorphism`].
pub fn is_split_monomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<bool, PathAlgebraModuleError> {
    if let Some(cached) = homomorphism.cached_is_split_monomorphism() {
        return Ok(cached);
    }
    let value = right_inverse_of_homomorphism(homomorphism)?.is_some();
    Ok(homomorphism.cache_is_split_monomorphism(value))
}

/// Translate QPA's `IsSplitEpimorphism(f)` property method.
///
/// # Errors
///
/// Propagates [`left_inverse_of_homomorphism`].
pub fn is_split_epimorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<bool, PathAlgebraModuleError> {
    if let Some(cached) = homomorphism.cached_is_split_epimorphism() {
        return Ok(cached);
    }
    let value = left_inverse_of_homomorphism(homomorphism)?.is_some();
    Ok(homomorphism.cache_is_split_epimorphism(value))
}

/// Translate QPA's one-step `MoreRightMinimalVersion(f)` operation.
///
/// Hom bases, the nested loop order, the Fitting exponent, and the returned
/// kernel/image restrictions follow `modulehom.gi` literally. `Unchanged`
/// contains the identical entered map, matching QPA's heterogeneous map-or-list
/// result.
///
/// # Errors
///
/// Propagates kernel, Hom-space, composition, power, and image construction
/// failures.
pub fn more_right_minimal_version<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<MoreMinimalVersionResult<F>, PathAlgebraModuleError> {
    let source = homomorphism.source();
    let kernel = kernel_inclusion(homomorphism)?;
    let kernel_module = kernel.source();
    let hom_source_kernel = hom_over_algebra(source, kernel_module)?
        .ok_or(PathAlgebraModuleError::InternalMinimalHomRejected)?;
    if hom_source_kernel.is_empty() {
        return Ok(MoreMinimalVersionResult::Unchanged(homomorphism.clone()));
    }
    let hom_kernel_kernel = hom_over_algebra(kernel_module, kernel_module)?
        .ok_or(PathAlgebraModuleError::InternalMinimalHomRejected)?;
    let exponent = kernel_module
        .dimension_vector()
        .iter()
        .chain(source.dimension_vector())
        .copied()
        .max()
        .unwrap_or(0);

    for source_to_kernel in &hom_source_kernel {
        for kernel_endomorphism in &hom_kernel_kernel {
            let kernel_cycle = kernel
                .compose(source_to_kernel)?
                .compose(kernel_endomorphism)?;
            let fitting_kernel = homomorphism_power(&kernel_cycle, exponent)?;
            if !fitting_kernel.is_zero() {
                let source_cycle = source_to_kernel
                    .compose(kernel_endomorphism)?
                    .compose(&kernel)?;
                let fitting_source = homomorphism_power(&source_cycle, exponent)?;
                return Ok(MoreMinimalVersionResult::Reduced {
                    minimal_part: kernel_inclusion(&fitting_source)?.compose(homomorphism)?,
                    removed_part: image_inclusion(&fitting_source)?.compose(homomorphism)?,
                });
            }
        }
    }
    Ok(MoreMinimalVersionResult::Unchanged(homomorphism.clone()))
}

/// Translate QPA's one-step `MoreLeftMinimalVersion(f)` operation.
///
/// This preserves QPA's cokernel Hom-basis order, its Fitting powers, and the
/// right inverse used to corestrict to the retained kernel summand.
///
/// # Errors
///
/// Propagates cokernel, Hom-space, composition, power, split-inverse, kernel,
/// and image construction failures. An absent internal Fitting retraction is
/// QPA's subsequent map-times-`false` error and is reported explicitly.
pub fn more_left_minimal_version<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<MoreMinimalVersionResult<F>, PathAlgebraModuleError> {
    let range = homomorphism.range();
    let cokernel = cokernel_projection(homomorphism)?;
    let cokernel_module = cokernel.range();
    let hom_cokernel_range = hom_over_algebra(cokernel_module, range)?
        .ok_or(PathAlgebraModuleError::InternalMinimalHomRejected)?;
    if hom_cokernel_range.is_empty() {
        return Ok(MoreMinimalVersionResult::Unchanged(homomorphism.clone()));
    }
    let hom_cokernel_cokernel = hom_over_algebra(cokernel_module, cokernel_module)?
        .ok_or(PathAlgebraModuleError::InternalMinimalHomRejected)?;
    let exponent = range
        .dimension_vector()
        .iter()
        .chain(cokernel_module.dimension_vector())
        .copied()
        .max()
        .unwrap_or(0);

    for cokernel_endomorphism in &hom_cokernel_cokernel {
        for cokernel_to_range in &hom_cokernel_range {
            let cokernel_cycle = cokernel_endomorphism
                .compose(cokernel_to_range)?
                .compose(&cokernel)?;
            let fitting_cokernel = homomorphism_power(&cokernel_cycle, exponent)?;
            if !fitting_cokernel.is_zero() {
                let range_cycle = cokernel
                    .compose(cokernel_endomorphism)?
                    .compose(cokernel_to_range)?;
                let fitting_range = homomorphism_power(&range_cycle, exponent)?;
                let retained_inclusion = kernel_inclusion(&fitting_range)?;
                let retraction = right_inverse_of_homomorphism(&retained_inclusion)?
                    .ok_or(PathAlgebraModuleError::InternalMinimalKernelNotSplit)?;
                return Ok(MoreMinimalVersionResult::Reduced {
                    minimal_part: homomorphism.compose(&retraction)?,
                    removed_part: homomorphism.compose(&image_projection(&fitting_range)?)?,
                });
            }
        }
    }
    Ok(MoreMinimalVersionResult::Unchanged(homomorphism.clone()))
}

/// Translate QPA's cached `RightMinimalVersion(f)` attribute.
///
/// One source summand is removed per successful one-step call and recorded in
/// order. The final map has QPA's `IsRightMinimal` property installed as true.
///
/// # Errors
///
/// Propagates [`more_right_minimal_version`].
pub fn right_minimal_version<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleMinimalVersion<F>, PathAlgebraModuleError> {
    if let Some(cached) = homomorphism.cached_right_minimal_version() {
        return Ok(cached);
    }
    let mut removed_modules = Vec::new();
    let mut minimal = homomorphism.clone();
    loop {
        match more_right_minimal_version(&minimal)? {
            MoreMinimalVersionResult::Unchanged(_) => break,
            MoreMinimalVersionResult::Reduced {
                minimal_part,
                removed_part,
            } => {
                removed_modules.push(removed_part.source().clone());
                minimal = minimal_part;
            }
        }
    }
    let _ = minimal.cache_is_right_minimal(true);
    Ok(homomorphism.cache_right_minimal_version(&minimal, &removed_modules))
}

/// Translate QPA's cached `LeftMinimalVersion(f)` attribute.
///
/// One range summand is removed per successful one-step call and recorded in
/// order. The final map has QPA's `IsLeftMinimal` property installed as true.
///
/// # Errors
///
/// Propagates [`more_left_minimal_version`].
pub fn left_minimal_version<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleMinimalVersion<F>, PathAlgebraModuleError> {
    if let Some(cached) = homomorphism.cached_left_minimal_version() {
        return Ok(cached);
    }
    let mut removed_modules = Vec::new();
    let mut minimal = homomorphism.clone();
    loop {
        match more_left_minimal_version(&minimal)? {
            MoreMinimalVersionResult::Unchanged(_) => break,
            MoreMinimalVersionResult::Reduced {
                minimal_part,
                removed_part,
            } => {
                removed_modules.push(removed_part.range().clone());
                minimal = minimal_part;
            }
        }
    }
    let _ = minimal.cache_is_left_minimal(true);
    Ok(homomorphism.cache_left_minimal_version(&minimal, &removed_modules))
}

fn homomorphism_power<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
    exponent: usize,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    if exponent == 0 {
        return identity_module_mapping(homomorphism.source().clone());
    }
    let mut power = homomorphism.clone();
    for _ in 1..exponent {
        power = power.compose(homomorphism)?;
    }
    Ok(power)
}

fn solve_map_linear_combination<F: Field>(
    basis: &[PathAlgebraModuleHomomorphism<F>],
    target: &PathAlgebraModuleHomomorphism<F>,
) -> Result<Option<Vec<F::Element>>, PathAlgebraModuleError> {
    let field = target.source().right_acting_algebra().field().clone();
    let row_length = target.flat_vertex_matrix_entries().len();
    let row_count = basis.len();
    let entry_count = row_count
        .checked_mul(row_length)
        .ok_or(PathAlgebraModuleError::HomEquationSizeOverflow)?;
    let mut entries = Vec::with_capacity(entry_count);
    for homomorphism in basis {
        entries.extend(homomorphism.flat_vertex_matrix_entries());
    }
    let rows = DenseMatrix::from_flat(row_count, row_length, entries, field)?;
    Ok(rows.solution_row(&target.flat_vertex_matrix_entries())?)
}

fn linear_combination_of_homomorphisms<F: Field>(
    basis: &[PathAlgebraModuleHomomorphism<F>],
    coefficients: &[F::Element],
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    if basis.len() != coefficients.len() {
        return Err(PathAlgebraModuleError::SplittingCoefficientCount {
            actual: coefficients.len(),
            expected: basis.len(),
        });
    }
    let first = basis
        .first()
        .ok_or(PathAlgebraModuleError::EmptySplittingHomBasis)?;
    let field = first.source().right_acting_algebra().field();
    let mut combination = zero_module_mapping(first.source().clone(), first.range().clone())?;
    for (homomorphism, coefficient) in basis.iter().zip(coefficients) {
        if !field.is_zero(coefficient) {
            combination = combination.add(&homomorphism.scale(coefficient)?)?;
        }
    }
    Ok(combination)
}

fn zero_module_for_algebra<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<PathAlgebraModule<F>, PathAlgebraModuleError> {
    right_module_over_path_algebra(
        algebra.clone(),
        vec![0; algebra.quiver().number_of_vertices()],
        &[],
    )?
    .ok_or(PathAlgebraModuleError::GeneratedModuleRejectedRelations)
}

/// Translate QPA's cached `KernelInclusion(f)` attribute.
///
/// Component kernels use the exact GAP `NullspaceMat` row convention. Their
/// arrow actions are reconstructed by solving in those ordered bases, and the
/// zero-kernel branch returns QPA's physical zero mapping.
///
/// # Errors
///
/// Propagates exact matrix operations and checked module/map construction.
pub fn kernel_inclusion<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    if let Some(cached) = homomorphism.cached_kernel_inclusion() {
        return Ok(cached);
    }

    let source = homomorphism.source();
    let algebra = source.right_acting_algebra().clone();
    let field = algebra.field().clone();
    let kernel_bases = source
        .dimension_vector()
        .iter()
        .enumerate()
        .map(|(vertex, &dimension)| {
            if dimension == 0 {
                DenseMatrix::from_flat(0, 0, Vec::new(), field.clone())
            } else {
                Ok(logical_vertex_map(homomorphism, vertex).nullspace_rows())
            }
        })
        .collect::<Result<Vec<_>, MatrixError>>()?;
    let dimensions = kernel_bases
        .iter()
        .map(DenseMatrix::nrows)
        .collect::<Vec<_>>();

    let inclusion = if dimensions.iter().all(|&dimension| dimension == 0) {
        homomorphism.set_is_injective(true);
        zero_module_mapping(zero_module_for_algebra(&algebra)?, source.clone())?
    } else {
        let mut arrow_matrices = Vec::new();
        for arrow in algebra.quiver().arrows() {
            let source_vertex = arrow.source().index().zero_based();
            let target_vertex = arrow.target().index().zero_based();
            if dimensions[source_vertex] == 0 || dimensions[target_vertex] == 0 {
                continue;
            }
            let source_action = logical_submatrix(
                &source.arrow_matrices()[arrow.index().zero_based()],
                source.dimension_vector()[source_vertex],
                source.dimension_vector()[target_vertex],
            )?;
            let acted = kernel_bases[source_vertex].multiply(&source_action)?;
            let mut entries =
                Vec::with_capacity(dimensions[source_vertex] * dimensions[target_vertex]);
            for row in 0..acted.nrows() {
                let coordinates = kernel_bases[target_vertex]
                    .solution_row(acted.row(row))?
                    .ok_or(PathAlgebraModuleError::InternalKernelCoordinates)?;
                entries.extend(coordinates);
            }
            arrow_matrices.push(NamedArrowMatrix::new(
                arrow.name(),
                DenseMatrix::from_flat(
                    dimensions[source_vertex],
                    dimensions[target_vertex],
                    entries,
                    field.clone(),
                )?,
            ));
        }
        let kernel = right_module_over_path_algebra(algebra, dimensions.clone(), &arrow_matrices)?
            .ok_or(PathAlgebraModuleError::GeneratedModuleRejectedRelations)?;
        let inclusion_matrices = kernel_bases
            .iter()
            .enumerate()
            .map(|(vertex, basis)| {
                if dimensions[vertex] == 0 {
                    DenseMatrix::zero(
                        1,
                        physical_dimension(source.dimension_vector()[vertex]),
                        field.clone(),
                    )
                } else {
                    basis.clone()
                }
            })
            .collect();
        right_module_hom_over_algebra(kernel, source.clone(), inclusion_matrices)?
    };
    inclusion.set_kernel_of_what(homomorphism);
    inclusion.set_is_injective(true);
    Ok(homomorphism.cache_kernel_inclusion(&inclusion))
}

/// Return the source of QPA's canonical kernel inclusion.
///
/// # Errors
///
/// Propagates [`kernel_inclusion`].
pub fn kernel_of_module_homomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModule<F>, PathAlgebraModuleError> {
    Ok(kernel_inclusion(homomorphism)?.source().clone())
}

/// Translate QPA's cached image projection/inclusion factorization.
///
/// At each vertex, QPA takes the fully reduced canonical basis of the row
/// image, reconstructs the induced range action, and solves each source-basis
/// image in that basis. The returned pair is projection first, inclusion
/// second, exactly as in `ImageProjectionInclusion`.
///
/// # Errors
///
/// Propagates exact row-space/coordinate operations and checked module/map
/// construction.
#[allow(clippy::too_many_lines)]
pub fn image_projection_inclusion<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<
    (
        PathAlgebraModuleHomomorphism<F>,
        PathAlgebraModuleHomomorphism<F>,
    ),
    PathAlgebraModuleError,
> {
    if let Some(cached) = homomorphism.cached_image_projection_inclusion() {
        return Ok(cached);
    }

    let source = homomorphism.source();
    let range = homomorphism.range();
    let algebra = source.right_acting_algebra().clone();
    let field = algebra.field().clone();
    let image_bases = source
        .dimension_vector()
        .iter()
        .enumerate()
        .map(|(vertex, _)| logical_vertex_map(homomorphism, vertex).canonical_row_space_basis())
        .collect::<Vec<_>>();
    let dimensions = image_bases
        .iter()
        .map(DenseMatrix::nrows)
        .collect::<Vec<_>>();

    if dimensions.iter().all(|&dimension| dimension == 0) {
        let zero = zero_module_for_algebra(&algebra)?;
        let projection = zero_module_mapping(source.clone(), zero.clone())?;
        let inclusion = zero_module_mapping(zero, range.clone())?;
        return Ok(homomorphism.cache_image_projection_inclusion(&projection, &inclusion));
    }

    let mut arrow_matrices = Vec::new();
    for arrow in algebra.quiver().arrows() {
        let source_vertex = arrow.source().index().zero_based();
        let target_vertex = arrow.target().index().zero_based();
        if dimensions[source_vertex] == 0 || dimensions[target_vertex] == 0 {
            continue;
        }
        let range_action = logical_submatrix(
            &range.arrow_matrices()[arrow.index().zero_based()],
            range.dimension_vector()[source_vertex],
            range.dimension_vector()[target_vertex],
        )?;
        let acted = image_bases[source_vertex].multiply(&range_action)?;
        let mut entries = Vec::with_capacity(dimensions[source_vertex] * dimensions[target_vertex]);
        for row in 0..acted.nrows() {
            let coordinates = image_bases[target_vertex]
                .solution_row(acted.row(row))?
                .ok_or(PathAlgebraModuleError::InternalImageCoordinates)?;
            entries.extend(coordinates);
        }
        arrow_matrices.push(NamedArrowMatrix::new(
            arrow.name(),
            DenseMatrix::from_flat(
                dimensions[source_vertex],
                dimensions[target_vertex],
                entries,
                field.clone(),
            )?,
        ));
    }
    let image = right_module_over_path_algebra(algebra, dimensions.clone(), &arrow_matrices)?
        .ok_or(PathAlgebraModuleError::GeneratedModuleRejectedRelations)?;

    let inclusion_matrices = image_bases
        .iter()
        .enumerate()
        .map(|(vertex, basis)| {
            if dimensions[vertex] == 0 {
                DenseMatrix::zero(
                    1,
                    physical_dimension(range.dimension_vector()[vertex]),
                    field.clone(),
                )
            } else {
                basis.clone()
            }
        })
        .collect();
    let inclusion =
        right_module_hom_over_algebra(image.clone(), range.clone(), inclusion_matrices)?;
    inclusion.set_image_of_what(homomorphism);
    inclusion.set_is_injective(true);

    let projection_matrices = source
        .dimension_vector()
        .iter()
        .enumerate()
        .map(
            |(vertex, &source_dimension)| -> Result<DenseMatrix<F>, PathAlgebraModuleError> {
                if source_dimension == 0 || dimensions[vertex] == 0 {
                    return Ok(DenseMatrix::zero(
                        physical_dimension(source_dimension),
                        physical_dimension(dimensions[vertex]),
                        field.clone(),
                    ));
                }
                let logical_map = logical_vertex_map(homomorphism, vertex);
                let mut entries = Vec::with_capacity(source_dimension * dimensions[vertex]);
                for row in 0..source_dimension {
                    let coordinates = image_bases[vertex]
                        .solution_row(logical_map.row(row))?
                        .ok_or(PathAlgebraModuleError::InternalImageCoordinates)?;
                    entries.extend(coordinates);
                }
                Ok(DenseMatrix::from_flat(
                    source_dimension,
                    dimensions[vertex],
                    entries,
                    field.clone(),
                )?)
            },
        )
        .collect::<Result<Vec<_>, _>>()?;
    let projection = right_module_hom_over_algebra(source.clone(), image, projection_matrices)?;
    projection.set_image_of_what(homomorphism);
    projection.set_is_surjective(true);

    Ok(homomorphism.cache_image_projection_inclusion(&projection, &inclusion))
}

/// Return QPA's canonical projection onto the image.
///
/// # Errors
///
/// Propagates [`image_projection_inclusion`].
pub fn image_projection<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    Ok(image_projection_inclusion(homomorphism)?.0)
}

/// Return QPA's canonical inclusion of the image.
///
/// # Errors
///
/// Propagates [`image_projection_inclusion`].
pub fn image_inclusion<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    Ok(image_projection_inclusion(homomorphism)?.1)
}

/// Return QPA's image module object.
///
/// # Errors
///
/// Propagates [`image_projection_inclusion`].
pub fn image_of_module_homomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModule<F>, PathAlgebraModuleError> {
    Ok(image_projection(homomorphism)?.range().clone())
}

/// Translate QPA's linear-equation `HomOverAlgebra(M, N)` method.
///
/// Unknown vertex-map entries are ordered vertex-major, then row-major. The
/// equation columns follow supported source vertices, their outgoing arrows,
/// source rows, and target coordinates exactly as in `modulehom.gi`.
/// `Ok(None)` is QPA's printed-then-`fail` result for modules over different
/// algebras.
///
/// # Errors
///
/// Reports address-size overflow and propagates exact matrix or checked-map
/// construction failures.
#[allow(clippy::too_many_lines)]
pub fn hom_over_algebra<F: Field>(
    source: &PathAlgebraModule<F>,
    range: &PathAlgebraModule<F>,
) -> Result<Option<Vec<PathAlgebraModuleHomomorphism<F>>>, PathAlgebraModuleError> {
    if source.right_acting_algebra() != range.right_acting_algebra() {
        return Ok(None);
    }
    let algebra = source.right_acting_algebra();
    let field = algebra.field().clone();
    let source_dimensions = source.dimension_vector();
    let range_dimensions = range.dimension_vector();
    let vertex_count = source_dimensions.len();
    let source_support = source_dimensions
        .iter()
        .map(|&dimension| dimension != 0)
        .collect::<Vec<_>>();
    let range_support = range_dimensions
        .iter()
        .map(|&dimension| dimension != 0)
        .collect::<Vec<_>>();

    let mut row_count = 0_usize;
    let mut block_rows = vec![None; vertex_count];
    for vertex in 0..vertex_count {
        if source_support[vertex] {
            block_rows[vertex] = Some(row_count);
            row_count = row_count
                .checked_add(
                    source_dimensions[vertex]
                        .checked_mul(range_dimensions[vertex])
                        .ok_or(PathAlgebraModuleError::HomEquationSizeOverflow)?,
                )
                .ok_or(PathAlgebraModuleError::HomEquationSizeOverflow)?;
        }
    }

    let arrows = algebra.quiver().arrows().collect::<Vec<_>>();
    let mut column_count = 0_usize;
    for vertex in 0..vertex_count {
        if !source_support[vertex] {
            continue;
        }
        for arrow in arrows
            .iter()
            .filter(|arrow| arrow.source().index().zero_based() == vertex)
        {
            let arrow_source = arrow.source().index().zero_based();
            let arrow_target = arrow.target().index().zero_based();
            if range_support[arrow_target]
                && (range_support[arrow_source] || source_support[arrow_target])
            {
                column_count = column_count
                    .checked_add(
                        source_dimensions[arrow_source]
                            .checked_mul(range_dimensions[arrow_target])
                            .ok_or(PathAlgebraModuleError::HomEquationSizeOverflow)?,
                    )
                    .ok_or(PathAlgebraModuleError::HomEquationSizeOverflow)?;
            }
        }
    }

    if row_count != 0 && column_count != 0 {
        let mut equations = DenseMatrix::zero(row_count, column_count, field.clone());
        let mut previous_column = 0_usize;
        for vertex in 0..vertex_count {
            if !source_support[vertex] {
                continue;
            }
            for arrow in arrows
                .iter()
                .filter(|arrow| arrow.source().index().zero_based() == vertex)
            {
                let arrow_source = arrow.source().index().zero_based();
                let arrow_target = arrow.target().index().zero_based();
                if !range_support[arrow_target]
                    || !(range_support[arrow_source] || source_support[arrow_target])
                {
                    continue;
                }
                let source_action = logical_submatrix(
                    &source.arrow_matrices()[arrow.index().zero_based()],
                    source_dimensions[arrow_source],
                    source_dimensions[arrow_target],
                )?;
                let range_action = logical_submatrix(
                    &range.arrow_matrices()[arrow.index().zero_based()],
                    range_dimensions[arrow_source],
                    range_dimensions[arrow_target],
                )?;
                for source_row in 0..source_dimensions[arrow_source] {
                    let column_start =
                        previous_column + source_row * range_dimensions[arrow_target];
                    if range_support[arrow_source] {
                        let row_start = block_rows[arrow_source].ok_or(
                            PathAlgebraModuleError::MissingHomUnknownBlock {
                                vertex: arrow_source + 1,
                            },
                        )? + source_row * range_dimensions[arrow_source];
                        for range_row in 0..range_dimensions[arrow_source] {
                            for target_column in 0..range_dimensions[arrow_target] {
                                equations[(row_start + range_row, column_start + target_column)] =
                                    range_action[(range_row, target_column)].clone();
                            }
                        }
                    }
                    if source_support[arrow_target] {
                        let target_row_start = block_rows[arrow_target].ok_or(
                            PathAlgebraModuleError::MissingHomUnknownBlock {
                                vertex: arrow_target + 1,
                            },
                        )?;
                        for target_row in 0..source_dimensions[arrow_target] {
                            let negative = field.neg(&source_action[(source_row, target_row)]);
                            for target_column in 0..range_dimensions[arrow_target] {
                                field.add_assign(
                                    &mut equations[(
                                        target_row_start
                                            + target_row * range_dimensions[arrow_target]
                                            + target_column,
                                        column_start + target_column,
                                    )],
                                    &negative,
                                );
                            }
                        }
                    }
                }
                previous_column += source_dimensions[arrow_source] * range_dimensions[arrow_target];
            }
        }

        let nullspace = equations.nullspace_rows();
        let mut homomorphisms = Vec::with_capacity(nullspace.nrows());
        for row in 0..nullspace.nrows() {
            homomorphisms.push(module_homomorphism_from_logical_flat(
                source,
                range,
                nullspace.row(row),
            )?);
        }
        return Ok(Some(homomorphisms));
    }

    if source.dimension() == 0 || range.dimension() == 0 {
        return Ok(Some(Vec::new()));
    }

    let mut homomorphisms = Vec::new();
    for vertex in 0..vertex_count {
        if !source_support[vertex] || !range_support[vertex] {
            continue;
        }
        let block_start = block_rows[vertex]
            .ok_or(PathAlgebraModuleError::MissingHomUnknownBlock { vertex: vertex + 1 })?;
        for row in 0..source_dimensions[vertex] {
            for column in 0..range_dimensions[vertex] {
                let mut entries = vec![field.zero(); row_count];
                entries[block_start + row * range_dimensions[vertex] + column] = field.one();
                homomorphisms.push(module_homomorphism_from_logical_flat(
                    source, range, &entries,
                )?);
            }
        }
    }
    Ok(Some(homomorphisms))
}

/// Translate QPA's `HomOverAlgebraWithBasisFunction` result pair.
///
/// `Ok(None)` is the underlying `HomOverAlgebra` failure for different acting
/// algebras.
///
/// # Errors
///
/// Propagates [`hom_over_algebra`] and exact flattened-basis construction.
pub fn hom_over_algebra_with_basis_function<F: Field>(
    source: &PathAlgebraModule<F>,
    range: &PathAlgebraModule<F>,
) -> Result<Option<PathAlgebraModuleHomomorphismBasis<F>>, PathAlgebraModuleError> {
    let Some(basis) = hom_over_algebra(source, range)? else {
        return Ok(None);
    };
    let flattened_basis = if let Some(first) = basis.first() {
        let column_count = flatten_module_homomorphism(first).len();
        Some(DenseMatrix::from_flat(
            basis.len(),
            column_count,
            basis.iter().flat_map(flatten_module_homomorphism).collect(),
            source.right_acting_algebra().field().clone(),
        )?)
    } else {
        None
    };
    Ok(Some(PathAlgebraModuleHomomorphismBasis {
        basis,
        flattened_basis,
        field: source.right_acting_algebra().field().clone(),
    }))
}

/// Translate QPA's cached `EndOverAlgebra(M)` matrix-algebra construction.
///
/// QPA embeds each vertex map as a diagonal block in the underlying
/// vertex-major vector space, forms the direct sum of full matrix algebras,
/// and then asks GAP for the basis of the generated subalgebra with one. The
/// latter is GAP's insertion-ordered semi-echelon basis of the flattened block
/// matrices.
///
/// # Errors
///
/// Propagates [`hom_over_algebra`], address-size overflow, and exact matrix
/// construction failures.
pub fn end_over_algebra<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<ModuleEndomorphismAlgebra<F>, PathAlgebraModuleError> {
    if let Some(cached) = module.cached_endomorphism_algebra() {
        return Ok(cached);
    }
    if module.dimension() == 0 {
        return Err(PathAlgebraModuleError::EmptyEndomorphismDirectSum);
    }
    let homomorphisms = hom_over_algebra(module, module)?
        .ok_or(PathAlgebraModuleError::EndomorphismAlgebraMismatch)?;
    let dimension = module.dimension();
    let flattened_dimension = dimension
        .checked_mul(dimension)
        .ok_or(PathAlgebraModuleError::HomEquationSizeOverflow)?;
    let field = module.right_acting_algebra().field().clone();
    let embedded = homomorphisms
        .iter()
        .map(from_hom_mm_to_end_m)
        .collect::<Result<Vec<_>, _>>()?;
    let generator_rows = DenseMatrix::from_flat(
        embedded.len(),
        flattened_dimension,
        embedded.iter().flat_map(flatten_matrix).collect(),
        field.clone(),
    )?;
    let flattened_basis = generator_rows.semi_echelon_transformation().vectors;
    let basis_matrices = (0..flattened_basis.nrows())
        .map(|row| {
            DenseMatrix::from_flat(
                dimension,
                dimension,
                flattened_basis.row(row).to_vec(),
                field.clone(),
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let algebra = ModuleEndomorphismAlgebra(Arc::new(ModuleEndomorphismAlgebraData {
        field: field.clone(),
        basis_matrices,
        flattened_basis,
        one: DenseMatrix::identity(dimension, field),
    }));
    Ok(module.cache_endomorphism_algebra(&algebra))
}

/// Translate QPA's `FromHomMMToEndM(f)` block-diagonal embedding.
///
/// # Errors
///
/// The source and range must be structurally equal modules; exact matrix
/// construction failures propagate.
pub fn from_hom_mm_to_end_m<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<DenseMatrix<F>, PathAlgebraModuleError> {
    if homomorphism.source() != homomorphism.range() {
        return Err(PathAlgebraModuleError::NonEndomorphismConversion);
    }
    let dimensions = homomorphism.source().dimension_vector();
    let total = homomorphism.source().dimension();
    let field = homomorphism.source().right_acting_algebra().field().clone();
    let offsets = dimension_offsets(dimensions);
    let mut matrix = DenseMatrix::zero(total, total, field);
    for (vertex, &dimension) in dimensions.iter().enumerate() {
        if dimension == 0 {
            continue;
        }
        let block = logical_vertex_map(homomorphism, vertex);
        for row in 0..dimension {
            for column in 0..dimension {
                matrix[(offsets[vertex] + row, offsets[vertex] + column)] =
                    block[(row, column)].clone();
            }
        }
    }
    Ok(matrix)
}

/// Translate QPA's `FromEndMToHomMM(M, mat)` diagonal-block extraction.
///
/// Entries outside the vertex-diagonal blocks are ignored, just as in QPA.
/// A matrix larger than the underlying module dimension is therefore
/// accepted, while a smaller or wrong-field matrix is rejected explicitly.
///
/// # Errors
///
/// Propagates field/size validation and QPA's checked module-map constructor.
pub fn from_end_m_to_hom_mm<F: Field>(
    module: &PathAlgebraModule<F>,
    matrix: &DenseMatrix<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    let field = module.right_acting_algebra().field().clone();
    if matrix.field() != &field {
        return Err(PathAlgebraModuleError::EndomorphismMatrixField);
    }
    let total = module.dimension();
    if matrix.nrows() < total || matrix.ncols() < total {
        return Err(PathAlgebraModuleError::EndomorphismMatrixTooSmall {
            actual: (matrix.nrows(), matrix.ncols()),
            minimum: (total, total),
        });
    }
    let offsets = dimension_offsets(module.dimension_vector());
    let matrices = module
        .dimension_vector()
        .iter()
        .enumerate()
        .map(|(vertex, &dimension)| {
            if dimension == 0 {
                return Ok(DenseMatrix::zero(1, 1, field.clone()));
            }
            let offset = offsets[vertex];
            let mut entries = Vec::with_capacity(dimension * dimension);
            for row in offset..offset + dimension {
                entries.extend_from_slice(&matrix.row(row)[offset..offset + dimension]);
            }
            Ok(DenseMatrix::from_flat(
                dimension,
                dimension,
                entries,
                field.clone(),
            )?)
        })
        .collect::<Result<Vec<_>, PathAlgebraModuleError>>()?;
    right_module_hom_over_algebra(module.clone(), module.clone(), matrices)
}

/// Construct QPA's natural quotient homomorphism by the two-sided ideal
/// generated by represented endomorphism matrices.
///
/// Generator order is retained through GAP's insertion-ordered semi-echelon
/// basis algorithm. The ideal is closed under left and right multiplication
/// by the parent basis before the natural row-space quotient is formed.
///
/// # Errors
///
/// Every generator and generated product must lie in the parent algebra;
/// exact multiplication, coordinate, ideal-closure, and quotient failures
/// propagate.
pub fn natural_endomorphism_quotient_by_ideal<F: Field>(
    parent: &ModuleEndomorphismAlgebra<F>,
    generators: &[DenseMatrix<F>],
) -> Result<ModuleEndomorphismQuotientHomomorphism<F>, PathAlgebraModuleError> {
    let ideal_coordinates = endomorphism_ideal_basis_coordinates(parent, generators)?;
    let ideal_basis_matrices = (0..ideal_coordinates.nrows())
        .map(|row| parent.matrix_from_coefficients(ideal_coordinates.row(row)))
        .collect::<Result<Vec<_>, _>>()?;
    let quotient = natural_quotient_by_row_space(&ideal_coordinates)?;
    let range =
        ModuleEndomorphismQuotientAlgebra(Arc::new(ModuleEndomorphismQuotientAlgebraData {
            parent: parent.clone(),
            ideal_basis_matrices,
            projection: quotient.projection,
            lifts: quotient.lifts,
        }));
    Ok(ModuleEndomorphismQuotientHomomorphism {
        source: parent.clone(),
        range,
    })
}

fn endomorphism_ideal_basis_coordinates<F: Field>(
    parent: &ModuleEndomorphismAlgebra<F>,
    generators: &[DenseMatrix<F>],
) -> Result<DenseMatrix<F>, PathAlgebraModuleError> {
    let dimension = parent.dimension();
    let capacity = generators
        .len()
        .checked_mul(dimension)
        .ok_or(PathAlgebraModuleError::HomEquationSizeOverflow)?;
    let mut generator_coordinates = Vec::with_capacity(capacity);
    for generator in generators {
        generator_coordinates.extend(
            parent
                .coefficients(generator)?
                .ok_or(PathAlgebraModuleError::EndomorphismIdealGeneratorOutsideAlgebra)?,
        );
    }
    let mut basis = DenseMatrix::from_flat(
        generators.len(),
        dimension,
        generator_coordinates,
        parent.field().clone(),
    )?
    .semi_echelon_transformation()
    .vectors;

    loop {
        let ideal_matrices = (0..basis.nrows())
            .map(|row| parent.matrix_from_coefficients(basis.row(row)))
            .collect::<Result<Vec<_>, _>>()?;
        let product_rows = ideal_matrices
            .len()
            .checked_mul(parent.basis_matrices().len())
            .and_then(|count| count.checked_mul(2))
            .ok_or(PathAlgebraModuleError::HomEquationSizeOverflow)?;
        let expanded_rows = basis
            .nrows()
            .checked_add(product_rows)
            .ok_or(PathAlgebraModuleError::HomEquationSizeOverflow)?;
        let expanded_capacity = expanded_rows
            .checked_mul(dimension)
            .ok_or(PathAlgebraModuleError::HomEquationSizeOverflow)?;
        let mut expanded_coordinates = Vec::with_capacity(expanded_capacity);
        expanded_coordinates.extend_from_slice(basis.as_flat());
        for ideal_element in &ideal_matrices {
            for algebra_element in parent.basis_matrices() {
                for product in [
                    parent.multiply(algebra_element, ideal_element)?,
                    parent.multiply(ideal_element, algebra_element)?,
                ] {
                    expanded_coordinates.extend(
                        parent.coefficients(&product)?.ok_or(
                            PathAlgebraModuleError::EndomorphismIdealProductOutsideAlgebra,
                        )?,
                    );
                }
            }
        }
        let expanded = DenseMatrix::from_flat(
            expanded_rows,
            dimension,
            expanded_coordinates,
            parent.field().clone(),
        )?
        .semi_echelon_transformation()
        .vectors;
        if expanded.nrows() == basis.nrows() {
            return Ok(basis);
        }
        basis = expanded;
    }
}

fn module_homomorphism_from_logical_flat<F: Field>(
    source: &PathAlgebraModule<F>,
    range: &PathAlgebraModule<F>,
    entries: &[F::Element],
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    let field = source.right_acting_algebra().field().clone();
    let mut offset = 0_usize;
    let matrices = source
        .dimension_vector()
        .iter()
        .zip(range.dimension_vector())
        .map(
            |(&source_dimension, &range_dimension)| -> Result<DenseMatrix<F>, MatrixError> {
                if source_dimension == 0 || range_dimension == 0 {
                    return Ok(DenseMatrix::zero(
                        physical_dimension(source_dimension),
                        physical_dimension(range_dimension),
                        field.clone(),
                    ));
                }
                let length = source_dimension * range_dimension;
                let matrix = DenseMatrix::from_flat(
                    source_dimension,
                    range_dimension,
                    entries[offset..offset + length].to_vec(),
                    field.clone(),
                )?;
                offset += length;
                Ok(matrix)
            },
        )
        .collect::<Result<Vec<_>, _>>()?;
    right_module_hom_over_algebra(source.clone(), range.clone(), matrices)
}

fn flatten_module_homomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Vec<F::Element> {
    homomorphism.flat_vertex_matrix_entries()
}

fn flatten_matrix<F: Field>(matrix: &DenseMatrix<F>) -> Vec<F::Element> {
    (0..matrix.nrows())
        .flat_map(|row| matrix.row(row).iter().cloned())
        .collect()
}

/// Translate QPA's `SubRepresentationInclusion(M, generators)` method.
///
/// QPA first makes every generator vertex-uniform, repeatedly acts by every
/// arrow until no nonzero elements remain, takes GAP's fully reduced canonical
/// row-space basis, and then reconstructs both the subrepresentation and its
/// inclusion.  This function follows that sequence and ordering.
///
/// # Errors
///
/// Rejects generators outside the identical module family and propagates
/// algebra, matrix, and checked module-construction failures.
#[allow(clippy::too_many_lines)]
pub fn subrepresentation_inclusion<F: Field>(
    module: PathAlgebraModule<F>,
    generators: &[PathAlgebraModuleElement<F>],
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    if generators
        .iter()
        .any(|generator| !module.is_identical(generator.module()))
    {
        return Err(PathAlgebraModuleError::ForeignSubrepresentationGenerator);
    }

    let algebra = module.right_acting_algebra().clone();
    let vertex_count = module.dimension_vector().len();
    if generators.is_empty() {
        let zero = right_module_over_path_algebra(algebra, vec![0; vertex_count], &[])?
            .ok_or(PathAlgebraModuleError::GeneratedModuleRejectedRelations)?;
        return zero_module_mapping(zero, module);
    }

    let mut uniform_generators = Vec::new();
    for generator in generators {
        for vertex in 0..vertex_count {
            let mut uniform = module.zero_element();
            uniform.coordinates[vertex].clone_from(&generator.coordinates[vertex]);
            if !uniform.is_zero() {
                uniform_generators.push(uniform);
            }
        }
    }

    let arrow_elements = algebra
        .quiver()
        .arrows()
        .map(|arrow| quiver_algebra_path_element(&algebra, arrow.path()))
        .collect::<Result<Vec<_>, _>>()?;
    let mut spanning_elements = Vec::new();
    let mut new = uniform_generators;
    while !new.is_empty() {
        let mut next = Vec::new();
        for element in &new {
            for arrow in &arrow_elements {
                let image = element.act(arrow)?;
                if !image.is_zero() {
                    next.push(image);
                }
            }
        }
        spanning_elements.extend(new);
        new = next;
    }

    let field = algebra.field().clone();
    let spanning_rows = DenseMatrix::from_flat(
        spanning_elements.len(),
        module.dimension(),
        spanning_elements
            .iter()
            .flat_map(PathAlgebraModuleElement::basis_coefficients)
            .collect(),
        field.clone(),
    )?;
    let canonical_basis = spanning_rows.canonical_row_space_basis();
    let basis_elements = (0..canonical_basis.nrows())
        .map(|row| {
            PathAlgebraModuleElement::from_logical_flat(module.clone(), canonical_basis.row(row))
        })
        .collect::<Result<Vec<_>, _>>()?;

    let mut sub_dimension = vec![0; vertex_count];
    for element in &basis_elements {
        for vertex in element.support_vertices() {
            sub_dimension[vertex - 1] += 1;
        }
    }
    if canonical_basis.nrows() != sub_dimension.iter().sum::<usize>() {
        return Err(PathAlgebraModuleError::NonuniformSubrepresentationBasis);
    }
    let sub_offsets = dimension_offsets(&sub_dimension);
    let module_offsets = dimension_offsets(module.dimension_vector());

    let mut arrow_matrices = Vec::new();
    for (arrow, arrow_element) in algebra.quiver().arrows().zip(&arrow_elements) {
        let source = arrow.source().index().zero_based();
        let target = arrow.target().index().zero_based();
        if sub_dimension[source] == 0 || sub_dimension[target] == 0 {
            continue;
        }
        let mut rows = Vec::with_capacity(sub_dimension[source] * sub_dimension[target]);
        for basis_element in basis_elements
            .iter()
            .skip(sub_offsets[source])
            .take(sub_dimension[source])
        {
            let image = basis_element.act(arrow_element)?;
            let coefficients = canonical_basis
                .solution_row(&image.basis_coefficients())?
                .ok_or(PathAlgebraModuleError::InternalSubrepresentationCoordinates)?;
            rows.extend_from_slice(
                &coefficients[sub_offsets[target]..sub_offsets[target] + sub_dimension[target]],
            );
        }
        arrow_matrices.push(NamedArrowMatrix::new(
            arrow.name(),
            DenseMatrix::from_flat(
                sub_dimension[source],
                sub_dimension[target],
                rows,
                field.clone(),
            )?,
        ));
    }

    let submodule =
        right_module_over_path_algebra(algebra, sub_dimension.clone(), &arrow_matrices)?
            .ok_or(PathAlgebraModuleError::GeneratedModuleRejectedRelations)?;

    let inclusion_matrices = (0..vertex_count)
        .map(|vertex| -> Result<DenseMatrix<F>, MatrixError> {
            if sub_dimension[vertex] == 0 || module.dimension_vector()[vertex] == 0 {
                return Ok(DenseMatrix::zero(
                    physical_dimension(sub_dimension[vertex]),
                    physical_dimension(module.dimension_vector()[vertex]),
                    field.clone(),
                ));
            }
            let mut entries =
                Vec::with_capacity(sub_dimension[vertex] * module.dimension_vector()[vertex]);
            for row in sub_offsets[vertex]..sub_offsets[vertex] + sub_dimension[vertex] {
                entries.extend_from_slice(
                    &canonical_basis.row(row)[module_offsets[vertex]
                        ..module_offsets[vertex] + module.dimension_vector()[vertex]],
                );
            }
            DenseMatrix::from_flat(
                sub_dimension[vertex],
                module.dimension_vector()[vertex],
                entries,
                field.clone(),
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    right_module_hom_over_algebra(submodule, module, inclusion_matrices)
}

/// Translate QPA's `SubRepresentation(M, generators)` convenience method.
///
/// # Errors
///
/// Propagates [`subrepresentation_inclusion`].
pub fn subrepresentation<F: Field>(
    module: PathAlgebraModule<F>,
    generators: &[PathAlgebraModuleElement<F>],
) -> Result<PathAlgebraModule<F>, PathAlgebraModuleError> {
    Ok(subrepresentation_inclusion(module, generators)?
        .source()
        .clone())
}

/// Translate QPA's `RadicalOfModuleInclusion` method for finite path algebras
/// and finite admissible quotients.
///
/// # Errors
///
/// Reports the same source-domain restriction and propagates generated
/// subrepresentation failures.
pub fn radical_of_module_inclusion<F: Field>(
    module: PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    let algebra = module.right_acting_algebra();
    let supported = match algebra {
        QuiverAlgebra::Path(algebra) => algebra.has_finite_canonical_basis(),
        QuiverAlgebra::Quotient(algebra) => {
            algebra.is_finite_dimensional()? && algebra.is_admissible()?
        }
    };
    if !supported {
        return Err(PathAlgebraModuleError::UnsupportedRadicalAlgebra);
    }

    let basis = module.basis_elements();
    let arrows = algebra
        .quiver()
        .arrows()
        .map(|arrow| quiver_algebra_path_element(algebra, arrow.path()))
        .collect::<Result<Vec<_>, _>>()?;
    let mut generators = Vec::new();
    for arrow in &arrows {
        for element in &basis {
            let image = element.act(arrow)?;
            if !image.is_zero() && !generators.contains(&image) {
                generators.push(image);
            }
        }
    }
    subrepresentation_inclusion(module, &generators)
}

/// Translate QPA's `RadicalOfModule(M)` convenience method.
///
/// # Errors
///
/// Propagates [`radical_of_module_inclusion`].
pub fn radical_of_module<F: Field>(
    module: PathAlgebraModule<F>,
) -> Result<PathAlgebraModule<F>, PathAlgebraModuleError> {
    Ok(radical_of_module_inclusion(module)?.source().clone())
}

/// Translate QPA's cached `IsSimpleQPAModule(M)` dimension test.
///
/// QPA defines this property to be true exactly in total dimension one. A
/// true result also installs its declared implications to
/// `IsIndecomposableModule` and `IsSemisimpleModule`.
#[must_use]
pub fn is_simple_qpa_module<F: Field>(module: &PathAlgebraModule<F>) -> bool {
    if let Some(cached) = module.cached_is_simple_qpa_module() {
        return cached;
    }
    module.cache_is_simple_qpa_module(module.dimension() == 1)
}

/// Translate QPA's cached `IsSemisimpleModule(M)` radical test.
///
/// # Errors
///
/// Propagates the exact radical-module construction.
pub fn is_semisimple_module<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<bool, PathAlgebraModuleError> {
    if let Some(cached) = module.cached_is_semisimple_module() {
        return Ok(cached);
    }
    let radical = radical_of_module(module.clone())?;
    Ok(module.cache_is_semisimple_module(radical.dimension() == 0))
}

/// Translate QPA's module `LoewyLength` method by repeated radicals.
///
/// # Errors
///
/// Propagates the source-domain and generated-radical failures of
/// [`radical_of_module`].
pub fn loewy_length_of_module<F: Field>(
    mut module: PathAlgebraModule<F>,
) -> Result<usize, PathAlgebraModuleError> {
    let mut length = 0;
    while module.dimension() != 0 {
        module = radical_of_module(module)?;
        length += 1;
    }
    Ok(length)
}

/// Translate QPA's `CoKernelProjection` for a path-algebra module map.
///
/// At each vertex this uses GAP's `NaturalHomomorphismBySubspace`: a
/// semi-echelon image basis is greedily extended by ambient standard vectors,
/// and those selected standard vectors are retained as canonical quotient
/// basis lifts.
///
/// # Errors
///
/// Propagates exact quotient-space, module-construction, and commutativity
/// failures.
pub fn cokernel_projection<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    if let Some(cached) = homomorphism.cached_cokernel_projection() {
        return Ok(cached);
    }
    let range = homomorphism.range();
    let algebra = range.right_acting_algebra().clone();
    let field = algebra.field().clone();
    let mut quotients = Vec::with_capacity(range.dimension_vector().len());
    for (vertex, &range_dimension) in range.dimension_vector().iter().enumerate() {
        let source_dimension = homomorphism.source().dimension_vector()[vertex];
        let mut image_entries = Vec::with_capacity(source_dimension * range_dimension);
        for row in 0..source_dimension {
            image_entries.extend_from_slice(
                &homomorphism.vertex_matrices()[vertex].row(row)[..range_dimension],
            );
        }
        let image_rows = DenseMatrix::from_flat(
            source_dimension,
            range_dimension,
            image_entries,
            field.clone(),
        )?;
        quotients.push(natural_quotient_by_row_space(&image_rows)?);
    }
    let quotient_dimension = quotients
        .iter()
        .map(|quotient| quotient.lifts.nrows())
        .collect::<Vec<_>>();

    let mut arrow_matrices = Vec::new();
    for arrow in algebra.quiver().arrows() {
        let source = arrow.source().index().zero_based();
        let target = arrow.target().index().zero_based();
        if quotient_dimension[source] == 0 || quotient_dimension[target] == 0 {
            continue;
        }
        let range_action = logical_submatrix(
            &range.arrow_matrices()[arrow.index().zero_based()],
            range.dimension_vector()[source],
            range.dimension_vector()[target],
        )?;
        let action = quotients[source]
            .lifts
            .multiply(&range_action)?
            .multiply(&quotients[target].projection)?;
        arrow_matrices.push(NamedArrowMatrix::new(arrow.name(), action));
    }

    let cokernel =
        right_module_over_path_algebra(algebra, quotient_dimension.clone(), &arrow_matrices)?
            .ok_or(PathAlgebraModuleError::GeneratedModuleRejectedRelations)?;
    let projection_matrices = range
        .dimension_vector()
        .iter()
        .enumerate()
        .map(|(vertex, &range_dimension)| {
            let quotient_dimension = quotient_dimension[vertex];
            if range_dimension == 0 || quotient_dimension == 0 {
                DenseMatrix::zero(
                    physical_dimension(range_dimension),
                    physical_dimension(quotient_dimension),
                    field.clone(),
                )
            } else {
                quotients[vertex].projection.clone()
            }
        })
        .collect();
    let projection = right_module_hom_over_algebra(range.clone(), cokernel, projection_matrices)?;
    projection.set_cokernel_of_what(homomorphism);
    projection.set_is_surjective(true);
    Ok(homomorphism.cache_cokernel_projection(&projection))
}

/// Translate QPA's cokernel-object convenience method as the range of its
/// canonical projection.
///
/// # Errors
///
/// Propagates [`cokernel_projection`].
pub fn cokernel_of_module_homomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModule<F>, PathAlgebraModuleError> {
    Ok(cokernel_projection(homomorphism)?.range().clone())
}

/// Translate QPA's `TopOfModuleProjection` attribute computation.
///
/// The total-zero branch returns the zero endomorphism of the identical input
/// module, exactly as the source method does.
///
/// # Errors
///
/// Propagates radical and cokernel construction failures.
#[allow(clippy::needless_pass_by_value)]
pub fn top_of_module_projection<F: Field>(
    module: PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    if let Some(cached) = module.cached_top_projection() {
        return Ok(cached);
    }
    let projection = if module.dimension() == 0 {
        zero_module_mapping(module.clone(), module.clone())?
    } else {
        let generated = cokernel_projection(&radical_of_module_inclusion(module.clone())?)?;
        let top = module.cache_top_module(generated.range());
        if top.is_identical(generated.range()) {
            generated
        } else {
            let projection = right_module_hom_over_algebra(
                module.clone(),
                top,
                generated.vertex_matrices().to_vec(),
            )?;
            projection.set_is_surjective(true);
            projection
        }
    };
    Ok(module.cache_top_projection(&projection))
}

/// Translate QPA's `TopOfModule(M)` convenience method.
///
/// # Errors
///
/// Propagates [`top_of_module_projection`].
pub fn top_of_module<F: Field>(
    module: PathAlgebraModule<F>,
) -> Result<PathAlgebraModule<F>, PathAlgebraModuleError> {
    if module.dimension() == 0 {
        return Ok(module);
    }
    if let Some(cached) = module.cached_top_module() {
        return Ok(cached);
    }
    Ok(top_of_module_projection(module)?.range().clone())
}

/// Translate QPA's `TopOfModule(f)` overload for a module homomorphism.
///
/// QPA lifts each vertex-major basis vector of the source top through its
/// canonical top projection, applies `f`, and projects to the range top. This
/// preserves the source method's chosen quotient representatives rather than
/// constructing an abstract induced linear map.
///
/// # Errors
///
/// Propagates top construction, preimage, image, and checked-map failures.
pub fn top_of_module_homomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    let source_projection = top_of_module_projection(homomorphism.source().clone())?;
    let range_projection = top_of_module_projection(homomorphism.range().clone())?;
    let source_top = source_projection.range();
    let range_top = range_projection.range();
    let field = source_top.right_acting_algebra().field().clone();
    let basis = source_top.basis_elements();
    let offsets = dimension_offsets(source_top.dimension_vector());
    let matrices = source_top
        .dimension_vector()
        .iter()
        .enumerate()
        .map(
            |(vertex, &source_dimension)| -> Result<DenseMatrix<F>, PathAlgebraModuleError> {
                let range_dimension = range_top.dimension_vector()[vertex];
                if source_dimension == 0 || range_dimension == 0 {
                    return Ok(DenseMatrix::zero(
                        physical_dimension(source_dimension),
                        physical_dimension(range_dimension),
                        field.clone(),
                    ));
                }
                let mut entries = Vec::with_capacity(source_dimension * range_dimension);
                for basis_element in basis.iter().skip(offsets[vertex]).take(source_dimension) {
                    let lift = source_projection
                        .preimage_representative(basis_element)?
                        .ok_or(PathAlgebraModuleError::MissingTopBasisPreimage)?;
                    let image = homomorphism.image(&lift)?;
                    let top_image = range_projection.image(&image)?;
                    entries.extend_from_slice(
                        &top_image.external_representation()[vertex][..range_dimension],
                    );
                }
                Ok(DenseMatrix::from_flat(
                    source_dimension,
                    range_dimension,
                    entries,
                    field.clone(),
                )?)
            },
        )
        .collect::<Result<Vec<_>, _>>()?;
    right_module_hom_over_algebra(source_top.clone(), range_top.clone(), matrices)
}

/// Translate QPA's cached `MinimalGeneratingSetOfModule(M)` construction.
///
/// QPA takes the vertex-major basis of the canonical top and lifts every
/// basis element through `TopOfModuleProjection` using its chosen natural
/// quotient representatives.
///
/// # Errors
///
/// Propagates top construction and exact preimage solving. A missing lift is
/// reported as an internal source-translation failure.
#[allow(clippy::needless_pass_by_value)]
pub fn minimal_generating_set_of_module<F: Field>(
    module: PathAlgebraModule<F>,
) -> Result<Vec<PathAlgebraModuleElement<F>>, PathAlgebraModuleError> {
    if let Some(cached) = module.cached_minimal_generating_set()? {
        return Ok(cached);
    }
    let projection = top_of_module_projection(module.clone())?;
    let generators = projection
        .range()
        .basis_elements()
        .iter()
        .map(|basis_element| {
            projection
                .preimage_representative(basis_element)?
                .ok_or(PathAlgebraModuleError::MissingTopBasisPreimage)
        })
        .collect::<Result<Vec<_>, _>>()?;
    module.cache_minimal_generating_set(&generators)
}

/// Translate QPA's arrow-kernel `SocleOfModuleInclusion(M)` method.
///
/// At each vertex with outgoing arrows, the source method concatenates their
/// action rows and takes GAP's left nullspace. Sink spaces are retained in
/// full, and zero spaces preserve QPA's physical `1 x 1` zero matrix.
///
/// # Errors
///
/// Propagates exact matrix arithmetic and checked module/map construction.
pub fn socle_of_module_inclusion<F: Field>(
    module: PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    if module.dimension() == 0 {
        return zero_module_mapping(module.clone(), module);
    }

    let algebra = module.right_acting_algebra();
    let field = algebra.field().clone();
    let arrows = algebra.quiver().arrows().collect::<Vec<_>>();
    let mut socle_dimensions = Vec::with_capacity(module.dimension_vector().len());
    let mut inclusion_matrices = Vec::with_capacity(module.dimension_vector().len());

    for (vertex, &dimension) in module.dimension_vector().iter().enumerate() {
        let outgoing = arrows
            .iter()
            .filter(|arrow| arrow.source().index().zero_based() == vertex)
            .collect::<Vec<_>>();
        if dimension == 0 {
            socle_dimensions.push(0);
            inclusion_matrices.push(DenseMatrix::zero(1, 1, field.clone()));
            continue;
        }
        if outgoing.is_empty() {
            socle_dimensions.push(dimension);
            inclusion_matrices.push(DenseMatrix::identity(dimension, field.clone()));
            continue;
        }

        let concatenated_columns = outgoing
            .iter()
            .map(|arrow| module.dimension_vector()[arrow.target().index().zero_based()])
            .sum::<usize>();
        let mut entries = Vec::with_capacity(dimension * concatenated_columns);
        for row in 0..dimension {
            for arrow in &outgoing {
                let target = arrow.target().index().zero_based();
                let action = logical_submatrix(
                    &module.arrow_matrices()[arrow.index().zero_based()],
                    dimension,
                    module.dimension_vector()[target],
                )?;
                entries.extend_from_slice(action.row(row));
            }
        }
        let concatenated =
            DenseMatrix::from_flat(dimension, concatenated_columns, entries, field.clone())?;
        let kernel = concatenated.nullspace_rows();
        socle_dimensions.push(kernel.nrows());
        if kernel.nrows() == 0 {
            inclusion_matrices.push(DenseMatrix::zero(1, dimension, field.clone()));
        } else {
            inclusion_matrices.push(kernel);
        }
    }

    let socle = right_module_over_path_algebra(
        module.right_acting_algebra().clone(),
        socle_dimensions,
        &[],
    )?
    .ok_or(PathAlgebraModuleError::GeneratedSocleViolatesRelations)?;
    right_module_hom_over_algebra(socle, module, inclusion_matrices)
}

/// Return QPA's `SocleOfModule(M)` object as the source of its inclusion.
///
/// # Errors
///
/// Propagates [`socle_of_module_inclusion`].
pub fn socle_of_module<F: Field>(
    module: PathAlgebraModule<F>,
) -> Result<PathAlgebraModule<F>, PathAlgebraModuleError> {
    Ok(socle_of_module_inclusion(module)?.source().clone())
}

/// Translate the earlier, runtime-shadowed QPA `RightFacApproximation`
/// method from `modulehom.gi:1735-1750`.
///
/// QPA later installs a higher-ranked method on the synonymous
/// `TraceOfModule` operation; user-facing aliases therefore call
/// [`trace_of_module`]. This helper retains the earlier source algorithm for
/// the operation ledger and method-level differential tests.
///
/// # Errors
///
/// Rejects different acting algebras and propagates Hom, image, and
/// subrepresentation construction failures.
#[doc(hidden)]
pub fn right_fac_approximation_via_all_basis_images<F: Field>(
    module: &PathAlgebraModule<F>,
    target: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    let homomorphisms = hom_over_algebra(module, target)?
        .ok_or(PathAlgebraModuleError::DifferentHomomorphismAlgebras)?;
    let basis = module.basis_elements();
    let mut generators = Vec::new();
    for homomorphism in homomorphisms {
        generators.extend(homomorphism.images_set(&basis)?);
    }
    subrepresentation_inclusion(target.clone(), &generators)
}

/// Translate QPA's selected `TraceOfModule(M, N)` method.
///
/// This is also the operation behind `RightFacMApproximation`,
/// `RightFacApproximation`, and `MinimalRightFacMApproximation`. QPA orders
/// the generators first by a minimal top lift of `M`, then by its ordered Hom
/// basis, before closing their images inside `N`.
///
/// # Errors
///
/// Rejects modules over different algebras and propagates top, Hom, image,
/// zero-module, and subrepresentation construction failures.
pub fn trace_of_module<F: Field>(
    module: &PathAlgebraModule<F>,
    target: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    if module.right_acting_algebra() != target.right_acting_algebra() {
        return Err(PathAlgebraModuleError::DifferentHomomorphismAlgebras);
    }
    if module.dimension() == 0 || target.dimension() == 0 {
        return zero_submodule_inclusion(target);
    }
    let homomorphisms = hom_over_algebra(module, target)?
        .ok_or(PathAlgebraModuleError::DifferentHomomorphismAlgebras)?;
    if homomorphisms.is_empty() {
        return zero_submodule_inclusion(target);
    }
    let generators = minimal_generating_set_of_module(module.clone())?;
    let mut trace_generators = Vec::with_capacity(generators.len() * homomorphisms.len());
    for generator in &generators {
        for homomorphism in &homomorphisms {
            trace_generators.push(homomorphism.image(generator)?);
        }
    }
    subrepresentation_inclusion(target.clone(), &trace_generators)
}

/// QPA synonym of [`trace_of_module`].
///
/// # Errors
///
/// Propagates [`trace_of_module`].
pub fn right_fac_approximation<F: Field>(
    module: &PathAlgebraModule<F>,
    target: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    trace_of_module(module, target)
}

/// QPA synonym of [`trace_of_module`].
///
/// # Errors
///
/// Propagates [`trace_of_module`].
pub fn right_fac_m_approximation<F: Field>(
    module: &PathAlgebraModule<F>,
    target: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    trace_of_module(module, target)
}

/// QPA synonym of [`trace_of_module`].
///
/// # Errors
///
/// Propagates [`trace_of_module`].
pub fn minimal_right_fac_m_approximation<F: Field>(
    module: &PathAlgebraModule<F>,
    target: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    trace_of_module(module, target)
}

fn zero_submodule_inclusion<F: Field>(
    target: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    let vertex_count = target.dimension_vector().len();
    let zero = right_module_over_path_algebra(
        target.right_acting_algebra().clone(),
        vec![0; vertex_count],
        &[],
    )?
    .ok_or(PathAlgebraModuleError::GeneratedZeroModuleViolatesRelations)?;
    zero_module_mapping(zero, target.clone())
}

/// Translate QPA's `DirectSumOfQPAModules` block-matrix construction.
///
/// Summands and every vertex block retain their entered order. The returned
/// module carries QPA's direct-sum filter together with lazily recoverable
/// projection and inclusion attributes. `Ok(None)` is QPA's `fail` for an
/// empty list.
///
/// # Errors
///
/// Rejects different acting algebras, reports dimension overflow, and
/// propagates checked module/map construction failures.
pub fn direct_sum_of_qpa_modules<F: Field>(
    modules: &[PathAlgebraModule<F>],
) -> Result<Option<PathAlgebraModule<F>>, PathAlgebraModuleError> {
    let Some(first) = modules.first() else {
        return Ok(None);
    };
    if modules
        .iter()
        .any(|module| module.right_acting_algebra() != first.right_acting_algebra())
    {
        return Err(PathAlgebraModuleError::DirectSumDifferentAlgebras);
    }

    let vertex_count = first.dimension_vector().len();
    let mut dimensions = vec![0_usize; vertex_count];
    for module in modules {
        for (total, &dimension) in dimensions.iter_mut().zip(module.dimension_vector()) {
            *total = total
                .checked_add(dimension)
                .ok_or(PathAlgebraModuleError::DimensionOverflow)?;
        }
    }

    let field = first.right_acting_algebra().field().clone();
    let mut specification = Vec::new();
    for arrow in first.right_acting_algebra().quiver().arrows() {
        let source = arrow.source().index().zero_based();
        let target = arrow.target().index().zero_based();
        if dimensions[source] == 0 || dimensions[target] == 0 {
            continue;
        }
        let mut matrix = DenseMatrix::zero(dimensions[source], dimensions[target], field.clone());
        let mut row_offset = 0;
        let mut column_offset = 0;
        for module in modules {
            let source_dimension = module.dimension_vector()[source];
            let target_dimension = module.dimension_vector()[target];
            if source_dimension != 0 && target_dimension != 0 {
                let block = logical_submatrix(
                    &module.arrow_matrices()[arrow.index().zero_based()],
                    source_dimension,
                    target_dimension,
                )?;
                for row in 0..source_dimension {
                    for column in 0..target_dimension {
                        matrix[(row_offset + row, column_offset + column)] =
                            block[(row, column)].clone();
                    }
                }
            }
            row_offset += source_dimension;
            column_offset += target_dimension;
        }
        specification.push(NamedArrowMatrix::new(arrow.name(), matrix));
    }

    let direct_sum = right_module_over_path_algebra(
        first.right_acting_algebra().clone(),
        dimensions,
        &specification,
    )?
    .ok_or(PathAlgebraModuleError::GeneratedDirectSumViolatesRelations)?;
    let projections = construct_direct_sum_component_maps(&direct_sum, modules, true)?;
    let inclusions = construct_direct_sum_component_maps(&direct_sum, modules, false)?;
    direct_sum
        .0
        .direct_sum
        .set(DirectSumMetadata {
            summands: modules.to_vec(),
            projections: RwLock::new(
                projections
                    .iter()
                    .map(PathAlgebraModuleHomomorphism::downgrade)
                    .collect(),
            ),
            inclusions: RwLock::new(
                inclusions
                    .iter()
                    .map(PathAlgebraModuleHomomorphism::downgrade)
                    .collect(),
            ),
        })
        .map_err(|_| PathAlgebraModuleError::DirectSumMetadataAlreadySet)?;
    Ok(Some(direct_sum))
}

/// Return QPA's cached direct-sum projections in summand order.
///
/// # Errors
///
/// Rejects a module not built by [`direct_sum_of_qpa_modules`] and propagates
/// checked map reconstruction failures.
pub fn direct_sum_projections<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<Vec<PathAlgebraModuleHomomorphism<F>>, PathAlgebraModuleError> {
    direct_sum_component_maps(module, true)
}

/// Return QPA's cached direct-sum inclusions in summand order.
///
/// # Errors
///
/// Rejects a module not built by [`direct_sum_of_qpa_modules`] and propagates
/// checked map reconstruction failures.
pub fn direct_sum_inclusions<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<Vec<PathAlgebraModuleHomomorphism<F>>, PathAlgebraModuleError> {
    direct_sum_component_maps(module, false)
}

fn direct_sum_component_maps<F: Field>(
    module: &PathAlgebraModule<F>,
    projections: bool,
) -> Result<Vec<PathAlgebraModuleHomomorphism<F>>, PathAlgebraModuleError> {
    let metadata = module
        .0
        .direct_sum
        .get()
        .ok_or(PathAlgebraModuleError::NotDirectSumModule)?;
    let cache = if projections {
        &metadata.projections
    } else {
        &metadata.inclusions
    };
    let mut slot = cache
        .write()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    let mut result = Vec::with_capacity(metadata.summands.len());
    for summand_index in 0..metadata.summands.len() {
        if let Some(cached) = slot[summand_index].upgrade() {
            result.push(cached);
            continue;
        }
        let generated = construct_direct_sum_component_map(
            module,
            &metadata.summands,
            summand_index,
            projections,
        )?;
        slot[summand_index] = generated.downgrade();
        result.push(generated);
    }
    Ok(result)
}

fn construct_direct_sum_component_maps<F: Field>(
    direct_sum: &PathAlgebraModule<F>,
    summands: &[PathAlgebraModule<F>],
    projections: bool,
) -> Result<Vec<PathAlgebraModuleHomomorphism<F>>, PathAlgebraModuleError> {
    let mut result = Vec::with_capacity(summands.len());
    for summand_index in 0..summands.len() {
        result.push(construct_direct_sum_component_map(
            direct_sum,
            summands,
            summand_index,
            projections,
        )?);
    }
    Ok(result)
}

fn construct_direct_sum_component_map<F: Field>(
    direct_sum: &PathAlgebraModule<F>,
    summands: &[PathAlgebraModule<F>],
    summand_index: usize,
    projections: bool,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    let summand = &summands[summand_index];
    let field = direct_sum.right_acting_algebra().field().clone();
    let matrices = direct_sum
        .dimension_vector()
        .iter()
        .zip(summand.dimension_vector())
        .enumerate()
        .map(|(vertex, (&total, &component))| {
            let offset = summands[..summand_index]
                .iter()
                .map(|previous| previous.dimension_vector()[vertex])
                .sum::<usize>();
            let (rows, columns) = if projections {
                (physical_dimension(total), physical_dimension(component))
            } else {
                (physical_dimension(component), physical_dimension(total))
            };
            let mut matrix = DenseMatrix::zero(rows, columns, field.clone());
            if total != 0 && component != 0 {
                for coordinate in 0..component {
                    if projections {
                        matrix[(offset + coordinate, coordinate)] = field.one();
                    } else {
                        matrix[(coordinate, offset + coordinate)] = field.one();
                    }
                }
            }
            matrix
        })
        .collect::<Vec<_>>();
    if projections {
        let map = right_module_hom_over_algebra(direct_sum.clone(), summand.clone(), matrices)?;
        map.set_is_surjective(true);
        Ok(map)
    } else {
        let map = right_module_hom_over_algebra(summand.clone(), direct_sum.clone(), matrices)?;
        map.set_is_injective(true);
        Ok(map)
    }
}

/// Translate QPA's `PushOut(f, g)` construction.
///
/// The returned maps are ordered from `Range(g)` and `Range(f)` into the
/// pushout, matching QPA's `[f', g']`. `Ok(None)` is QPA's printed `fail`
/// result when the maps do not have structurally equal sources.
///
/// # Errors
///
/// Propagates direct-sum, map arithmetic, and cokernel construction failures.
pub fn push_out<F: Field>(
    first: &PathAlgebraModuleHomomorphism<F>,
    second: &PathAlgebraModuleHomomorphism<F>,
) -> Result<Option<PathAlgebraModuleHomomorphismPair<F>>, PathAlgebraModuleError> {
    if first.source() != second.source() {
        return Ok(None);
    }
    let direct_sum = direct_sum_of_qpa_modules(&[second.range().clone(), first.range().clone()])?
        .ok_or(PathAlgebraModuleError::InternalNonemptyDirectSumFailed)?;
    let inclusions = direct_sum_inclusions(&direct_sum)?;
    let difference = first
        .compose(&inclusions[1])?
        .add(&second.compose(&inclusions[0])?.negate()?)?;
    let cokernel = cokernel_projection(&difference)?;
    Ok(Some((
        inclusions[0].compose(&cokernel)?,
        inclusions[1].compose(&cokernel)?,
    )))
}

/// Translate QPA's `PullBack(f, g)` construction.
///
/// The returned maps are ordered into `Source(g)` and `Source(f)`, matching
/// QPA's `[f', g']`. `Ok(None)` is QPA's printed `fail` result when the maps
/// do not have structurally equal ranges.
///
/// # Errors
///
/// Propagates direct-sum, map arithmetic, and kernel construction failures.
pub fn pull_back<F: Field>(
    first: &PathAlgebraModuleHomomorphism<F>,
    second: &PathAlgebraModuleHomomorphism<F>,
) -> Result<Option<PathAlgebraModuleHomomorphismPair<F>>, PathAlgebraModuleError> {
    if first.range() != second.range() {
        return Ok(None);
    }
    let direct_sum = direct_sum_of_qpa_modules(&[first.source().clone(), second.source().clone()])?
        .ok_or(PathAlgebraModuleError::InternalNonemptyDirectSumFailed)?;
    let projections = direct_sum_projections(&direct_sum)?;
    let difference = projections[0]
        .compose(first)?
        .add(&projections[1].compose(second)?.negate()?)?;
    let kernel = kernel_inclusion(&difference)?;
    Ok(Some((
        kernel.compose(&projections[1])?,
        kernel.compose(&projections[0])?,
    )))
}

/// Translate QPA's `MorphismOnKernel(f, g, beta, alpha)`.
///
/// `Ok(None)` is QPA's `fail` for a noncommutative square.
///
/// # Errors
///
/// Rejects malformed diagram compositions and propagates exact kernel,
/// preimage, image, and checked-map construction failures.
pub fn morphism_on_kernel<F: Field>(
    first: &PathAlgebraModuleHomomorphism<F>,
    second: &PathAlgebraModuleHomomorphism<F>,
    beta: &PathAlgebraModuleHomomorphism<F>,
    alpha: &PathAlgebraModuleHomomorphism<F>,
) -> Result<Option<PathAlgebraModuleHomomorphism<F>>, PathAlgebraModuleError> {
    if first.compose(alpha)? != beta.compose(second)? {
        return Ok(None);
    }
    let first_kernel = kernel_inclusion(first)?;
    let second_kernel = kernel_inclusion(second)?;
    let route = first_kernel.compose(beta)?;
    let mut images = Vec::with_capacity(first_kernel.source().dimension());
    for basis_element in first_kernel.source().basis_elements() {
        let image = route.image(&basis_element)?;
        images.push(
            second_kernel
                .preimage_representative(&image)?
                .ok_or(PathAlgebraModuleError::MissingInducedKernelPreimage)?,
        );
    }
    Ok(Some(homomorphism_from_images(
        first_kernel.source(),
        second_kernel.source(),
        &images,
    )?))
}

/// Translate QPA's `MorphismOnImage(f, g, beta, alpha)`.
///
/// `Ok(None)` is QPA's `fail` for a noncommutative square.
///
/// # Errors
///
/// Rejects malformed diagram compositions and propagates exact image,
/// preimage, and checked-map construction failures.
pub fn morphism_on_image<F: Field>(
    first: &PathAlgebraModuleHomomorphism<F>,
    second: &PathAlgebraModuleHomomorphism<F>,
    beta: &PathAlgebraModuleHomomorphism<F>,
    alpha: &PathAlgebraModuleHomomorphism<F>,
) -> Result<Option<PathAlgebraModuleHomomorphism<F>>, PathAlgebraModuleError> {
    if first.compose(alpha)? != beta.compose(second)? {
        return Ok(None);
    }
    let first_image = image_projection(first)?;
    let second_image = image_projection(second)?;
    let route = beta.compose(&second_image)?;
    let mut images = Vec::with_capacity(first_image.range().dimension());
    for basis_element in first_image.range().basis_elements() {
        let lift = first_image
            .preimage_representative(&basis_element)?
            .ok_or(PathAlgebraModuleError::MissingInducedImagePreimage)?;
        images.push(route.image(&lift)?);
    }
    Ok(Some(homomorphism_from_images(
        first_image.range(),
        second_image.range(),
        &images,
    )?))
}

/// Translate QPA's `MorphismOnCoKernel(f, g, beta, alpha)`.
///
/// `Ok(None)` is QPA's `fail` for a noncommutative square.
///
/// # Errors
///
/// Rejects malformed diagram compositions and propagates exact cokernel,
/// preimage, image, and checked-map construction failures.
pub fn morphism_on_cokernel<F: Field>(
    first: &PathAlgebraModuleHomomorphism<F>,
    second: &PathAlgebraModuleHomomorphism<F>,
    beta: &PathAlgebraModuleHomomorphism<F>,
    alpha: &PathAlgebraModuleHomomorphism<F>,
) -> Result<Option<PathAlgebraModuleHomomorphism<F>>, PathAlgebraModuleError> {
    if first.compose(alpha)? != beta.compose(second)? {
        return Ok(None);
    }
    let first_cokernel = cokernel_projection(first)?;
    let second_cokernel = cokernel_projection(second)?;
    let route = alpha.compose(&second_cokernel)?;
    let mut images = Vec::with_capacity(first_cokernel.range().dimension());
    for basis_element in first_cokernel.range().basis_elements() {
        let lift = first_cokernel
            .preimage_representative(&basis_element)?
            .ok_or(PathAlgebraModuleError::MissingInducedCokernelPreimage)?;
        images.push(route.image(&lift)?);
    }
    Ok(Some(homomorphism_from_images(
        first_cokernel.range(),
        second_cokernel.range(),
        &images,
    )?))
}

/// Translate QPA's `LiftingInclusionMorphisms(f, g)`.
///
/// QPA first constructs the cokernel of `f`, then returns `fail` unless both
/// maps are injective with a common structural range and the image of `g`
/// lies in the image of `f`. `Ok(None)` represents that `fail` result.
///
/// # Errors
///
/// Propagates cokernel, image, preimage, and checked-map construction
/// failures. A missing preimage after the containment test is reported as an
/// internal source-translation failure.
pub fn lifting_inclusion_morphisms<F: Field>(
    first: &PathAlgebraModuleHomomorphism<F>,
    second: &PathAlgebraModuleHomomorphism<F>,
) -> Result<Option<PathAlgebraModuleHomomorphism<F>>, PathAlgebraModuleError> {
    let cokernel = cokernel_projection(first)?;
    if !first.is_injective()
        || !second.is_injective()
        || first.range() != second.range()
        || !second.compose(&cokernel)?.is_zero()
    {
        return Ok(None);
    }
    let mut images = Vec::with_capacity(second.source().dimension());
    for basis_element in second.source().basis_elements() {
        let image = second.image(&basis_element)?;
        images.push(
            first
                .preimage_representative(&image)?
                .ok_or(PathAlgebraModuleError::MissingInclusionLiftPreimage)?,
        );
    }
    Ok(Some(homomorphism_from_images(
        second.source(),
        first.source(),
        &images,
    )?))
}

/// Translate QPA's `HomomorphismFromImages(M, N, genImages)` operation.
///
/// Images correspond to QPA's standard basis in vertex-major order. A
/// nonzero image must be supported at the same vertex as its basis vector;
/// additional images are ignored exactly as in the source loop.
///
/// # Errors
///
/// Reports a short image list or wrong support and propagates exact matrix and
/// checked module-map construction failures.
pub fn homomorphism_from_images<F: Field>(
    source: &PathAlgebraModule<F>,
    range: &PathAlgebraModule<F>,
    images: &[PathAlgebraModuleElement<F>],
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    if images.len() < source.dimension() {
        return Err(PathAlgebraModuleError::BasisImageCount {
            actual: images.len(),
            expected: source.dimension(),
        });
    }
    let mut image_position = 0;
    for (vertex, &dimension) in source.dimension_vector().iter().enumerate() {
        for _ in 0..dimension {
            let image = &images[image_position];
            let support = image.support_vertices();
            if !support.is_empty() && support != [vertex + 1] {
                return Err(PathAlgebraModuleError::WrongBasisImageSupport {
                    position: image_position + 1,
                });
            }
            image_position += 1;
        }
    }
    let field = source.right_acting_algebra().field().clone();
    let mut offset = 0_usize;
    let matrices = source
        .dimension_vector()
        .iter()
        .zip(range.dimension_vector())
        .enumerate()
        .map(|(vertex, (&source_dimension, &range_dimension))| {
            if source_dimension == 0 {
                return Ok(DenseMatrix::zero(
                    1,
                    physical_dimension(range_dimension),
                    field.clone(),
                ));
            }
            let mut entries =
                Vec::with_capacity(source_dimension * physical_dimension(range_dimension));
            for image in &images[offset..offset + source_dimension] {
                entries.extend_from_slice(&image.external_representation()[vertex]);
            }
            offset += source_dimension;
            Ok(DenseMatrix::from_flat(
                source_dimension,
                physical_dimension(range_dimension),
                entries,
                field.clone(),
            )?)
        })
        .collect::<Result<Vec<_>, PathAlgebraModuleError>>()?;
    right_module_hom_over_algebra(source.clone(), range.clone(), matrices)
}

/// Translate QPA's `MultiplyListsOfMaps(projections, matrix, inclusions)`.
///
/// The computation is the literal row-by-matrix-by-column product in QPA's
/// left-to-right map-composition convention. Extra matrix rows are ignored by
/// the source loop.
///
/// # Errors
///
/// The intended QPA domain requires nonempty outer lists, at least one matrix
/// row per inclusion, and one entry per projection in each consumed row.
/// Composition, addition, and checked map failures propagate.
pub fn multiply_lists_of_maps<F: Field>(
    projections: &[PathAlgebraModuleHomomorphism<F>],
    matrix: &[Vec<PathAlgebraModuleHomomorphism<F>>],
    inclusions: &[PathAlgebraModuleHomomorphism<F>],
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    if projections.is_empty() {
        return Err(PathAlgebraModuleError::EmptyMapProjectionList);
    }
    if inclusions.is_empty() {
        return Err(PathAlgebraModuleError::EmptyMapInclusionList);
    }
    if matrix.len() < inclusions.len() {
        return Err(PathAlgebraModuleError::MapMatrixRowCount {
            actual: matrix.len(),
            expected: inclusions.len(),
        });
    }
    let mut outer_terms = Vec::with_capacity(inclusions.len());
    for (row_position, (row, inclusion)) in matrix.iter().zip(inclusions).enumerate() {
        if row.len() != projections.len() {
            return Err(PathAlgebraModuleError::MapMatrixColumnCount {
                row: row_position + 1,
                actual: row.len(),
                expected: projections.len(),
            });
        }
        let mut row_terms = projections
            .iter()
            .zip(row)
            .map(|(projection, entry)| projection.compose(entry));
        let Some(first_term) = row_terms.next() else {
            return Err(PathAlgebraModuleError::EmptyMapProjectionList);
        };
        let mut row_sum = first_term?;
        for term in row_terms {
            row_sum = row_sum.add(&term?)?;
        }
        outer_terms.push(row_sum.compose(inclusion)?);
    }
    let mut outer_terms = outer_terms.into_iter();
    let Some(mut sum) = outer_terms.next() else {
        return Err(PathAlgebraModuleError::EmptyMapInclusionList);
    };
    for term in outer_terms {
        sum = sum.add(&term)?;
    }
    Ok(sum)
}

/// Translate QPA's `UnderlyingLinearMap(f)` block-diagonal construction.
///
/// Vertex spaces occur in quiver order and the stored row-action matrix at
/// each vertex is copied into the corresponding diagonal block. QPA retains
/// one physical zero row or column for every zero-dimensional vertex space;
/// consequently the returned dimensions are the sums of the physical vertex
/// dimensions rather than the logical module dimensions.
///
/// # Errors
///
/// Reports address-size overflow while summing physical dimensions.
pub fn underlying_linear_map<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<DenseMatrix<F>, PathAlgebraModuleError> {
    let source_dimension =
        homomorphism
            .vertex_matrices()
            .iter()
            .try_fold(0_usize, |total, matrix| {
                total
                    .checked_add(matrix.nrows())
                    .ok_or(PathAlgebraModuleError::DimensionOverflow)
            })?;
    let range_dimension =
        homomorphism
            .vertex_matrices()
            .iter()
            .try_fold(0_usize, |total, matrix| {
                total
                    .checked_add(matrix.ncols())
                    .ok_or(PathAlgebraModuleError::DimensionOverflow)
            })?;
    let field = homomorphism.source().right_acting_algebra().field().clone();
    let mut result = DenseMatrix::zero(source_dimension, range_dimension, field);
    let mut source_offset = 0;
    let mut range_offset = 0;
    for matrix in homomorphism.vertex_matrices() {
        for row in 0..matrix.nrows() {
            for column in 0..matrix.ncols() {
                result[(source_offset + row, range_offset + column)] =
                    matrix[(row, column)].clone();
            }
        }
        source_offset += matrix.nrows();
        range_offset += matrix.ncols();
    }
    Ok(result)
}

/// Translate QPA's recursive `IsomorphismOfModules(M, N)` operation.
///
/// `Ok(None)` is QPA's `false`. After the same algebra, dimension-vector,
/// zero, and structural-equality branches, this follows QPA's ordered Hom-basis
/// search for a common Fitting summand and recursively identifies the two
/// complements.
///
/// # Errors
///
/// Propagates exact Hom, composition, kernel/image, split-inverse, preimage,
/// basis-image, and checked map construction failures. Failures of the
/// theoretically split internal inclusions are reported explicitly.
#[allow(clippy::too_many_lines)]
pub fn isomorphism_of_modules<F: Field>(
    source: &PathAlgebraModule<F>,
    range: &PathAlgebraModule<F>,
) -> Result<Option<PathAlgebraModuleHomomorphism<F>>, PathAlgebraModuleError> {
    if source.right_acting_algebra() != range.right_acting_algebra()
        || source.dimension_vector() != range.dimension_vector()
    {
        return Ok(None);
    }
    if source.dimension() == 0 && range.dimension() == 0 {
        return Ok(Some(zero_module_mapping(source.clone(), range.clone())?));
    }
    if source == range {
        let field = source.right_acting_algebra().field().clone();
        let matrices = source
            .dimension_vector()
            .iter()
            .map(|&dimension| {
                if dimension == 0 {
                    DenseMatrix::zero(1, 1, field.clone())
                } else {
                    DenseMatrix::identity(dimension, field.clone())
                }
            })
            .collect();
        return Ok(Some(right_module_hom_over_algebra(
            source.clone(),
            range.clone(),
            matrices,
        )?));
    }

    let Some((source_endomorphism, range_endomorphism, source_to_range)) =
        common_splitting_maps(source, range)?
    else {
        return Ok(None);
    };
    let source_summand_inclusion = image_inclusion(&source_endomorphism)?;
    let source_summand_retraction = right_inverse_of_homomorphism(&source_summand_inclusion)?
        .ok_or(PathAlgebraModuleError::InternalIsomorphismSummandNotSplit)?;
    let range_summand_inclusion = image_inclusion(&range_endomorphism)?;
    let range_summand_retraction = right_inverse_of_homomorphism(&range_summand_inclusion)?
        .ok_or(PathAlgebraModuleError::InternalIsomorphismSummandNotSplit)?;

    let source_idempotent = source_summand_retraction.compose(&source_summand_inclusion)?;
    let source_common_inclusion = image_inclusion(&source_idempotent)?;
    let source_complement_inclusion = kernel_inclusion(&source_idempotent)?;
    let range_idempotent = range_summand_retraction.compose(&range_summand_inclusion)?;
    let range_common_inclusion = image_inclusion(&range_idempotent)?;
    let range_complement_inclusion = kernel_inclusion(&range_idempotent)?;

    if source_complement_inclusion.source().dimension() == 0
        && range_complement_inclusion.source().dimension() == 0
    {
        return Ok(Some(source_endomorphism.compose(&source_to_range)?));
    }

    let common_route = source_common_inclusion
        .compose(&source_endomorphism)?
        .compose(&source_to_range)?;
    let mut common_images = Vec::with_capacity(common_route.source().dimension());
    for basis_element in common_route.source().basis_elements() {
        let image = common_route.image(&basis_element)?;
        common_images.push(
            range_common_inclusion
                .preimage_representative(&image)?
                .ok_or(PathAlgebraModuleError::MissingIsomorphismSummandPreimage)?,
        );
    }
    let common_isomorphism = homomorphism_from_images(
        source_common_inclusion.source(),
        range_common_inclusion.source(),
        &common_images,
    )?;
    let Some(complement_isomorphism) = isomorphism_of_modules(
        source_complement_inclusion.source(),
        range_complement_inclusion.source(),
    )?
    else {
        return Ok(None);
    };

    let source_common_retraction = right_inverse_of_homomorphism(&source_common_inclusion)?
        .ok_or(PathAlgebraModuleError::InternalIsomorphismSummandNotSplit)?;
    let source_complement_retraction = right_inverse_of_homomorphism(&source_complement_inclusion)?
        .ok_or(PathAlgebraModuleError::InternalIsomorphismSummandNotSplit)?;
    let common_map = source_common_retraction
        .compose(&common_isomorphism)?
        .compose(&range_common_inclusion)?;
    let complement_map = source_complement_retraction
        .compose(&complement_isomorphism)?
        .compose(&range_complement_inclusion)?;
    Ok(Some(common_map.add(&complement_map)?))
}

fn common_splitting_maps<F: Field>(
    source: &PathAlgebraModule<F>,
    range: &PathAlgebraModule<F>,
) -> Result<Option<PathAlgebraModuleHomomorphismTriple<F>>, PathAlgebraModuleError> {
    let source_to_range = hom_over_algebra(source, range)?
        .ok_or(PathAlgebraModuleError::InternalIsomorphismHomRejected)?;
    let range_to_source = hom_over_algebra(range, source)?
        .ok_or(PathAlgebraModuleError::InternalIsomorphismHomRejected)?;
    if source_to_range.is_empty() || range_to_source.is_empty() {
        return Ok(None);
    }
    let source_maximum = source.dimension_vector().iter().copied().max().unwrap_or(0);
    let range_maximum = range.dimension_vector().iter().copied().max().unwrap_or(0);
    let fitting_length = if source_maximum == range_maximum {
        range_maximum
    } else {
        source_maximum.min(range_maximum) + 1
    };

    for range_map in &range_to_source {
        for source_map in &source_to_range {
            let fitting_map = if fitting_length > 1 {
                homomorphism_power(&source_map.compose(range_map)?, fitting_length - 1)?
                    .compose(source_map)?
            } else {
                source_map.clone()
            };
            let source_endomorphism = fitting_map.compose(range_map)?;
            if !source_endomorphism.is_zero() {
                return Ok(Some((
                    source_endomorphism,
                    range_map.compose(&fitting_map)?,
                    source_map.clone(),
                )));
            }
        }
    }
    Ok(None)
}

/// Translate QPA's binary `IntersectionOfSubmodules(f, g)`.
///
/// The tuple contains the maps into the second and first source followed by
/// the canonical inclusion of their intersection into the common range.
///
/// # Errors
///
/// Requires two monomorphisms with structurally equal ranges and propagates
/// QPA's pullback construction failures.
pub fn intersection_of_submodules<F: Field>(
    first: &PathAlgebraModuleHomomorphism<F>,
    second: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphismTriple<F>, PathAlgebraModuleError> {
    validate_submodule_pair(first, second)?;
    let (to_second, to_first) = pull_back(first, second)?
        .ok_or(PathAlgebraModuleError::InternalCommonRangePullbackFailed)?;
    let inclusion = to_first.compose(first)?;
    Ok((to_second, to_first, inclusion))
}

/// Translate QPA's list `IntersectionOfSubmodules` method.
///
/// A singleton list returns the identical entered map, while longer lists
/// fold binary pullbacks in entered order.
///
/// # Errors
///
/// Preserves QPA's validation order for an empty list, source algebras,
/// injectivity, and common ranges, then propagates binary construction errors.
pub fn intersection_of_submodules_list<F: Field>(
    homomorphisms: &[PathAlgebraModuleHomomorphism<F>],
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    validate_submodule_list(homomorphisms, true)?;
    let mut intersection = homomorphisms[0].clone();
    for homomorphism in &homomorphisms[1..] {
        intersection = intersection_of_submodules(&intersection, homomorphism)?.2;
    }
    Ok(intersection)
}

/// Translate QPA's binary `SumOfSubmodules(f, g)`.
///
/// The tuple is the inclusion of the sum into the common range followed by
/// the maps from the first and second source into that sum.
///
/// # Errors
///
/// Requires two monomorphisms with structurally equal ranges and propagates
/// QPA's direct-sum and image-factorization construction failures.
pub fn sum_of_submodules<F: Field>(
    first: &PathAlgebraModuleHomomorphism<F>,
    second: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphismTriple<F>, PathAlgebraModuleError> {
    validate_submodule_pair(first, second)?;
    let direct_sum = direct_sum_of_qpa_modules(&[first.source().clone(), second.source().clone()])?
        .ok_or(PathAlgebraModuleError::InternalNonemptyDirectSumFailed)?;
    let inclusions = direct_sum_inclusions(&direct_sum)?;
    let projections = direct_sum_projections(&direct_sum)?;
    let combined = projections[0]
        .compose(first)?
        .add(&projections[1].compose(second)?)?;
    let (projection, inclusion) = image_projection_inclusion(&combined)?;
    Ok((
        inclusion,
        inclusions[0].compose(&projection)?,
        inclusions[1].compose(&projection)?,
    ))
}

/// Translate QPA's list `SumOfSubmodules` method.
///
/// A singleton list returns the identical entered map, while longer lists
/// fold binary image sums in entered order.
///
/// # Errors
///
/// Preserves QPA's validation order for an empty list, source algebras,
/// common ranges, and injectivity, then propagates binary construction errors.
pub fn sum_of_submodules_list<F: Field>(
    homomorphisms: &[PathAlgebraModuleHomomorphism<F>],
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    validate_submodule_list(homomorphisms, false)?;
    let mut sum = homomorphisms[0].clone();
    for homomorphism in &homomorphisms[1..] {
        sum = sum_of_submodules(&sum, homomorphism)?.0;
    }
    Ok(sum)
}

fn validate_submodule_pair<F: Field>(
    first: &PathAlgebraModuleHomomorphism<F>,
    second: &PathAlgebraModuleHomomorphism<F>,
) -> Result<(), PathAlgebraModuleError> {
    if !first.is_injective() {
        return Err(PathAlgebraModuleError::FirstSubmoduleMapNotInjective);
    }
    if !second.is_injective() {
        return Err(PathAlgebraModuleError::SecondSubmoduleMapNotInjective);
    }
    if first.range() != second.range() {
        return Err(PathAlgebraModuleError::SubmoduleRangesDiffer);
    }
    Ok(())
}

fn validate_submodule_list<F: Field>(
    homomorphisms: &[PathAlgebraModuleHomomorphism<F>],
    injectivity_before_range: bool,
) -> Result<(), PathAlgebraModuleError> {
    let Some(first) = homomorphisms.first() else {
        return Err(PathAlgebraModuleError::EmptySubmoduleList);
    };
    if homomorphisms.iter().any(|homomorphism| {
        homomorphism.source().right_acting_algebra() != first.source().right_acting_algebra()
    }) {
        return Err(PathAlgebraModuleError::SubmoduleListDifferentAlgebras);
    }
    let injective = || {
        homomorphisms
            .iter()
            .all(PathAlgebraModuleHomomorphism::is_injective)
    };
    let common_range = || {
        homomorphisms
            .iter()
            .all(|homomorphism| homomorphism.range() == first.range())
    };
    if injectivity_before_range {
        if !injective() {
            return Err(PathAlgebraModuleError::SomeSubmoduleMapNotInjective);
        }
        if !common_range() {
            return Err(PathAlgebraModuleError::SubmoduleRangesDiffer);
        }
    } else {
        if !common_range() {
            return Err(PathAlgebraModuleError::SubmoduleRangesDiffer);
        }
        if !injective() {
            return Err(PathAlgebraModuleError::SomeSubmoduleMapNotInjective);
        }
    }
    Ok(())
}

/// Translate QPA's `RejectOfModule(N, M)` intersection-of-kernels method.
///
/// # Errors
///
/// Requires both modules over the same algebra and propagates exact Hom,
/// kernel, and submodule-intersection failures.
pub fn reject_of_module<F: Field>(
    module: &PathAlgebraModule<F>,
    rejected_by: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    if module.right_acting_algebra() != rejected_by.right_acting_algebra() {
        return Err(PathAlgebraModuleError::RejectModulesOverDifferentAlgebras);
    }
    let homomorphisms = hom_over_algebra(module, rejected_by)?
        .ok_or(PathAlgebraModuleError::InternalCommonAlgebraHomFailed)?;
    match homomorphisms.as_slice() {
        [] => identity_module_mapping(module.clone()),
        [homomorphism] => kernel_inclusion(homomorphism),
        _ => {
            let kernels = homomorphisms
                .iter()
                .map(kernel_inclusion)
                .collect::<Result<Vec<_>, _>>()?;
            intersection_of_submodules_list(&kernels)
        }
    }
}

#[derive(Clone, Debug)]
struct NaturalRowQuotient<F: Field> {
    projection: DenseMatrix<F>,
    lifts: DenseMatrix<F>,
}

fn natural_quotient_by_row_space<F: Field>(
    image_rows: &DenseMatrix<F>,
) -> Result<NaturalRowQuotient<F>, MatrixError> {
    let field = image_rows.field().clone();
    let dimension = image_rows.ncols();
    let kernel_basis = image_rows.semi_echelon_transformation().vectors;
    let mut mutable_basis = kernel_basis.clone();
    let mut complement = Vec::new();
    for coordinate in 0..dimension {
        let mut standard = vec![field.zero(); dimension];
        standard[coordinate] = field.one();
        let mut candidate = mutable_basis.as_flat().to_vec();
        candidate.extend_from_slice(&standard);
        let candidate = DenseMatrix::from_flat(
            mutable_basis.nrows() + 1,
            dimension,
            candidate,
            field.clone(),
        )?;
        let semi = candidate.semi_echelon_transformation().vectors;
        if semi.nrows() > mutable_basis.nrows() {
            complement.push(standard);
            mutable_basis = semi;
        }
    }

    let quotient_dimension = complement.len();
    let lifts = DenseMatrix::from_flat(
        quotient_dimension,
        dimension,
        complement.iter().flatten().cloned().collect(),
        field.clone(),
    )?;
    let mut basis_entries = kernel_basis.as_flat().to_vec();
    basis_entries.extend(complement.iter().flatten().cloned());
    let full_basis = DenseMatrix::from_flat(dimension, dimension, basis_entries, field.clone())?;
    let mut projection_entries = Vec::with_capacity(dimension * quotient_dimension);
    for coordinate in 0..dimension {
        let mut standard = vec![field.zero(); dimension];
        standard[coordinate] = field.one();
        let coefficients = full_basis
            .solution_row(&standard)?
            .expect("a kernel basis extended by independent vectors spans the ambient space");
        projection_entries
            .extend_from_slice(&coefficients[dimension - quotient_dimension..dimension]);
    }
    Ok(NaturalRowQuotient {
        projection: DenseMatrix::from_flat(
            dimension,
            quotient_dimension,
            projection_entries,
            field,
        )?,
        lifts,
    })
}

fn logical_submatrix<F: Field>(
    matrix: &DenseMatrix<F>,
    rows: usize,
    columns: usize,
) -> Result<DenseMatrix<F>, MatrixError> {
    let mut entries = Vec::with_capacity(rows * columns);
    for row in 0..rows {
        entries.extend_from_slice(&matrix.row(row)[..columns]);
    }
    DenseMatrix::from_flat(rows, columns, entries, matrix.field().clone())
}

fn logical_vertex_map<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
    vertex: usize,
) -> DenseMatrix<F> {
    logical_submatrix(
        &homomorphism.vertex_matrices()[vertex],
        homomorphism.source().dimension_vector()[vertex],
        homomorphism.range().dimension_vector()[vertex],
    )
    .expect("a checked module homomorphism retains every logical vertex submatrix")
}

fn dimension_offsets(dimensions: &[usize]) -> Vec<usize> {
    let mut offset = 0;
    dimensions
        .iter()
        .map(|&dimension| {
            let result = offset;
            offset += dimension;
            result
        })
        .collect()
}

fn quiver_algebra_path_element<F: Field>(
    algebra: &QuiverAlgebra<F>,
    path: Path,
) -> Result<QuiverAlgebraElement<F>, PathAlgebraModuleError> {
    match algebra {
        QuiverAlgebra::Path(algebra) => Ok(QuiverAlgebraElement::Path(algebra.element(path)?)),
        QuiverAlgebra::Quotient(algebra) => Ok(QuiverAlgebraElement::Quotient(algebra.path(path)?)),
    }
}

/// Translate QPA's checked three-argument
/// `RightModuleOverPathAlgebra(A, dimension_vector, generators)`.
///
/// `Ok(None)` is QPA's `fail` result when the entered actions do not satisfy
/// the quotient relations. Malformed input is returned as an error, matching
/// the source method's GAP errors.
///
/// # Errors
///
/// Rejects a malformed dimension vector, matrices over another field, unknown
/// labels, nonzero maps incident with a zero space, wrong shapes, and matrix
/// arithmetic failures during relation evaluation.
pub fn right_module_over_path_algebra<F: Field>(
    algebra: QuiverAlgebra<F>,
    dimension_vector: Vec<usize>,
    generators: &[NamedArrowMatrix<F>],
) -> Result<Option<PathAlgebraModule<F>>, PathAlgebraModuleError> {
    right_module_over_path_algebra_impl(algebra, dimension_vector, generators, None)
}

/// Construct a checked module while retaining the tensor object underlying
/// its acting algebra. This is the cross-crate hook used by QPA's enveloping
/// algebra attributes.
#[doc(hidden)]
pub fn right_module_over_path_algebra_with_tensor<F: Field>(
    algebra: QuiverAlgebra<F>,
    dimension_vector: Vec<usize>,
    generators: &[NamedArrowMatrix<F>],
    tensor: TensorProductAlgebra<F>,
) -> Result<Option<PathAlgebraModule<F>>, PathAlgebraModuleError> {
    if tensor.algebra() != &algebra {
        return Err(PathAlgebraModuleError::RetainedTensorAlgebraMismatch);
    }
    right_module_over_path_algebra_impl(algebra, dimension_vector, generators, Some(tensor))
}

fn right_module_over_path_algebra_impl<F: Field>(
    algebra: QuiverAlgebra<F>,
    dimension_vector: Vec<usize>,
    generators: &[NamedArrowMatrix<F>],
    retained_tensor: Option<TensorProductAlgebra<F>>,
) -> Result<Option<PathAlgebraModule<F>>, PathAlgebraModuleError> {
    let quiver = algebra.quiver();
    if dimension_vector.len() != quiver.number_of_vertices() {
        return Err(PathAlgebraModuleError::DimensionVectorLength {
            actual: dimension_vector.len(),
            expected: quiver.number_of_vertices(),
        });
    }
    let dimension = dimension_vector.iter().try_fold(0_usize, |sum, &entry| {
        sum.checked_add(entry)
            .ok_or(PathAlgebraModuleError::DimensionOverflow)
    })?;

    for generator in generators {
        if generator.matrix().field() != algebra.field() {
            return Err(PathAlgebraModuleError::DifferentCoefficientField {
                label: generator.label().into(),
            });
        }
    }
    let missing_labels = generators
        .iter()
        .filter(|generator| {
            !quiver
                .arrows()
                .any(|arrow| arrow.name() == generator.label())
        })
        .map(|generator| Box::<str>::from(generator.label()))
        .collect::<Vec<_>>();
    if !missing_labels.is_empty() {
        return Err(PathAlgebraModuleError::UnknownArrowLabels {
            labels: missing_labels,
        });
    }

    let field = algebra.field().clone();
    let mut arrow_matrices = quiver
        .arrows()
        .map(|arrow| {
            let source = dimension_vector[arrow.source().index().zero_based()];
            let target = dimension_vector[arrow.target().index().zero_based()];
            DenseMatrix::zero(
                physical_dimension(source),
                physical_dimension(target),
                field.clone(),
            )
        })
        .collect::<Vec<_>>();

    for generator in generators {
        let Generator::Arrow(arrow) = quiver
            .generator_by_name(generator.label())
            .expect("membership in the arrow-label list was checked above")
        else {
            unreachable!("an arrow name takes precedence in QPA's quiver record")
        };
        let source = dimension_vector[arrow.source().index().zero_based()];
        let target = dimension_vector[arrow.target().index().zero_based()];
        if source == 0 || target == 0 {
            if !generator.matrix().is_zero() {
                return Err(PathAlgebraModuleError::NonzeroMapAtZeroSpace {
                    label: generator.label().into(),
                });
            }
            return Err(PathAlgebraModuleError::ZeroSpaceMapMustBeOmitted {
                label: generator.label().into(),
            });
        }
        let actual = (generator.matrix().nrows(), generator.matrix().ncols());
        let expected = (source, target);
        if actual != expected {
            return Err(PathAlgebraModuleError::WrongArrowMatrixShape {
                label: generator.label().into(),
                actual,
                expected,
            });
        }
        arrow_matrices[arrow.index().zero_based()] = generator.matrix().clone();
    }

    if dimension != 0
        && algebra.as_quotient().is_some()
        && !relations_are_satisfied(&algebra, &dimension_vector, &arrow_matrices, &field)?
    {
        return Ok(None);
    }

    Ok(Some(PathAlgebraModule(Arc::new(
        PathAlgebraModuleData::new(
            algebra,
            dimension_vector,
            dimension,
            arrow_matrices,
            retained_tensor,
        ),
    ))))
}

fn physical_dimension(dimension: usize) -> usize {
    dimension.max(1)
}

fn relations_are_satisfied<F: Field>(
    algebra: &QuiverAlgebra<F>,
    dimension_vector: &[usize],
    arrow_matrices: &[DenseMatrix<F>],
    field: &F,
) -> Result<bool, MatrixError> {
    for relation in algebra.relations() {
        let mut value: Option<DenseMatrix<F>> = None;
        for term in relation.terms() {
            let term_value =
                path_action_matrix(term.path(), dimension_vector, arrow_matrices, field.clone())?
                    .scale_left(term.coefficient());
            value = Some(match value {
                Some(sum) => sum.add(&term_value)?,
                None => term_value,
            });
        }
        if value.is_some_and(|matrix| !matrix.is_zero()) {
            return Ok(false);
        }
    }
    Ok(true)
}

fn path_action_matrix<F: Field>(
    path: &Path,
    dimension_vector: &[usize],
    arrow_matrices: &[DenseMatrix<F>],
    field: F,
) -> Result<DenseMatrix<F>, MatrixError> {
    let Some(walk) = path.walk() else {
        unreachable!("normalized path-algebra relations contain no zero-path term")
    };
    let mut walk = walk.iter();
    let Some(first) = walk.next() else {
        let PathEndpoint::Vertex(vertex) = path.source() else {
            unreachable!("a vertex path has a vertex source")
        };
        return Ok(DenseMatrix::identity(
            physical_dimension(dimension_vector[vertex.index().zero_based()]),
            field,
        ));
    };
    let mut result = arrow_matrices[first.zero_based()].clone();
    for arrow in walk {
        result = result.multiply(&arrow_matrices[arrow.zero_based()])?;
    }
    Ok(result)
}

/// Failure in QPA's checked right-module constructor.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum PathAlgebraModuleError {
    /// The dimension vector must have one entry per quiver vertex.
    #[error("dimension vector has length {actual}, expected {expected}")]
    DimensionVectorLength {
        /// Entered length.
        actual: usize,
        /// Required vertex count.
        expected: usize,
    },
    /// Total module dimension exceeded Rust's addressable size.
    #[error("the total module dimension overflows usize")]
    DimensionOverflow,
    /// QPA's Hom-space equation matrix exceeded Rust's addressable size.
    #[error("the HomOverAlgebra equation matrix dimensions overflow usize")]
    HomEquationSizeOverflow,
    /// A supported vertex must have a row block in QPA's Hom unknown vector.
    #[error("supported vertex {vertex} has no HomOverAlgebra unknown block")]
    MissingHomUnknownBlock {
        /// One-based vertex position.
        vertex: usize,
    },
    /// An endomorphism-space query cannot fail its identical algebra check.
    #[error("HomOverAlgebra unexpectedly rejected an endomorphism pair")]
    EndomorphismAlgebraMismatch,
    /// A reverse Hom query for one module map necessarily has a common
    /// acting algebra.
    #[error("HomOverAlgebra unexpectedly rejected a splitting pair")]
    InternalSplittingHomRejected,
    /// Kernel/cokernel-derived modules necessarily share the entered map's
    /// acting algebra during a minimal-version computation.
    #[error("HomOverAlgebra unexpectedly rejected a minimal-version pair")]
    InternalMinimalHomRejected,
    /// QPA's Fitting decomposition requires this canonical kernel inclusion
    /// to have the retraction used by `MoreLeftMinimalVersion`.
    #[error("the internal Fitting kernel inclusion has no right inverse")]
    InternalMinimalKernelNotSplit,
    /// `InverseOfIsomorphism` rejects a map before inspecting its matrices.
    #[error("entered map is not an isomorphism")]
    EnteredMapNotIsomorphism,
    /// Every nonzero component matrix of a module isomorphism is invertible.
    #[error("an isomorphism has a singular nonzero vertex matrix")]
    InternalIsomorphismMatrixSingular,
    /// The two modules in QPA's recursive isomorphism search have already
    /// passed their common-algebra check.
    #[error("HomOverAlgebra unexpectedly rejected an isomorphism-search pair")]
    InternalIsomorphismHomRejected,
    /// QPA's Fitting summands and complements must have canonical retractions.
    #[error("an internal isomorphism-search summand inclusion has no right inverse")]
    InternalIsomorphismSummandNotSplit,
    /// QPA's common-summand route must land in the chosen range summand.
    #[error("a common-summand image has no preimage in the range summand")]
    MissingIsomorphismSummandPreimage,
    /// GAP's splitting solver returns one coefficient per reverse-Hom basis
    /// vector.
    #[error("splitting combination has {actual} coefficients, expected {expected}")]
    SplittingCoefficientCount {
        /// Entered coefficient count.
        actual: usize,
        /// Reverse-Hom basis length.
        expected: usize,
    },
    /// Internal linear-combination construction follows a nonempty-basis
    /// source branch.
    #[error("an internal splitting combination has an empty Hom basis")]
    EmptySplittingHomBasis,
    /// GAP's `DirectSumOfAlgebras` rejects the empty list for the zero module.
    #[error("EndOverAlgebra cannot form a direct sum for a zero module")]
    EmptyEndomorphismDirectSum,
    /// `FromHomMMToEndM` is defined for an endomorphism of one module.
    #[error("FromHomMMToEndM requires structurally equal source and range modules")]
    NonEndomorphismConversion,
    /// `FromEndMToHomMM` requires the module's coefficient field.
    #[error("the endomorphism matrix is over a different coefficient field")]
    EndomorphismMatrixField,
    /// QPA's diagonal-block extraction indexes the full underlying module.
    #[error("endomorphism matrix has shape {actual:?}, but must be at least {minimum:?}")]
    EndomorphismMatrixTooSmall {
        /// Entered matrix shape.
        actual: (usize, usize),
        /// Minimum row and column counts.
        minimum: (usize, usize),
    },
    /// A natural quotient accepts only matrices represented in its source.
    #[error("the matrix is not an element of the represented endomorphism algebra")]
    MatrixOutsideEndomorphismAlgebra,
    /// Every generator of an endomorphism ideal belongs to the parent algebra.
    #[error("an ideal generator is outside the represented endomorphism algebra")]
    EndomorphismIdealGeneratorOutsideAlgebra,
    /// Products of parent and ideal elements must remain in the parent.
    #[error("an endomorphism-ideal product left the represented parent algebra")]
    EndomorphismIdealProductOutsideAlgebra,
    /// Quotient multiplication requires one identical quotient family.
    #[error("the quotient elements belong to different endomorphism quotient families")]
    ForeignEndomorphismQuotientElement,
    /// Quotient elements carry one coefficient per quotient basis element.
    #[error("endomorphism quotient element has {actual} coordinates, expected {expected}")]
    EndomorphismQuotientElementCoordinateCount {
        /// Entered coordinate count.
        actual: usize,
        /// Quotient dimension.
        expected: usize,
    },
    /// The arrow-kernel construction must define a module satisfying the relations.
    #[error("the generated socle module does not satisfy the algebra relations")]
    GeneratedSocleViolatesRelations,
    /// The all-zero representation must satisfy every supported algebra relation.
    #[error("the generated zero module does not satisfy the algebra relations")]
    GeneratedZeroModuleViolatesRelations,
    /// Every direct summand must act over the same algebra.
    #[error("all direct summands must be modules over the same algebra")]
    DirectSumDifferentAlgebras,
    /// The source block-diagonal action must satisfy the common relations.
    #[error("the generated direct sum does not satisfy the algebra relations")]
    GeneratedDirectSumViolatesRelations,
    /// Direct-sum metadata can be installed only on its fresh module object.
    #[error("direct-sum metadata was already installed")]
    DirectSumMetadataAlreadySet,
    /// Projection and inclusion attributes require a QPA direct sum.
    #[error("the entered module was not built by DirectSumOfQPAModules")]
    NotDirectSumModule,
    /// An explicitly installed minimal generator must belong to the module.
    #[error("a minimal generator does not belong to the entered module")]
    MinimalGeneratorInWrongModule,
    /// Installing a minimal-generator attribute must make it readable.
    #[error("the minimal-generating-set cache could not be initialized")]
    MinimalGeneratingSetCacheFailed,
    /// Every internal use of a direct sum supplies at least one module.
    #[error("an internal nonempty DirectSumOfQPAModules call returned fail")]
    InternalNonemptyDirectSumFailed,
    /// A commutative square must carry the first kernel into the second.
    #[error("the induced kernel image unexpectedly has no preimage")]
    MissingInducedKernelPreimage,
    /// The canonical image projection must lift every image-basis vector.
    #[error("an image basis vector unexpectedly has no preimage")]
    MissingInducedImagePreimage,
    /// The canonical cokernel projection must lift every cokernel-basis vector.
    #[error("a cokernel basis vector unexpectedly has no preimage")]
    MissingInducedCokernelPreimage,
    /// A successful inclusion-containment test must provide every lift.
    #[error("an included basis image unexpectedly has no preimage")]
    MissingInclusionLiftPreimage,
    /// A basis-image map needs at least one image per source basis vector.
    #[error("basis-image list has length {actual}, expected at least {expected}")]
    BasisImageCount {
        /// Entered image count.
        actual: usize,
        /// Required source dimension.
        expected: usize,
    },
    /// QPA requires a nonzero basis image at the source basis vector's vertex.
    #[error("basis image {position} has the wrong vertex support")]
    WrongBasisImageSupport {
        /// One-based position in the entered image list.
        position: usize,
    },
    /// The map-matrix utility cannot infer a map family from no projections.
    #[error("MultiplyListsOfMaps requires a nonempty projection list")]
    EmptyMapProjectionList,
    /// The map-matrix utility cannot return a map from no inclusions.
    #[error("MultiplyListsOfMaps requires a nonempty inclusion list")]
    EmptyMapInclusionList,
    /// QPA indexes one matrix row for every inclusion.
    #[error("map matrix has {actual} rows, expected at least {expected}")]
    MapMatrixRowCount {
        /// Entered row count.
        actual: usize,
        /// Required row count.
        expected: usize,
    },
    /// Each consumed matrix row is dotted with the whole projection list.
    #[error("map matrix row {row} has {actual} columns, expected {expected}")]
    MapMatrixColumnCount {
        /// One-based row position.
        row: usize,
        /// Entered column count.
        actual: usize,
        /// Required projection count.
        expected: usize,
    },
    /// The first binary submodule argument must be a monomorphism.
    #[error("the first submodule argument is not a monomorphism")]
    FirstSubmoduleMapNotInjective,
    /// The second binary submodule argument must be a monomorphism.
    #[error("the second submodule argument is not a monomorphism")]
    SecondSubmoduleMapNotInjective,
    /// Submodule inclusions must have the same structural range.
    #[error("submodule homomorphisms must have the same range")]
    SubmoduleRangesDiffer,
    /// QPA requires a nonempty list of submodule inclusions.
    #[error("the submodule homomorphism list must be nonempty")]
    EmptySubmoduleList,
    /// Every source in a submodule list must act over the same algebra.
    #[error("all submodule-list sources must be modules over the same algebra")]
    SubmoduleListDifferentAlgebras,
    /// Every map in a submodule list must be a monomorphism.
    #[error("not all submodule-list arguments are monomorphisms")]
    SomeSubmoduleMapNotInjective,
    /// Equal ranges make QPA's internal pullback branch infallible.
    #[error("an internal pullback of maps with a common range returned fail")]
    InternalCommonRangePullbackFailed,
    /// `RejectOfModule` requires a common acting algebra.
    #[error("the reject modules are not modules over the same algebra")]
    RejectModulesOverDifferentAlgebras,
    /// A Hom query cannot reject modules after their algebras were checked.
    #[error("an internal HomOverAlgebra call rejected a common algebra")]
    InternalCommonAlgebraHomFailed,
    /// A module element needs one coordinate row per quiver vertex.
    #[error("module element has {actual} vertex rows, expected {expected}")]
    ElementVertexCount {
        /// Entered row count.
        actual: usize,
        /// Required vertex count.
        expected: usize,
    },
    /// A logical coordinate row must match its vertex-space dimension.
    #[error("module element row at vertex {vertex} has length {actual}, expected {expected}")]
    ElementCoordinateCount {
        /// One-based vertex position.
        vertex: usize,
        /// Entered coordinate count.
        actual: usize,
        /// Required logical dimension.
        expected: usize,
    },
    /// A flat logical coordinate list must have total module dimension.
    #[error("module element has {actual} flat coordinates, expected {expected}")]
    ElementFlatCoordinateCount {
        /// Entered coordinate count.
        actual: usize,
        /// Required total dimension.
        expected: usize,
    },
    /// GAP module-element arithmetic requires an identical element family.
    #[error("module elements belong to different element families")]
    DifferentElementFamilies,
    /// Right action requires an element of the identical acting algebra.
    #[error("the acting element belongs to a different algebra family")]
    ForeignActingElement,
    /// A normalized nonzero algebra element cannot contain QPA's zero path.
    #[error("the acting algebra element contains the zero path")]
    ZeroPathAction,
    /// Module homomorphisms require one common structurally equal algebra.
    #[error("the two modules are not over the same algebra")]
    DifferentHomomorphismAlgebras,
    /// A module homomorphism stores one matrix per vertex.
    #[error("received {actual} module-map matrices, expected {expected}")]
    HomomorphismMatrixCount {
        /// Entered matrix count.
        actual: usize,
        /// Required vertex count.
        expected: usize,
    },
    /// Every vertex map must use the acting algebra's coefficient field.
    #[error("module-map matrix at vertex {vertex} is over a different field")]
    HomomorphismMatrixField {
        /// One-based vertex position.
        vertex: usize,
    },
    /// A vertex map has the wrong physical shape or a nonzero zero-space map.
    #[error(
        "module-map matrix at vertex {vertex} has shape {actual:?}, expected zero matrix of shape {expected:?}"
    )]
    HomomorphismMatrixShape {
        /// One-based vertex position.
        vertex: usize,
        /// Entered shape.
        actual: (usize, usize),
        /// Required physical shape.
        expected: (usize, usize),
    },
    /// The entered vertex maps fail an arrow commutativity square.
    #[error("entered map does not commute with arrow {arrow:?}")]
    NoncommutingHomomorphismSquare {
        /// First failing arrow label.
        arrow: Box<str>,
    },
    /// Addition requires structurally equal source and range modules.
    #[error("the two maps do not belong to the same homomorphism set")]
    DifferentHomomorphismSets,
    /// Composition requires the first range to equal the second source.
    #[error("the first map's range is not the second map's source")]
    NoncomposableHomomorphisms,
    /// Mapping application requires an element of the identical source.
    #[error("the module element does not belong to the mapping source")]
    ForeignMapSourceElement,
    /// Preimage computation requires an element of the identical range.
    #[error("the module element does not belong to the mapping range")]
    ForeignMapRangeElement,
    /// Every generator of a subrepresentation must belong to the input module.
    #[error("a subrepresentation generator is not an element of the input module")]
    ForeignSubrepresentationGenerator,
    /// QPA's canonical subspace basis must remain vertex-uniform.
    #[error("the canonical generated-submodule basis is not vertex-uniform")]
    NonuniformSubrepresentationBasis,
    /// Coordinates in the generated submodule's canonical basis must exist.
    #[error("failed to express a generated-submodule element in its canonical basis")]
    InternalSubrepresentationCoordinates,
    /// A kernel action must remain in the next component kernel.
    #[error("failed to express a kernel action in its component kernel basis")]
    InternalKernelCoordinates,
    /// A range action and every source image must remain in the image basis.
    #[error("failed to express a module-map image in its canonical image basis")]
    InternalImageCoordinates,
    /// Every basis vector of a constructed top must have QPA's chosen lift.
    #[error("a source-top basis vector has no preimage under its top projection")]
    MissingTopBasisPreimage,
    /// A source-generated module unexpectedly failed its defining relations.
    #[error("a module reconstructed by a QPA source method failed the acting-algebra relations")]
    GeneratedModuleRejectedRelations,
    /// QPA's radical method applies only in its finite/admissible source domain.
    #[error(
        "RadicalOfModuleInclusion requires a finite path algebra or finite admissible quotient"
    )]
    UnsupportedRadicalAlgebra,
    /// An entered matrix does not use the acting algebra's field instance.
    #[error("matrix for arrow {label:?} is over a different coefficient field")]
    DifferentCoefficientField {
        /// Arrow label attached to the matrix.
        label: Box<str>,
    },
    /// At least one entered label is absent from the quiver's arrow list.
    #[error("arrow labels {labels:?} do not exist in the acting quiver")]
    UnknownArrowLabels {
        /// Missing labels in entered order, including repetitions.
        labels: Vec<Box<str>>,
    },
    /// QPA forbids a nonzero action incident with a zero-dimensional space.
    #[error("a nonzero matrix was entered for arrow {label:?} at a zero-dimensional space")]
    NonzeroMapAtZeroSpace {
        /// Rejected arrow label.
        label: Box<str>,
    },
    /// GAP cannot retain a zero-row or zero-column entered matrix; such an
    /// arrow must be omitted so the constructor creates its physical stand-in.
    #[error("the zero matrix for arrow {label:?} at a zero-dimensional space must be omitted")]
    ZeroSpaceMapMustBeOmitted {
        /// Rejected arrow label.
        label: Box<str>,
    },
    /// An entered arrow map has dimensions different from its source and target spaces.
    #[error("matrix for arrow {label:?} has shape {actual:?}, expected {expected:?}")]
    WrongArrowMatrixShape {
        /// Rejected arrow label.
        label: Box<str>,
        /// Entered shape.
        actual: (usize, usize),
        /// Source-by-target dimension.
        expected: (usize, usize),
    },
    /// The retained tensor metadata must wrap the identical acting algebra.
    #[error("retained tensor product does not wrap the module's acting algebra")]
    RetainedTensorAlgebraMismatch,
    /// Concurrent construction lost QPA's reciprocal dual-module attribute.
    #[error("concurrent dual-module cache initialization failed")]
    DualCacheRace,
    /// Relation evaluation encountered incompatible explicit matrices.
    #[error(transparent)]
    Matrix(#[from] MatrixError),
    /// Quiver-algebra element construction or arithmetic failed.
    #[error(transparent)]
    AlgebraElement(#[from] QuiverAlgebraElementError),
    /// Path-algebra element construction failed.
    #[error(transparent)]
    PathAlgebra(#[from] crate::PathAlgebraError),
    /// Quotient computation or normal-form reduction failed.
    #[error(transparent)]
    Quotient(#[from] QuotientError),
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use numerica::domains::{
        Field, Ring,
        finite_field::{Z2, Zp},
        integer::Integer,
        rational::{Q, Rational},
    };
    use serde_json::{Value, json};

    use crate::{
        ArrowSpec, DenseMatrix, PathAlgebra, Quiver, QuiverAlgebra, QuotientPathAlgebra,
        VertexIndex,
    };

    use super::{
        MoreMinimalVersionResult, NamedArrowMatrix, PathAlgebraModuleError,
        PathAlgebraModuleHomomorphism, PathAlgebraModuleMatrix, cokernel_projection,
        direct_sum_inclusions, direct_sum_of_qpa_modules, direct_sum_projections, end_over_algebra,
        hom_over_algebra, hom_over_algebra_with_basis_function, homomorphism_from_images,
        identity_module_mapping, image_projection_inclusion, intersection_of_submodules,
        intersection_of_submodules_list, inverse_of_isomorphism, is_split_epimorphism,
        is_split_monomorphism, isomorphism_of_modules, kernel_inclusion,
        left_inverse_of_homomorphism, left_minimal_version, minimal_generating_set_of_module,
        more_left_minimal_version, more_right_minimal_version, morphism_on_cokernel,
        morphism_on_image, morphism_on_kernel, multiply_lists_of_maps,
        natural_endomorphism_quotient_by_ideal, pull_back, push_out, radical_of_module_inclusion,
        reject_of_module, right_fac_approximation, right_fac_approximation_via_all_basis_images,
        right_inverse_of_homomorphism, right_minimal_version, right_module_hom_over_algebra,
        right_module_over_path_algebra, socle_of_module, socle_of_module_inclusion,
        sum_of_submodules, sum_of_submodules_list, top_of_module, top_of_module_homomorphism,
        top_of_module_projection, trace_of_module, zero_module_mapping,
    };

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn q(value: i64) -> Rational {
        Q.nth(value.into())
    }

    fn matrix_json<F: Field>(
        matrix: &DenseMatrix<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        json!({
            "shape": [matrix.nrows(), matrix.ncols()],
            "rows": (0..matrix.nrows())
                .map(|row| matrix.row(row).iter().map(coefficient).collect::<Vec<_>>())
                .collect::<Vec<_>>(),
        })
    }

    fn module_map_json<F: Field>(
        homomorphism: &PathAlgebraModuleHomomorphism<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        json!({
            "source_dimension": homomorphism.source().dimension_vector(),
            "range_dimension": homomorphism.range().dimension_vector(),
            "matrices": homomorphism.vertex_matrices().iter()
                .map(|matrix| matrix_json(matrix, coefficient)).collect::<Vec<_>>(),
        })
    }

    fn module_map_list_json<F: Field>(
        homomorphisms: &[&PathAlgebraModuleHomomorphism<F>],
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        json!(
            homomorphisms
                .iter()
                .map(|homomorphism| module_map_json(homomorphism, coefficient))
                .collect::<Vec<_>>()
        )
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn modulehom_split_case<F: Field>(
        field: F,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let zero = field.zero();
        let one = field.one();
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
        let layer = right_module_over_path_algebra(algebra.clone(), vec![1, 1], &[])
            .unwrap()
            .unwrap();
        let simple_one = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
            .unwrap()
            .unwrap();
        let simple_two = right_module_over_path_algebra(algebra.clone(), vec![0, 1], &[])
            .unwrap()
            .unwrap();
        let zero_module = right_module_over_path_algebra(algebra.clone(), vec![0, 0], &[])
            .unwrap()
            .unwrap();
        let projective_one = right_module_over_path_algebra(
            algebra,
            vec![1, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![one.clone()]], field.clone()).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let inclusion = right_module_hom_over_algebra(
            simple_one.clone(),
            layer.clone(),
            vec![
                DenseMatrix::from_rows(vec![vec![one.clone()]], field.clone()).unwrap(),
                DenseMatrix::from_rows(vec![vec![zero.clone()]], field.clone()).unwrap(),
            ],
        )
        .unwrap();
        let projection = right_module_hom_over_algebra(
            layer.clone(),
            simple_one.clone(),
            vec![
                DenseMatrix::from_rows(vec![vec![one.clone()]], field.clone()).unwrap(),
                DenseMatrix::from_rows(vec![vec![zero.clone()]], field.clone()).unwrap(),
            ],
        )
        .unwrap();
        let right_inverse = right_inverse_of_homomorphism(&inclusion).unwrap().unwrap();
        let left_inverse = left_inverse_of_homomorphism(&projection).unwrap().unwrap();
        let radical_inclusion = right_module_hom_over_algebra(
            simple_two,
            projective_one.clone(),
            vec![
                DenseMatrix::from_rows(vec![vec![zero.clone()]], field.clone()).unwrap(),
                DenseMatrix::from_rows(vec![vec![one.clone()]], field.clone()).unwrap(),
            ],
        )
        .unwrap();
        let top_projection = right_module_hom_over_algebra(
            projective_one,
            simple_one.clone(),
            vec![
                DenseMatrix::from_rows(vec![vec![one]], field.clone()).unwrap(),
                DenseMatrix::from_rows(vec![vec![zero]], field.clone()).unwrap(),
            ],
        )
        .unwrap();
        let zero_inclusion = zero_module_mapping(zero_module.clone(), layer.clone()).unwrap();
        let zero_projection = zero_module_mapping(layer.clone(), zero_module).unwrap();
        let identity = identity_module_mapping(layer.clone()).unwrap();
        let identity_right = right_inverse_of_homomorphism(&identity).unwrap().unwrap();
        let identity_left = left_inverse_of_homomorphism(&identity).unwrap().unwrap();
        let property_inclusion = right_module_hom_over_algebra(
            simple_one.clone(),
            layer.clone(),
            vec![
                DenseMatrix::identity(1, field.clone()),
                DenseMatrix::zero(1, 1, field.clone()),
            ],
        )
        .unwrap();
        let property_projection = right_module_hom_over_algebra(
            layer,
            simple_one.clone(),
            vec![
                DenseMatrix::identity(1, field.clone()),
                DenseMatrix::zero(1, 1, field),
            ],
        )
        .unwrap();
        let property_values = [
            is_split_monomorphism(&property_inclusion).unwrap(),
            is_split_epimorphism(&property_projection).unwrap(),
        ];
        let property_right = right_inverse_of_homomorphism(&property_inclusion)
            .unwrap()
            .unwrap();
        let property_left = left_inverse_of_homomorphism(&property_projection)
            .unwrap()
            .unwrap();

        json!({
            "right_inverse": module_map_json(&right_inverse, coefficient),
            "left_inverse": module_map_json(&left_inverse, coefficient),
            "triangles": [
                inclusion.compose(&right_inverse).unwrap()
                    == identity_module_mapping(simple_one.clone()).unwrap(),
                left_inverse.compose(&projection).unwrap()
                    == identity_module_mapping(simple_one.clone()).unwrap(),
            ],
            "cache_identity": [
                right_inverse_of_homomorphism(&inclusion).unwrap().unwrap()
                    .is_identical(&right_inverse),
                left_inverse_of_homomorphism(&projection).unwrap().unwrap()
                    .is_identical(&left_inverse),
                left_inverse_of_homomorphism(&right_inverse).unwrap().unwrap()
                    .is_identical(&inclusion),
                right_inverse_of_homomorphism(&left_inverse).unwrap().unwrap()
                    .is_identical(&projection),
            ],
            "split_properties": [
                is_split_monomorphism(&inclusion).unwrap(),
                is_split_epimorphism(&right_inverse).unwrap(),
                is_split_epimorphism(&projection).unwrap(),
                is_split_monomorphism(&left_inverse).unwrap(),
            ],
            "nonsplit": [
                right_inverse_of_homomorphism(&radical_inclusion).unwrap().is_none(),
                is_split_monomorphism(&radical_inclusion).unwrap(),
                left_inverse_of_homomorphism(&top_projection).unwrap().is_none(),
                is_split_epimorphism(&top_projection).unwrap(),
            ],
            "zero_space_bug": [
                right_inverse_of_homomorphism(&zero_inclusion).unwrap().is_none(),
                is_split_monomorphism(&zero_inclusion).unwrap(),
                left_inverse_of_homomorphism(&zero_projection).unwrap().is_none(),
                is_split_epimorphism(&zero_projection).unwrap(),
            ],
            "identity_right": module_map_json(&identity_right, coefficient),
            "identity_left": module_map_json(&identity_left, coefficient),
            "identity_cache_same": identity_right.is_identical(&identity_left),
            "property_values": property_values,
            "property_right": module_map_json(&property_right, coefficient),
            "property_left": module_map_json(&property_left, coefficient),
            "property_triangles": [
                property_inclusion.compose(&property_right).unwrap()
                    == identity_module_mapping(simple_one.clone()).unwrap(),
                property_left.compose(&property_projection).unwrap()
                    == identity_module_mapping(simple_one).unwrap(),
            ],
        })
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn modulehom_map_utils_case<F: Field>(
        field: F,
        integer: impl Fn(i64) -> F::Element + Copy,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let quiver = Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
        let two = right_module_over_path_algebra(algebra.clone(), vec![2], &[])
            .unwrap()
            .unwrap();
        let automorphism = right_module_hom_over_algebra(
            two.clone(),
            two.clone(),
            vec![
                DenseMatrix::from_rows(
                    vec![vec![integer(1), integer(1)], vec![integer(0), integer(1)]],
                    field.clone(),
                )
                .unwrap(),
            ],
        )
        .unwrap();
        let inverse = inverse_of_isomorphism(&automorphism).unwrap();
        let zero_module = right_module_over_path_algebra(algebra, vec![0], &[])
            .unwrap()
            .unwrap();
        let zero_inverse =
            inverse_of_isomorphism(&zero_module_mapping(zero_module.clone(), zero_module).unwrap())
                .unwrap();

        let quiver_a2 =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra_a2 = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver_a2).unwrap());
        let representation = right_module_over_path_algebra(
            algebra_a2,
            vec![2, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![integer(1)], vec![integer(0)]], field.clone())
                    .unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let source_map = right_module_hom_over_algebra(
            representation.clone(),
            representation.clone(),
            vec![
                DenseMatrix::from_rows(
                    vec![vec![integer(1), integer(1)], vec![integer(0), integer(1)]],
                    field.clone(),
                )
                .unwrap(),
                DenseMatrix::identity(1, field.clone()),
            ],
        )
        .unwrap();
        let basis = representation.basis_elements();
        let images = basis
            .iter()
            .map(|element| source_map.image(element).unwrap())
            .collect::<Vec<_>>();
        let reconstructed =
            homomorphism_from_images(&representation, &representation, &images).unwrap();
        let mut images_with_extra = images;
        images_with_extra.push(representation.zero_element());
        let reconstructed_extra =
            homomorphism_from_images(&representation, &representation, &images_with_extra).unwrap();

        let make_map = |rows: [[i64; 2]; 2]| {
            right_module_hom_over_algebra(
                two.clone(),
                two.clone(),
                vec![
                    DenseMatrix::from_rows(
                        rows.into_iter()
                            .map(|row| row.into_iter().map(integer).collect())
                            .collect(),
                        field.clone(),
                    )
                    .unwrap(),
                ],
            )
            .unwrap()
        };
        let p1 = make_map([[1, 1], [0, 1]]);
        let p2 = make_map([[1, 0], [1, 1]]);
        let map_matrix = vec![
            vec![make_map([[1, 0], [0, 0]]), make_map([[1, 0], [0, 1]])],
            vec![make_map([[0, 1], [1, 0]]), make_map([[1, 1], [1, 0]])],
        ];
        let i1 = make_map([[1, 0], [1, 1]]);
        let i2 = make_map([[0, 1], [1, 0]]);
        let product = multiply_lists_of_maps(
            &[p1.clone(), p2.clone()],
            &map_matrix,
            &[i1.clone(), i2.clone()],
        )
        .unwrap();
        let mut map_matrix_extra = map_matrix;
        map_matrix_extra.push(vec![p1.clone(), p2.clone()]);
        let product_extra =
            multiply_lists_of_maps(&[p1, p2], &map_matrix_extra, &[i1, i2]).unwrap();

        json!({
            "inverse": module_map_json(&inverse, coefficient),
            "inverse_triangles": [
                automorphism.compose(&inverse).unwrap()
                    == identity_module_mapping(two.clone()).unwrap(),
                inverse.compose(&automorphism).unwrap()
                    == identity_module_mapping(two).unwrap(),
            ],
            "zero_inverse": module_map_json(&zero_inverse, coefficient),
            "from_images": module_map_json(&reconstructed, coefficient),
            "from_images_equal": reconstructed == source_map,
            "extra_image_ignored": reconstructed_extra == source_map,
            "list_product": module_map_json(&product, coefficient),
            "extra_row_ignored": product_extra == product,
        })
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn modulehom_isomorphism_case<F: Field>(
        field: F,
        integer: impl Fn(i64) -> F::Element + Copy,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let quiver = Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
        let one = right_module_over_path_algebra(algebra.clone(), vec![1], &[])
            .unwrap()
            .unwrap();
        let two = right_module_over_path_algebra(algebra.clone(), vec![2], &[])
            .unwrap()
            .unwrap();
        let equal_one = right_module_over_path_algebra(algebra.clone(), vec![2], &[])
            .unwrap()
            .unwrap();
        let equal_two = right_module_over_path_algebra(algebra.clone(), vec![2], &[])
            .unwrap()
            .unwrap();
        let zero_module = right_module_over_path_algebra(algebra, vec![0], &[])
            .unwrap()
            .unwrap();
        let zero_isomorphism = isomorphism_of_modules(&zero_module, &zero_module)
            .unwrap()
            .unwrap();
        let equal_isomorphism = isomorphism_of_modules(&equal_one, &equal_two)
            .unwrap()
            .unwrap();
        let different_dimension = isomorphism_of_modules(&one, &two).unwrap();

        let other_quiver = Quiver::from_indexed(2, Vec::<ArrowSpec>::new()).unwrap();
        let other_algebra =
            QuiverAlgebra::Path(PathAlgebra::new(field.clone(), other_quiver).unwrap());
        let other_module = right_module_over_path_algebra(other_algebra, vec![1, 0], &[])
            .unwrap()
            .unwrap();
        let different_algebra = isomorphism_of_modules(&one, &other_module).unwrap();

        let quiver_a2 =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra_a2 = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver_a2).unwrap());
        let make_representation = |rows: [[i64; 2]; 2]| {
            right_module_over_path_algebra(
                algebra_a2.clone(),
                vec![2, 2],
                &[NamedArrowMatrix::new(
                    "a",
                    DenseMatrix::from_rows(
                        rows.into_iter()
                            .map(|row| row.into_iter().map(integer).collect())
                            .collect(),
                        field.clone(),
                    )
                    .unwrap(),
                )],
            )
            .unwrap()
            .unwrap()
        };
        let rank_one = make_representation([[1, 0], [0, 0]]);
        let conjugate_rank_one = make_representation([[1, 1], [1, 1]]);
        let rank_two = make_representation([[1, 0], [0, 1]]);
        let recursive = isomorphism_of_modules(&rank_one, &conjugate_rank_one)
            .unwrap()
            .unwrap();
        let reverse = isomorphism_of_modules(&conjugate_rank_one, &rank_one)
            .unwrap()
            .unwrap();
        let recursive_inverse = inverse_of_isomorphism(&recursive).unwrap();
        let nonisomorphic = isomorphism_of_modules(&rank_one, &rank_two).unwrap();

        json!({
            "zero": module_map_json(&zero_isomorphism, coefficient),
            "equal": module_map_json(&equal_isomorphism, coefficient),
            "equal_endpoints_identical": [
                equal_isomorphism.source().is_identical(&equal_one),
                equal_isomorphism.range().is_identical(&equal_two),
            ],
            "different_dimension": different_dimension.is_none(),
            "different_algebra": different_algebra.is_none(),
            "recursive": module_map_json(&recursive, coefficient),
            "recursive_isomorphism": recursive.is_isomorphism(),
            "recursive_triangles": [
                recursive.compose(&recursive_inverse).unwrap()
                    == identity_module_mapping(rank_one.clone()).unwrap(),
                recursive_inverse.compose(&recursive).unwrap()
                    == identity_module_mapping(conjugate_rank_one.clone()).unwrap(),
            ],
            "reverse": module_map_json(&reverse, coefficient),
            "nonisomorphic": nonisomorphic.is_none(),
        })
    }

    fn more_minimal_json<F: Field>(
        result: &MoreMinimalVersionResult<F>,
        original: &PathAlgebraModuleHomomorphism<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        match result {
            MoreMinimalVersionResult::Unchanged(homomorphism) => json!({
                "changed": false,
                "same_object": homomorphism.is_identical(original),
                "map": module_map_json(homomorphism, coefficient),
            }),
            MoreMinimalVersionResult::Reduced {
                minimal_part,
                removed_part,
            } => json!({
                "changed": true,
                "minimal_part": module_map_json(minimal_part, coefficient),
                "removed_part": module_map_json(removed_part, coefficient),
            }),
        }
    }

    fn minimal_version_json<F: Field>(
        result: &(
            PathAlgebraModuleHomomorphism<F>,
            Vec<super::PathAlgebraModule<F>>,
        ),
        original: &PathAlgebraModuleHomomorphism<F>,
        right: bool,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let property = if right {
            result.0.cached_is_right_minimal()
        } else {
            result.0.cached_is_left_minimal()
        };
        json!({
            "map": module_map_json(&result.0, coefficient),
            "removed_dimensions": result.1.iter()
                .map(super::PathAlgebraModule::dimension_vector).collect::<Vec<_>>(),
            "map_is_original": result.0.is_identical(original),
            "property_bound": property.is_some(),
            "property_value": property.unwrap_or(false),
        })
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn modulehom_minimal_case<F: Field>(
        field: F,
        integer: impl Fn(i64) -> F::Element + Copy,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let quiver = Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
        let one = right_module_over_path_algebra(algebra.clone(), vec![1], &[])
            .unwrap()
            .unwrap();
        let three = right_module_over_path_algebra(algebra, vec![3], &[])
            .unwrap()
            .unwrap();
        let right_matrix = || {
            DenseMatrix::from_rows(
                vec![vec![integer(1)], vec![integer(0)], vec![integer(0)]],
                field.clone(),
            )
            .unwrap()
        };
        let left_matrix = || {
            DenseMatrix::from_rows(
                vec![vec![integer(1), integer(0), integer(0)]],
                field.clone(),
            )
            .unwrap()
        };
        let right_more_map =
            right_module_hom_over_algebra(three.clone(), one.clone(), vec![right_matrix()])
                .unwrap();
        let right_more = more_right_minimal_version(&right_more_map).unwrap();
        let left_more_map =
            right_module_hom_over_algebra(one.clone(), three.clone(), vec![left_matrix()]).unwrap();
        let left_more = more_left_minimal_version(&left_more_map).unwrap();

        let right_map =
            right_module_hom_over_algebra(three.clone(), one.clone(), vec![right_matrix()])
                .unwrap();
        let right_result = right_minimal_version(&right_map).unwrap();
        let right_again = right_minimal_version(&right_map).unwrap();
        let left_map =
            right_module_hom_over_algebra(one.clone(), three, vec![left_matrix()]).unwrap();
        let left_result = left_minimal_version(&left_map).unwrap();
        let left_again = left_minimal_version(&left_map).unwrap();
        let identity = identity_module_mapping(one.clone()).unwrap();
        let right_more_unchanged = more_right_minimal_version(&identity).unwrap();
        let left_more_unchanged = more_left_minimal_version(&identity).unwrap();
        let right_identity = right_minimal_version(&identity).unwrap();
        let left_identity = left_minimal_version(&identity).unwrap();
        let zero_map = zero_module_mapping(one.clone(), one).unwrap();
        let right_zero = right_minimal_version(&zero_map).unwrap();

        let quiver_a2 =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra_a2 = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver_a2).unwrap());
        let projective_one = right_module_over_path_algebra(
            algebra_a2.clone(),
            vec![1, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![integer(1)]], field.clone()).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let simple_two = right_module_over_path_algebra(algebra_a2, vec![0, 1], &[])
            .unwrap()
            .unwrap();
        let direct_sum = direct_sum_of_qpa_modules(&[projective_one, simple_two])
            .unwrap()
            .unwrap();
        let projection = direct_sum_projections(&direct_sum).unwrap()[0].clone();
        let inclusion = direct_sum_inclusions(&direct_sum).unwrap()[0].clone();
        let a2_right = right_minimal_version(&projection).unwrap();
        let a2_left = left_minimal_version(&inclusion).unwrap();

        json!({
            "right_more": more_minimal_json(&right_more, &right_more_map, coefficient),
            "left_more": more_minimal_json(&left_more, &left_more_map, coefficient),
            "right_more_unchanged": more_minimal_json(
                &right_more_unchanged, &identity, coefficient),
            "left_more_unchanged": more_minimal_json(
                &left_more_unchanged, &identity, coefficient),
            "right_minimal": minimal_version_json(&right_result, &right_map, true, coefficient),
            "left_minimal": minimal_version_json(&left_result, &left_map, false, coefficient),
            "attribute_cache_identity": [
                right_result.0.is_identical(&right_again.0),
                right_result.1.iter().zip(&right_again.1)
                    .all(|(left, right)| left.is_identical(right)),
                left_result.0.is_identical(&left_again.0),
                left_result.1.iter().zip(&left_again.1)
                    .all(|(left, right)| left.is_identical(right)),
            ],
            "right_identity": minimal_version_json(&right_identity, &identity, true, coefficient),
            "left_identity": minimal_version_json(&left_identity, &identity, false, coefficient),
            "right_zero": minimal_version_json(&right_zero, &zero_map, true, coefficient),
            "a2_right": minimal_version_json(&a2_right, &projection, true, coefficient),
            "a2_left": minimal_version_json(&a2_left, &inclusion, false, coefficient),
        })
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn modulehom_categorical_case<F: Field>(
        field: F,
        integer: impl Fn(i64) -> F::Element + Copy,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
        let module = right_module_over_path_algebra(
            algebra.clone(),
            vec![2, 2],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(
                    vec![vec![integer(1), integer(1)], vec![integer(0), integer(0)]],
                    field.clone(),
                )
                .unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let f = right_module_hom_over_algebra(
            module.clone(),
            module.clone(),
            vec![
                DenseMatrix::from_rows(
                    vec![vec![integer(1), integer(0)], vec![integer(0), integer(0)]],
                    field.clone(),
                )
                .unwrap(),
                DenseMatrix::identity(2, field.clone()),
            ],
        )
        .unwrap();
        let t = right_module_hom_over_algebra(
            module.clone(),
            module.clone(),
            vec![
                DenseMatrix::from_rows(
                    vec![vec![integer(2), integer(0)], vec![integer(0), integer(3)]],
                    field.clone(),
                )
                .unwrap(),
                DenseMatrix::identity(2, field.clone()).scale_left(&integer(2)),
            ],
        )
        .unwrap();
        let identity = identity_module_mapping(module.clone()).unwrap();
        let layer = right_module_over_path_algebra(algebra.clone(), vec![1, 1], &[])
            .unwrap()
            .unwrap();
        let first = right_module_hom_over_algebra(
            layer.clone(),
            module.clone(),
            vec![
                DenseMatrix::from_rows(vec![vec![integer(0), integer(1)]], field.clone()).unwrap(),
                DenseMatrix::from_rows(vec![vec![integer(1), integer(0)]], field.clone()).unwrap(),
            ],
        )
        .unwrap();
        let second = right_module_hom_over_algebra(
            layer.clone(),
            module.clone(),
            vec![
                DenseMatrix::from_rows(vec![vec![integer(0), integer(1)]], field.clone()).unwrap(),
                DenseMatrix::from_rows(vec![vec![integer(0), integer(1)]], field.clone()).unwrap(),
            ],
        )
        .unwrap();
        let kernel_module = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
            .unwrap()
            .unwrap();
        let kernel_to_module = right_module_hom_over_algebra(
            kernel_module.clone(),
            module.clone(),
            vec![
                DenseMatrix::from_rows(vec![vec![integer(0), integer(1)]], field.clone()).unwrap(),
                DenseMatrix::zero(1, 2, field.clone()),
            ],
        )
        .unwrap();
        let simple_one = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
            .unwrap()
            .unwrap();
        let simple_two = right_module_over_path_algebra(algebra, vec![0, 1], &[])
            .unwrap()
            .unwrap();

        let direct_sum =
            direct_sum_of_qpa_modules(&[layer.clone(), kernel_module, simple_two.clone()])
                .unwrap()
                .unwrap();
        let projections = direct_sum_projections(&direct_sum).unwrap();
        let inclusions = direct_sum_inclusions(&direct_sum).unwrap();
        let singleton = direct_sum_of_qpa_modules(std::slice::from_ref(&layer))
            .unwrap()
            .unwrap();
        let induced_kernel = morphism_on_kernel(&f, &f, &t, &t).unwrap().unwrap();
        let induced_image = morphism_on_image(&f, &f, &t, &t).unwrap().unwrap();
        let induced_cokernel = morphism_on_cokernel(&f, &f, &t, &t).unwrap().unwrap();
        let intersection = intersection_of_submodules(&first, &second).unwrap();
        let intersection_list = intersection_of_submodules_list(&[
            first.clone(),
            second.clone(),
            kernel_to_module.clone(),
        ])
        .unwrap();
        let sum = sum_of_submodules(&first, &second).unwrap();
        let sum_list =
            sum_of_submodules_list(&[first.clone(), second.clone(), kernel_to_module]).unwrap();
        let pullback = pull_back(&first, &second).unwrap().unwrap();
        let pushout = push_out(&intersection.0, &intersection.1).unwrap().unwrap();
        let reject_zero = reject_of_module(&simple_one, &simple_two).unwrap();
        let reject_one = reject_of_module(&module, &simple_two).unwrap();
        let reject_many = reject_of_module(&module, &simple_one).unwrap();

        json!({
            "direct_sum": {
                "empty_is_fail": direct_sum_of_qpa_modules::<F>(&[]).unwrap().is_none(),
                "filter": direct_sum.is_direct_sum_of_modules(),
                "dimension": direct_sum.dimension_vector(),
                "matrices": direct_sum.arrow_matrices().iter()
                    .map(|matrix| matrix_json(matrix, coefficient)).collect::<Vec<_>>(),
                "projections": projections.iter()
                    .map(|map| module_map_json(map, coefficient)).collect::<Vec<_>>(),
                "inclusions": inclusions.iter()
                    .map(|map| module_map_json(map, coefficient)).collect::<Vec<_>>(),
                "projection_cache_identity": projections[0]
                    .is_identical(&direct_sum_projections(&direct_sum).unwrap()[0]),
                "inclusion_cache_identity": inclusions[0]
                    .is_identical(&direct_sum_inclusions(&direct_sum).unwrap()[0]),
                "singleton_equal": singleton == layer,
                "singleton_identical": singleton.is_identical(&layer),
            },
            "induced": {
                "kernel": module_map_json(&induced_kernel, coefficient),
                "image": module_map_json(&induced_image, coefficient),
                "cokernel": module_map_json(&induced_cokernel, coefficient),
                "noncommuting_fails": [
                    morphism_on_kernel(&f, &f, &identity, &t).unwrap().is_none(),
                    morphism_on_image(&f, &f, &identity, &t).unwrap().is_none(),
                    morphism_on_cokernel(&f, &f, &identity, &t).unwrap().is_none(),
                ],
            },
            "pullback": module_map_list_json(&[&pullback.0, &pullback.1], coefficient),
            "pushout": module_map_list_json(&[&pushout.0, &pushout.1], coefficient),
            "intersection": {
                "pair": module_map_list_json(
                    &[&intersection.0, &intersection.1, &intersection.2], coefficient),
                "list": module_map_json(&intersection_list, coefficient),
                "singleton_identity": intersection_of_submodules_list(
                    std::slice::from_ref(&first)).unwrap().is_identical(&first),
            },
            "sum": {
                "pair": module_map_list_json(&[&sum.0, &sum.1, &sum.2], coefficient),
                "list": module_map_json(&sum_list, coefficient),
                "singleton_identity": sum_of_submodules_list(
                    std::slice::from_ref(&first)).unwrap().is_identical(&first),
            },
            "reject": {
                "hom_zero": module_map_json(&reject_zero, coefficient),
                "hom_one": module_map_json(&reject_one, coefficient),
                "hom_many": module_map_json(&reject_many, coefficient),
            },
        })
    }

    #[test]
    fn checked_constructor_sets_omitted_zero_maps_in_arrow_order() {
        let quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(2), vertex(3), "b"),
            ],
        )
        .unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
        let a = DenseMatrix::from_rows(vec![vec![q(1), q(2)]], Q).unwrap();
        let module = right_module_over_path_algebra(
            algebra.clone(),
            vec![1, 2, 0],
            &[NamedArrowMatrix::new("a", a.clone())],
        )
        .unwrap()
        .unwrap();

        assert_eq!(module.right_acting_algebra(), &algebra);
        assert_eq!(module.dimension_vector(), &[1, 2, 0]);
        assert_eq!(module.dimension(), 3);
        assert_eq!(module.arrow_matrices()[0], a);
        assert_eq!(module.arrow_matrices()[1], DenseMatrix::zero(2, 1, Q));
        assert_eq!(
            module.matrices_of_path_algebra_module(),
            vec![
                PathAlgebraModuleMatrix::Matrix(module.arrow_matrices()[0].clone()),
                PathAlgebraModuleMatrix::Matrix(module.arrow_matrices()[1].clone()),
            ]
        );
    }

    #[test]
    fn zero_module_reports_gap_scalar_zero_entries() {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
        let module = right_module_over_path_algebra(algebra, vec![0, 0], &[])
            .unwrap()
            .unwrap();

        assert_eq!(module.arrow_matrices()[0], DenseMatrix::zero(1, 1, Q));
        assert_eq!(
            module.matrices_of_path_algebra_module(),
            vec![PathAlgebraModuleMatrix::ZeroScalar]
        );
    }

    #[test]
    fn quotient_relations_return_fail_when_not_satisfied() {
        let quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let parent = PathAlgebra::new(Q, quiver).unwrap();
        let x = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let quotient =
            QuotientPathAlgebra::by_relations_high_level(parent, &[x.multiply(&x).unwrap()])
                .unwrap();
        let algebra = QuiverAlgebra::Quotient(quotient);
        let nilpotent =
            DenseMatrix::from_rows(vec![vec![q(0), q(1)], vec![q(0), q(0)]], Q).unwrap();
        let identity = DenseMatrix::identity(2, Q);

        assert!(
            right_module_over_path_algebra(
                algebra.clone(),
                vec![2],
                &[NamedArrowMatrix::new("x", nilpotent)],
            )
            .unwrap()
            .is_some()
        );
        assert!(
            right_module_over_path_algebra(
                algebra,
                vec![2],
                &[NamedArrowMatrix::new("x", identity)],
            )
            .unwrap()
            .is_none()
        );
    }

    #[test]
    fn checked_constructor_preserves_qpa_input_errors() {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());

        assert_eq!(
            right_module_over_path_algebra(algebra.clone(), vec![1], &[]).unwrap_err(),
            PathAlgebraModuleError::DimensionVectorLength {
                actual: 1,
                expected: 2,
            }
        );
        assert_eq!(
            right_module_over_path_algebra(
                algebra.clone(),
                vec![1, 1],
                &[NamedArrowMatrix::new("missing", DenseMatrix::zero(1, 1, Q))],
            )
            .unwrap_err(),
            PathAlgebraModuleError::UnknownArrowLabels {
                labels: vec![Box::from("missing")],
            }
        );
        assert!(matches!(
            right_module_over_path_algebra(
                algebra.clone(),
                vec![0, 1],
                &[NamedArrowMatrix::new(
                    "a",
                    DenseMatrix::from_rows(vec![vec![q(1)]], Q).unwrap(),
                )],
            ),
            Err(PathAlgebraModuleError::NonzeroMapAtZeroSpace { .. })
        ));
        assert!(matches!(
            right_module_over_path_algebra(
                algebra.clone(),
                vec![0, 1],
                &[NamedArrowMatrix::new(
                    "a",
                    DenseMatrix::from_flat(0, 1, Vec::new(), Q).unwrap(),
                )],
            ),
            Err(PathAlgebraModuleError::ZeroSpaceMapMustBeOmitted { .. })
        ));
        assert!(matches!(
            right_module_over_path_algebra(
                algebra.clone(),
                vec![1, 2],
                &[NamedArrowMatrix::new("a", DenseMatrix::zero(1, 1, Q))],
            ),
            Err(PathAlgebraModuleError::WrongArrowMatrixShape { .. })
        ));

        let prime_three = Zp::new(3);
        let prime_five = Zp::new(5);
        let prime_quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let prime_algebra =
            QuiverAlgebra::Path(PathAlgebra::new(prime_three, prime_quiver).unwrap());
        assert!(matches!(
            right_module_over_path_algebra(
                prime_algebra,
                vec![1, 1],
                &[NamedArrowMatrix::new(
                    "a",
                    DenseMatrix::zero(1, 1, prime_five),
                )],
            ),
            Err(PathAlgebraModuleError::DifferentCoefficientField { .. })
        ));
    }

    #[test]
    fn clones_are_identical_but_reconstruction_is_only_equal() {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
        let first = right_module_over_path_algebra(algebra.clone(), vec![1, 1], &[])
            .unwrap()
            .unwrap();
        let clone = first.clone();
        let second = right_module_over_path_algebra(algebra, vec![1, 1], &[])
            .unwrap()
            .unwrap();

        assert!(first.is_identical(&clone));
        assert!(!first.is_identical(&second));
        assert_eq!(first, second);
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn direct_sum_uses_qpa_block_order_component_maps_and_live_identity_cache() {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
        let first = right_module_over_path_algebra(
            algebra.clone(),
            vec![1, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![q(1)]], Q).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let second = right_module_over_path_algebra(
            algebra.clone(),
            vec![2, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![q(1)], vec![q(0)]], Q).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let third = right_module_over_path_algebra(algebra.clone(), vec![0, 1], &[])
            .unwrap()
            .unwrap();

        assert!(
            direct_sum_of_qpa_modules::<numerica::domains::rational::RationalField>(&[])
                .unwrap()
                .is_none()
        );
        let direct_sum = direct_sum_of_qpa_modules(&[first.clone(), second.clone(), third.clone()])
            .unwrap()
            .unwrap();
        assert!(direct_sum.is_direct_sum_of_modules());
        assert_eq!(direct_sum.dimension_vector(), &[3, 3]);
        assert_eq!(
            direct_sum.arrow_matrices()[0],
            DenseMatrix::from_rows(
                vec![
                    vec![q(1), q(0), q(0)],
                    vec![q(0), q(1), q(0)],
                    vec![q(0), q(0), q(0)],
                ],
                Q,
            )
            .unwrap()
        );
        let summands = direct_sum.direct_sum_summands().unwrap();
        assert!(summands[0].is_identical(&first));
        assert!(summands[1].is_identical(&second));
        assert!(summands[2].is_identical(&third));

        let projections = direct_sum_projections(&direct_sum).unwrap();
        let inclusions = direct_sum_inclusions(&direct_sum).unwrap();
        assert_eq!(
            projections[1].vertex_matrices(),
            &[
                DenseMatrix::from_rows(
                    vec![vec![q(0), q(0)], vec![q(1), q(0)], vec![q(0), q(1)]],
                    Q,
                )
                .unwrap(),
                DenseMatrix::from_rows(vec![vec![q(0)], vec![q(1)], vec![q(0)]], Q).unwrap(),
            ]
        );
        assert_eq!(
            inclusions[2].vertex_matrices(),
            &[
                DenseMatrix::zero(1, 3, Q),
                DenseMatrix::from_rows(vec![vec![q(0), q(0), q(1)]], Q).unwrap(),
            ]
        );
        for index in 0..3 {
            assert!(projections[index].is_surjective());
            assert!(inclusions[index].is_injective());
            assert_eq!(
                inclusions[index].compose(&projections[index]).unwrap(),
                identity_module_mapping(summands[index].clone()).unwrap()
            );
            for (other, projection) in projections.iter().enumerate() {
                if index != other {
                    assert!(inclusions[index].compose(projection).unwrap().is_zero());
                }
            }
        }
        let identity_decomposition = projections
            .iter()
            .zip(&inclusions)
            .map(|(projection, inclusion)| projection.compose(inclusion).unwrap())
            .reduce(|left, right| left.add(&right).unwrap())
            .unwrap();
        assert_eq!(
            identity_decomposition,
            identity_module_mapping(direct_sum.clone()).unwrap()
        );

        let retained_projection = projections[0].clone();
        drop(projections);
        let regenerated = direct_sum_projections(&direct_sum).unwrap();
        assert!(retained_projection.is_identical(&regenerated[0]));
        let repeated = direct_sum_projections(&direct_sum).unwrap();
        assert!(
            regenerated
                .iter()
                .zip(repeated)
                .all(|(left, right)| left.is_identical(&right))
        );

        let singleton = direct_sum_of_qpa_modules(std::slice::from_ref(&first))
            .unwrap()
            .unwrap();
        assert!(!singleton.is_identical(&first));
        assert_eq!(singleton, first);
        assert!(direct_sum_projections(&first).is_err());
    }

    #[test]
    fn characteristic_two_direct_sum_uses_the_same_qpa_block_algorithm() {
        let e = |value: i64| Z2.nth(Integer::from(value));
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, quiver).unwrap());
        let first = right_module_over_path_algebra(
            algebra.clone(),
            vec![1, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![e(1)]], Z2).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let second = right_module_over_path_algebra(algebra, vec![1, 0], &[])
            .unwrap()
            .unwrap();
        let direct_sum = direct_sum_of_qpa_modules(&[first, second])
            .unwrap()
            .unwrap();

        assert_eq!(direct_sum.dimension_vector(), &[2, 1]);
        assert_eq!(
            direct_sum.arrow_matrices()[0],
            DenseMatrix::from_rows(vec![vec![e(1)], vec![e(0)]], Z2).unwrap()
        );
        let projections = direct_sum_projections(&direct_sum).unwrap();
        let inclusions = direct_sum_inclusions(&direct_sum).unwrap();
        assert!(
            inclusions[0]
                .compose(&projections[0])
                .unwrap()
                .is_isomorphism()
        );
        assert!(
            inclusions[1]
                .compose(&projections[1])
                .unwrap()
                .is_isomorphism()
        );
    }

    #[test]
    fn pinned_qpa_modulehom_categorical_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehom_categorical.schema-v1.json"
        ))
        .unwrap();
        let rational =
            modulehom_categorical_case(Q, |value| Q.nth(Integer::from(value)), ToString::to_string);
        let binary = modulehom_categorical_case(
            Z2,
            |value| Z2.nth(Integer::from(value)),
            |entry| {
                if Z2.is_zero(entry) {
                    "0*Z(2)".to_owned()
                } else {
                    "Z(2)^0".to_owned()
                }
            },
        );

        assert_eq!(rational, expected["rational"]);
        assert_eq!(binary, expected["binary"]);
    }

    #[test]
    fn pinned_qpa_modulehom_split_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehom_split.schema-v1.json"
        ))
        .unwrap();
        let rational = modulehom_split_case(Q, ToString::to_string);
        let binary = modulehom_split_case(Z2, |entry| {
            if Z2.is_zero(entry) {
                "0*Z(2)".to_owned()
            } else {
                "Z(2)^0".to_owned()
            }
        });

        assert_eq!(rational, expected["rational"]);
        assert_eq!(binary, expected["binary"]);
    }

    #[test]
    fn pinned_qpa_modulehom_map_utils_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehom_map_utils.schema-v1.json"
        ))
        .unwrap();
        let rational =
            modulehom_map_utils_case(Q, |value| Q.nth(Integer::from(value)), ToString::to_string);
        let binary = modulehom_map_utils_case(
            Z2,
            |value| Z2.nth(Integer::from(value)),
            |entry| {
                if Z2.is_zero(entry) {
                    "0*Z(2)".to_owned()
                } else {
                    "Z(2)^0".to_owned()
                }
            },
        );

        assert_eq!(rational, expected["rational"]);
        assert_eq!(binary, expected["binary"]);
    }

    #[test]
    fn pinned_qpa_modulehom_isomorphism_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehom_isomorphism.schema-v1.json"
        ))
        .unwrap();
        let rational =
            modulehom_isomorphism_case(Q, |value| Q.nth(Integer::from(value)), ToString::to_string);
        let binary = modulehom_isomorphism_case(
            Z2,
            |value| Z2.nth(Integer::from(value)),
            |entry| {
                if Z2.is_zero(entry) {
                    "0*Z(2)".to_owned()
                } else {
                    "Z(2)^0".to_owned()
                }
            },
        );

        assert_eq!(rational, expected["rational"]);
        assert_eq!(binary, expected["binary"]);
    }

    #[test]
    fn map_utility_failures_preserve_qpa_validation_boundaries() {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
        let module = right_module_over_path_algebra(
            algebra,
            vec![2, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![q(1)], vec![q(0)]], Q).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let zero = zero_module_mapping(module.clone(), module.clone()).unwrap();
        assert_eq!(
            inverse_of_isomorphism(&zero).unwrap_err(),
            PathAlgebraModuleError::EnteredMapNotIsomorphism
        );

        let basis = module.basis_elements();
        let mut wrong_support = basis.clone();
        wrong_support[0] = basis[2].clone();
        assert_eq!(
            homomorphism_from_images(&module, &module, &wrong_support).unwrap_err(),
            PathAlgebraModuleError::WrongBasisImageSupport { position: 1 }
        );
        assert_eq!(
            homomorphism_from_images(&module, &module, &basis[..2]).unwrap_err(),
            PathAlgebraModuleError::BasisImageCount {
                actual: 2,
                expected: 3,
            }
        );
        assert_eq!(
            multiply_lists_of_maps::<numerica::domains::rational::Q>(&[], &[], &[]).unwrap_err(),
            PathAlgebraModuleError::EmptyMapProjectionList
        );
    }

    #[test]
    fn pinned_qpa_modulehom_minimal_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehom_minimal.schema-v1.json"
        ))
        .unwrap();
        let rational =
            modulehom_minimal_case(Q, |value| Q.nth(Integer::from(value)), ToString::to_string);
        let binary = modulehom_minimal_case(
            Z2,
            |value| Z2.nth(Integer::from(value)),
            |entry| {
                if Z2.is_zero(entry) {
                    "0*Z(2)".to_owned()
                } else {
                    "Z(2)^0".to_owned()
                }
            },
        );

        assert_eq!(rational, expected["rational"]);
        assert_eq!(binary, expected["binary"]);
    }

    #[test]
    fn left_minimal_zero_endomorphism_preserves_qpa_f_times_false_failure() {
        let quiver = Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
        let one = right_module_over_path_algebra(algebra, vec![1], &[])
            .unwrap()
            .unwrap();
        let zero = zero_module_mapping(one.clone(), one).unwrap();

        assert_eq!(
            more_left_minimal_version(&zero).unwrap_err(),
            PathAlgebraModuleError::InternalMinimalKernelNotSplit
        );
        assert_eq!(
            left_minimal_version(&zero).unwrap_err(),
            PathAlgebraModuleError::InternalMinimalKernelNotSplit
        );
    }

    #[test]
    fn categorical_validation_preserves_qpa_fail_and_error_precedence() {
        let quiver = Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
        let one = right_module_over_path_algebra(algebra.clone(), vec![1], &[])
            .unwrap()
            .unwrap();
        let two = right_module_over_path_algebra(algebra.clone(), vec![2], &[])
            .unwrap()
            .unwrap();
        let three = right_module_over_path_algebra(algebra.clone(), vec![3], &[])
            .unwrap()
            .unwrap();
        let inclusion = right_module_hom_over_algebra(
            one.clone(),
            two.clone(),
            vec![DenseMatrix::from_rows(vec![vec![q(1), q(0)]], Q).unwrap()],
        )
        .unwrap();
        let zero = zero_module_mapping(one.clone(), two.clone()).unwrap();
        let different_range = right_module_hom_over_algebra(
            one.clone(),
            three,
            vec![DenseMatrix::from_rows(vec![vec![q(1), q(0), q(0)]], Q).unwrap()],
        )
        .unwrap();

        assert!(
            push_out(&inclusion, &identity_module_mapping(two.clone()).unwrap())
                .unwrap()
                .is_none()
        );
        assert!(pull_back(&inclusion, &different_range).unwrap().is_none());
        assert_eq!(
            intersection_of_submodules(&zero, &inclusion).unwrap_err(),
            PathAlgebraModuleError::FirstSubmoduleMapNotInjective
        );
        assert_eq!(
            sum_of_submodules(&inclusion, &zero).unwrap_err(),
            PathAlgebraModuleError::SecondSubmoduleMapNotInjective
        );
        assert_eq!(
            intersection_of_submodules(&inclusion, &different_range).unwrap_err(),
            PathAlgebraModuleError::SubmoduleRangesDiffer
        );
        assert_eq!(
            intersection_of_submodules_list::<numerica::domains::rational::RationalField>(&[])
                .unwrap_err(),
            PathAlgebraModuleError::EmptySubmoduleList
        );
        assert_eq!(
            sum_of_submodules_list::<numerica::domains::rational::RationalField>(&[]).unwrap_err(),
            PathAlgebraModuleError::EmptySubmoduleList
        );
        assert_eq!(
            intersection_of_submodules_list(&[
                inclusion.clone(),
                zero.clone(),
                different_range.clone(),
            ])
            .unwrap_err(),
            PathAlgebraModuleError::SomeSubmoduleMapNotInjective
        );
        assert_eq!(
            sum_of_submodules_list(&[inclusion.clone(), zero, different_range]).unwrap_err(),
            PathAlgebraModuleError::SubmoduleRangesDiffer
        );

        let foreign_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let foreign_algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, foreign_quiver).unwrap());
        let foreign = right_module_over_path_algebra(
            foreign_algebra,
            vec![1],
            &[NamedArrowMatrix::new("x", DenseMatrix::zero(1, 1, Q))],
        )
        .unwrap()
        .unwrap();
        let foreign_identity = identity_module_mapping(foreign.clone()).unwrap();
        assert_eq!(
            direct_sum_of_qpa_modules(&[one.clone(), foreign.clone()]).unwrap_err(),
            PathAlgebraModuleError::DirectSumDifferentAlgebras
        );
        assert_eq!(
            intersection_of_submodules_list(&[inclusion, foreign_identity]).unwrap_err(),
            PathAlgebraModuleError::SubmoduleListDifferentAlgebras
        );
        assert_eq!(
            reject_of_module(&one, &foreign).unwrap_err(),
            PathAlgebraModuleError::RejectModulesOverDifferentAlgebras
        );
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn foundational_module_map_factorizations_follow_qpa_row_bases() {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
        let action = DenseMatrix::from_rows(vec![vec![q(1), q(1)], vec![q(0), q(0)]], Q).unwrap();
        let module = right_module_over_path_algebra(
            algebra,
            vec![2, 2],
            &[NamedArrowMatrix::new("a", action)],
        )
        .unwrap()
        .unwrap();
        let map = right_module_hom_over_algebra(
            module.clone(),
            module.clone(),
            vec![
                DenseMatrix::from_rows(vec![vec![q(1), q(0)], vec![q(0), q(0)]], Q).unwrap(),
                DenseMatrix::identity(2, Q),
            ],
        )
        .unwrap();

        assert!(!map.is_injective());
        assert!(!map.is_surjective());
        assert!(!map.is_isomorphism());
        assert!(!map.is_zero());

        let kernel = kernel_inclusion(&map).unwrap();
        assert_eq!(kernel.source().dimension_vector(), &[1, 0]);
        assert_eq!(
            kernel.vertex_matrices(),
            &[
                DenseMatrix::from_rows(vec![vec![q(0), q(1)]], Q).unwrap(),
                DenseMatrix::zero(1, 2, Q),
            ]
        );
        assert!(kernel.is_injective());
        assert!(kernel.kernel_of_what().unwrap().is_identical(&map));
        assert!(kernel_inclusion(&map).unwrap().is_identical(&kernel));
        assert!(kernel.compose(&map).unwrap().is_zero());

        let (projection, inclusion) = image_projection_inclusion(&map).unwrap();
        assert_eq!(projection.range().dimension_vector(), &[1, 2]);
        assert_eq!(
            projection.vertex_matrices(),
            &[
                DenseMatrix::from_rows(vec![vec![q(1)], vec![q(0)]], Q).unwrap(),
                DenseMatrix::identity(2, Q),
            ]
        );
        assert_eq!(
            inclusion.vertex_matrices(),
            &[
                DenseMatrix::from_rows(vec![vec![q(1), q(0)]], Q).unwrap(),
                DenseMatrix::identity(2, Q),
            ]
        );
        assert!(projection.is_surjective());
        assert!(inclusion.is_injective());
        assert!(projection.image_of_what().unwrap().is_identical(&map));
        assert!(inclusion.image_of_what().unwrap().is_identical(&map));
        assert!(
            image_projection_inclusion(&map)
                .unwrap()
                .0
                .is_identical(&projection)
        );
        assert_eq!(projection.compose(&inclusion).unwrap(), map);

        let cokernel = cokernel_projection(&map).unwrap();
        assert_eq!(cokernel.range().dimension_vector(), &[1, 0]);
        assert_eq!(
            cokernel.vertex_matrices(),
            &[
                DenseMatrix::from_rows(vec![vec![q(0)], vec![q(1)]], Q).unwrap(),
                DenseMatrix::zero(2, 1, Q),
            ]
        );
        assert!(cokernel.is_surjective());
        assert!(cokernel.cokernel_of_what().unwrap().is_identical(&map));
        assert!(cokernel_projection(&map).unwrap().is_identical(&cokernel));
        assert!(map.compose(&cokernel).unwrap().is_zero());

        let zero = zero_module_mapping(module.clone(), module.clone()).unwrap();
        let identity = identity_module_mapping(module.clone()).unwrap();
        assert_eq!(map.add(&zero).unwrap(), map);
        assert!(map.add(&map.negate().unwrap()).unwrap().is_zero());
        assert!(identity.is_isomorphism());
        assert_eq!(identity.compose(&map).unwrap(), map);
        assert_eq!(map.compose(&identity).unwrap(), map);
        assert!(zero.zero().unwrap().is_zero());

        let top_map = top_of_module_homomorphism(&map).unwrap();
        assert_eq!(top_map.source().dimension_vector(), &[2, 1]);
        assert_eq!(top_map.range().dimension_vector(), &[2, 1]);
        assert_eq!(
            top_map.vertex_matrices(),
            &[
                DenseMatrix::from_rows(vec![vec![q(1), q(0)], vec![q(0), q(0)]], Q).unwrap(),
                DenseMatrix::identity(1, Q),
            ]
        );

        let hom_basis = hom_over_algebra(&module, &module).unwrap().unwrap();
        assert_eq!(hom_basis.len(), 5);
        let z = DenseMatrix::zero(2, 2, Q);
        let expected = vec![
            vec![
                DenseMatrix::from_rows(vec![vec![q(0), q(1)], vec![q(0), q(0)]], Q).unwrap(),
                z.clone(),
            ],
            vec![
                DenseMatrix::from_rows(vec![vec![q(0), q(0)], vec![q(0), q(1)]], Q).unwrap(),
                z.clone(),
            ],
            vec![
                DenseMatrix::from_rows(vec![vec![q(1), q(0)], vec![q(0), q(0)]], Q).unwrap(),
                DenseMatrix::from_rows(vec![vec![q(1), q(1)], vec![q(0), q(0)]], Q).unwrap(),
            ],
            vec![
                z,
                DenseMatrix::from_rows(vec![vec![q(-1), q(0)], vec![q(1), q(0)]], Q).unwrap(),
            ],
            map.vertex_matrices().to_vec(),
        ];
        assert_eq!(
            hom_basis
                .iter()
                .map(|homomorphism| homomorphism.vertex_matrices().to_vec())
                .collect::<Vec<_>>(),
            expected
        );
        let hom_with_basis = hom_over_algebra_with_basis_function(&module, &module)
            .unwrap()
            .unwrap();
        assert_eq!(hom_with_basis.basis(), hom_basis);
        assert_eq!(
            hom_with_basis.coefficients(&map).unwrap().unwrap(),
            vec![q(0), q(0), q(0), q(0), q(1)]
        );

        let endomorphism_algebra = end_over_algebra(&module).unwrap();
        assert_eq!(endomorphism_algebra.dimension(), 5);
        assert!(
            endomorphism_algebra.is_identical(&end_over_algebra(&module).unwrap()),
            "QPA declares EndOverAlgebra as a cached attribute"
        );
        assert_eq!(
            endomorphism_algebra.basis_matrices(),
            &[
                DenseMatrix::from_rows(
                    vec![
                        vec![q(0), q(1), q(0), q(0)],
                        vec![q(0), q(0), q(0), q(0)],
                        vec![q(0), q(0), q(0), q(0)],
                        vec![q(0), q(0), q(0), q(0)],
                    ],
                    Q,
                )
                .unwrap(),
                DenseMatrix::from_rows(
                    vec![
                        vec![q(0), q(0), q(0), q(0)],
                        vec![q(0), q(1), q(0), q(0)],
                        vec![q(0), q(0), q(0), q(0)],
                        vec![q(0), q(0), q(0), q(0)],
                    ],
                    Q,
                )
                .unwrap(),
                DenseMatrix::from_rows(
                    vec![
                        vec![q(1), q(0), q(0), q(0)],
                        vec![q(0), q(0), q(0), q(0)],
                        vec![q(0), q(0), q(1), q(1)],
                        vec![q(0), q(0), q(0), q(0)],
                    ],
                    Q,
                )
                .unwrap(),
                DenseMatrix::from_rows(
                    vec![
                        vec![q(0), q(0), q(0), q(0)],
                        vec![q(0), q(0), q(0), q(0)],
                        vec![q(0), q(0), q(1), q(0)],
                        vec![q(0), q(0), q(-1), q(0)],
                    ],
                    Q,
                )
                .unwrap(),
                DenseMatrix::from_rows(
                    vec![
                        vec![q(0), q(0), q(0), q(0)],
                        vec![q(0), q(0), q(0), q(0)],
                        vec![q(0), q(0), q(0), q(1)],
                        vec![q(0), q(0), q(0), q(-1)],
                    ],
                    Q,
                )
                .unwrap(),
            ]
        );
        assert_eq!(*endomorphism_algebra.one(), DenseMatrix::identity(4, Q));
        assert_eq!(
            endomorphism_algebra
                .coefficients(endomorphism_algebra.one())
                .unwrap()
                .unwrap(),
            vec![q(0), q(1), q(1), q(0), q(-1)]
        );
        assert_eq!(
            endomorphism_algebra
                .multiply(
                    &endomorphism_algebra.basis_matrices()[0],
                    endomorphism_algebra.one(),
                )
                .unwrap(),
            endomorphism_algebra.basis_matrices()[0]
        );
        let cached_data = Arc::as_ptr(&endomorphism_algebra.0);
        drop(endomorphism_algebra);
        assert_eq!(
            Arc::as_ptr(&end_over_algebra(&module).unwrap().0),
            cached_data,
            "the GAP attribute cache survives after the caller drops its handle"
        );

        let basis = module.basis_elements();
        let nonzero_images = map.images_set(&basis).unwrap();
        assert_eq!(nonzero_images.len(), 3);
        assert!(nonzero_images.iter().all(|image| !image.is_zero()));
    }

    #[test]
    fn characteristic_two_module_map_factorizations_match_the_same_qpa_fixture() {
        let e = |value: i64| Z2.nth(Integer::from(value));
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, quiver).unwrap());
        let module = right_module_over_path_algebra(
            algebra,
            vec![2, 2],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![e(1), e(1)], vec![e(0), e(0)]], Z2).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let map = right_module_hom_over_algebra(
            module.clone(),
            module,
            vec![
                DenseMatrix::from_rows(vec![vec![e(1), e(0)], vec![e(0), e(0)]], Z2).unwrap(),
                DenseMatrix::identity(2, Z2),
            ],
        )
        .unwrap();

        let kernel = kernel_inclusion(&map).unwrap();
        assert_eq!(kernel.source().dimension_vector(), &[1, 0]);
        assert_eq!(
            kernel.vertex_matrices()[0],
            DenseMatrix::from_rows(vec![vec![e(0), e(1)]], Z2).unwrap()
        );
        let (projection, inclusion) = image_projection_inclusion(&map).unwrap();
        assert_eq!(projection.range().dimension_vector(), &[1, 2]);
        assert_eq!(projection.compose(&inclusion).unwrap(), map);
        let cokernel = cokernel_projection(&map).unwrap();
        assert_eq!(cokernel.range().dimension_vector(), &[1, 0]);
        assert_eq!(
            cokernel.vertex_matrices()[0],
            DenseMatrix::from_rows(vec![vec![e(0)], vec![e(1)]], Z2).unwrap()
        );
        let hom_basis = hom_over_algebra(map.source(), map.range())
            .unwrap()
            .unwrap();
        assert_eq!(hom_basis.len(), 5);
        assert_eq!(
            hom_basis[3].vertex_matrices()[1],
            DenseMatrix::from_rows(vec![vec![e(1), e(0)], vec![e(1), e(0)]], Z2,).unwrap()
        );
        assert_eq!(
            hom_over_algebra_with_basis_function(map.source(), map.range())
                .unwrap()
                .unwrap()
                .coefficients(&map)
                .unwrap()
                .unwrap(),
            vec![e(0), e(0), e(0), e(0), e(1)]
        );
        let endomorphism_algebra = end_over_algebra(map.source()).unwrap();
        assert_eq!(endomorphism_algebra.dimension(), 5);
        assert_eq!(
            endomorphism_algebra.basis_matrices()[3],
            DenseMatrix::from_rows(
                vec![
                    vec![e(0), e(0), e(0), e(0)],
                    vec![e(0), e(0), e(0), e(0)],
                    vec![e(0), e(0), e(1), e(0)],
                    vec![e(0), e(0), e(1), e(0)],
                ],
                Z2,
            )
            .unwrap()
        );
        assert_eq!(
            endomorphism_algebra
                .coefficients(endomorphism_algebra.one())
                .unwrap()
                .unwrap(),
            vec![e(0), e(1), e(1), e(0), e(1)]
        );
    }

    #[test]
    fn unconstrained_and_disjoint_hom_branches_match_qpa() {
        let free_quiver = Quiver::from_indexed(1, vec![]).unwrap();
        let free_algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, free_quiver).unwrap());
        let free_source = right_module_over_path_algebra(free_algebra.clone(), vec![1], &[])
            .unwrap()
            .unwrap();
        let free_range = right_module_over_path_algebra(free_algebra, vec![2], &[])
            .unwrap()
            .unwrap();
        let free_basis = hom_over_algebra(&free_source, &free_range)
            .unwrap()
            .unwrap();
        assert_eq!(free_basis.len(), 2);
        assert_eq!(
            free_basis
                .iter()
                .map(|homomorphism| homomorphism.vertex_matrices()[0].clone())
                .collect::<Vec<_>>(),
            vec![
                DenseMatrix::from_rows(vec![vec![q(1), q(0)]], Q).unwrap(),
                DenseMatrix::from_rows(vec![vec![q(0), q(1)]], Q).unwrap(),
            ]
        );

        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
        let disjoint_source = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
            .unwrap()
            .unwrap();
        let disjoint_range = right_module_over_path_algebra(algebra, vec![0, 1], &[])
            .unwrap()
            .unwrap();
        assert!(
            hom_over_algebra(&disjoint_source, &disjoint_range)
                .unwrap()
                .unwrap()
                .is_empty()
        );
        let zero = zero_module_mapping(disjoint_source.clone(), disjoint_range.clone()).unwrap();
        assert_eq!(
            hom_over_algebra_with_basis_function(&disjoint_source, &disjoint_range)
                .unwrap()
                .unwrap()
                .coefficients(&zero)
                .unwrap()
                .unwrap(),
            vec![q(0)]
        );

        let zero_module = right_module_over_path_algebra(
            disjoint_source.right_acting_algebra().clone(),
            vec![0, 0],
            &[],
        )
        .unwrap()
        .unwrap();
        assert_eq!(
            end_over_algebra(&zero_module).unwrap_err(),
            PathAlgebraModuleError::EmptyEndomorphismDirectSum
        );
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn rational_minimal_generators_socle_and_trace_match_qpa() {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
        let module = right_module_over_path_algebra(
            algebra.clone(),
            vec![2, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![q(1)], vec![q(0)]], Q).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();

        let minimal = minimal_generating_set_of_module(module.clone()).unwrap();
        assert_eq!(minimal.len(), 2);
        assert_eq!(
            minimal[0].external_representation(),
            &[vec![q(1), q(0)], vec![q(0)]]
        );
        assert_eq!(
            minimal[1].external_representation(),
            &[vec![q(0), q(1)], vec![q(0)]]
        );
        let top_projection = top_of_module_projection(module.clone()).unwrap();
        assert!(top_projection.is_identical(&top_of_module_projection(module.clone()).unwrap()));
        let top = top_of_module(module.clone()).unwrap();
        assert!(top.is_identical(&top_of_module(module.clone()).unwrap()));
        assert!(top.is_identical(top_projection.range()));

        let socle = socle_of_module_inclusion(module.clone()).unwrap();
        assert_eq!(socle.source().dimension_vector(), &[1, 1]);
        assert_eq!(socle.range(), &module);
        assert_eq!(
            socle.vertex_matrices(),
            &[
                DenseMatrix::from_rows(vec![vec![q(0), q(1)]], Q).unwrap(),
                DenseMatrix::identity(1, Q),
            ]
        );
        assert_eq!(
            socle_of_module(module.clone()).unwrap().dimension_vector(),
            &[1, 1]
        );

        let zero_module = right_module_over_path_algebra(algebra.clone(), vec![0, 0], &[])
            .unwrap()
            .unwrap();
        let zero_socle = socle_of_module_inclusion(zero_module.clone()).unwrap();
        assert!(zero_socle.source().is_identical(&zero_module));
        assert!(zero_socle.range().is_identical(&zero_module));
        assert!(
            minimal_generating_set_of_module(zero_module.clone())
                .unwrap()
                .is_empty()
        );

        let target = right_module_over_path_algebra(
            algebra.clone(),
            vec![1, 2],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![q(1), q(1)]], Q).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let simple_one = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
            .unwrap()
            .unwrap();
        let simple_two = right_module_over_path_algebra(algebra, vec![0, 1], &[])
            .unwrap()
            .unwrap();

        let zero_trace = trace_of_module(&simple_one, &target).unwrap();
        assert_eq!(zero_trace.source().dimension_vector(), &[0, 0]);
        assert_eq!(zero_trace.range(), &target);
        assert_eq!(
            zero_trace.vertex_matrices(),
            &[DenseMatrix::zero(1, 1, Q), DenseMatrix::zero(1, 2, Q),]
        );
        assert_eq!(trace_of_module(&zero_module, &target).unwrap(), zero_trace);
        let zero_target_trace = trace_of_module(&simple_two, &zero_module).unwrap();
        assert_eq!(zero_target_trace.source().dimension_vector(), &[0, 0]);
        assert!(zero_target_trace.range().is_identical(&zero_module));

        let trace = trace_of_module(&simple_two, &target).unwrap();
        assert_eq!(trace.source().dimension_vector(), &[0, 2]);
        assert_eq!(
            trace.vertex_matrices(),
            &[DenseMatrix::zero(1, 1, Q), DenseMatrix::identity(2, Q)]
        );
        assert_eq!(
            right_fac_approximation(&simple_two, &target).unwrap(),
            trace
        );
        assert_eq!(
            right_fac_approximation_via_all_basis_images(&simple_two, &target).unwrap(),
            trace
        );

        let foreign = QuiverAlgebra::Path(
            PathAlgebra::new(Q, Quiver::from_indexed(1, vec![]).unwrap()).unwrap(),
        );
        let foreign_module = right_module_over_path_algebra(foreign, vec![1], &[])
            .unwrap()
            .unwrap();
        assert_eq!(
            trace_of_module(&foreign_module, &target).unwrap_err(),
            PathAlgebraModuleError::DifferentHomomorphismAlgebras
        );
    }

    #[test]
    fn binary_socle_and_trace_use_the_same_qpa_nullspaces() {
        let e = |value: i64| Z2.nth(Integer::from(value));
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, quiver).unwrap());
        let module = right_module_over_path_algebra(
            algebra.clone(),
            vec![2, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![e(1)], vec![e(0)]], Z2).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let socle = socle_of_module_inclusion(module).unwrap();
        assert_eq!(socle.source().dimension_vector(), &[1, 1]);
        assert_eq!(
            socle.vertex_matrices()[0],
            DenseMatrix::from_rows(vec![vec![e(0), e(1)]], Z2).unwrap()
        );

        let target = right_module_over_path_algebra(
            algebra.clone(),
            vec![1, 2],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![e(1), e(1)]], Z2).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let simple_two = right_module_over_path_algebra(algebra, vec![0, 1], &[])
            .unwrap()
            .unwrap();
        let trace = trace_of_module(&simple_two, &target).unwrap();
        assert_eq!(trace.source().dimension_vector(), &[0, 2]);
        assert_eq!(trace.vertex_matrices()[1], DenseMatrix::identity(2, Z2));
    }

    #[test]
    fn rational_radical_and_top_preserve_gap_quotient_lifts() {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
        let action = DenseMatrix::from_rows(vec![vec![q(1), q(1)], vec![q(0), q(0)]], Q).unwrap();
        let module = right_module_over_path_algebra(
            algebra,
            vec![2, 2],
            &[NamedArrowMatrix::new("a", action)],
        )
        .unwrap()
        .unwrap();

        let radical = radical_of_module_inclusion(module.clone()).unwrap();
        assert_eq!(radical.source().dimension_vector(), &[0, 1]);
        assert_eq!(radical.vertex_matrices()[0], DenseMatrix::zero(1, 2, Q));
        assert_eq!(
            radical.vertex_matrices()[1],
            DenseMatrix::from_rows(vec![vec![q(1), q(1)]], Q).unwrap()
        );

        let top = top_of_module_projection(module).unwrap();
        assert_eq!(top.range().dimension_vector(), &[2, 1]);
        assert_eq!(top.range().arrow_matrices()[0], DenseMatrix::zero(2, 1, Q));
        assert_eq!(top.vertex_matrices()[0], DenseMatrix::identity(2, Q));
        assert_eq!(
            top.vertex_matrices()[1],
            DenseMatrix::from_rows(vec![vec![q(1)], vec![q(-1)]], Q).unwrap()
        );

        let lifts = top
            .range()
            .basis_elements()
            .iter()
            .map(|basis| {
                top.preimage_representative(basis)
                    .unwrap()
                    .unwrap()
                    .external_representation()
                    .to_vec()
            })
            .collect::<Vec<_>>();
        assert_eq!(
            lifts,
            vec![
                vec![vec![q(1), q(0)], vec![q(0), q(0)]],
                vec![vec![q(0), q(1)], vec![q(0), q(0)]],
                vec![vec![q(0), q(0)], vec![q(1), q(0)]],
            ]
        );
    }

    #[test]
    fn characteristic_two_top_uses_the_same_qpa_source_method() {
        let e = |value: i64| Z2.nth(Integer::from(value));
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, quiver).unwrap());
        let action = DenseMatrix::from_rows(vec![vec![e(1), e(1)], vec![e(0), e(0)]], Z2).unwrap();
        let module = right_module_over_path_algebra(
            algebra,
            vec![2, 2],
            &[NamedArrowMatrix::new("a", action)],
        )
        .unwrap()
        .unwrap();

        let top = top_of_module_projection(module).unwrap();
        assert_eq!(top.range().dimension_vector(), &[2, 1]);
        assert_eq!(
            top.vertex_matrices()[1],
            DenseMatrix::from_rows(vec![vec![e(1)], vec![e(1)]], Z2).unwrap()
        );
    }

    #[test]
    fn zero_top_is_qpa_zero_endomorphism_of_the_identical_module() {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
        let module = right_module_over_path_algebra(algebra, vec![0, 0], &[])
            .unwrap()
            .unwrap();

        let top = top_of_module_projection(module.clone()).unwrap();
        assert!(top.source().is_identical(&module));
        assert!(top.range().is_identical(&module));
        assert_eq!(
            top.vertex_matrices(),
            &[DenseMatrix::zero(1, 1, Q), DenseMatrix::zero(1, 1, Q)]
        );
    }

    #[test]
    fn endomorphism_ideal_closure_is_two_sided_before_quotienting() {
        let quiver = Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
        let module = right_module_over_path_algebra(algebra, vec![2], &[])
            .unwrap()
            .unwrap();
        let endomorphisms = end_over_algebra(&module).unwrap();
        let matrix_unit =
            DenseMatrix::from_rows(vec![vec![q(0), q(1)], vec![q(0), q(0)]], Q).unwrap();
        let natural =
            natural_endomorphism_quotient_by_ideal(&endomorphisms, &[matrix_unit]).unwrap();

        assert_eq!(endomorphisms.dimension(), 4);
        assert_eq!(natural.kernel_basis_matrices().len(), 4);
        assert_eq!(natural.range().dimension(), 0);
        assert!(natural.range().one().unwrap().is_none());
        assert!(
            endomorphisms
                .basis_matrices()
                .iter()
                .all(|matrix| natural.image(matrix).unwrap().is_zero())
        );
    }
}
