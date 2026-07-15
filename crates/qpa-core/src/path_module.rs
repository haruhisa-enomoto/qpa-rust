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

/// A QPA homomorphism between right path-algebra matrix modules.
///
/// The matrix at vertex `i` acts on row vectors from the source vertex space
/// to the range vertex space.  Physical one-row/one-column zero-space
/// stand-ins are retained exactly as in `RightModuleHomOverAlgebra`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PathAlgebraModuleHomomorphism<F: Field> {
    source: PathAlgebraModule<F>,
    range: PathAlgebraModule<F>,
    vertex_matrices: Vec<DenseMatrix<F>>,
}

impl<F: Field> PathAlgebraModuleHomomorphism<F> {
    /// Return QPA's `Source` module.
    #[must_use]
    pub const fn source(&self) -> &PathAlgebraModule<F> {
        &self.source
    }

    /// Return QPA's `Range` module.
    #[must_use]
    pub const fn range(&self) -> &PathAlgebraModule<F> {
        &self.range
    }

    /// Translate `MatricesOfPathAlgebraMatModuleHomomorphism`.
    #[must_use]
    pub fn vertex_matrices(&self) -> &[DenseMatrix<F>] {
        &self.vertex_matrices
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
        if !self.source.is_identical(element.module()) {
            return Err(PathAlgebraModuleError::ForeignMapSourceElement);
        }
        let mut image = self.range.zero_element();
        for (vertex, matrix) in self.vertex_matrices.iter().enumerate() {
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
        if !self.range.is_identical(element.module()) {
            return Err(PathAlgebraModuleError::ForeignMapRangeElement);
        }
        let mut preimage = self.source.zero_element();
        for (vertex, matrix) in self.vertex_matrices.iter().enumerate() {
            let Some(row) = matrix.solution_row(&element.coordinates[vertex])? else {
                return Ok(None);
            };
            preimage.coordinates[vertex] = row;
        }
        Ok(Some(preimage))
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

    Ok(PathAlgebraModuleHomomorphism {
        source,
        range,
        vertex_matrices,
    })
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
    right_module_hom_over_algebra(range.clone(), cokernel, projection_matrices)
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
pub fn top_of_module_projection<F: Field>(
    module: PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    if module.dimension() == 0 {
        zero_module_mapping(module.clone(), module)
    } else {
        cokernel_projection(&radical_of_module_inclusion(module)?)
    }
}

/// Translate QPA's `TopOfModule(M)` convenience method.
///
/// # Errors
///
/// Propagates [`top_of_module_projection`].
pub fn top_of_module<F: Field>(
    module: PathAlgebraModule<F>,
) -> Result<PathAlgebraModule<F>, PathAlgebraModuleError> {
    Ok(top_of_module_projection(module)?.range().clone())
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

    Ok(Some(PathAlgebraModule(Arc::new(PathAlgebraModuleData {
        acting_algebra: algebra,
        dimension_vector,
        dimension,
        arrow_matrices,
        retained_tensor,
        dual_weak: RwLock::new(None),
        dual_strong: OnceLock::new(),
    }))))
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
    use numerica::domains::{
        Ring,
        finite_field::{Z2, Zp},
        integer::Integer,
        rational::{Q, Rational},
    };

    use crate::{
        ArrowSpec, DenseMatrix, PathAlgebra, Quiver, QuiverAlgebra, QuotientPathAlgebra,
        VertexIndex,
    };

    use super::{
        NamedArrowMatrix, PathAlgebraModuleError, PathAlgebraModuleMatrix,
        radical_of_module_inclusion, right_module_over_path_algebra, top_of_module_projection,
    };

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn q(value: i64) -> Rational {
        Q.nth(value.into())
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
}
