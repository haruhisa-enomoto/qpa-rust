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

/// One arrow action accepted by QPA's dimension-inferring module constructor.
///
/// GAP matrices cannot retain a zero row or column count, so QPA uses the
/// integer pair `[source_dimension, target_dimension]` when either endpoint
/// is the zero space. The enum keeps that input distinction explicit.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum PathAlgebraModuleArrowAction<F: Field> {
    /// An ordinary positive-by-positive arrow matrix.
    Matrix(DenseMatrix<F>),
    /// QPA's special zero-space dimension pair.
    ZeroSpaceDimensions {
        /// Logical dimension at the arrow source.
        source: usize,
        /// Logical dimension at the arrow target.
        target: usize,
    },
}

impl<F: Field> PathAlgebraModuleArrowAction<F> {
    /// Wrap an ordinary arrow matrix.
    #[must_use]
    pub const fn matrix(matrix: DenseMatrix<F>) -> Self {
        Self::Matrix(matrix)
    }

    /// Construct QPA's `[source_dimension, target_dimension]` notation.
    #[must_use]
    pub const fn zero_space_dimensions(source: usize, target: usize) -> Self {
        Self::ZeroSpaceDimensions { source, target }
    }
}

/// A string-labeled entry in QPA's two-argument module constructor.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct NamedPathAlgebraModuleArrowAction<F: Field> {
    label: Box<str>,
    action: PathAlgebraModuleArrowAction<F>,
}

impl<F: Field> NamedPathAlgebraModuleArrowAction<F> {
    /// Pair an arrow label with its proposed action.
    #[must_use]
    pub fn new(label: impl Into<Box<str>>, action: PathAlgebraModuleArrowAction<F>) -> Self {
        Self {
            label: label.into(),
            action,
        }
    }

    /// Return the entered arrow label.
    #[must_use]
    pub fn label(&self) -> &str {
        &self.label
    }

    /// Return the entered action.
    #[must_use]
    pub const fn action(&self) -> &PathAlgebraModuleArrowAction<F> {
        &self.action
    }
}

/// An arrow-object entry in QPA's two-argument module constructor.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct IndexedPathAlgebraModuleArrowAction<F: Field> {
    arrow: ArrowIndex,
    action: PathAlgebraModuleArrowAction<F>,
}

impl<F: Field> IndexedPathAlgebraModuleArrowAction<F> {
    /// Pair a one-based QPA arrow index with its proposed action.
    #[must_use]
    pub const fn new(arrow: ArrowIndex, action: PathAlgebraModuleArrowAction<F>) -> Self {
        Self { arrow, action }
    }

    /// Return the entered one-based arrow index.
    #[must_use]
    pub const fn arrow(&self) -> ArrowIndex {
        self.arrow
    }

    /// Return the entered action.
    #[must_use]
    pub const fn action(&self) -> &PathAlgebraModuleArrowAction<F> {
        &self.action
    }
}

/// The three homogeneous list formats dispatched by QPA's two-argument
/// `RightModuleOverPathAlgebra` methods.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum PathAlgebraModuleArrowActions<F: Field> {
    /// Entries of the form `["arrow_name", action]`.
    Named(Vec<NamedPathAlgebraModuleArrowAction<F>>),
    /// Bare matrices in quiver-arrow order.
    Positional(Vec<DenseMatrix<F>>),
    /// Entries of the form `[arrow_object, action]`.
    Indexed(Vec<IndexedPathAlgebraModuleArrowAction<F>>),
}

/// Translate QPA's `IsInFullMatrixRing(M, R)` for a typed Numerica matrix.
///
/// QPA flattens the matrix and checks every entry for membership in `R`.
/// [`DenseMatrix`] already stores entries of one coefficient type, so the
/// remaining runtime check is equality of its retained field and `ring`.
#[must_use]
pub fn is_in_full_matrix_ring<F: Field>(matrix: &DenseMatrix<F>, ring: &F) -> bool {
    matrix.field() == ring
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

/// QPA's dynamically shaped result from `RadicalSeries` and `SocleSeries`.
///
/// For a zero module QPA returns its dimension vector directly. For a
/// nonzero module it returns a list of layer dimension vectors. The enum
/// preserves that distinction without an untyped nested list.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ModuleDimensionSeries {
    /// The direct dimension-vector return used only for a zero module.
    ZeroDimensionVector(Vec<usize>),
    /// Ordered radical or socle layer dimension vectors for a nonzero module.
    Layers(Vec<Vec<usize>>),
}

#[derive(Clone, Debug)]
pub(crate) enum ModuleDecompositionMultiplicitiesCache<F: Field> {
    Fail,
    Value {
        basic_summands: Vec<PathAlgebraModule<F>>,
        multiplicities: Vec<usize>,
    },
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
    transpose_of_module: OnceLock<PathAlgebraModule<F>>,
    d_tr: OnceLock<PathAlgebraModule<F>>,
    tr_d: OnceLock<PathAlgebraModule<F>>,
    star_of_module: OnceLock<PathAlgebraModule<F>>,
    nakayama_functor: OnceLock<PathAlgebraModule<F>>,
    opposite_nakayama_functor: OnceLock<PathAlgebraModule<F>>,
    endomorphism_algebra: OnceLock<ModuleEndomorphismAlgebra<F>>,
    decomposition: OnceLock<Vec<PathAlgebraModule<F>>>,
    decomposition_inclusions: WeakPathAlgebraModuleHomomorphismListCache<F>,
    decomposition_multiplicities: OnceLock<ModuleDecompositionMultiplicitiesCache<F>>,
    top_module: OnceLock<PathAlgebraModule<F>>,
    top_projection_weak: RwLock<Option<WeakPathAlgebraModuleHomomorphism<F>>>,
    direct_sum: OnceLock<DirectSumMetadata<F>>,
    minimal_generators: OnceLock<Vec<Vec<Vec<F::Element>>>>,
    is_projective: OnceLock<bool>,
    is_injective: OnceLock<bool>,
    projective_dimension: OnceLock<usize>,
    injective_dimension: OnceLock<usize>,
    faithful_dimension: OnceLock<Option<usize>>,
    is_indecomposable: OnceLock<bool>,
    is_simple_qpa: OnceLock<bool>,
    is_semisimple: OnceLock<bool>,
    is_rigid: OnceLock<bool>,
    is_exceptional: OnceLock<bool>,
    is_tilting: OnceLock<bool>,
    is_cotilting: OnceLock<bool>,
    is_tau_rigid: OnceLock<bool>,
    projective_cover_weak: RwLock<Option<WeakPathAlgebraModuleHomomorphism<F>>>,
    injective_envelope_weak: RwLock<Option<WeakPathAlgebraModuleHomomorphism<F>>>,
    identity_to_double_star_weak: RwLock<Option<WeakPathAlgebraModuleHomomorphism<F>>>,
    almost_split_sequence: OnceLock<Option<PathAlgebraModuleHomomorphismPair<F>>>,
    irreducible_ending: WeakPathAlgebraModuleHomomorphismListCache<F>,
    irreducible_starting: WeakPathAlgebraModuleHomomorphismListCache<F>,
    first_syzygy: OnceLock<PathAlgebraModule<F>>,
    annihilator: OnceLock<Vec<QuiverAlgebraElement<F>>>,
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
            transpose_of_module: OnceLock::new(),
            d_tr: OnceLock::new(),
            tr_d: OnceLock::new(),
            star_of_module: OnceLock::new(),
            nakayama_functor: OnceLock::new(),
            opposite_nakayama_functor: OnceLock::new(),
            endomorphism_algebra: OnceLock::new(),
            decomposition: OnceLock::new(),
            decomposition_inclusions: WeakPathAlgebraModuleHomomorphismListCache::new(),
            decomposition_multiplicities: OnceLock::new(),
            top_module: OnceLock::new(),
            top_projection_weak: RwLock::new(None),
            direct_sum: OnceLock::new(),
            minimal_generators: OnceLock::new(),
            is_projective: OnceLock::new(),
            is_injective: OnceLock::new(),
            projective_dimension: OnceLock::new(),
            injective_dimension: OnceLock::new(),
            faithful_dimension: OnceLock::new(),
            is_indecomposable: OnceLock::new(),
            is_simple_qpa: OnceLock::new(),
            is_semisimple: OnceLock::new(),
            is_rigid: OnceLock::new(),
            is_exceptional: OnceLock::new(),
            is_tilting: OnceLock::new(),
            is_cotilting: OnceLock::new(),
            is_tau_rigid: OnceLock::new(),
            projective_cover_weak: RwLock::new(None),
            injective_envelope_weak: RwLock::new(None),
            identity_to_double_star_weak: RwLock::new(None),
            almost_split_sequence: OnceLock::new(),
            irreducible_ending: WeakPathAlgebraModuleHomomorphismListCache::new(),
            irreducible_starting: WeakPathAlgebraModuleHomomorphismListCache::new(),
            first_syzygy: OnceLock::new(),
            annihilator: OnceLock::new(),
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

struct PathAlgebraModuleSubmoduleData<F: Field> {
    parent: PathAlgebraModule<F>,
    generators: Box<[PathAlgebraModuleElement<F>]>,
}

/// QPA's `SubAlgebraModule` object inside a path-algebra matrix module.
///
/// QPA stores the nonzero vertex-uniform components of the entered generators
/// and retains the ambient module as `Parent`. Clones preserve object
/// identity.
#[derive(Clone)]
pub struct PathAlgebraModuleSubmodule<F: Field>(Arc<PathAlgebraModuleSubmoduleData<F>>);

struct PathAlgebraModuleBasisMappingData<F: Field> {
    source: PathAlgebraModule<F>,
    range: PathAlgebraModule<F>,
    basis_images: Box<[PathAlgebraModuleElement<F>]>,
    surjective: bool,
}

/// A QPA `LeftModuleHomomorphismByImagesNC` represented by ordered basis
/// images.
///
/// This is deliberately distinct from [`PathAlgebraModuleHomomorphism`]: the
/// mapping returned by QPA's `NaturalHomomorphismBySubAlgebraModule` can have
/// basis images that do not preserve vertex support because of a source-list
/// aliasing bug.
#[derive(Clone)]
pub struct PathAlgebraModuleBasisMapping<F: Field>(Arc<PathAlgebraModuleBasisMappingData<F>>);

/// An object accepted by QPA's generic `Hom(R, M, N)` method.
///
/// QPA uses the same `IsPathAlgebraMatModule` filter for ambient matrix
/// modules and for `SubAlgebraModule` objects carrying a `Parent` attribute.
/// The enum keeps that dispatch distinction explicit in Rust.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum PathAlgebraMatModuleObject<F: Field> {
    /// A matrix module without a `Parent` attribute.
    Module(PathAlgebraModule<F>),
    /// A `SubAlgebraModule` with its ambient module retained as `Parent`.
    Submodule(PathAlgebraModuleSubmodule<F>),
}

/// A finite-dimensional based right module over a quiver algebra.
///
/// This is the Rust boundary for QPA's generic GAP
/// `IsRightAlgebraModuleElementCollection` input. Each matrix describes right
/// action on row coordinates in one common entered basis. Vertex actions are
/// ordered by quiver vertex and arrow actions by quiver arrow.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BasedRightAlgebraModule<F: Field> {
    algebra: QuiverAlgebra<F>,
    dimension: usize,
    vertex_actions: Vec<DenseMatrix<F>>,
    arrow_actions: Vec<DenseMatrix<F>>,
}

struct GenericPathAlgebraModuleHomomorphismData<F: Field> {
    source: PathAlgebraMatModuleObject<F>,
    range: PathAlgebraMatModuleObject<F>,
    source_basis: Box<[PathAlgebraModuleElement<F>]>,
    range_basis: Box<[PathAlgebraModuleElement<F>]>,
    matrix: DenseMatrix<F>,
}

/// One generator returned by QPA's generic `Hom(R, M, N)` method.
///
/// QPA represents these as `LeftModuleHomomorphismByMatrix` mappings and
/// installs `IsAlgebraModuleHomomorphism`. The matrix is expressed in the
/// ordinary module basis or, for a submodule, its ordered canonical basis.
#[derive(Clone)]
pub struct GenericPathAlgebraModuleHomomorphism<F: Field>(
    Arc<GenericPathAlgebraModuleHomomorphismData<F>>,
);

/// QPA's generic `Hom(R, M, N)` vector space.
///
/// The `generators` retain the exact `HOPF_MatrixModuleHom` nullspace order.
/// GAP's `GeneratorsOfLeftModule` exposes this list even in cases where its
/// later `BasisVectors(Basis(...))` dispatch fails.
#[derive(Clone, Debug)]
pub struct GenericPathAlgebraModuleHomSpace<F: Field> {
    field: F,
    source: PathAlgebraMatModuleObject<F>,
    range: PathAlgebraMatModuleObject<F>,
    generators: Vec<GenericPathAlgebraModuleHomomorphism<F>>,
    flattened_generators: DenseMatrix<F>,
}

/// QPA's generic `End(R, M)` algebra-with-one wrapper.
///
/// QPA passes the ordered `HOPF_MatrixModuleHom` list directly to
/// `AlgebraWithOne(..., "basis")`; no semi-echelon rebasing is performed.
#[derive(Clone, Debug)]
pub struct GenericPathAlgebraModuleEndomorphismAlgebra<F: Field> {
    field: F,
    module: PathAlgebraMatModuleObject<F>,
    generators: Vec<GenericPathAlgebraModuleHomomorphism<F>>,
    flattened_generators: DenseMatrix<F>,
    one: GenericPathAlgebraModuleHomomorphism<F>,
}

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

    /// Translate QPA's `IsZero(M)` method for path-algebra matrix modules.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.dimension() == 0
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

    /// Translate QPA's raw `PathModuleElem(family, external_representation)`
    /// constructor.
    ///
    /// Unlike [`Self::element_from_vertex_coordinates`], this accepts the
    /// physical external representation. QPA deliberately permits any row
    /// length at a zero-dimensional vertex; retaining that behavior is useful
    /// for exact `ExtRepOfObj` round trips, including an empty row.
    ///
    /// # Errors
    ///
    /// Rejects a wrong vertex count or a nonzero-dimensional row with the
    /// wrong length.
    pub fn element_from_external_representation(
        &self,
        coordinates: Vec<Vec<F::Element>>,
    ) -> Result<PathAlgebraModuleElement<F>, PathAlgebraModuleError> {
        if coordinates.len() != self.dimension_vector().len() {
            return Err(PathAlgebraModuleError::ElementVertexCount {
                actual: coordinates.len(),
                expected: self.dimension_vector().len(),
            });
        }
        for (vertex, (row, &dimension)) in
            coordinates.iter().zip(self.dimension_vector()).enumerate()
        {
            if dimension != 0 && row.len() != dimension {
                return Err(PathAlgebraModuleError::ElementCoordinateCount {
                    vertex: vertex + 1,
                    actual: row.len(),
                    expected: dimension,
                });
            }
        }
        Ok(PathAlgebraModuleElement {
            module: self.clone(),
            coordinates,
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

    /// Return QPA's cached `IsRigidModule` value, if installed.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_is_rigid_module(&self) -> Option<bool> {
        self.0.is_rigid.get().copied()
    }

    /// Install QPA's `IsRigidModule` property value.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_is_rigid_module(&self, value: bool) -> bool {
        *self.0.is_rigid.get_or_init(|| value)
    }

    /// Return QPA's cached `IsExceptionalModule` value, if installed.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_is_exceptional_module(&self) -> Option<bool> {
        self.0.is_exceptional.get().copied()
    }

    /// Install QPA's `IsExceptionalModule` property value.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_is_exceptional_module(&self, value: bool) -> bool {
        *self.0.is_exceptional.get_or_init(|| value)
    }

    /// Return QPA's installed `IsTiltingModule` attribute, if bound.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_is_tilting_module(&self) -> Option<bool> {
        self.0.is_tilting.get().copied()
    }

    /// Install QPA's successful `IsTiltingModule` attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_is_tilting_module(&self, value: bool) -> bool {
        *self.0.is_tilting.get_or_init(|| value)
    }

    /// Return QPA's installed `IsCotiltingModule` attribute, if bound.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_is_cotilting_module(&self) -> Option<bool> {
        self.0.is_cotilting.get().copied()
    }

    /// Install QPA's `IsCotiltingModule` attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_is_cotilting_module(&self, value: bool) -> bool {
        *self.0.is_cotilting.get_or_init(|| value)
    }

    /// Return QPA's cached `TransposeOfModule` attribute, if installed.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_transpose_of_module(&self) -> Option<Self> {
        self.0.transpose_of_module.get().cloned()
    }

    /// Install QPA's `TransposeOfModule` attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_transpose_of_module(&self, value: &Self) -> Self {
        self.0
            .transpose_of_module
            .get_or_init(|| value.clone())
            .clone()
    }

    /// Return QPA's cached module `DTr` attribute, if installed.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_d_tr(&self) -> Option<Self> {
        self.0.d_tr.get().cloned()
    }

    /// Install QPA's module `DTr` attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_d_tr(&self, value: &Self) -> Self {
        self.0.d_tr.get_or_init(|| value.clone()).clone()
    }

    /// Return QPA's cached module `TrD` attribute, if installed.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_tr_d(&self) -> Option<Self> {
        self.0.tr_d.get().cloned()
    }

    /// Install QPA's module `TrD` attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_tr_d(&self, value: &Self) -> Self {
        self.0.tr_d.get_or_init(|| value.clone()).clone()
    }

    /// Return QPA's cached `IsTauRigidModule` property, if installed.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_is_tau_rigid_module(&self) -> Option<bool> {
        self.0.is_tau_rigid.get().copied()
    }

    /// Install QPA's `IsTauRigidModule` property.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_is_tau_rigid_module(&self, value: bool) -> bool {
        *self.0.is_tau_rigid.get_or_init(|| value)
    }

    /// Return QPA's cached `StarOfModule` attribute, if installed.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_star_of_module(&self) -> Option<Self> {
        self.0.star_of_module.get().cloned()
    }

    /// Install QPA's `StarOfModule` attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_star_of_module(&self, value: &Self) -> Self {
        self.0.star_of_module.get_or_init(|| value.clone()).clone()
    }

    /// Return QPA's cached `NakayamaFunctorOfModule` attribute, if installed.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_nakayama_functor_of_module(&self) -> Option<Self> {
        self.0.nakayama_functor.get().cloned()
    }

    /// Install QPA's `NakayamaFunctorOfModule` attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_nakayama_functor_of_module(&self, value: &Self) -> Self {
        self.0
            .nakayama_functor
            .get_or_init(|| value.clone())
            .clone()
    }

    /// Return QPA's cached `OppositeNakayamaFunctorOfModule` attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_opposite_nakayama_functor_of_module(&self) -> Option<Self> {
        self.0.opposite_nakayama_functor.get().cloned()
    }

    /// Install QPA's `OppositeNakayamaFunctorOfModule` attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_opposite_nakayama_functor_of_module(&self, value: &Self) -> Self {
        self.0
            .opposite_nakayama_functor
            .get_or_init(|| value.clone())
            .clone()
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

    /// Return QPA's installed `ProjDimension` attribute, if bound.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_projective_dimension(&self) -> Option<usize> {
        self.0.projective_dimension.get().copied()
    }

    /// Install QPA's finite `ProjDimension` attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_projective_dimension(&self, value: usize) -> usize {
        *self.0.projective_dimension.get_or_init(|| value)
    }

    /// Return QPA's installed `InjDimension` attribute, if bound.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_injective_dimension(&self) -> Option<usize> {
        self.0.injective_dimension.get().copied()
    }

    /// Install QPA's finite `InjDimension` attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_injective_dimension(&self, value: usize) -> usize {
        *self.0.injective_dimension.get_or_init(|| value)
    }

    /// Return QPA's installed `FaithfulDimension` attribute, if bound.
    ///
    /// The inner `None` is QPA's `infinity` value.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_faithful_dimension(&self) -> Option<Option<usize>> {
        self.0.faithful_dimension.get().copied()
    }

    /// Install QPA's `FaithfulDimension` attribute.
    ///
    /// `None` represents QPA's `infinity` value.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_faithful_dimension(&self, value: Option<usize>) -> Option<usize> {
        *self.0.faithful_dimension.get_or_init(|| value)
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

    /// Return an already-live cached `FromIdentityToDoubleStarHomomorphism` attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_identity_to_double_star(&self) -> Option<PathAlgebraModuleHomomorphism<F>> {
        self.0
            .identity_to_double_star_weak
            .read()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .as_ref()
            .and_then(WeakPathAlgebraModuleHomomorphism::upgrade)
    }

    /// Install `FromIdentityToDoubleStarHomomorphism` without an `Arc` cycle.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_identity_to_double_star(
        &self,
        homomorphism: &PathAlgebraModuleHomomorphism<F>,
    ) -> PathAlgebraModuleHomomorphism<F> {
        cache_weak_homomorphism(&self.0.identity_to_double_star_weak, homomorphism)
    }

    /// Return QPA's cached `AlmostSplitSequence` attribute state.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_almost_split_sequence(
        &self,
    ) -> CachedAttribute<PathAlgebraModuleHomomorphismPair<F>> {
        match self.0.almost_split_sequence.get() {
            None => CachedAttribute::Unbound,
            Some(None) => CachedAttribute::Fail,
            Some(Some(sequence)) => CachedAttribute::Value(sequence.clone()),
        }
    }

    /// Install QPA's `AlmostSplitSequence` value or `fail` sentinel.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_almost_split_sequence(
        &self,
        sequence: Option<&PathAlgebraModuleHomomorphismPair<F>>,
    ) -> Option<PathAlgebraModuleHomomorphismPair<F>> {
        self.0
            .almost_split_sequence
            .get_or_init(|| sequence.cloned())
            .clone()
    }

    /// Return already-live cached `IrreducibleMorphismsEndingIn` maps.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_irreducible_morphisms_ending_in(
        &self,
    ) -> Option<Vec<PathAlgebraModuleHomomorphism<F>>> {
        self.0.irreducible_ending.cached()
    }

    /// Install `IrreducibleMorphismsEndingIn` without retaining an `Arc` cycle.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_irreducible_morphisms_ending_in(
        &self,
        homomorphisms: &[PathAlgebraModuleHomomorphism<F>],
    ) -> Vec<PathAlgebraModuleHomomorphism<F>> {
        self.0.irreducible_ending.merge(homomorphisms)
    }

    /// Return already-live cached `IrreducibleMorphismsStartingIn` maps.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_irreducible_morphisms_starting_in(
        &self,
    ) -> Option<Vec<PathAlgebraModuleHomomorphism<F>>> {
        self.0.irreducible_starting.cached()
    }

    /// Install `IrreducibleMorphismsStartingIn` without retaining an `Arc` cycle.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_irreducible_morphisms_starting_in(
        &self,
        homomorphisms: &[PathAlgebraModuleHomomorphism<F>],
    ) -> Vec<PathAlgebraModuleHomomorphism<F>> {
        self.0.irreducible_starting.merge(homomorphisms)
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

    /// Return QPA's cached `AnnihilatorOfModule` basis, if installed.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_annihilator_of_module(&self) -> Option<Vec<QuiverAlgebraElement<F>>> {
        self.0.annihilator.get().cloned()
    }

    /// Install and return QPA's `AnnihilatorOfModule` basis attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_annihilator_of_module(
        &self,
        basis: &[QuiverAlgebraElement<F>],
    ) -> Vec<QuiverAlgebraElement<F>> {
        self.0.annihilator.get_or_init(|| basis.to_vec()).clone()
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

    pub(crate) fn cached_decomposition(&self) -> Option<Vec<Self>> {
        self.0.decomposition.get().cloned()
    }

    pub(crate) fn cache_decomposition(&self, decomposition: &[Self]) -> Vec<Self> {
        self.0
            .decomposition
            .get_or_init(|| decomposition.to_vec())
            .clone()
    }

    pub(crate) fn cached_decomposition_inclusions(
        &self,
    ) -> Option<Vec<PathAlgebraModuleHomomorphism<F>>> {
        self.0.decomposition_inclusions.cached()
    }

    pub(crate) fn cache_decomposition_inclusions(
        &self,
        inclusions: &[PathAlgebraModuleHomomorphism<F>],
    ) -> Vec<PathAlgebraModuleHomomorphism<F>> {
        self.0.decomposition_inclusions.merge(inclusions)
    }

    pub(crate) fn cached_decomposition_multiplicities(
        &self,
    ) -> Option<ModuleDecompositionMultiplicitiesCache<F>> {
        self.0.decomposition_multiplicities.get().cloned()
    }

    pub(crate) fn cache_decomposition_multiplicities(
        &self,
        decomposition: ModuleDecompositionMultiplicitiesCache<F>,
    ) -> ModuleDecompositionMultiplicitiesCache<F> {
        self.0
            .decomposition_multiplicities
            .get_or_init(|| decomposition)
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

    /// Return QPA's `PrintObj` text for the raw path-module element.
    ///
    /// The callback supplies GAP's scalar `String` representation.
    #[must_use]
    pub fn qpa_print_string_with(
        &self,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> String {
        qpa_list(
            self.coordinates
                .iter()
                .map(|row| qpa_list(row.iter().map(coefficient)))
                .collect::<Vec<_>>(),
        )
    }

    /// Add elements from the identical QPA element family.
    ///
    /// GAP vector addition retains unmatched trailing entries. This matters
    /// only for the arbitrary physical rows accepted at zero-dimensional
    /// vertices, and is preserved in both operand orders.
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
                (0..left.len().max(right.len()))
                    .map(|position| match (left.get(position), right.get(position)) {
                        (Some(left), Some(right)) => field.add(left, right),
                        (Some(left), None) => left.clone(),
                        (None, Some(right)) => right.clone(),
                        (None, None) => unreachable!("position is below the maximum row length"),
                    })
                    .collect()
            })
            .collect();
        Ok(Self {
            module: self.module.clone(),
            coordinates,
        })
    }

    /// Subtract elements from the identical QPA element family.
    ///
    /// # Errors
    ///
    /// Rejects different module families.
    pub fn subtract(&self, right: &Self) -> Result<Self, PathAlgebraModuleError> {
        self.add(&right.additive_inverse())
    }

    /// Return QPA's `AdditiveInverseOp` result.
    #[must_use]
    pub fn additive_inverse(&self) -> Self {
        let field = self.module.right_acting_algebra().field();
        Self {
            module: self.module.clone(),
            coordinates: self
                .coordinates
                .iter()
                .map(|row| row.iter().map(|entry| field.neg(entry)).collect())
                .collect(),
        }
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

    /// Translate QPA's scalar-times-module-element method.
    #[must_use]
    pub fn scale_left(&self, scalar: &F::Element) -> Self {
        let field = self.module.right_acting_algebra().field();
        Self {
            module: self.module.clone(),
            coordinates: self
                .coordinates
                .iter()
                .map(|row| row.iter().map(|entry| field.mul(scalar, entry)).collect())
                .collect(),
        }
    }

    /// Translate QPA's module-element-times-scalar method.
    #[must_use]
    pub fn scale_right(&self, scalar: &F::Element) -> Self {
        self.scale(scalar)
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

/// A QPA free left-module space whose elements are raw path-module elements.
///
/// This models the `LeftModuleByGenerators` objects consumed by QPA's
/// `NewBasis`, `BasisNC`, `Basis`, and `BasisOfDomain` methods while retaining
/// the identical ambient path-module element family.
#[derive(Clone, Debug)]
pub struct PathAlgebraModuleElementSpace<F: Field> {
    module: PathAlgebraModule<F>,
    generators: Vec<PathAlgebraModuleElement<F>>,
}

impl<F: Field> PathAlgebraModuleElementSpace<F> {
    /// Construct QPA's raw element space from entered generators.
    ///
    /// # Errors
    ///
    /// Rejects a generator from another path-module element family.
    pub fn by_generators(
        module: PathAlgebraModule<F>,
        generators: Vec<PathAlgebraModuleElement<F>>,
    ) -> Result<Self, PathAlgebraModuleError> {
        if generators
            .iter()
            .any(|generator| !module.is_identical(generator.module()))
        {
            return Err(PathAlgebraModuleError::ForeignElementSpaceGenerator);
        }
        Ok(Self { module, generators })
    }

    /// Return the identical ambient path module.
    #[must_use]
    pub const fn module(&self) -> &PathAlgebraModule<F> {
        &self.module
    }

    /// Return `GeneratorsOfLeftModule(V)` in entered order.
    #[must_use]
    pub fn generators(&self) -> &[PathAlgebraModuleElement<F>] {
        &self.generators
    }

    /// Translate QPA's structural `NewBasis(V, vectors)` operation.
    ///
    /// This intentionally does not attach the delegate flattened-vector
    /// basis used by `Coefficients`.
    ///
    /// # Errors
    ///
    /// Rejects vectors outside the ambient element family.
    pub fn new_basis(
        &self,
        vectors: Vec<PathAlgebraModuleElement<F>>,
    ) -> Result<PathAlgebraModuleElementBasis<F>, PathAlgebraModuleError> {
        self.ensure_ambient_family(&vectors)?;
        Ok(PathAlgebraModuleElementBasis {
            space: self.clone(),
            vectors,
            delegate_basis: None,
            canonical: false,
            empty: false,
        })
    }

    /// Translate QPA's unchecked `BasisNC(V, vectors)` method.
    ///
    /// Zero vectors are removed in entered order before the flattened
    /// delegate basis is attached.
    ///
    /// # Errors
    ///
    /// Rejects foreign vectors or ragged physical external representations.
    pub fn basis_nc(
        &self,
        vectors: Vec<PathAlgebraModuleElement<F>>,
    ) -> Result<PathAlgebraModuleElementBasis<F>, PathAlgebraModuleError> {
        self.ensure_ambient_family(&vectors)?;
        let vectors = vectors
            .into_iter()
            .filter(|vector| !vector.is_zero())
            .collect::<Vec<_>>();
        let delegate_basis = element_basis_matrix(&self.module, &vectors)?;
        Ok(PathAlgebraModuleElementBasis {
            space: self.clone(),
            vectors,
            delegate_basis: Some(delegate_basis),
            canonical: false,
            empty: false,
        })
    }

    /// Translate QPA's checked `Basis(V, vectors)` overload.
    ///
    /// `Ok(None)` is QPA's `fail` when an entered vector lies outside `V` or
    /// the nonzero entered vectors are linearly dependent.
    ///
    /// # Errors
    ///
    /// Rejects foreign families and propagates exact flattened-row solving.
    pub fn basis_from_vectors(
        &self,
        vectors: Vec<PathAlgebraModuleElement<F>>,
    ) -> Result<Option<PathAlgebraModuleElementBasis<F>>, PathAlgebraModuleError> {
        self.ensure_ambient_family(&vectors)?;
        for vector in &vectors {
            if !self.contains(vector)? {
                return Ok(None);
            }
        }
        let basis = self.basis_nc(vectors)?;
        let Some(delegate) = basis.delegate_basis.as_ref() else {
            return Err(PathAlgebraModuleError::MissingElementBasisDelegate);
        };
        if delegate.canonical_row_space_basis().nrows() != basis.vectors.len() {
            return Ok(None);
        }
        Ok(Some(basis))
    }

    /// Translate QPA's `BasisOfDomain(V)` method.
    ///
    /// # Errors
    ///
    /// Propagates unchecked-basis construction from the stored generators.
    pub fn basis_of_domain(
        &self,
    ) -> Result<PathAlgebraModuleElementBasis<F>, PathAlgebraModuleError> {
        self.basis_nc(self.generators.clone())
    }

    /// Translate QPA's default `Basis(V)` method.
    ///
    /// # Errors
    ///
    /// Propagates unchecked-basis construction from the stored generators.
    pub fn basis(&self) -> Result<PathAlgebraModuleElementBasis<F>, PathAlgebraModuleError> {
        self.basis_of_domain()
    }

    /// Translate QPA's membership test `element in V`.
    ///
    /// QPA computes `Coefficients(BasisOfDomain(V), element)` and returns
    /// whether that operation succeeds. An element from another ambient
    /// module family is therefore not a member.
    ///
    /// # Errors
    ///
    /// Propagates exact flattened-row basis construction and solving.
    pub fn contains(
        &self,
        element: &PathAlgebraModuleElement<F>,
    ) -> Result<bool, PathAlgebraModuleError> {
        if !self.module.is_identical(element.module()) {
            return Ok(false);
        }
        if self.generators.is_empty() {
            return Ok(element.is_zero());
        }
        let rows = element_basis_matrix(&self.module, &self.generators)?;
        Ok(rows
            .solution_row(&element.flat_external_representation())?
            .is_some())
    }

    fn ensure_ambient_family(
        &self,
        vectors: &[PathAlgebraModuleElement<F>],
    ) -> Result<(), PathAlgebraModuleError> {
        if vectors
            .iter()
            .any(|vector| !self.module.is_identical(vector.module()))
        {
            Err(PathAlgebraModuleError::ForeignElementBasisVector)
        } else {
            Ok(())
        }
    }
}

/// A QPA basis of a raw path-module element space.
#[derive(Clone, Debug)]
pub struct PathAlgebraModuleElementBasis<F: Field> {
    space: PathAlgebraModuleElementSpace<F>,
    vectors: Vec<PathAlgebraModuleElement<F>>,
    delegate_basis: Option<DenseMatrix<F>>,
    canonical: bool,
    empty: bool,
}

impl<F: Field> PathAlgebraModuleElementBasis<F> {
    /// Return QPA's `UnderlyingLeftModule` value.
    #[must_use]
    pub const fn underlying_space(&self) -> &PathAlgebraModuleElementSpace<F> {
        &self.space
    }

    /// Return QPA's ordered `BasisVectors` value.
    #[must_use]
    pub fn basis_vectors(&self) -> &[PathAlgebraModuleElement<F>] {
        &self.vectors
    }

    /// Return QPA's installed `IsCanonicalBasis` property.
    #[must_use]
    pub const fn is_canonical_basis(&self) -> bool {
        self.canonical
    }

    /// Return QPA's explicitly stored `IsEmpty` property.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.empty
    }

    /// Translate QPA's `Coefficients(B, v)` delegate-basis call.
    ///
    /// `Ok(None)` is GAP's `fail` for a vector outside the basis span.
    ///
    /// # Errors
    ///
    /// Rejects a foreign element family, a structural `NewBasis` without a
    /// delegate, or malformed physical coordinates.
    pub fn coefficients(
        &self,
        element: &PathAlgebraModuleElement<F>,
    ) -> Result<Option<Vec<F::Element>>, PathAlgebraModuleError> {
        if !self.space.module.is_identical(element.module()) {
            return Err(PathAlgebraModuleError::ForeignElementBasisVector);
        }
        let delegate = self
            .delegate_basis
            .as_ref()
            .ok_or(PathAlgebraModuleError::MissingElementBasisDelegate)?;
        Ok(delegate.solution_row(&element.flat_external_representation())?)
    }

    /// Translate QPA's `CanonicalBasis` method for a path-algebra matrix
    /// module.
    ///
    /// # Errors
    ///
    /// Preserves QPA's zero-module failure as a typed error and propagates raw
    /// element and delegate-basis construction.
    pub fn canonical_for_module(
        module: &PathAlgebraModule<F>,
    ) -> Result<Self, PathAlgebraModuleError> {
        if module.dimension() == 0 {
            return Err(PathAlgebraModuleError::ZeroModuleCanonicalBasisFailure);
        }
        let field = module.right_acting_algebra().field();
        let generators = module.basis_elements();
        let mut vertex_rows = vec![Vec::new(); module.dimension_vector().len()];
        for generator in &generators {
            for (vertex, row) in generator.external_representation().iter().enumerate() {
                if row.iter().any(|entry| !field.is_zero(entry)) {
                    vertex_rows[vertex].push(row.clone());
                }
            }
        }
        let mut vectors = Vec::new();
        for (vertex, rows) in vertex_rows.into_iter().enumerate() {
            for row in rows {
                let mut coordinates = module.zero_element().external_representation().to_vec();
                coordinates[vertex] = row;
                vectors.push(module.element_from_external_representation(coordinates)?);
            }
        }
        let space = PathAlgebraModuleElementSpace::by_generators(module.clone(), generators)?;
        let mut basis = space.basis_nc(vectors)?;
        basis.canonical = true;
        Ok(basis)
    }
}

fn element_basis_matrix<F: Field>(
    module: &PathAlgebraModule<F>,
    vectors: &[PathAlgebraModuleElement<F>],
) -> Result<DenseMatrix<F>, MatrixError> {
    if vectors.is_empty() {
        let physical_dimension = module
            .dimension_vector()
            .iter()
            .map(|&dimension| dimension.max(1))
            .sum();
        return Ok(DenseMatrix::zero(
            0,
            physical_dimension,
            module.right_acting_algebra().field().clone(),
        ));
    }
    DenseMatrix::from_rows(
        vectors
            .iter()
            .map(PathAlgebraModuleElement::flat_external_representation)
            .collect(),
        module.right_acting_algebra().field().clone(),
    )
}

impl<F: Field> PartialEq for PathAlgebraModuleElement<F> {
    fn eq(&self, other: &Self) -> bool {
        self.module.is_identical(&other.module) && self.coordinates == other.coordinates
    }
}

impl<F: Field> Eq for PathAlgebraModuleElement<F> {}

impl<F> fmt::Display for PathAlgebraModuleElement<F>
where
    F: Field,
    F::Element: fmt::Display,
{
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.qpa_print_string_with(ToString::to_string))
    }
}

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
        write!(formatter, "<{}>", qpa_usize_list(self.dimension_vector()))
    }
}

impl<F: QpaCoefficientField> PathAlgebraModule<F> {
    /// Return QPA's reconstructive `PrintObj(M)` text.
    #[must_use]
    pub fn qpa_print_string(&self) -> String {
        format!(
            "<Module over {} with dimension vector {}>",
            qpa_quiver_algebra_view(self.right_acting_algebra()),
            qpa_usize_list(self.dimension_vector())
        )
    }

    /// Return QPA's multiline `Display(M)` text.
    ///
    /// The coefficient callback supplies GAP's scalar `String` representation;
    /// Rust performs the exact default (non-compact) `PrintArray` alignment.
    ///
    /// # Errors
    ///
    /// QPA itself errors when `Display` reaches the scalar-zero arrow entry of
    /// a zero module, rather than a dense matrix. Rust reports that branch as
    /// a typed error.
    pub fn qpa_display_string_with(
        &self,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Result<String, PathAlgebraModuleError> {
        let mut output = format!("{} and linear maps given by\n", self.qpa_print_string());
        for (arrow, matrix) in self
            .right_acting_algebra()
            .quiver()
            .arrows()
            .zip(self.matrices_of_path_algebra_module())
        {
            writeln!(&mut output, "for arrow {}:", arrow.name())
                .expect("writing to a String cannot fail");
            let PathAlgebraModuleMatrix::Matrix(matrix) = matrix else {
                return Err(PathAlgebraModuleError::DisplayScalarZeroArrowMatrix);
            };
            output.push_str(&qpa_print_array_with(&matrix, coefficient));
        }
        Ok(output)
    }
}

impl<F: Field> PartialEq for PathAlgebraModule<F> {
    fn eq(&self, other: &Self) -> bool {
        if self.is_identical(other) {
            return true;
        }
        self.right_acting_algebra() == other.right_acting_algebra()
            && self.dimension_vector() == other.dimension_vector()
            && if self.dimension() == 0 {
                // MatricesOfPathAlgebraModule exposes only scalar zeros here.
                self.arrow_matrices().len() == other.arrow_matrices().len()
            } else {
                self.arrow_matrices() == other.arrow_matrices()
            }
    }
}

impl<F: Field> Eq for PathAlgebraModule<F> {}

impl<F: Field> PathAlgebraModuleSubmodule<F> {
    /// Return QPA's ambient `Parent` module.
    #[must_use]
    pub fn parent(&self) -> &PathAlgebraModule<F> {
        &self.0.parent
    }

    /// Return QPA's nonzero vertex-uniform `GeneratorsOfAlgebraModule` list.
    #[must_use]
    pub fn generators(&self) -> &[PathAlgebraModuleElement<F>] {
        &self.0.generators
    }

    /// Return whether two handles denote the identical QPA submodule object.
    #[must_use]
    pub fn is_identical(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }

    /// Return QPA's ordered `BasisVectors(CanonicalBasis(N))` as elements of
    /// the parent module.
    ///
    /// # Errors
    ///
    /// Returns the empty basis for QPA's generic zero-submodule dispatch,
    /// preserves its dependent-generator failure, and propagates exact action
    /// and row-space computations.
    pub fn canonical_basis_elements(
        &self,
    ) -> Result<Vec<PathAlgebraModuleElement<F>>, PathAlgebraModuleError> {
        canonical_basis_of_subalgebra_module(self)
    }

    /// Translate QPA's `SubmoduleAsModule(N)` conversion.
    ///
    /// # Errors
    ///
    /// Propagates QPA canonical-basis and checked representation construction.
    pub fn as_module(&self) -> Result<PathAlgebraModule<F>, PathAlgebraModuleError> {
        submodule_as_module(self)
    }
}

impl<F: Field> fmt::Debug for PathAlgebraModuleSubmodule<F> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("PathAlgebraModuleSubmodule")
            .field("parent", self.parent())
            .field("generators", &self.generators())
            .finish()
    }
}

impl<F: Field> PartialEq for PathAlgebraModuleSubmodule<F> {
    fn eq(&self, other: &Self) -> bool {
        self.is_identical(other)
    }
}

impl<F: Field> Eq for PathAlgebraModuleSubmodule<F> {}

impl<F: Field> PathAlgebraModuleBasisMapping<F> {
    /// Return QPA's source module.
    #[must_use]
    pub fn source(&self) -> &PathAlgebraModule<F> {
        &self.0.source
    }

    /// Return QPA's range module.
    #[must_use]
    pub fn range(&self) -> &PathAlgebraModule<F> {
        &self.0.range
    }

    /// Return the images corresponding to `BasisVectors(Basis(Source(map)))`.
    #[must_use]
    pub fn basis_images(&self) -> &[PathAlgebraModuleElement<F>] {
        &self.0.basis_images
    }

    /// Return QPA's installed `IsSurjective` value.
    #[must_use]
    pub fn is_surjective(&self) -> bool {
        self.0.surjective
    }

    /// Return whether two handles denote the identical QPA mapping object.
    #[must_use]
    pub fn is_identical(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }

    /// Apply the basis-image mapping by QPA's underlying left-module linear
    /// extension.
    ///
    /// # Errors
    ///
    /// Rejects a foreign source family and propagates exact element
    /// arithmetic.
    pub fn image(
        &self,
        element: &PathAlgebraModuleElement<F>,
    ) -> Result<PathAlgebraModuleElement<F>, PathAlgebraModuleError> {
        if !self.source().is_identical(element.module()) {
            return Err(PathAlgebraModuleError::ForeignBasisMappingSourceElement);
        }
        let mut image = self.range().zero_element();
        for (coefficient, basis_image) in
            element.basis_coefficients().iter().zip(self.basis_images())
        {
            image = image.add(&basis_image.scale_left(coefficient))?;
        }
        Ok(image)
    }
}

impl<F: Field> fmt::Debug for PathAlgebraModuleBasisMapping<F> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("PathAlgebraModuleBasisMapping")
            .field("source", self.source())
            .field("range", self.range())
            .field("basis_images", &self.basis_images())
            .field("surjective", &self.is_surjective())
            .finish()
    }
}

impl<F: Field> PartialEq for PathAlgebraModuleBasisMapping<F> {
    fn eq(&self, other: &Self) -> bool {
        self.source() == other.source()
            && self.range() == other.range()
            && self.basis_images() == other.basis_images()
            && self.is_surjective() == other.is_surjective()
    }
}

impl<F: Field> Eq for PathAlgebraModuleBasisMapping<F> {}

impl<F: Field> From<PathAlgebraModule<F>> for PathAlgebraMatModuleObject<F> {
    fn from(module: PathAlgebraModule<F>) -> Self {
        Self::Module(module)
    }
}

impl<F: Field> From<PathAlgebraModuleSubmodule<F>> for PathAlgebraMatModuleObject<F> {
    fn from(submodule: PathAlgebraModuleSubmodule<F>) -> Self {
        Self::Submodule(submodule)
    }
}

impl<F: Field> BasedRightAlgebraModule<F> {
    /// Construct a based finite-dimensional right module from generator
    /// action matrices.
    ///
    /// The caller contract is that these matrices come from an actual module
    /// action. As in QPA's source conversion, this constructor checks only the
    /// representation boundary (counts, field, and common square shape), not
    /// the module axioms or algebra relations.
    ///
    /// # Errors
    ///
    /// Rejects a missing or extra generator action, a foreign coefficient
    /// field, or matrices without one common square shape.
    pub fn new(
        algebra: QuiverAlgebra<F>,
        vertex_actions: Vec<DenseMatrix<F>>,
        arrow_actions: Vec<DenseMatrix<F>>,
    ) -> Result<Self, PathAlgebraModuleError> {
        let expected_vertices = algebra.quiver().number_of_vertices();
        if vertex_actions.len() != expected_vertices {
            return Err(PathAlgebraModuleError::BasedModuleVertexActionCount {
                actual: vertex_actions.len(),
                expected: expected_vertices,
            });
        }
        let expected_arrows = algebra.quiver().number_of_arrows();
        if arrow_actions.len() != expected_arrows {
            return Err(PathAlgebraModuleError::BasedModuleArrowActionCount {
                actual: arrow_actions.len(),
                expected: expected_arrows,
            });
        }
        let dimension = vertex_actions[0].nrows();
        let expected_shape = (dimension, dimension);
        for (label, action) in algebra
            .quiver()
            .vertices()
            .map(|vertex| format!("vertex {}", vertex.name()))
            .chain(
                algebra
                    .quiver()
                    .arrows()
                    .map(|arrow| format!("arrow {}", arrow.name())),
            )
            .zip(vertex_actions.iter().chain(&arrow_actions))
        {
            if action.field() != algebra.field() {
                return Err(PathAlgebraModuleError::BasedModuleActionField {
                    generator: label.into_boxed_str(),
                });
            }
            let actual = (action.nrows(), action.ncols());
            if actual != expected_shape {
                return Err(PathAlgebraModuleError::BasedModuleActionShape {
                    generator: label.into_boxed_str(),
                    actual,
                    expected: expected_shape,
                });
            }
        }
        Ok(Self {
            algebra,
            dimension,
            vertex_actions,
            arrow_actions,
        })
    }

    /// Expose a QPA path-algebra matrix module through the generic based-module
    /// input boundary.
    ///
    /// # Errors
    ///
    /// Propagates the common based-action shape validation.
    pub fn from_path_algebra_module(
        module: &PathAlgebraModule<F>,
    ) -> Result<Self, PathAlgebraModuleError> {
        let algebra = module.right_acting_algebra().clone();
        let field = algebra.field().clone();
        let dimension = module.dimension();
        let mut offsets = Vec::with_capacity(module.dimension_vector().len());
        let mut offset = 0_usize;
        for &vertex_dimension in module.dimension_vector() {
            offsets.push(offset);
            offset += vertex_dimension;
        }
        let vertex_actions = module
            .dimension_vector()
            .iter()
            .enumerate()
            .map(|(vertex, &vertex_dimension)| {
                let mut action = DenseMatrix::zero(dimension, dimension, field.clone());
                for row in 0..vertex_dimension {
                    action[(offsets[vertex] + row, offsets[vertex] + row)] = field.one();
                }
                action
            })
            .collect::<Vec<_>>();
        let mut arrow_actions = Vec::with_capacity(algebra.quiver().number_of_arrows());
        for arrow in algebra.quiver().arrows() {
            let source = arrow.source().index().zero_based();
            let target = arrow.target().index().zero_based();
            let source_dimension = module.dimension_vector()[source];
            let target_dimension = module.dimension_vector()[target];
            let logical = logical_submatrix(
                &module.arrow_matrices()[arrow.index().zero_based()],
                source_dimension,
                target_dimension,
            )?;
            let mut action = DenseMatrix::zero(dimension, dimension, field.clone());
            for row in 0..source_dimension {
                for column in 0..target_dimension {
                    action[(offsets[source] + row, offsets[target] + column)] =
                        logical[(row, column)].clone();
                }
            }
            arrow_actions.push(action);
        }
        Self::new(algebra, vertex_actions, arrow_actions)
    }

    /// Return the right acting quiver algebra.
    #[must_use]
    pub const fn algebra(&self) -> &QuiverAlgebra<F> {
        &self.algebra
    }

    /// Return the entered common basis dimension.
    #[must_use]
    pub const fn dimension(&self) -> usize {
        self.dimension
    }

    /// Return the vertex-idempotent right-action matrices.
    #[must_use]
    pub fn vertex_actions(&self) -> &[DenseMatrix<F>] {
        &self.vertex_actions
    }

    /// Return the arrow right-action matrices.
    #[must_use]
    pub fn arrow_actions(&self) -> &[DenseMatrix<F>] {
        &self.arrow_actions
    }
}

impl<F: Field> PathAlgebraMatModuleObject<F> {
    /// Return QPA's `ActingAlgebra` value for this matrix-module object.
    #[must_use]
    pub fn right_acting_algebra(&self) -> &QuiverAlgebra<F> {
        match self {
            Self::Module(module) => module.right_acting_algebra(),
            Self::Submodule(submodule) => submodule.parent().right_acting_algebra(),
        }
    }

    /// Return the module whose element family represents this object.
    #[must_use]
    pub fn ambient_module(&self) -> &PathAlgebraModule<F> {
        match self {
            Self::Module(module) => module,
            Self::Submodule(submodule) => submodule.parent(),
        }
    }

    /// Return whether two handles denote the identical QPA object.
    #[must_use]
    pub fn is_identical(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Module(left), Self::Module(right)) => left.is_identical(right),
            (Self::Submodule(left), Self::Submodule(right)) => left.is_identical(right),
            _ => false,
        }
    }

    /// Return the basis used by QPA's generic matrix-homomorphism wrapper.
    ///
    /// # Errors
    ///
    /// Preserves failures of a submodule's ordered canonical basis.
    pub fn canonical_basis_elements(
        &self,
    ) -> Result<Vec<PathAlgebraModuleElement<F>>, PathAlgebraModuleError> {
        match self {
            Self::Module(module) => Ok(module.basis_elements()),
            Self::Submodule(submodule) => submodule.canonical_basis_elements(),
        }
    }
}

impl<F: Field> GenericPathAlgebraModuleHomomorphism<F> {
    /// Return QPA's generic mapping source.
    #[must_use]
    pub fn source(&self) -> &PathAlgebraMatModuleObject<F> {
        &self.0.source
    }

    /// Return QPA's generic mapping range.
    #[must_use]
    pub fn range(&self) -> &PathAlgebraMatModuleObject<F> {
        &self.0.range
    }

    /// Return the full row-action matrix in the retained source/range bases.
    #[must_use]
    pub fn matrix(&self) -> &DenseMatrix<F> {
        &self.0.matrix
    }

    /// Return QPA's explicitly installed mapping filter.
    #[must_use]
    pub const fn is_algebra_module_homomorphism(&self) -> bool {
        true
    }

    /// Return whether two handles denote the identical GAP mapping object.
    #[must_use]
    pub fn is_identical(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }

    /// Apply QPA's basis-matrix mapping to an ambient module element.
    ///
    /// # Errors
    ///
    /// Rejects a foreign element family or an element outside a submodule
    /// source, and propagates exact row-coordinate arithmetic.
    pub fn image(
        &self,
        element: &PathAlgebraModuleElement<F>,
    ) -> Result<PathAlgebraModuleElement<F>, PathAlgebraModuleError> {
        if !self
            .source()
            .ambient_module()
            .is_identical(element.module())
        {
            return Err(PathAlgebraModuleError::ForeignGenericHomomorphismSourceElement);
        }
        let field = self.source().right_acting_algebra().field().clone();
        let source_basis = DenseMatrix::from_rows(
            self.0
                .source_basis
                .iter()
                .map(PathAlgebraModuleElement::basis_coefficients)
                .collect(),
            field.clone(),
        )?;
        let coefficients = source_basis
            .solution_row(&element.basis_coefficients())?
            .ok_or(PathAlgebraModuleError::ElementOutsideGenericHomomorphismSource)?;
        let coefficient_row = DenseMatrix::from_flat(1, coefficients.len(), coefficients, field)?;
        let image_coefficients = coefficient_row.multiply(self.matrix())?.row(0).to_vec();
        let mut image = self.range().ambient_module().zero_element();
        for (coefficient, basis_element) in image_coefficients.iter().zip(self.0.range_basis.iter())
        {
            image = image.add(&basis_element.scale_left(coefficient))?;
        }
        Ok(image)
    }
}

impl<F: Field> fmt::Debug for GenericPathAlgebraModuleHomomorphism<F> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("GenericPathAlgebraModuleHomomorphism")
            .field("source", self.source())
            .field("range", self.range())
            .field("matrix", self.matrix())
            .finish()
    }
}

impl<F: Field> PartialEq for GenericPathAlgebraModuleHomomorphism<F> {
    fn eq(&self, other: &Self) -> bool {
        self.source() == other.source()
            && self.range() == other.range()
            && self.matrix() == other.matrix()
    }
}

impl<F: Field> Eq for GenericPathAlgebraModuleHomomorphism<F> {}

impl<F: Field> GenericPathAlgebraModuleHomSpace<F> {
    /// Return the common coefficient field.
    #[must_use]
    pub const fn field(&self) -> &F {
        &self.field
    }

    /// Return QPA's entered source object.
    #[must_use]
    pub const fn source(&self) -> &PathAlgebraMatModuleObject<F> {
        &self.source
    }

    /// Return QPA's entered range object.
    #[must_use]
    pub const fn range(&self) -> &PathAlgebraMatModuleObject<F> {
        &self.range
    }

    /// Return `Dimension(Hom(R, M, N))`.
    #[must_use]
    pub fn dimension(&self) -> usize {
        self.generators.len()
    }

    /// Return QPA's `GeneratorsOfLeftModule` list.
    #[must_use]
    pub fn generators(&self) -> &[GenericPathAlgebraModuleHomomorphism<F>] {
        &self.generators
    }

    /// Express a compatible full matrix in QPA's retained generator basis.
    ///
    /// `Ok(None)` means that the matrix is outside the represented Hom space.
    ///
    /// # Errors
    ///
    /// Propagates a different field or matrix shape.
    pub fn coefficients(
        &self,
        matrix: &DenseMatrix<F>,
    ) -> Result<Option<Vec<F::Element>>, MatrixError> {
        if matrix.field() != &self.field {
            return Err(MatrixError::DifferentFields);
        }
        self.flattened_generators
            .solution_row(&flatten_matrix(matrix))
    }
}

impl<F: Field> GenericPathAlgebraModuleEndomorphismAlgebra<F> {
    /// Return the common coefficient field.
    #[must_use]
    pub const fn field(&self) -> &F {
        &self.field
    }

    /// Return the entered QPA matrix-module object.
    #[must_use]
    pub const fn module(&self) -> &PathAlgebraMatModuleObject<F> {
        &self.module
    }

    /// Return `Dimension(End(R, M))`.
    #[must_use]
    pub fn dimension(&self) -> usize {
        self.generators.len()
    }

    /// Return the basis list passed to QPA's `AlgebraWithOne` constructor.
    #[must_use]
    pub fn generators(&self) -> &[GenericPathAlgebraModuleHomomorphism<F>] {
        &self.generators
    }

    /// Return the identity endomorphism on the retained object basis.
    #[must_use]
    pub const fn one(&self) -> &GenericPathAlgebraModuleHomomorphism<F> {
        &self.one
    }

    /// Return the zero endomorphism in the retained object basis.
    #[must_use]
    pub fn zero(&self) -> GenericPathAlgebraModuleHomomorphism<F> {
        self.mapping_from_matrix(DenseMatrix::zero(
            self.one.matrix().nrows(),
            self.one.matrix().ncols(),
            self.field.clone(),
        ))
    }

    /// Multiply two represented generic endomorphisms in row-action order.
    ///
    /// # Errors
    ///
    /// Rejects mappings outside this algebra object and propagates exact
    /// matrix multiplication.
    pub fn multiply(
        &self,
        left: &GenericPathAlgebraModuleHomomorphism<F>,
        right: &GenericPathAlgebraModuleHomomorphism<F>,
    ) -> Result<GenericPathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
        if !left.source().is_identical(&self.module)
            || !left.range().is_identical(&self.module)
            || !right.source().is_identical(&self.module)
            || !right.range().is_identical(&self.module)
        {
            return Err(PathAlgebraModuleError::ForeignGenericEndomorphism);
        }
        Ok(self.mapping_from_matrix(left.matrix().multiply(right.matrix())?))
    }

    /// Express an endomorphism matrix in QPA's supplied algebra basis.
    ///
    /// # Errors
    ///
    /// Propagates a different field or matrix shape.
    pub fn coefficients(
        &self,
        matrix: &DenseMatrix<F>,
    ) -> Result<Option<Vec<F::Element>>, MatrixError> {
        if matrix.field() != &self.field {
            return Err(MatrixError::DifferentFields);
        }
        self.flattened_generators
            .solution_row(&flatten_matrix(matrix))
    }

    fn mapping_from_matrix(
        &self,
        matrix: DenseMatrix<F>,
    ) -> GenericPathAlgebraModuleHomomorphism<F> {
        GenericPathAlgebraModuleHomomorphism(Arc::new(GenericPathAlgebraModuleHomomorphismData {
            source: self.module.clone(),
            range: self.module.clone(),
            source_basis: self.one.0.source_basis.clone(),
            range_basis: self.one.0.range_basis.clone(),
            matrix,
        }))
    }
}

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
    star_homomorphism: OnceLock<PathAlgebraModuleHomomorphism<F>>,
    transpose_homomorphism: OnceLock<PathAlgebraModuleHomomorphism<F>>,
    d_tr_homomorphism: OnceLock<PathAlgebraModuleHomomorphism<F>>,
    tr_d_homomorphism: OnceLock<PathAlgebraModuleHomomorphism<F>>,
    nakayama_homomorphism: OnceLock<PathAlgebraModuleHomomorphism<F>>,
    opposite_nakayama_homomorphism: OnceLock<PathAlgebraModuleHomomorphism<F>>,
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
