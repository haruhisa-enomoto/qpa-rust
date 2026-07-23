/// A QPA homomorphism between right path-algebra matrix modules.
///
/// The matrix at vertex `i` acts on row vectors from the source vertex space
/// to the range vertex space. Physical one-row/one-column zero-space
/// stand-ins are retained exactly as in `RightModuleHomOverAlgebra`. Clones
/// retain mapping identity so QPA's cached kernel, image, and cokernel
/// attributes can return the identical map object.
#[derive(Clone, Debug)]
pub struct PathAlgebraModuleHomomorphism<F: Field>(Arc<PathAlgebraModuleHomomorphismData<F>>);

impl<F: Field> fmt::Display for PathAlgebraModuleHomomorphism<F> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "<{} ---> {}>", self.source(), self.range())
    }
}
impl<F: QpaCoefficientField> PathAlgebraModuleHomomorphism<F> {
    /// Return QPA's reconstructive `PrintObj(f)` text.
    #[must_use]
    pub fn qpa_print_string(&self) -> String {
        format!(
            "<{} ---> {}>",
            self.source().qpa_print_string(),
            self.range().qpa_print_string()
        )
    }

    /// Return QPA's multiline `Display(f)` text.
    ///
    /// The coefficient callback supplies GAP's scalar `String`
    /// representation used by the same `PrintArray` translation as module
    /// display.
    #[must_use]
    pub fn qpa_display_string_with(
        &self,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> String {
        let mut output = format!("{}\nwith ", self.qpa_print_string());
        for (position, matrix) in self.vertex_matrices().iter().enumerate() {
            writeln!(
                &mut output,
                "linear map for vertex number {}:",
                position + 1
            )
            .expect("writing to a String cannot fail");
            output.push_str(&qpa_print_array_with(matrix, coefficient));
        }
        output
    }
}

fn qpa_usize_list(values: &[usize]) -> String {
    qpa_list(values.iter().map(ToString::to_string))
}

fn qpa_list(values: impl IntoIterator<Item = String>) -> String {
    format!("[ {} ]", values.into_iter().collect::<Vec<_>>().join(", "))
}

fn qpa_quiver_algebra_view<F: QpaCoefficientField>(algebra: &QuiverAlgebra<F>) -> String {
    let field = match algebra.field().qpa_characteristic() {
        QpaFieldCharacteristic::Zero => "Rationals".to_owned(),
        QpaFieldCharacteristic::Prime(prime) => algebra
            .field()
            .qpa_order()
            .map_or_else(|| format!("GF({prime})"), |order| format!("GF({order})")),
    };
    let quiver = algebra.quiver();
    let parent = format!(
        "<{field}[<quiver with {} vertices and {} arrows>]>",
        quiver.number_of_vertices(),
        quiver.number_of_arrows()
    );
    match algebra {
        QuiverAlgebra::Path(_) => parent,
        QuiverAlgebra::Quotient(quotient) => {
            let count = quotient.ideal().generators().len();
            let suffix = if count == 1 { "" } else { "s" };
            format!(
                "<{field}[<quiver with {} vertices and {} arrows>]/<two-sided ideal in {parent}, ({count} generator{suffix})>>",
                quiver.number_of_vertices(),
                quiver.number_of_arrows()
            )
        }
    }
}

fn qpa_print_array_with<F: Field>(
    matrix: &DenseMatrix<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> String {
    if matrix.nrows() == 0 {
        return "[ ]\n".to_owned();
    }
    if matrix.nrows() == 1 && matrix.ncols() == 0 {
        return "[ [ ] ]\n".to_owned();
    }
    let rows = (0..matrix.nrows())
        .map(|row| matrix.row(row).iter().map(coefficient).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let maximum = rows
        .iter()
        .flat_map(|row| row.iter())
        .map(String::len)
        .max()
        .unwrap_or(1);
    let width = maximum + 1;
    let mut output = "[ ".to_owned();
    for (row_index, row) in rows.iter().enumerate() {
        if row_index != 0 {
            output.push_str("  ");
        }
        output.push_str("[ ");
        if row.is_empty() {
            output.push_str("  ]");
        } else {
            for (column, entry) in row.iter().enumerate() {
                write!(&mut output, "{entry:>width$}").expect("writing to a String cannot fail");
                if column + 1 == row.len() {
                    output.push_str(" ]");
                } else {
                    output.push_str(", ");
                }
            }
        }
        if row_index + 1 == rows.len() {
            output.push_str(" ]\n");
        } else {
            output.push_str(",\n");
        }
    }
    output
}

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

/// QPA's `NumberOfNonIsoDirSummands` result.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct NonisomorphicDirectSummandData {
    count: usize,
    block_dimensions: Vec<usize>,
}

impl NonisomorphicDirectSummandData {
    /// Number of nonisomorphic indecomposable direct-summand types.
    #[must_use]
    pub const fn count(&self) -> usize {
        self.count
    }

    /// Dimensions of the simple blocks of `End(M) / rad End(M)` in QPA order.
    #[must_use]
    pub fn block_dimensions(&self) -> &[usize] {
        &self.block_dimensions
    }
}

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

/// QPA's successful `CommonDirectSummand(M, N)` result.
///
/// The two common modules are isomorphic images in `M` and `N`; the other two
/// fields are the corresponding Fitting kernels retained as complements.
#[derive(Clone, Debug)]
pub struct CommonDirectSummand<F: Field> {
    source_common: PathAlgebraModule<F>,
    source_complement: PathAlgebraModule<F>,
    range_common: PathAlgebraModule<F>,
    range_complement: PathAlgebraModule<F>,
}

impl<F: Field> CommonDirectSummand<F> {
    /// Return the common summand found inside the first entered module.
    #[must_use]
    pub const fn source_common(&self) -> &PathAlgebraModule<F> {
        &self.source_common
    }

    /// Return its Fitting complement inside the first entered module.
    #[must_use]
    pub const fn source_complement(&self) -> &PathAlgebraModule<F> {
        &self.source_complement
    }

    /// Return the isomorphic common summand inside the second entered module.
    #[must_use]
    pub const fn range_common(&self) -> &PathAlgebraModule<F> {
        &self.range_common
    }

    /// Return its Fitting complement inside the second entered module.
    #[must_use]
    pub const fn range_complement(&self) -> &PathAlgebraModule<F> {
        &self.range_complement
    }
}

/// QPA's heterogeneous `CommonDirectSummand` return contract.
#[derive(Clone, Debug)]
pub enum CommonDirectSummandResult<F: Field> {
    /// QPA prints a diagnostic and returns `fail` for different acting algebras.
    DifferentAlgebras,
    /// QPA returns `false` after exhausting the ordered Hom-basis search.
    NoCommonSummand,
    /// QPA returns the images and kernels of the two Fitting endomorphisms.
    Found(CommonDirectSummand<F>),
}

/// QPA's successful `MaximalCommonDirectSummand(M, N)` result.
#[derive(Clone, Debug)]
pub struct MaximalCommonDirectSummand<F: Field> {
    common_summands: Vec<PathAlgebraModule<F>>,
    source_remainder: PathAlgebraModule<F>,
    range_remainder: PathAlgebraModule<F>,
}

impl<F: Field> MaximalCommonDirectSummand<F> {
    /// Return the common summands in QPA discovery order.
    #[must_use]
    pub fn common_summands(&self) -> &[PathAlgebraModule<F>] {
        &self.common_summands
    }

    /// Return the remainder of the first entered module.
    #[must_use]
    pub const fn source_remainder(&self) -> &PathAlgebraModule<F> {
        &self.source_remainder
    }

    /// Return the remainder of the second entered module.
    #[must_use]
    pub const fn range_remainder(&self) -> &PathAlgebraModule<F> {
        &self.range_remainder
    }
}

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

#[derive(Debug)]
struct WeakPathAlgebraModuleHomomorphismListCache<F: Field> {
    homomorphisms: RwLock<Vec<WeakPathAlgebraModuleHomomorphism<F>>>,
    bound: OnceLock<()>,
}

impl<F: Field> WeakPathAlgebraModuleHomomorphismListCache<F> {
    const fn new() -> Self {
        Self {
            homomorphisms: RwLock::new(Vec::new()),
            bound: OnceLock::new(),
        }
    }

    fn cached(&self) -> Option<Vec<PathAlgebraModuleHomomorphism<F>>> {
        self.bound.get()?;
        self.homomorphisms
            .read()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .iter()
            .map(WeakPathAlgebraModuleHomomorphism::upgrade)
            .collect()
    }

    fn merge(
        &self,
        generated: &[PathAlgebraModuleHomomorphism<F>],
    ) -> Vec<PathAlgebraModuleHomomorphism<F>> {
        let mut slot = self
            .homomorphisms
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if self.bound.get().is_none() || slot.len() != generated.len() {
            *slot = generated
                .iter()
                .map(PathAlgebraModuleHomomorphism::downgrade)
                .collect();
            let _ = self.bound.set(());
            return generated.to_vec();
        }
        generated
            .iter()
            .enumerate()
            .map(|(index, homomorphism)| {
                if let Some(cached) = slot[index].upgrade() {
                    cached
                } else {
                    slot[index] = homomorphism.downgrade();
                    homomorphism.clone()
                }
            })
            .collect()
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

    /// Return the exact common coefficient field of the represented Hom space.
    #[must_use]
    pub fn field(&self) -> &F {
        &self.field
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

impl<F: QpaCoefficientField> ModuleEndomorphismAlgebra<F> {
    /// Convert QPA's ordered `EndOverAlgebra` matrix basis to structure
    /// constants without changing its multiplication order.
    ///
    /// # Errors
    ///
    /// Every matrix product must remain in the represented basis. Exact
    /// multiplication, coordinate, and structure-constant construction
    /// failures propagate.
    pub fn structure_constants_algebra(
        &self,
    ) -> Result<StructureConstantsAlgebra<F>, PathAlgebraModuleError> {
        let dimension = self.dimension();
        let capacity = dimension
            .checked_mul(dimension)
            .and_then(|square| square.checked_mul(dimension))
            .ok_or(PathAlgebraModuleError::HomEquationSizeOverflow)?;
        let mut products = Vec::with_capacity(capacity);
        for left in self.basis_matrices() {
            for right in self.basis_matrices() {
                let product = self.multiply(left, right)?;
                products.extend(
                    self.coefficients(&product)?
                        .ok_or(PathAlgebraModuleError::EndomorphismAlgebraProductOutsideBasis)?,
                );
            }
        }
        Ok(StructureConstantsAlgebra::new(
            self.field().clone(),
            dimension,
            products,
        )?)
    }

    /// Construct the structure-constant model generated by the transposes of
    /// QPA's ordered `EndOverAlgebra` basis matrices.
    ///
    /// This is the exact `AlgebraWithOne(K, genmats)` input used by
    /// `DecomposeModule`. The returned matrices are the transposes of GAP's
    /// canonical `AlgebraWithOne` basis and therefore convert its coordinates
    /// back to `EndOverAlgebra` matrices.
    ///
    /// # Errors
    ///
    /// Every matrix product must remain in GAP's generated canonical basis.
    /// Exact multiplication, coordinate, and structure-constant construction
    /// failures propagate.
    pub(crate) fn transposed_structure_constants_model(
        &self,
    ) -> Result<(StructureConstantsAlgebra<F>, Vec<DenseMatrix<F>>), PathAlgebraModuleError> {
        let matrix_dimension = self.0.one.nrows();
        let flattened_dimension = matrix_dimension
            .checked_mul(matrix_dimension)
            .ok_or(PathAlgebraModuleError::HomEquationSizeOverflow)?;
        let generator_count = self
            .dimension()
            .checked_add(1)
            .ok_or(PathAlgebraModuleError::HomEquationSizeOverflow)?;
        let generator_capacity = generator_count
            .checked_mul(flattened_dimension)
            .ok_or(PathAlgebraModuleError::HomEquationSizeOverflow)?;
        let mut generators = Vec::with_capacity(generator_capacity);
        generators.extend(flatten_matrix(self.one()));
        for basis in self.basis_matrices() {
            generators.extend(flatten_matrix(&basis.transpose()));
        }
        let flattened_basis = DenseMatrix::from_flat(
            generator_count,
            flattened_dimension,
            generators,
            self.field().clone(),
        )?
        .semi_echelon_transformation()
        .vectors;
        let transposed_basis = (0..flattened_basis.nrows())
            .map(|row| {
                DenseMatrix::from_flat(
                    matrix_dimension,
                    matrix_dimension,
                    flattened_basis.row(row).to_vec(),
                    self.field().clone(),
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        let dimension = transposed_basis.len();
        let capacity = dimension
            .checked_mul(dimension)
            .and_then(|square| square.checked_mul(dimension))
            .ok_or(PathAlgebraModuleError::HomEquationSizeOverflow)?;
        let mut products = Vec::with_capacity(capacity);
        for left in &transposed_basis {
            for right in &transposed_basis {
                let product = left.multiply(right)?;
                products.extend(
                    flattened_basis
                        .solution_row(&flatten_matrix(&product))?
                        .ok_or(PathAlgebraModuleError::EndomorphismAlgebraProductOutsideBasis)?,
                );
            }
        }
        let algebra = StructureConstantsAlgebra::new(self.field().clone(), dimension, products)?;
        let endomorphism_basis = transposed_basis
            .into_iter()
            .map(|matrix| matrix.transpose())
            .collect();
        Ok((algebra, endomorphism_basis))
    }

    /// Translate QPA's `RadicalOfAlgebra` conversion for an algebra represented
    /// by linear mappings.
    ///
    /// QPA transposes the mapping matrices, computes the radical in the
    /// resulting ordinary matrix algebra, and transposes the radical basis
    /// back. The represented endomorphism algebra already stores the logical
    /// ordinary matrices, so Rust constructs the identical multiplication
    /// table directly and lifts the radical coordinates to those matrices.
    ///
    /// # Errors
    ///
    /// Propagates matrix multiplication/coordinate failures and the exact
    /// structure-constant radical computation.
    pub fn radical_basis_matrices(&self) -> Result<Vec<DenseMatrix<F>>, PathAlgebraModuleError> {
        let algebra = self.structure_constants_algebra()?;
        let radical = algebra.radical()?;
        (0..radical.basis().nrows())
            .map(|row| {
                self.matrix_from_coefficients(radical.basis().row(row))
                    .map_err(PathAlgebraModuleError::from)
            })
            .collect()
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

    /// Return the exact coefficient field inherited from `End(M)`.
    #[must_use]
    pub fn field(&self) -> &F {
        self.0.parent.field()
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

    /// Return QPA's cached `StarOfModuleHomomorphism` attribute, if set.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_star_homomorphism(&self) -> Option<Self> {
        self.0.star_homomorphism.get().cloned()
    }

    /// Install QPA's `StarOfModuleHomomorphism` attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_star_homomorphism(&self, value: &Self) -> Self {
        self.0
            .star_homomorphism
            .get_or_init(|| value.clone())
            .clone()
    }

    /// Return QPA's cached `TransposeOfModuleHomomorphism` attribute, if set.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_transpose_homomorphism(&self) -> Option<Self> {
        self.0.transpose_homomorphism.get().cloned()
    }

    /// Install QPA's one-way `TransposeOfModuleHomomorphism` attribute cache.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_transpose_homomorphism(&self, value: &Self) -> Self {
        self.0
            .transpose_homomorphism
            .get_or_init(|| value.clone())
            .clone()
    }

    /// Return QPA's cached `DTr` homomorphism attribute, if set.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_d_tr_homomorphism(&self) -> Option<Self> {
        self.0.d_tr_homomorphism.get().cloned()
    }

    /// Install QPA's one-way `DTr` homomorphism attribute cache.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_d_tr_homomorphism(&self, value: &Self) -> Self {
        self.0
            .d_tr_homomorphism
            .get_or_init(|| value.clone())
            .clone()
    }

    /// Return QPA's cached `TrD` homomorphism attribute, if set.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_tr_d_homomorphism(&self) -> Option<Self> {
        self.0.tr_d_homomorphism.get().cloned()
    }

    /// Install QPA's one-way `TrD` homomorphism attribute cache.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_tr_d_homomorphism(&self, value: &Self) -> Self {
        self.0
            .tr_d_homomorphism
            .get_or_init(|| value.clone())
            .clone()
    }

    /// Return QPA's cached `NakayamaFunctorOfModuleHomomorphism` attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_nakayama_functor_homomorphism(&self) -> Option<Self> {
        self.0.nakayama_homomorphism.get().cloned()
    }

    /// Install QPA's `NakayamaFunctorOfModuleHomomorphism` attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_nakayama_functor_homomorphism(&self, value: &Self) -> Self {
        self.0
            .nakayama_homomorphism
            .get_or_init(|| value.clone())
            .clone()
    }

    /// Return QPA's cached opposite-Nakayama homomorphism attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_opposite_nakayama_functor_homomorphism(&self) -> Option<Self> {
        self.0.opposite_nakayama_homomorphism.get().cloned()
    }

    /// Install QPA's opposite-Nakayama homomorphism attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_opposite_nakayama_functor_homomorphism(&self, value: &Self) -> Self {
        self.0
            .opposite_nakayama_homomorphism
            .get_or_init(|| value.clone())
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
