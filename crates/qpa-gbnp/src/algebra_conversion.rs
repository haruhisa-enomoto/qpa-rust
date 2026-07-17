//! QPA's conversion of a finite-dimensional algebra to a quiver algebra.
//!
//! The implementation follows `AlgebraAsQuiverAlgebra` in `lib/algebra.gi`.
//! The first field contract is the split elementary prime-field branch: this
//! is exactly the branch used by QPA once its constructed field `K` equals the
//! entered prime field. Extension-field scalar restriction remains explicit.

use numerica::domains::Field;
use qpa_core::{
    DenseMatrix, MatrixError, ModuleEndomorphismAlgebra, PathAlgebra, PathAlgebraElement,
    PathAlgebraError, PathAlgebraModule, PathAlgebraModuleError, PathEndpoint, QpaCoefficientField,
    Quiver, QuiverAlgebra, QuiverError, QuotientError, StructureConstantsAlgebra,
    StructureConstantsError, end_over_algebra,
};
use thiserror::Error;

use crate::{PathalgGbnpError, quotient_by_relations};

/// QPA's `[KQ/I, images]` return value from `AlgebraAsQuiverAlgebra`.
#[derive(Clone, Debug)]
pub struct AlgebraAsQuiverAlgebraResult<F: Field> {
    algebra: QuiverAlgebra<F>,
    generator_images: Vec<Vec<F::Element>>,
    loewy_length: usize,
}

impl<F: Field> AlgebraAsQuiverAlgebraResult<F> {
    /// Return the path algebra or quotient path algebra constructed by QPA.
    #[must_use]
    pub const fn algebra(&self) -> &QuiverAlgebra<F> {
        &self.algebra
    }

    /// Return images of vertices followed by arrows in the entered algebra's
    /// standard coordinates.
    #[must_use]
    pub fn generator_images(&self) -> &[Vec<F::Element>] {
        &self.generator_images
    }

    /// Return QPA's `Length(RadicalSeriesOfAlgebra(A)) - 1` exponent.
    #[must_use]
    pub const fn loewy_length(&self) -> usize {
        self.loewy_length
    }

    /// Consume the result and return QPA's two public list components.
    #[must_use]
    pub fn into_parts(self) -> (QuiverAlgebra<F>, Vec<Vec<F::Element>>) {
        (self.algebra, self.generator_images)
    }
}

/// QPA's successful `EndOfModuleAsQuiverAlgebra` return value.
#[derive(Clone, Debug)]
pub struct EndOfModuleAsQuiverAlgebraResult<F: Field> {
    endomorphism_algebra: ModuleEndomorphismAlgebra<F>,
    quiver_algebra: QuiverAlgebra<F>,
}

impl<F: Field> EndOfModuleAsQuiverAlgebraResult<F> {
    /// Return the entered module's QPA endomorphism algebra.
    #[must_use]
    pub const fn endomorphism_algebra(&self) -> &ModuleEndomorphismAlgebra<F> {
        &self.endomorphism_algebra
    }

    /// Return the adjacency matrix of the converted endomorphism quiver.
    #[must_use]
    pub fn adjacency_matrix(&self) -> &[Box<[usize]>] {
        self.quiver_algebra.quiver().adjacency_matrix()
    }

    /// Return the endomorphism ring as QPA's quiver algebra.
    #[must_use]
    pub const fn quiver_algebra(&self) -> &QuiverAlgebra<F> {
        &self.quiver_algebra
    }

    /// Consume QPA's three-component result in the same order.
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        ModuleEndomorphismAlgebra<F>,
        Vec<Vec<usize>>,
        QuiverAlgebra<F>,
    ) {
        let adjacency = self
            .quiver_algebra
            .quiver()
            .adjacency_matrix()
            .iter()
            .map(|row| row.to_vec())
            .collect();
        (self.endomorphism_algebra, adjacency, self.quiver_algebra)
    }
}

/// Translate QPA's `EndOfModuleAsQuiverAlgebra(M)` adapter.
///
/// QPA constructs `EndOverAlgebra(M)`, delegates the actual computation to
/// `AlgebraAsQuiverAlgebra`, and returns the endomorphism algebra, converted
/// quiver adjacency matrix, and converted quiver algebra. Rust preserves that
/// exact dependency chain rather than reconstructing a quiver from a module
/// decomposition.
///
/// # Errors
///
/// Propagates endomorphism structure-constant conversion and QPA's finite,
/// indecomposable, basic, elementary prime-field conversion boundaries.
pub fn end_of_module_as_quiver_algebra<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
) -> Result<EndOfModuleAsQuiverAlgebraResult<F>, EndOfModuleAsQuiverAlgebraError> {
    let endomorphism_algebra = end_over_algebra(module)?;
    let structure_constants = endomorphism_algebra.structure_constants_algebra()?;
    let converted = algebra_as_quiver_algebra(&structure_constants)?;
    let (quiver_algebra, _) = converted.into_parts();
    Ok(EndOfModuleAsQuiverAlgebraResult {
        endomorphism_algebra,
        quiver_algebra,
    })
}

/// Translate QPA's `AlgebraAsQuiverAlgebra` on its split prime-field branch.
///
/// QPA checks indecomposability, basicness, and elementarity; lifts primitive
/// idempotents across the radical; obtains arrows from
/// `e_i rad(A)/rad(A)^2 e_j`; constructs `KQ/J^t`; and finds the final
/// relations as the left nullspace of the resulting surjection to `A`. This
/// function preserves that sequence and uses the translated default GBNP
/// backend for both path-algebra quotients.
///
/// The current Numerica contract supports the case in which QPA's common
/// elementary component field is the entered prime field. A common component
/// of larger dimension would make QPA construct `GF(p^r)` and is reported
/// explicitly until an extension-field backend exists.
///
/// # Errors
///
/// Reproduces QPA's finite-field, indecomposable, basic, and elementary
/// rejections; reports the deferred extension-field branch; and propagates
/// exact structure-constant, matrix, quiver, path, quotient, and GBNP errors.
#[allow(clippy::too_many_lines)]
pub fn algebra_as_quiver_algebra<F: QpaCoefficientField>(
    algebra: &StructureConstantsAlgebra<F>,
) -> Result<AlgebraAsQuiverAlgebraResult<F>, AlgebraAsQuiverAlgebraError> {
    if algebra.field().qpa_elements()?.is_none() {
        return Err(AlgebraAsQuiverAlgebraError::InfiniteCoefficientField);
    }
    let identity = algebra
        .find_multiplicative_identity()?
        .ok_or(AlgebraAsQuiverAlgebraError::MissingMultiplicativeIdentity)?;

    // QPA first decomposes the center and rejects more than one block.
    let central_blocks = algebra.central_idempotents_split_prime()?;
    if central_blocks.len() > 1 {
        return Err(AlgebraAsQuiverAlgebraError::DecomposableAlgebra);
    }

    let radical = algebra.radical()?;
    let semisimple = algebra.quotient_by_ideal_basis(radical.basis().clone())?;
    let primitive_central = semisimple.algebra().central_idempotents_split_prime()?;
    if primitive_central.is_empty() {
        return Err(AlgebraAsQuiverAlgebraError::MissingSemisimpleComponent);
    }

    // `factoralgebradecomp`: c * BasisVectors(Basis(A/rad A)).
    let mut component_dimensions = Vec::with_capacity(primitive_central.len());
    for idempotent in &primitive_central {
        let mut generators = Vec::with_capacity(semisimple.algebra().dimension());
        for basis_index in 0..semisimple.algebra().dimension() {
            let mut basis = vec![algebra.field().zero(); semisimple.algebra().dimension()];
            basis[basis_index] = algebra.field().one();
            let product = semisimple
                .algebra()
                .multiply_coordinates(idempotent, &basis)?;
            if !is_zero_row(algebra.field(), &product) {
                generators.push(product);
            }
        }
        let component = span_basis(
            algebra.field(),
            semisimple.algebra().dimension(),
            generators,
        )?;
        if !subspace_is_commutative(semisimple.algebra(), &component)? {
            return Err(AlgebraAsQuiverAlgebraError::NonBasicAlgebra);
        }
        component_dimensions.push(component.nrows());
    }
    let component_dimension = component_dimensions[0];
    if component_dimensions
        .iter()
        .any(|&dimension| dimension != component_dimension)
    {
        return Err(AlgebraAsQuiverAlgebraError::NonElementaryAlgebra);
    }
    if component_dimension != 1 {
        return Err(AlgebraAsQuiverAlgebraError::ExtensionFieldRequired {
            degree: component_dimension,
        });
    }

    let vertices = semisimple.lift_complete_set_of_orthogonal_idempotents(&primitive_central)?;

    // ProductSpace(rad A, rad A), followed by A/rad(A)^2.
    let radical_square = product_space_basis(algebra, radical.basis(), radical.basis())?;
    let modulo_square = algebra.quotient_by_ideal_basis(radical_square)?;
    let projected_radical = (0..radical.basis().nrows())
        .map(|row| modulo_square.project_coordinates(radical.basis().row(row)))
        .collect::<Result<Vec<_>, _>>()?;
    let radical_modulo_square = span_basis(
        algebra.field(),
        modulo_square.algebra().dimension(),
        projected_radical,
    )?;
    let quotient_vertices = vertices
        .iter()
        .map(|vertex| modulo_square.project_coordinates(vertex))
        .collect::<Result<Vec<_>, _>>()?;

    // QPA's row-major arrows[i][j] list and adjacency matrix.
    let vertex_count = vertices.len();
    let mut adjacency = vec![vec![0usize; vertex_count]; vertex_count];
    let mut arrows = Vec::new();
    for left_vertex in 0..vertex_count {
        for right_vertex in 0..vertex_count {
            let mut candidates = Vec::new();
            for row in 0..radical_modulo_square.nrows() {
                let left_product = modulo_square.algebra().multiply_coordinates(
                    &quotient_vertices[left_vertex],
                    radical_modulo_square.row(row),
                )?;
                let component = modulo_square
                    .algebra()
                    .multiply_coordinates(&left_product, &quotient_vertices[right_vertex])?;
                if !is_zero_row(algebra.field(), &component) {
                    candidates.push(component);
                }
            }
            let arrow_space = span_basis(
                algebra.field(),
                modulo_square.algebra().dimension(),
                candidates,
            )?;
            adjacency[left_vertex][right_vertex] = arrow_space.nrows();
            for row in 0..arrow_space.nrows() {
                let representative = modulo_square.lift_coordinates(arrow_space.row(row))?;
                let left_product =
                    algebra.multiply_coordinates(&vertices[left_vertex], &representative)?;
                arrows.push(algebra.multiply_coordinates(&left_product, &vertices[right_vertex])?);
            }
        }
    }

    let loewy_length = radical_nilpotency_index(algebra, radical.basis())?;
    let exponent = i64::try_from(loewy_length)
        .map_err(|_| AlgebraAsQuiverAlgebraError::LoewyLengthOverflow)?;
    let quiver = Quiver::from_adjacency(adjacency)?;
    let path_algebra = PathAlgebra::new(algebra.field().clone(), quiver)?;
    let arrow_power = path_algebra.nth_power_of_arrow_ideal(exponent)?;
    let mut generator_images = vertices;
    generator_images.extend(arrows);

    // QPA computes a basis of AA = KQ/J^t, applies the generator map to each
    // basis vector, and takes the left nullspace of that image matrix.
    let (basis_representatives, initial_algebra) = if arrow_power.is_empty() {
        let basis = path_algebra
            .canonical_basis()
            .ok_or(AlgebraAsQuiverAlgebraError::MissingTruncatedBasis)?;
        (basis, QuiverAlgebra::Path(path_algebra.clone()))
    } else {
        let truncated = quotient_by_relations(path_algebra.clone(), &arrow_power)?;
        let basis = truncated
            .canonical_basis()?
            .ok_or(AlgebraAsQuiverAlgebraError::MissingTruncatedBasis)?;
        (
            basis
                .into_iter()
                .map(|element| element.representative().clone())
                .collect(),
            QuiverAlgebra::Quotient(truncated),
        )
    };
    let image_rows = basis_representatives
        .iter()
        .map(|basis| {
            image_of_path_element(
                algebra,
                identity.coordinates(),
                basis,
                vertex_count,
                &generator_images,
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let image_matrix = DenseMatrix::from_flat(
        image_rows.len(),
        algebra.dimension(),
        image_rows.into_iter().flatten().collect(),
        algebra.field().clone(),
    )?;
    let kernel = image_matrix.nullspace_rows();
    let mut ideal_generators = Vec::with_capacity(kernel.nrows() + arrow_power.len());
    for row in 0..kernel.nrows() {
        ideal_generators.push(linear_combination_of_path_basis(
            &path_algebra,
            &basis_representatives,
            kernel.row(row),
        )?);
    }
    if !arrow_power.is_empty() {
        ideal_generators.extend(arrow_power);
    }

    let result_algebra = if ideal_generators.is_empty() {
        initial_algebra
    } else {
        QuiverAlgebra::Quotient(quotient_by_relations(path_algebra, &ideal_generators)?)
    };
    Ok(AlgebraAsQuiverAlgebraResult {
        algebra: result_algebra,
        generator_images,
        loewy_length,
    })
}

fn span_basis<F: Field>(
    field: &F,
    width: usize,
    rows: Vec<Vec<F::Element>>,
) -> Result<DenseMatrix<F>, MatrixError> {
    DenseMatrix::from_flat(
        rows.len(),
        width,
        rows.into_iter().flatten().collect(),
        field.clone(),
    )
    .map(|matrix| matrix.semi_echelon_transformation().vectors)
}

fn product_space_basis<F: Field>(
    algebra: &StructureConstantsAlgebra<F>,
    left: &DenseMatrix<F>,
    right: &DenseMatrix<F>,
) -> Result<DenseMatrix<F>, AlgebraAsQuiverAlgebraError> {
    let mut products = Vec::with_capacity(left.nrows().saturating_mul(right.nrows()));
    for left_row in 0..left.nrows() {
        for right_row in 0..right.nrows() {
            products.push(algebra.multiply_coordinates(left.row(left_row), right.row(right_row))?);
        }
    }
    Ok(span_basis(algebra.field(), algebra.dimension(), products)?)
}

fn subspace_is_commutative<F: Field>(
    algebra: &StructureConstantsAlgebra<F>,
    basis: &DenseMatrix<F>,
) -> Result<bool, StructureConstantsError> {
    for left in 0..basis.nrows() {
        for right in left + 1..basis.nrows() {
            if algebra.multiply_coordinates(basis.row(left), basis.row(right))?
                != algebra.multiply_coordinates(basis.row(right), basis.row(left))?
            {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

fn radical_nilpotency_index<F: Field>(
    algebra: &StructureConstantsAlgebra<F>,
    radical: &DenseMatrix<F>,
) -> Result<usize, AlgebraAsQuiverAlgebraError> {
    let mut power = radical.clone();
    let mut exponent = 1usize;
    while power.nrows() > 0 {
        power = product_space_basis(algebra, &power, radical)?;
        exponent = exponent
            .checked_add(1)
            .ok_or(AlgebraAsQuiverAlgebraError::LoewyLengthOverflow)?;
        if exponent > algebra.dimension().saturating_add(1) {
            return Err(AlgebraAsQuiverAlgebraError::RadicalNotNilpotent);
        }
    }
    Ok(exponent)
}

fn image_of_path_element<F: Field>(
    algebra: &StructureConstantsAlgebra<F>,
    identity: &[F::Element],
    element: &PathAlgebraElement<F>,
    vertex_count: usize,
    images: &[Vec<F::Element>],
) -> Result<Vec<F::Element>, AlgebraAsQuiverAlgebraError> {
    let mut total = vec![algebra.field().zero(); algebra.dimension()];
    for term in element.terms() {
        let walk = term
            .path()
            .walk()
            .ok_or(AlgebraAsQuiverAlgebraError::ZeroPathInTruncatedBasis)?;
        let mut image = if walk.is_empty() {
            let PathEndpoint::Vertex(vertex) = term.path().source() else {
                return Err(AlgebraAsQuiverAlgebraError::InvalidVertexBasisPath);
            };
            images
                .get(vertex.index().get() - 1)
                .cloned()
                .ok_or(AlgebraAsQuiverAlgebraError::MissingGeneratorImage)?
        } else {
            let mut product = identity.to_vec();
            for arrow in walk {
                let arrow_image = images
                    .get(vertex_count + arrow.get() - 1)
                    .ok_or(AlgebraAsQuiverAlgebraError::MissingGeneratorImage)?;
                product = algebra.multiply_coordinates(&product, arrow_image)?;
            }
            product
        };
        for coefficient in &mut image {
            algebra.field().mul_assign(coefficient, term.coefficient());
        }
        add_row_assign(algebra.field(), &mut total, &image);
    }
    Ok(total)
}

fn linear_combination_of_path_basis<F: Field>(
    algebra: &PathAlgebra<F>,
    basis: &[PathAlgebraElement<F>],
    coefficients: &[F::Element],
) -> Result<PathAlgebraElement<F>, PathAlgebraError> {
    algebra.from_terms(
        basis
            .iter()
            .zip(coefficients)
            .flat_map(|(element, scalar)| {
                element.terms().iter().map(move |term| {
                    (
                        term.path().clone(),
                        algebra.field().mul(scalar, term.coefficient()),
                    )
                })
            }),
    )
}

fn is_zero_row<F: Field>(field: &F, row: &[F::Element]) -> bool {
    row.iter().all(|entry| field.is_zero(entry))
}

fn add_row_assign<F: Field>(field: &F, target: &mut [F::Element], source: &[F::Element]) {
    for (target, source) in target.iter_mut().zip(source) {
        field.add_assign(target, source);
    }
}

/// Failure in QPA's finite-dimensional-algebra conversion.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum AlgebraAsQuiverAlgebraError {
    /// QPA's method is restricted to finite coefficient fields.
    #[error("the entered algebra is not an algebra over a finite field")]
    InfiniteCoefficientField,
    /// QPA could not find an identity in the entered algebra.
    #[error("the entered algebra does not have a multiplicative identity")]
    MissingMultiplicativeIdentity,
    /// QPA rejects more than one central block.
    #[error("the entered algebra is not indecomposable")]
    DecomposableAlgebra,
    /// QPA's semisimple factor unexpectedly supplied no central component.
    #[error("the semisimple factor has no central component")]
    MissingSemisimpleComponent,
    /// One simple factor algebra is noncommutative.
    #[error("the entered algebra is not basic")]
    NonBasicAlgebra,
    /// The simple factor fields do not all have equal dimension.
    #[error("the entered algebra is not elementary, hence not a quotient of a path algebra")]
    NonElementaryAlgebra,
    /// QPA would construct an extension field not yet supplied by Numerica.
    #[error("the conversion requires a degree-{degree} extension of the prime coefficient field")]
    ExtensionFieldRequired {
        /// Common dimension of the elementary division-field components.
        degree: usize,
    },
    /// `Length(RadicalSeriesOfAlgebra(A)) - 1` did not fit QPA's exponent API.
    #[error("the radical-series length overflows the supported path exponent")]
    LoewyLengthOverflow,
    /// A purported Jacobson radical failed to become zero.
    #[error("the computed radical is not nilpotent")]
    RadicalNotNilpotent,
    /// QPA's `KQ/J^t` basis was unavailable.
    #[error("the truncated path algebra has no finite canonical basis")]
    MissingTruncatedBasis,
    /// A canonical basis must never contain QPA's distinguished zero path.
    #[error("the truncated canonical basis contains the zero path")]
    ZeroPathInTruncatedBasis,
    /// A length-zero canonical path did not identify a vertex.
    #[error("a truncated canonical basis path has no vertex endpoint")]
    InvalidVertexBasisPath,
    /// A basis path referred outside the QPA generator-image list.
    #[error("a truncated canonical basis path has no generator image")]
    MissingGeneratorImage,
    /// Structure-constant arithmetic or decomposition failed.
    #[error(transparent)]
    Structure(#[from] StructureConstantsError),
    /// Exact vector-space construction failed.
    #[error(transparent)]
    Matrix(#[from] MatrixError),
    /// The adjacency matrix did not define a QPA quiver.
    #[error(transparent)]
    Quiver(#[from] QuiverError),
    /// Parent path-algebra construction or arithmetic failed.
    #[error(transparent)]
    PathAlgebra(#[from] PathAlgebraError),
    /// Quotient canonical-basis construction failed.
    #[error(transparent)]
    Quotient(#[from] QuotientError),
    /// QPA's default GBNP-backed quotient construction failed.
    #[error(transparent)]
    Pathalg(#[from] PathalgGbnpError),
    /// Coefficient-field enumeration failed.
    #[error(transparent)]
    Coefficient(#[from] qpa_core::CoefficientDomainError),
}

/// Failure while applying QPA's module-endomorphism conversion adapter.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum EndOfModuleAsQuiverAlgebraError {
    /// Constructing or converting `EndOverAlgebra(M)` failed.
    #[error(transparent)]
    Module(#[from] PathAlgebraModuleError),
    /// QPA's delegated `AlgebraAsQuiverAlgebra` computation failed.
    #[error(transparent)]
    Conversion(#[from] AlgebraAsQuiverAlgebraError),
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        Field,
        finite_field::{FiniteFieldCore, Z2, Zp},
        integer::Integer,
        rational::Q,
    };
    use qpa_core::{
        ArrowSpec, DenseMatrix, NamedArrowMatrix, PathAlgebra, QpaCoefficientField, Quiver,
        QuiverAlgebra, StructureConstantsAlgebra, VertexIndex, right_module_over_path_algebra,
    };
    use serde_json::{Value, json};

    use crate::quotient_by_relations;

    use super::{
        AlgebraAsQuiverAlgebraError, EndOfModuleAsQuiverAlgebraError, algebra_as_quiver_algebra,
        end_of_module_as_quiver_algebra,
    };

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn set_basis_product<F: Field>(
        tensor: &mut [F::Element],
        field: &F,
        dimension: usize,
        left: usize,
        right: usize,
        output: usize,
    ) {
        tensor[(left * dimension + right) * dimension + output] = field.one();
    }

    fn upper_triangular_two_by_two<F: Field>(field: F) -> StructureConstantsAlgebra<F> {
        let dimension = 3;
        let mut tensor = vec![field.zero(); dimension * dimension * dimension];
        set_basis_product(&mut tensor, &field, dimension, 0, 0, 0);
        set_basis_product(&mut tensor, &field, dimension, 0, 1, 1);
        set_basis_product(&mut tensor, &field, dimension, 1, 2, 1);
        set_basis_product(&mut tensor, &field, dimension, 2, 2, 2);
        StructureConstantsAlgebra::new(field, dimension, tensor).unwrap()
    }

    fn dual_numbers<F: Field>(field: F) -> StructureConstantsAlgebra<F> {
        let dimension = 2;
        let mut tensor = vec![field.zero(); dimension * dimension * dimension];
        set_basis_product(&mut tensor, &field, dimension, 0, 0, 0);
        set_basis_product(&mut tensor, &field, dimension, 0, 1, 1);
        set_basis_product(&mut tensor, &field, dimension, 1, 0, 1);
        StructureConstantsAlgebra::new(field, dimension, tensor).unwrap()
    }

    fn product_of_two_fields<F: Field>(field: F) -> StructureConstantsAlgebra<F> {
        let dimension = 2;
        let mut tensor = vec![field.zero(); dimension * dimension * dimension];
        set_basis_product(&mut tensor, &field, dimension, 0, 0, 0);
        set_basis_product(&mut tensor, &field, dimension, 1, 1, 1);
        StructureConstantsAlgebra::new(field, dimension, tensor).unwrap()
    }

    fn two_by_two_matrix_algebra<F: Field>(field: F) -> StructureConstantsAlgebra<F> {
        // e11, e12, e21, e22 in row-major matrix-unit order.
        let dimension = 4;
        let mut tensor = vec![field.zero(); dimension * dimension * dimension];
        let units = [(0usize, 0usize), (0, 1), (1, 0), (1, 1)];
        for (left, &(i, j)) in units.iter().enumerate() {
            for (right, &(k, l)) in units.iter().enumerate() {
                if j == k {
                    let output = units
                        .iter()
                        .position(|&(row, column)| row == i && column == l)
                        .unwrap();
                    set_basis_product(&mut tensor, &field, dimension, left, right, output);
                }
            }
        }
        StructureConstantsAlgebra::new(field, dimension, tensor).unwrap()
    }

    fn conversion_json<F: qpa_core::QpaCoefficientField>(
        input: &StructureConstantsAlgebra<F>,
        coefficient_json: &impl Fn(&F::Element) -> Value,
    ) -> Value {
        let converted = algebra_as_quiver_algebra(input).unwrap();
        let result = converted.algebra();
        let dimension = match result {
            QuiverAlgebra::Path(algebra) => algebra.canonical_basis().unwrap().len(),
            QuiverAlgebra::Quotient(algebra) => algebra.canonical_basis().unwrap().unwrap().len(),
        };
        let arrows = result
            .quiver()
            .arrows()
            .map(|arrow| json!([arrow.source().index().get(), arrow.target().index().get()]))
            .collect::<Vec<_>>();
        let images = converted
            .generator_images()
            .iter()
            .map(|row| {
                row.iter()
                    .enumerate()
                    .filter(|(_, coefficient)| !input.field().is_zero(coefficient))
                    .map(|(coordinate, coefficient)| {
                        json!([coordinate + 1, coefficient_json(coefficient)])
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let relations = result
            .relations()
            .iter()
            .map(|relation| {
                relation
                    .terms()
                    .iter()
                    .map(|term| {
                        json!({
                            "walk": term.path().walk().unwrap().iter().map(|arrow| arrow.get()).collect::<Vec<_>>(),
                            "coefficient": coefficient_json(term.coefficient())
                        })
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        json!({
            "dimension": dimension,
            "vertices": result.quiver().number_of_vertices(),
            "arrows": arrows,
            "generator_images": images,
            "loewy_length": converted.loewy_length(),
            "relations": relations
        })
    }

    fn small_integer(value: &Integer) -> i64 {
        value.to_string().parse().unwrap()
    }

    fn end_of_module_json<F: QpaCoefficientField>(field: F) -> Value {
        let quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let parent = PathAlgebra::new(field.clone(), quiver).unwrap();
        let x = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let algebra = QuiverAlgebra::Quotient(
            quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap(),
        );
        let regular = right_module_over_path_algebra(
            algebra,
            vec![2],
            &[NamedArrowMatrix::new(
                "x",
                DenseMatrix::from_rows(
                    vec![
                        vec![field.zero(), field.one()],
                        vec![field.zero(), field.zero()],
                    ],
                    field,
                )
                .unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let converted = end_of_module_as_quiver_algebra(&regular).unwrap();
        let quiver_dimension = match converted.quiver_algebra() {
            QuiverAlgebra::Path(algebra) => algebra.canonical_basis().unwrap().len(),
            QuiverAlgebra::Quotient(algebra) => algebra.canonical_basis().unwrap().unwrap().len(),
        };
        json!({
            "endomorphism_dimension": converted.endomorphism_algebra().dimension(),
            "adjacency": converted.adjacency_matrix().iter()
                .map(|row| row.to_vec()).collect::<Vec<_>>(),
            "quiver_algebra_dimension": quiver_dimension,
        })
    }

    #[test]
    fn pinned_qpa_conversion_branches_match() {
        let fixture: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_module_preprojective_structure.schema-v1.json"
        ))
        .unwrap();
        assert_eq!(
            conversion_json(&upper_triangular_two_by_two(Z2), &|coefficient| {
                small_integer(&Z2.to_integer(coefficient)).into()
            }),
            fixture["conversion_cases"]["upper_triangular_f2"]
        );
        let field = Zp::new(5);
        assert_eq!(
            conversion_json(&dual_numbers(field.clone()), &|coefficient| {
                small_integer(&field.to_integer(coefficient)).into()
            }),
            fixture["conversion_cases"]["dual_numbers_f5"]
        );
    }

    #[test]
    fn qpa_conversion_rejects_infinite_decomposable_and_nonbasic_inputs() {
        assert!(matches!(
            algebra_as_quiver_algebra(&dual_numbers(Q)),
            Err(AlgebraAsQuiverAlgebraError::InfiniteCoefficientField)
        ));
        assert!(matches!(
            algebra_as_quiver_algebra(&product_of_two_fields(Z2)),
            Err(AlgebraAsQuiverAlgebraError::DecomposableAlgebra)
        ));
        assert!(matches!(
            algebra_as_quiver_algebra(&two_by_two_matrix_algebra(Z2)),
            Err(AlgebraAsQuiverAlgebraError::NonBasicAlgebra)
        ));
    }

    #[test]
    fn pinned_qpa_end_of_module_as_quiver_algebra_matches() {
        let fixture: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehom_end_quiver.schema-v1.json"
        ))
        .unwrap();
        assert_eq!(
            end_of_module_json(Z2),
            fixture["binary"]["regular_dual_numbers"]
        );
        assert_eq!(
            end_of_module_json(Zp::new(5)),
            fixture["quinary"]["regular_dual_numbers"]
        );

        let quiver = Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, quiver).unwrap());
        let repeated = right_module_over_path_algebra(algebra, vec![2], &[])
            .unwrap()
            .unwrap();
        assert!(matches!(
            end_of_module_as_quiver_algebra(&repeated),
            Err(EndOfModuleAsQuiverAlgebraError::Conversion(
                AlgebraAsQuiverAlgebraError::NonBasicAlgebra
            ))
        ));

        // QPA reaches the same zero radical here, but its multiplication of
        // the empty radical-basis list aborts before it can return K.  The
        // translated converter handles that zero space explicitly.
        let quiver = Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, quiver).unwrap());
        let simple = right_module_over_path_algebra(algebra, vec![1], &[])
            .unwrap()
            .unwrap();
        let converted = end_of_module_as_quiver_algebra(&simple).unwrap();
        assert_eq!(converted.endomorphism_algebra().dimension(), 1);
        assert_eq!(converted.adjacency_matrix(), &[vec![0].into_boxed_slice()]);
        assert!(matches!(
            converted.quiver_algebra(),
            QuiverAlgebra::Path(algebra) if algebra.canonical_basis().unwrap().len() == 1
        ));
    }
}
