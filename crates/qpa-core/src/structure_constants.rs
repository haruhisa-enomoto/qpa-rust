//! Finite-dimensional algebras given by structure constants.
//!
//! QPA constructs these algebras through GAP's `AlgebraByStructureConstants`.
//! The ordered standard basis and row-vector coefficient conventions are
//! observable in QPA's higher preprojective-algebra implementation.

use std::{
    cmp::Ordering,
    collections::BTreeMap,
    fmt,
    sync::{Arc, OnceLock, RwLock},
};

use numerica::domains::{Field, integer::Integer};
use thiserror::Error;

use crate::{
    AlgebraGlobalDimension, CoefficientDomainError, DenseMatrix, MatrixError, QpaCoefficientField,
    QpaFieldCharacteristic,
};

#[derive(Debug)]
struct StructureConstantsAlgebraData<F: Field> {
    field: F,
    dimension: usize,
    products: Vec<F::Element>,
    named_elements: RwLock<BTreeMap<Box<str>, Vec<F::Element>>>,
    semisimple: OnceLock<bool>,
    global_dimension: OnceLock<AlgebraGlobalDimension>,
    primitive_idempotents: OnceLock<Vec<Vec<F::Element>>>,
}

/// A finite-dimensional associative algebra in its ordered standard basis.
///
/// The product at `(left, right)` is stored as a full coordinate row in the
/// same format produced by GAP's `AlgebraByStructureConstants`.
#[derive(Clone, Debug)]
pub struct StructureConstantsAlgebra<F: Field>(Arc<StructureConstantsAlgebraData<F>>);

#[doc(hidden)]
#[derive(Debug)]
pub struct OppositeStructureConstantsAlgebraData<F: Field> {
    underlying: StructureConstantsAlgebra<F>,
    reversed: StructureConstantsAlgebra<F>,
}

/// QPA's generic `OppositeAlgebra(A)` result for a structure-constant algebra.
///
/// The `Original` variant preserves QPA's identical-object return when
/// commutativity is known. The `Opposite` variant is a fresh algebra object
/// whose elements wrap elements of the retained underlying algebra.
#[derive(Clone, Debug)]
pub enum OppositeStructureConstantsAlgebra<F: Field> {
    /// QPA returned the entered commutative algebra unchanged.
    Original(StructureConstantsAlgebra<F>),
    /// QPA constructed a noncommutative opposite algebra.
    Opposite(Arc<OppositeStructureConstantsAlgebraData<F>>),
}

/// An element of an [`OppositeStructureConstantsAlgebra`] together with the
/// exact underlying algebra element wrapped by QPA.
#[derive(Clone, Debug)]
pub struct OppositeStructureConstantsAlgebraElement<F: Field> {
    algebra: OppositeStructureConstantsAlgebra<F>,
    underlying: StructureConstantsAlgebraElement<F>,
}

impl<F: Field> StructureConstantsAlgebra<F> {
    /// Construct an algebra from a flattened `dimension × dimension ×
    /// dimension` multiplication tensor.
    ///
    /// # Errors
    ///
    /// Rejects an overflowing tensor size or a coefficient count different
    /// from the requested cube.
    pub fn new(
        field: F,
        dimension: usize,
        products: Vec<F::Element>,
    ) -> Result<Self, StructureConstantsError> {
        let expected = dimension
            .checked_mul(dimension)
            .and_then(|square| square.checked_mul(dimension))
            .ok_or(StructureConstantsError::TensorShapeOverflow { dimension })?;
        if products.len() != expected {
            return Err(StructureConstantsError::TensorLength {
                actual: products.len(),
                expected,
                dimension,
            });
        }
        Ok(Self(Arc::new(StructureConstantsAlgebraData {
            field,
            dimension,
            products,
            named_elements: RwLock::new(BTreeMap::new()),
            semisimple: OnceLock::new(),
            global_dimension: OnceLock::new(),
            primitive_idempotents: OnceLock::new(),
        })))
    }

    /// Return the exact coefficient field.
    #[must_use]
    pub fn field(&self) -> &F {
        &self.0.field
    }

    /// Return the dimension of the ordered standard basis.
    #[must_use]
    pub fn dimension(&self) -> usize {
        self.0.dimension
    }

    /// Return QPA's cached generic `IsSemisimpleAlgebra` property, when known.
    #[must_use]
    pub fn cached_is_semisimple_algebra(&self) -> Option<bool> {
        self.0.semisimple.get().copied()
    }

    /// Return QPA's cached `GlobalDimension` attribute, when known.
    #[must_use]
    pub fn cached_global_dimension(&self) -> Option<AlgebraGlobalDimension> {
        self.0.global_dimension.get().copied()
    }

    /// Return QPA's cached `PrimitiveIdempotents` attribute coordinates, when
    /// the attribute has already been evaluated.
    #[must_use]
    pub fn cached_primitive_idempotents(&self) -> Option<&[Vec<F::Element>]> {
        self.0.primitive_idempotents.get().map(Vec::as_slice)
    }

    pub(crate) fn cache_primitive_idempotents(
        &self,
        coordinates: Vec<Vec<F::Element>>,
    ) -> &[Vec<F::Element>] {
        self.0
            .primitive_idempotents
            .get_or_init(|| coordinates)
            .as_slice()
    }

    /// Invoke QPA's immediate `GlobalDimension` method when semisimplicity is
    /// already known true.
    ///
    /// `None` preserves GAP's missing-applicable-method state when the
    /// semisimple property is unknown or known false.
    #[must_use]
    pub fn global_dimension_from_known_properties(&self) -> Option<AlgebraGlobalDimension> {
        if let Some(cached) = self.cached_global_dimension() {
            return Some(cached);
        }
        (self.cached_is_semisimple_algebra() == Some(true)).then(|| {
            *self
                .0
                .global_dimension
                .get_or_init(|| AlgebraGlobalDimension::Finite(0))
        })
    }

    /// Return QPA's true-method implication from known semisimplicity to
    /// `IsFiniteGlobalDimensionAlgebra`.
    ///
    /// `None` means the source installs no conclusion from the currently
    /// known filters; it is not mathematical falsity.
    #[must_use]
    pub fn is_finite_global_dimension_from_known_properties(&self) -> Option<bool> {
        (self.cached_is_semisimple_algebra() == Some(true)).then_some(true)
    }

    /// Return the full coordinate row of a standard-basis product.
    ///
    /// # Panics
    ///
    /// Panics if either basis index is outside the algebra dimension.
    #[must_use]
    pub fn basis_product(&self, left: usize, right: usize) -> &[F::Element] {
        assert!(left < self.dimension(), "left basis index out of bounds");
        assert!(right < self.dimension(), "right basis index out of bounds");
        let start = (left * self.dimension() + right) * self.dimension();
        &self.0.products[start..start + self.dimension()]
    }

    /// Return the additive zero element.
    #[must_use]
    pub fn zero(&self) -> StructureConstantsAlgebraElement<F> {
        StructureConstantsAlgebraElement {
            algebra: self.clone(),
            coordinates: vec![self.field().zero(); self.dimension()],
        }
    }

    /// Return `BasisVectors(Basis(A))` in GAP's standard-basis order.
    #[must_use]
    pub fn basis_elements(&self) -> Vec<StructureConstantsAlgebraElement<F>> {
        (0..self.dimension())
            .map(|index| {
                let mut element = self.zero();
                element.coordinates[index] = self.field().one();
                element
            })
            .collect()
    }

    /// Construct an element from standard-basis coordinates.
    ///
    /// # Errors
    ///
    /// Rejects a coordinate row whose length differs from the dimension.
    pub fn element(
        &self,
        coordinates: Vec<F::Element>,
    ) -> Result<StructureConstantsAlgebraElement<F>, StructureConstantsError> {
        if coordinates.len() != self.dimension() {
            return Err(StructureConstantsError::ElementLength {
                actual: coordinates.len(),
                expected: self.dimension(),
            });
        }
        Ok(StructureConstantsAlgebraElement {
            algebra: self.clone(),
            coordinates,
        })
    }

    /// Bind a named algebra element for QPA's record-component lookup.
    ///
    /// GAP algebra objects may expose generators as dynamic components such as
    /// `A.x`. Rust records that optional surface explicitly and uses it for the
    /// generic opposite-algebra `A.name` method.
    ///
    /// # Errors
    ///
    /// Rejects an element from a different algebra object.
    pub fn bind_named_element(
        &self,
        name: impl Into<Box<str>>,
        element: &StructureConstantsAlgebraElement<F>,
    ) -> Result<(), StructureConstantsError> {
        if element.algebra() != self {
            return Err(StructureConstantsError::DifferentParents);
        }
        self.0
            .named_elements
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .insert(name.into(), element.coordinates().to_vec());
        Ok(())
    }

    /// Return a dynamically named element, when QPA's component is bound.
    ///
    /// # Errors
    ///
    /// Propagates checked reconstruction of the stored coordinate row.
    pub fn element_by_name(
        &self,
        name: &str,
    ) -> Result<Option<StructureConstantsAlgebraElement<F>>, StructureConstantsError> {
        self.0
            .named_elements
            .read()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .get(name)
            .cloned()
            .map(|coordinates| self.element(coordinates))
            .transpose()
    }

    /// Translate QPA's generic `OppositeAlgebra(A)` construction.
    ///
    /// QPA returns `A` itself when commutativity is known. Otherwise it wraps
    /// the entered generators and defines `op(a) * op(b) = op(b * a)`. The
    /// induced structure-constant table therefore swaps its first two tensor
    /// indices while retaining standard-basis order.
    ///
    /// # Errors
    ///
    /// Propagates checked construction of the reversed multiplication tensor.
    pub fn opposite(
        &self,
    ) -> Result<OppositeStructureConstantsAlgebra<F>, StructureConstantsError> {
        if is_commutative(self) {
            return Ok(OppositeStructureConstantsAlgebra::Original(self.clone()));
        }
        let dimension = self.dimension();
        let mut products = Vec::with_capacity(self.0.products.len());
        for left in 0..dimension {
            for right in 0..dimension {
                products.extend_from_slice(self.basis_product(right, left));
            }
        }
        let reversed = StructureConstantsAlgebra::new(self.field().clone(), dimension, products)?;
        Ok(OppositeStructureConstantsAlgebra::Opposite(Arc::new(
            OppositeStructureConstantsAlgebraData {
                underlying: self.clone(),
                reversed,
            },
        )))
    }

    /// Translate QPA's `FindMultiplicativeIdentity(A)` method.
    ///
    /// QPA concatenates the matrices of right multiplication by every basis
    /// vector and solves one row system. Consequently this returns the same
    /// first GAP-compatible left identity, or `None` for GAP's `fail`.
    ///
    /// # Errors
    ///
    /// Propagates checked matrix construction and solution failures.
    pub fn find_multiplicative_identity(
        &self,
    ) -> Result<Option<StructureConstantsAlgebraElement<F>>, StructureConstantsError> {
        let dimension = self.dimension();
        if dimension == 0 {
            return Err(StructureConstantsError::EmptyIdentityBlockMatrix);
        }
        let columns = dimension
            .checked_mul(dimension)
            .ok_or(StructureConstantsError::IdentityMatrixShapeOverflow { dimension })?;
        let mut matrix = Vec::with_capacity(
            dimension
                .checked_mul(columns)
                .ok_or(StructureConstantsError::IdentityMatrixShapeOverflow { dimension })?,
        );
        for unknown in 0..dimension {
            for right in 0..dimension {
                matrix.extend_from_slice(self.basis_product(unknown, right));
            }
        }
        let total_matrix =
            DenseMatrix::from_flat(dimension, columns, matrix, self.field().clone())?;
        let mut target = Vec::with_capacity(columns);
        for basis_index in 0..dimension {
            for coordinate in 0..dimension {
                target.push(if basis_index == coordinate {
                    self.field().one()
                } else {
                    self.field().zero()
                });
            }
        }
        total_matrix
            .solution_row(&target)?
            .map(|coordinates| self.element(coordinates))
            .transpose()
    }

    /// Translate GAP's `CentreFromSCTable` basis computation.
    ///
    /// The returned subalgebra basis consists of the exact left-nullspace rows
    /// of the structure-constant commutator system.
    ///
    /// # Errors
    ///
    /// Propagates exact matrix and subalgebra-coordinate construction.
    pub fn centre(&self) -> Result<StructureConstantsSubalgebra<F>, StructureConstantsError> {
        let dimension = self.dimension();
        let columns = dimension
            .checked_mul(dimension)
            .ok_or(StructureConstantsError::SubspaceMatrixShapeOverflow { dimension })?;
        let equation_length = dimension
            .checked_mul(columns)
            .ok_or(StructureConstantsError::SubspaceMatrixShapeOverflow { dimension })?;
        let mut equations = vec![self.field().zero(); equation_length];
        for unknown in 0..dimension {
            for commuting_with in 0..dimension {
                let offset = unknown * columns + commuting_with * dimension;
                for coordinate in 0..dimension {
                    equations[offset + coordinate] = self.field().sub(
                        &self.basis_product(unknown, commuting_with)[coordinate],
                        &self.basis_product(commuting_with, unknown)[coordinate],
                    );
                }
            }
        }
        let basis = DenseMatrix::from_flat(dimension, columns, equations, self.field().clone())?
            .nullspace_rows();
        StructureConstantsSubalgebra::from_basis(self.clone(), basis)
    }

    pub(crate) fn coordinate_product(
        &self,
        left: &[F::Element],
        right: &[F::Element],
    ) -> Vec<F::Element> {
        debug_assert_eq!(left.len(), self.dimension());
        debug_assert_eq!(right.len(), self.dimension());
        let mut coordinates = vec![self.field().zero(); self.dimension()];
        for (left_index, left_coefficient) in left.iter().enumerate() {
            if self.field().is_zero(left_coefficient) {
                continue;
            }
            for (right_index, right_coefficient) in right.iter().enumerate() {
                if self.field().is_zero(right_coefficient) {
                    continue;
                }
                let scalar = self.field().mul(left_coefficient, right_coefficient);
                for (coordinate, product_coefficient) in coordinates
                    .iter_mut()
                    .zip(self.basis_product(left_index, right_index))
                {
                    self.field()
                        .add_mul_assign(coordinate, &scalar, product_coefficient);
                }
            }
        }
        coordinates
    }

    /// Multiply two entered coordinate rows in the ordered standard basis.
    ///
    /// This is the coordinate-level operation used by QPA algorithms that
    /// temporarily regard ideals and quotient spaces as GAP vector spaces.
    ///
    /// # Errors
    ///
    /// Rejects either row when its width differs from the algebra dimension.
    pub fn multiply_coordinates(
        &self,
        left: &[F::Element],
        right: &[F::Element],
    ) -> Result<Vec<F::Element>, StructureConstantsError> {
        if left.len() != self.dimension() {
            return Err(StructureConstantsError::ElementLength {
                actual: left.len(),
                expected: self.dimension(),
            });
        }
        if right.len() != self.dimension() {
            return Err(StructureConstantsError::ElementLength {
                actual: right.len(),
                expected: self.dimension(),
            });
        }
        Ok(self.coordinate_product(left, right))
    }

    fn left_regular_matrix(
        &self,
        coordinates: &[F::Element],
    ) -> Result<DenseMatrix<F>, StructureConstantsError> {
        if coordinates.len() != self.dimension() {
            return Err(StructureConstantsError::ElementLength {
                actual: coordinates.len(),
                expected: self.dimension(),
            });
        }
        let mut matrix =
            DenseMatrix::zero(self.dimension(), self.dimension(), self.field().clone());
        for right in 0..self.dimension() {
            let mut standard = vec![self.field().zero(); self.dimension()];
            standard[right] = self.field().one();
            let product = self.coordinate_product(coordinates, &standard);
            for (row, coefficient) in product.into_iter().enumerate() {
                matrix[(row, right)] = coefficient;
            }
        }
        Ok(matrix)
    }

    /// Construct GAP's natural quotient by an entered ideal basis.
    ///
    /// QPA/GAP extends the ideal basis by successive standard basis vectors;
    /// those accepted vectors are retained as quotient-basis lifts.
    ///
    /// # Errors
    ///
    /// Rejects an incorrectly shaped or dependent ideal basis and propagates
    /// exact quotient multiplication and coordinate construction.
    pub fn quotient_by_ideal_basis(
        &self,
        ideal_basis: DenseMatrix<F>,
    ) -> Result<StructureConstantsQuotient<F>, StructureConstantsError> {
        StructureConstantsQuotient::new(self.clone(), ideal_basis)
    }
}

impl<F: QpaCoefficientField> StructureConstantsAlgebra<F> {
    /// Translate QPA's generic finite-dimensional `IsSemisimpleAlgebra`
    /// method as the zero-dimensionality test for `RadicalOfAlgebra(A)`.
    ///
    /// The result is cached as a GAP property value on the identical algebra
    /// object.
    ///
    /// # Errors
    ///
    /// Propagates the source-translated radical computation.
    pub fn is_semisimple_algebra(&self) -> Result<bool, StructureConstantsError> {
        if let Some(cached) = self.cached_is_semisimple_algebra() {
            return Ok(cached);
        }
        let result = self.radical()?.basis().nrows() == 0;
        let result = *self.0.semisimple.get_or_init(|| result);
        if result {
            let _ = self
                .0
                .global_dimension
                .get_or_init(|| AlgebraGlobalDimension::Finite(0));
        }
        Ok(result)
    }

    /// Translate QPA's `PowerSubalgebra(A, n)` on the currently supported
    /// prime-field contract.
    ///
    /// QPA first passes a noncommutative associative algebra to its centre.
    /// For a commutative algebra it forms the rows `v - v^(p^n)` in canonical
    /// basis order and returns the left-nullspace combinations. Numerica's
    /// initial QPA coefficient contract has prime-field degree one, so the
    /// source divisibility test accepts exactly `n = 1`.
    ///
    /// # Errors
    ///
    /// Requires a positive divisor of the supported field degree, a finite
    /// prime coefficient field, and the multiplicative identity used for
    /// exponentiation. Propagates exact centre and matrix failures.
    pub fn power_subalgebra(
        &self,
        degree: usize,
    ) -> Result<StructureConstantsSubalgebra<F>, StructureConstantsError> {
        const PRIME_FIELD_DEGREE: usize = 1;
        if degree == 0 {
            return Err(StructureConstantsError::PowerSubalgebraDegreeNotPositive);
        }
        if !PRIME_FIELD_DEGREE.is_multiple_of(degree) {
            return Err(StructureConstantsError::PowerSubalgebraDegreeNotDivisor {
                degree,
                field_degree: PRIME_FIELD_DEGREE,
            });
        }
        let prime = match self.field().qpa_characteristic() {
            QpaFieldCharacteristic::Prime(prime) => prime,
            QpaFieldCharacteristic::Zero => {
                return Err(StructureConstantsError::FiniteFieldRequired);
            }
        };

        // Calling `Centre` unconditionally is equivalent on a commutative
        // structure-constant algebra: its GAP canonical basis is the standard
        // basis, which is also the zero-commutator nullspace basis retained by
        // `centre()`.
        let centre = self.centre()?;
        let centre_algebra = centre.algebra();
        let identity = centre_algebra
            .find_multiplicative_identity()?
            .ok_or(StructureConstantsError::MissingMultiplicativeIdentity)?;
        let dimension = centre_algebra.dimension();
        let mut gamma = Vec::with_capacity(dimension.saturating_mul(dimension));
        for basis_index in 0..dimension {
            let mut basis_vector = vec![centre_algebra.field().zero(); dimension];
            basis_vector[basis_index] = centre_algebra.field().one();
            let power = algebra_power(centre_algebra, &basis_vector, prime, identity.coordinates());
            gamma.extend(subtract_rows(centre_algebra.field(), &basis_vector, &power));
        }
        let fixed_coordinates =
            DenseMatrix::from_flat(dimension, dimension, gamma, centre_algebra.field().clone())?
                .nullspace_rows();
        let ambient_basis = fixed_coordinates.multiply(centre.basis())?;
        StructureConstantsSubalgebra::from_basis(self.clone(), ambient_basis)
    }

    /// Translate GAP's `RadicalOfAlgebra` through its faithful left-regular
    /// matrix representation.
    ///
    /// In characteristic zero this is GAP's trace-form nullspace. In positive
    /// characteristic it follows the Cohen--Ivanyos--Wales ideal sequence from
    /// `algmat.gi`, including GAP's spinning characteristic polynomials and
    /// left-nullspace basis order. Prime-field Frobenius is the identity.
    ///
    /// # Errors
    ///
    /// Requires the multiplicative identity used by QPA's faithful adjoint
    /// representation and propagates exact matrix/subalgebra failures.
    pub fn radical(&self) -> Result<StructureConstantsSubalgebra<F>, StructureConstantsError> {
        if self.dimension() == 0 {
            return StructureConstantsSubalgebra::from_basis(
                self.clone(),
                DenseMatrix::zero(0, 0, self.field().clone()),
            );
        }
        let identity = self
            .find_multiplicative_identity()?
            .ok_or(StructureConstantsError::MissingMultiplicativeIdentity)?;
        match self.field().qpa_characteristic() {
            QpaFieldCharacteristic::Zero => self.radical_characteristic_zero(),
            QpaFieldCharacteristic::Prime(prime) => {
                self.radical_positive_characteristic(identity.coordinates(), prime)
            }
        }
    }

    /// Translate QPA's `RadicalSeriesOfAlgebra` basis sequence.
    ///
    /// The first matrix is the standard basis of `A`, followed by bases of
    /// `rad(A)`, `rad(A)^2`, and so on, including the final zero-dimensional
    /// space. Each product basis follows GAP's `ProductSpace` generator order
    /// and generic semi-echelon basis selection.
    ///
    /// # Errors
    ///
    /// Propagates radical, structure multiplication, and checked matrix
    /// construction failures.
    pub fn radical_series(&self) -> Result<Vec<DenseMatrix<F>>, StructureConstantsError> {
        let radical = self.radical()?;
        let mut series = vec![DenseMatrix::identity(
            self.dimension(),
            self.field().clone(),
        )];
        let mut power = radical.basis().clone();
        while power.nrows() != 0 {
            series.push(power.clone());
            power = self.product_space_basis(&power, radical.basis())?;
        }
        series.push(power);
        Ok(series)
    }

    /// Translate QPA's `IsRadicalSquareZeroAlgebra` method.
    ///
    /// # Errors
    ///
    /// Propagates the exact radical and `ProductSpace` computations.
    pub fn is_radical_square_zero(&self) -> Result<bool, StructureConstantsError> {
        let radical = self.radical()?;
        Ok(self
            .product_space_basis(radical.basis(), radical.basis())?
            .nrows()
            == 0)
    }

    /// Translate QPA's finite-field `IsBasicAlgebra` method on the supported
    /// split-prime-field contract.
    ///
    /// `None` is QPA's `TryNextMethod()` on an infinite coefficient field.
    /// Over a finite field, each central simple block of `A/rad(A)` is basic
    /// exactly when it is commutative, since finite division rings are fields.
    ///
    /// # Errors
    ///
    /// Propagates radical, quotient, central-idempotent, and exact subalgebra
    /// construction failures.
    pub fn is_basic_algebra(&self) -> Result<Option<bool>, StructureConstantsError> {
        if self.field().qpa_elements()?.is_none() {
            return Ok(None);
        }
        let components = self.semisimple_central_components()?;
        Ok(Some(
            components
                .iter()
                .all(|component| is_commutative(component.algebra())),
        ))
    }

    /// Translate QPA's finite-field `IsElementaryAlgebra` method on the
    /// supported split-prime-field contract.
    ///
    /// `None` preserves either `TryNextMethod()` branch: an infinite field or
    /// a nonbasic algebra. Otherwise QPA compares the dimensions of all simple
    /// division-field blocks of `A/rad(A)`.
    ///
    /// # Errors
    ///
    /// Propagates the basicness and semisimple-component computations.
    pub fn is_elementary_algebra(&self) -> Result<Option<bool>, StructureConstantsError> {
        if self.is_basic_algebra()? != Some(true) {
            return Ok(None);
        }
        let components = self.semisimple_central_components()?;
        let Some(first) = components.first() else {
            return Ok(Some(true));
        };
        Ok(Some(components.iter().all(|component| {
            component.algebra().dimension() == first.algebra().dimension()
        })))
    }

    /// Translate QPA's `BlocksOfAlgebra` on the supported finite prime-field
    /// contract.
    ///
    /// For each central idempotent in QPA order, the block basis is obtained
    /// from `BasisVectors(Basis(A)) * c` and the resulting GAP semi-echelon
    /// span. The induced algebra retains its ambient basis coordinates.
    ///
    /// # Errors
    ///
    /// Propagates central-idempotent, exact product, and subalgebra failures.
    pub fn blocks(&self) -> Result<Vec<StructureConstantsSubalgebra<F>>, StructureConstantsError> {
        let idempotents = self.central_idempotents_split_prime()?;
        idempotents
            .iter()
            .map(|idempotent| {
                let mut generators = Vec::with_capacity(self.dimension() + 1);
                // After constructing the FLMLOR, QPA sets its `One` before
                // asking for a basis. GAP's associative-with-one closure
                // therefore seeds its mutable basis with this idempotent.
                generators.push(idempotent.clone());
                generators.extend(
                    (0..self.dimension())
                        .map(|basis_index| {
                            let mut basis = vec![self.field().zero(); self.dimension()];
                            basis[basis_index] = self.field().one();
                            self.multiply_coordinates(&basis, idempotent)
                        })
                        .collect::<Result<Vec<_>, _>>()?,
                );
                let basis = span_rows(self.field(), self.dimension(), generators)?;
                StructureConstantsSubalgebra::from_basis(self.clone(), basis)
            })
            .collect()
    }

    fn semisimple_central_components(
        &self,
    ) -> Result<Vec<StructureConstantsSubalgebra<F>>, StructureConstantsError> {
        let radical = self.radical()?;
        let quotient = self.quotient_by_ideal_basis(radical.basis().clone())?;
        let idempotents = quotient.algebra().central_idempotents_split_prime()?;
        idempotents
            .iter()
            .map(|idempotent| {
                let generators = (0..quotient.algebra().dimension())
                    .map(|basis_index| {
                        let mut basis = vec![self.field().zero(); quotient.algebra().dimension()];
                        basis[basis_index] = self.field().one();
                        quotient.algebra().multiply_coordinates(idempotent, &basis)
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                let basis = span_rows(self.field(), quotient.algebra().dimension(), generators)?;
                StructureConstantsSubalgebra::from_basis(quotient.algebra().clone(), basis)
            })
            .collect()
    }

    fn product_space_basis(
        &self,
        left: &DenseMatrix<F>,
        right: &DenseMatrix<F>,
    ) -> Result<DenseMatrix<F>, StructureConstantsError> {
        if left.ncols() != self.dimension() {
            return Err(StructureConstantsError::SubalgebraBasisWidth {
                actual: left.ncols(),
                expected: self.dimension(),
            });
        }
        if right.ncols() != self.dimension() {
            return Err(StructureConstantsError::SubalgebraBasisWidth {
                actual: right.ncols(),
                expected: self.dimension(),
            });
        }
        let mut products = Vec::with_capacity(left.nrows().saturating_mul(right.nrows()));
        for left_row in 0..left.nrows() {
            for right_row in 0..right.nrows() {
                products.push(self.multiply_coordinates(left.row(left_row), right.row(right_row))?);
            }
        }
        Ok(span_rows(self.field(), self.dimension(), products)?)
    }

    fn radical_characteristic_zero(
        &self,
    ) -> Result<StructureConstantsSubalgebra<F>, StructureConstantsError> {
        let dimension = self.dimension();
        let basis = self.basis_elements();
        let matrices = basis
            .iter()
            .map(|element| self.left_regular_matrix(element.coordinates()))
            .collect::<Result<Vec<_>, _>>()?;
        let mut gram = DenseMatrix::zero(dimension, dimension, self.field().clone());
        for left in 0..dimension {
            for right in left..dimension {
                let product = matrices[left].multiply(&matrices[right])?;
                let mut trace = self.field().zero();
                for diagonal in 0..dimension {
                    self.field()
                        .add_assign(&mut trace, &product[(diagonal, diagonal)]);
                }
                gram[(left, right)] = trace.clone();
                gram[(right, left)] = trace;
            }
        }
        StructureConstantsSubalgebra::from_basis(self.clone(), gram.nullspace_rows())
    }

    fn radical_positive_characteristic(
        &self,
        identity: &[F::Element],
        prime: u64,
    ) -> Result<StructureConstantsSubalgebra<F>, StructureConstantsError> {
        let dimension = self.dimension();
        let mut ideal_basis = DenseMatrix::identity(dimension, self.field().clone());
        let identity_matrix = DenseMatrix::identity(dimension, self.field().clone());
        let mut characteristic_polynomials: Vec<Vec<Vec<F::Element>>> = Vec::new();
        let mut changed = true;
        let mut prime_power = 1_u64;
        let dimension_u64 = u64::try_from(dimension)
            .map_err(|_| StructureConstantsError::SubspaceMatrixShapeOverflow { dimension })?;

        while prime_power <= dimension_u64 {
            let mut extended_basis = ideal_basis.clone();
            if ideal_basis.solution_row(identity)?.is_none() {
                let mut entries = ideal_basis.as_flat().to_vec();
                entries.extend_from_slice(identity);
                extended_basis = DenseMatrix::from_flat(
                    ideal_basis.nrows() + 1,
                    dimension,
                    entries,
                    self.field().clone(),
                )?;
            }
            if changed {
                let ideal_matrices = (0..ideal_basis.nrows())
                    .map(|row| self.left_regular_matrix(ideal_basis.row(row)))
                    .collect::<Result<Vec<_>, _>>()?;
                let extended_matrices = (0..extended_basis.nrows())
                    .map(|row| {
                        if extended_basis.row(row) == identity {
                            Ok(identity_matrix.clone())
                        } else {
                            self.left_regular_matrix(extended_basis.row(row))
                        }
                    })
                    .collect::<Result<Vec<_>, StructureConstantsError>>()?;
                characteristic_polynomials = ideal_matrices
                    .iter()
                    .map(|left| {
                        extended_matrices
                            .iter()
                            .map(|right| {
                                Ok(left
                                    .multiply(right)?
                                    .characteristic_polynomial_coefficients()?)
                            })
                            .collect::<Result<Vec<_>, StructureConstantsError>>()
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                changed = false;
            }

            let coefficient_index = dimension
                .checked_sub(usize::try_from(prime_power).map_err(|_| {
                    StructureConstantsError::SubspaceMatrixShapeOverflow { dimension }
                })?)
                .ok_or(StructureConstantsError::RadicalCoefficientIndex {
                    dimension,
                    prime_power,
                })?;
            let sign = if prime_power.is_multiple_of(2) {
                self.field().one()
            } else {
                self.field().neg(&self.field().one())
            };
            let equations = DenseMatrix::from_flat(
                ideal_basis.nrows(),
                extended_basis.nrows(),
                characteristic_polynomials
                    .iter()
                    .flat_map(|row| {
                        row.iter().map(|polynomial| {
                            self.field().mul(&sign, &polynomial[coefficient_index])
                        })
                    })
                    .collect(),
                self.field().clone(),
            )?;
            let nullspace = equations.nullspace_rows();
            if nullspace.nrows() == 0 {
                return StructureConstantsSubalgebra::from_basis(
                    self.clone(),
                    DenseMatrix::zero(0, dimension, self.field().clone()),
                );
            }
            if nullspace.nrows() != ideal_basis.nrows() {
                ideal_basis = nullspace.multiply(&ideal_basis)?;
                changed = true;
            }
            let Some(next) = prime_power.checked_mul(prime) else {
                break;
            };
            prime_power = next;
        }
        StructureConstantsSubalgebra::from_basis(self.clone(), ideal_basis)
    }

    /// Translate GAP's `CentralIdempotentsOfAlgebra` on the split prime-field
    /// domain used by QPA's `AlgebraAsQuiverAlgebra` method.
    ///
    /// This follows GAP's centre/radical quotient, Frobenius-fixed small-field
    /// splitting, large-field splitting-element branch, polynomial CRT
    /// idempotents, radical lifting, and final sorted-list order. The current
    /// field contract requires every splitting minimal polynomial to factor
    /// linearly; nonsplit extension-field components are reported explicitly.
    ///
    /// # Errors
    ///
    /// Requires a finite prime field and a multiplicative identity, and
    /// propagates centre, radical, quotient, polynomial, and lifting failures.
    pub fn central_idempotents_split_prime(
        &self,
    ) -> Result<Vec<Vec<F::Element>>, StructureConstantsError> {
        let elements = self
            .field()
            .qpa_elements()?
            .ok_or(StructureConstantsError::FiniteFieldRequired)?;
        let centre = self.centre()?;
        let centre_radical = centre.algebra().radical()?;
        let quotient = centre
            .algebra()
            .quotient_by_ideal_basis(centre_radical.basis().clone())?;
        let quotient_identity = quotient
            .algebra()
            .find_multiplicative_identity()?
            .ok_or(StructureConstantsError::MissingMultiplicativeIdentity)?;
        let mut idempotents = split_commutative_semisimple_algebra(
            quotient.algebra(),
            quotient_identity.coordinates(),
            &elements,
        )?;

        for idempotent in &mut idempotents {
            *idempotent = quotient.lift_coordinates(idempotent)?;
        }
        lift_orthogonal_idempotents_across_radical(
            centre.algebra(),
            centre_radical.basis().nrows(),
            &mut idempotents,
        );
        let element_order = elements;
        idempotents
            .sort_by(|left, right| compare_field_rows_by_gap_order(left, right, &element_order));
        idempotents
            .into_iter()
            .map(|idempotent| centre.lift_coordinates(&idempotent))
            .collect()
    }
}

#[derive(Clone, Debug)]
struct CommutativeIdeal<F: Field> {
    basis: DenseMatrix<F>,
    identity: Vec<F::Element>,
}

fn split_commutative_semisimple_algebra<F: QpaCoefficientField>(
    algebra: &StructureConstantsAlgebra<F>,
    identity: &[F::Element],
    elements: &[F::Element],
) -> Result<Vec<Vec<F::Element>>, StructureConstantsError> {
    let dimension = algebra.dimension();
    if dimension == 0 {
        return Ok(Vec::new());
    }
    if dimension == 1 {
        return Ok(vec![identity.to_vec()]);
    }
    let threshold = dimension
        .checked_mul(dimension)
        .and_then(|square| square.checked_mul(2))
        .ok_or(StructureConstantsError::SubspaceMatrixShapeOverflow { dimension })?;
    if elements.len() > threshold {
        let splitting = deterministic_gap_splitting_element(algebra, identity, threshold)?;
        return lagrange_idempotents(algebra, identity, &splitting, elements);
    }

    let field = algebra.field();
    let mut ideals = vec![CommutativeIdeal {
        basis: DenseMatrix::identity(dimension, field.clone()),
        identity: identity.to_vec(),
    }];
    let mut idempotents = vec![identity.to_vec()];
    let mut position = 0usize;
    let prime = match field.qpa_characteristic() {
        QpaFieldCharacteristic::Prime(prime) => prime,
        QpaFieldCharacteristic::Zero => return Err(StructureConstantsError::FiniteFieldRequired),
    };

    while position < ideals.len() {
        if ideals[position].basis.nrows() == 1 {
            position += 1;
            continue;
        }
        let fixed = frobenius_fixed_subspace(algebra, &ideals[position], prime)?;
        if fixed.nrows() == 1 {
            position += 1;
            continue;
        }
        let first = fixed.row(0).to_vec();
        let identity_span = DenseMatrix::from_flat(
            1,
            dimension,
            ideals[position].identity.clone(),
            field.clone(),
        )?;
        let splitting = if identity_span.solution_row(&first)?.is_some() {
            fixed.row(1).to_vec()
        } else {
            first
        };
        let pieces =
            lagrange_idempotents(algebra, &ideals[position].identity, &splitting, elements)?;
        if pieces.len() < 2 {
            return Err(StructureConstantsError::FrobeniusElementDidNotSplit);
        }
        let old_basis = ideals[position].basis.clone();
        ideals.remove(position);
        idempotents.remove(position);
        for piece in pieces {
            let basis = idempotent_ideal_basis(algebra, &old_basis, &piece)?;
            ideals.push(CommutativeIdeal {
                basis,
                identity: piece.clone(),
            });
            idempotents.push(piece);
        }
    }
    Ok(idempotents)
}

fn deterministic_gap_splitting_element<F: QpaCoefficientField>(
    algebra: &StructureConstantsAlgebra<F>,
    identity: &[F::Element],
    threshold: usize,
) -> Result<Vec<F::Element>, StructureConstantsError> {
    let coefficient_count = threshold
        .checked_add(1)
        .ok_or(StructureConstantsError::SplittingSearchOverflow)?;
    let coefficients = (0..coefficient_count)
        .map(|value| {
            let integer = i64::try_from(value)
                .map(Integer::from)
                .map_err(|_| StructureConstantsError::SplittingSearchOverflow)?;
            Ok(algebra.field().nth(integer))
        })
        .collect::<Result<Vec<_>, StructureConstantsError>>()?;
    let mut digits = vec![0usize; algebra.dimension()];
    loop {
        let candidate = digits
            .iter()
            .map(|&digit| coefficients[digit].clone())
            .collect::<Vec<_>>();
        if minimum_polynomial(algebra, identity, &candidate)?.len() - 1 == algebra.dimension() {
            return Ok(candidate);
        }
        let mut coordinate = 0usize;
        while coordinate < digits.len() {
            digits[coordinate] += 1;
            if digits[coordinate] < coefficients.len() {
                break;
            }
            digits[coordinate] = 0;
            coordinate += 1;
        }
        if coordinate == digits.len() {
            return Err(StructureConstantsError::MissingSplittingElement);
        }
    }
}

fn frobenius_fixed_subspace<F: Field>(
    algebra: &StructureConstantsAlgebra<F>,
    ideal: &CommutativeIdeal<F>,
    field_size: u64,
) -> Result<DenseMatrix<F>, StructureConstantsError> {
    let mut equations =
        Vec::with_capacity(ideal.basis.nrows().checked_mul(ideal.basis.nrows()).ok_or(
            StructureConstantsError::SubspaceMatrixShapeOverflow {
                dimension: algebra.dimension(),
            },
        )?);
    for row in 0..ideal.basis.nrows() {
        let basis_vector = ideal.basis.row(row);
        let power = algebra_power(algebra, basis_vector, field_size, &ideal.identity);
        let difference = subtract_rows(algebra.field(), &power, basis_vector);
        let coordinates = ideal
            .basis
            .solution_row(&difference)?
            .ok_or(StructureConstantsError::IdealBasisCoordinateFailure)?;
        equations.extend(coordinates);
    }
    let fixed_coordinates = DenseMatrix::from_flat(
        ideal.basis.nrows(),
        ideal.basis.nrows(),
        equations,
        algebra.field().clone(),
    )?
    .nullspace_rows();
    fixed_coordinates.multiply(&ideal.basis).map_err(Into::into)
}

fn idempotent_ideal_basis<F: Field>(
    algebra: &StructureConstantsAlgebra<F>,
    old_basis: &DenseMatrix<F>,
    idempotent: &[F::Element],
) -> Result<DenseMatrix<F>, StructureConstantsError> {
    let mut generators = Vec::new();
    for row in 0..old_basis.nrows() {
        let left = algebra.coordinate_product(idempotent, old_basis.row(row));
        let product = algebra.coordinate_product(&left, idempotent);
        if product
            .iter()
            .any(|coefficient| !algebra.field().is_zero(coefficient))
        {
            generators.extend(product);
        }
    }
    if generators.is_empty() {
        return Ok(DenseMatrix::zero(
            0,
            algebra.dimension(),
            algebra.field().clone(),
        ));
    }
    let rows = generators.len() / algebra.dimension();
    Ok(DenseMatrix::from_flat(
        rows,
        algebra.dimension(),
        generators,
        algebra.field().clone(),
    )?
    .semi_echelon_transformation()
    .vectors)
}

pub(crate) fn minimum_polynomial<F: Field>(
    algebra: &StructureConstantsAlgebra<F>,
    identity: &[F::Element],
    element: &[F::Element],
) -> Result<Vec<F::Element>, StructureConstantsError> {
    let mut powers = vec![identity.to_vec()];
    let mut power = element.to_vec();
    loop {
        let basis = DenseMatrix::from_flat(
            powers.len(),
            algebra.dimension(),
            powers.iter().flatten().cloned().collect(),
            algebra.field().clone(),
        )?;
        if let Some(coordinates) = basis.solution_row(&power)? {
            let mut polynomial = coordinates
                .into_iter()
                .map(|coefficient| algebra.field().neg(&coefficient))
                .collect::<Vec<_>>();
            polynomial.push(algebra.field().one());
            return Ok(polynomial);
        }
        powers.push(power.clone());
        power = algebra.coordinate_product(&power, element);
    }
}

fn lagrange_idempotents<F: Field>(
    algebra: &StructureConstantsAlgebra<F>,
    identity: &[F::Element],
    element: &[F::Element],
    field_elements: &[F::Element],
) -> Result<Vec<Vec<F::Element>>, StructureConstantsError> {
    let polynomial = minimum_polynomial(algebra, identity, element)?;
    let degree = polynomial.len() - 1;
    let roots = field_elements
        .iter()
        .filter(|root| {
            algebra.field().is_zero(&evaluate_scalar_polynomial(
                algebra.field(),
                &polynomial,
                root,
            ))
        })
        .cloned()
        .collect::<Vec<_>>();
    if roots.len() != degree {
        return Err(StructureConstantsError::NonsplitMinimalPolynomial {
            degree,
            linear_factors: roots.len(),
        });
    }
    let mut idempotents = Vec::with_capacity(roots.len());
    for (position, root) in roots.iter().enumerate() {
        let mut indicator = vec![algebra.field().one()];
        let mut denominator = algebra.field().one();
        for (other_position, other) in roots.iter().enumerate() {
            if position == other_position {
                continue;
            }
            indicator = multiply_scalar_polynomials(
                algebra.field(),
                &indicator,
                &[algebra.field().neg(other), algebra.field().one()],
            );
            algebra
                .field()
                .mul_assign(&mut denominator, &algebra.field().sub(root, other));
        }
        let inverse = algebra.field().inv(&denominator);
        for coefficient in &mut indicator {
            algebra.field().mul_assign(coefficient, &inverse);
        }
        idempotents.push(evaluate_algebra_polynomial(
            algebra, &indicator, element, identity,
        ));
    }
    Ok(idempotents)
}

fn evaluate_scalar_polynomial<F: Field>(
    field: &F,
    polynomial: &[F::Element],
    value: &F::Element,
) -> F::Element {
    let mut result = field.zero();
    for coefficient in polynomial.iter().rev() {
        result = field.add(&field.mul(&result, value), coefficient);
    }
    result
}

pub(crate) fn evaluate_algebra_polynomial<F: Field>(
    algebra: &StructureConstantsAlgebra<F>,
    polynomial: &[F::Element],
    element: &[F::Element],
    identity: &[F::Element],
) -> Vec<F::Element> {
    let mut result = vec![algebra.field().zero(); algebra.dimension()];
    for coefficient in polynomial.iter().rev() {
        result = algebra.coordinate_product(&result, element);
        for (entry, identity_entry) in result.iter_mut().zip(identity) {
            algebra
                .field()
                .add_mul_assign(entry, coefficient, identity_entry);
        }
    }
    result
}

fn multiply_scalar_polynomials<F: Field>(
    field: &F,
    left: &[F::Element],
    right: &[F::Element],
) -> Vec<F::Element> {
    let mut result = vec![field.zero(); left.len() + right.len() - 1];
    for (left_degree, left_coefficient) in left.iter().enumerate() {
        for (right_degree, right_coefficient) in right.iter().enumerate() {
            field.add_mul_assign(
                &mut result[left_degree + right_degree],
                left_coefficient,
                right_coefficient,
            );
        }
    }
    result
}

fn algebra_power<F: Field>(
    algebra: &StructureConstantsAlgebra<F>,
    element: &[F::Element],
    exponent: u64,
    identity: &[F::Element],
) -> Vec<F::Element> {
    let mut result = identity.to_vec();
    let mut base = element.to_vec();
    let mut exponent = exponent;
    while exponent > 0 {
        if !exponent.is_multiple_of(2) {
            result = algebra.coordinate_product(&result, &base);
        }
        exponent /= 2;
        if exponent > 0 {
            base = algebra.coordinate_product(&base, &base);
        }
    }
    result
}

pub(crate) fn subtract_rows<F: Field>(
    field: &F,
    left: &[F::Element],
    right: &[F::Element],
) -> Vec<F::Element> {
    left.iter()
        .zip(right)
        .map(|(left, right)| field.sub(left, right))
        .collect()
}

pub(crate) fn add_rows<F: Field>(
    field: &F,
    left: &[F::Element],
    right: &[F::Element],
) -> Vec<F::Element> {
    left.iter()
        .zip(right)
        .map(|(left, right)| field.add(left, right))
        .collect()
}

fn lift_orthogonal_idempotents_across_radical<F: Field>(
    algebra: &StructureConstantsAlgebra<F>,
    radical_dimension: usize,
    idempotents: &mut [Vec<F::Element>],
) {
    if idempotents.is_empty() {
        return;
    }
    let mut correction_steps = 0usize;
    let mut power = 1usize;
    while power <= radical_dimension {
        correction_steps += 1;
        let Some(next) = power.checked_mul(2) else {
            break;
        };
        power = next;
    }
    let mut accumulated = vec![algebra.field().zero(); algebra.dimension()];
    let two = algebra.field().nth(Integer::from(2));
    for idempotent in idempotents {
        let accumulated_left = algebra.coordinate_product(&accumulated, idempotent);
        let accumulated_right = algebra.coordinate_product(idempotent, &accumulated);
        let sandwich = algebra.coordinate_product(&accumulated_left, &accumulated);
        let mut lifted = subtract_rows(algebra.field(), idempotent, &accumulated_left);
        lifted = subtract_rows(algebra.field(), &lifted, &accumulated_right);
        lifted = add_rows(algebra.field(), &lifted, &sandwich);
        let square = algebra.coordinate_product(&lifted, &lifted);
        let mut defect = subtract_rows(algebra.field(), &square, &lifted);
        for _ in 0..correction_steps {
            let lifted_defect = algebra.coordinate_product(&lifted, &defect);
            let mut next = add_rows(algebra.field(), &lifted, &defect);
            for (entry, correction) in next.iter_mut().zip(lifted_defect) {
                algebra.field().sub_mul_assign(entry, &two, &correction);
            }
            lifted = next;
            let square = algebra.coordinate_product(&lifted, &lifted);
            defect = subtract_rows(algebra.field(), &square, &lifted);
        }
        idempotent.clone_from(&lifted);
        accumulated = add_rows(algebra.field(), &accumulated, &lifted);
    }
}

fn compare_field_rows_by_gap_order<E: Eq>(
    left: &[E],
    right: &[E],
    element_order: &[E],
) -> Ordering {
    for (left, right) in left.iter().zip(right) {
        let left_position = element_order.iter().position(|element| element == left);
        let right_position = element_order.iter().position(|element| element == right);
        match left_position.cmp(&right_position) {
            Ordering::Equal => {}
            ordering => return ordering,
        }
    }
    left.len().cmp(&right.len())
}

/// A subalgebra with an explicitly retained ordered basis in its ambient
/// structure-constant algebra.
#[derive(Clone, Debug)]
pub struct StructureConstantsSubalgebra<F: Field> {
    ambient: StructureConstantsAlgebra<F>,
    basis: DenseMatrix<F>,
    algebra: StructureConstantsAlgebra<F>,
}

/// GAP's natural quotient of a structure-constant algebra by an ideal basis.
#[derive(Clone, Debug)]
pub struct StructureConstantsQuotient<F: Field> {
    original: StructureConstantsAlgebra<F>,
    ideal_basis: DenseMatrix<F>,
    projection: DenseMatrix<F>,
    lifts: DenseMatrix<F>,
    algebra: StructureConstantsAlgebra<F>,
}

impl<F: Field> StructureConstantsQuotient<F> {
    fn new(
        original: StructureConstantsAlgebra<F>,
        ideal_basis: DenseMatrix<F>,
    ) -> Result<Self, StructureConstantsError> {
        let dimension = original.dimension();
        if ideal_basis.ncols() != dimension {
            return Err(StructureConstantsError::IdealBasisWidth {
                actual: ideal_basis.ncols(),
                expected: dimension,
            });
        }
        let kernel_basis = ideal_basis.semi_echelon_transformation().vectors;
        if kernel_basis.nrows() != ideal_basis.nrows() {
            return Err(StructureConstantsError::DependentIdealBasis);
        }
        let mut mutable_basis = kernel_basis.clone();
        let mut complement = Vec::new();
        for coordinate in 0..dimension {
            let mut standard = vec![original.field().zero(); dimension];
            standard[coordinate] = original.field().one();
            let mut candidate = mutable_basis.as_flat().to_vec();
            candidate.extend_from_slice(&standard);
            let candidate = DenseMatrix::from_flat(
                mutable_basis.nrows() + 1,
                dimension,
                candidate,
                original.field().clone(),
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
            original.field().clone(),
        )?;
        let mut full_basis_entries = kernel_basis.as_flat().to_vec();
        full_basis_entries.extend(complement.iter().flatten().cloned());
        let full_basis = DenseMatrix::from_flat(
            dimension,
            dimension,
            full_basis_entries,
            original.field().clone(),
        )?;
        let mut projection_entries = Vec::with_capacity(
            dimension
                .checked_mul(quotient_dimension)
                .ok_or(StructureConstantsError::SubspaceMatrixShapeOverflow { dimension })?,
        );
        for coordinate in 0..dimension {
            let mut standard = vec![original.field().zero(); dimension];
            standard[coordinate] = original.field().one();
            let coefficients = full_basis
                .solution_row(&standard)?
                .ok_or(StructureConstantsError::IncompleteNaturalQuotientBasis)?;
            projection_entries
                .extend_from_slice(&coefficients[dimension - quotient_dimension..dimension]);
        }
        let projection = DenseMatrix::from_flat(
            dimension,
            quotient_dimension,
            projection_entries,
            original.field().clone(),
        )?;

        let tensor_length = quotient_dimension
            .checked_mul(quotient_dimension)
            .and_then(|square| square.checked_mul(quotient_dimension))
            .ok_or(StructureConstantsError::TensorShapeOverflow {
                dimension: quotient_dimension,
            })?;
        let mut products = Vec::with_capacity(tensor_length);
        for left in 0..quotient_dimension {
            for right in 0..quotient_dimension {
                let ambient_product =
                    original.coordinate_product(lifts.row(left), lifts.row(right));
                products.extend(multiply_row_matrix(&ambient_product, &projection));
            }
        }
        let algebra =
            StructureConstantsAlgebra::new(original.field().clone(), quotient_dimension, products)?;
        Ok(Self {
            original,
            ideal_basis,
            projection,
            lifts,
            algebra,
        })
    }

    /// Return the original algebra.
    #[must_use]
    pub const fn original(&self) -> &StructureConstantsAlgebra<F> {
        &self.original
    }

    /// Return the entered ideal basis in original coordinates.
    #[must_use]
    pub const fn ideal_basis(&self) -> &DenseMatrix<F> {
        &self.ideal_basis
    }

    /// Return the quotient's induced structure-constant algebra.
    #[must_use]
    pub const fn algebra(&self) -> &StructureConstantsAlgebra<F> {
        &self.algebra
    }

    /// Return standard quotient-basis lifts in original coordinates.
    #[must_use]
    pub const fn lifts(&self) -> &DenseMatrix<F> {
        &self.lifts
    }

    /// Project an original coordinate row to the natural quotient.
    ///
    /// # Errors
    ///
    /// Rejects a row with the wrong original dimension.
    pub fn project_coordinates(
        &self,
        coordinates: &[F::Element],
    ) -> Result<Vec<F::Element>, StructureConstantsError> {
        if coordinates.len() != self.original.dimension() {
            return Err(StructureConstantsError::ElementLength {
                actual: coordinates.len(),
                expected: self.original.dimension(),
            });
        }
        Ok(multiply_row_matrix(coordinates, &self.projection))
    }

    /// Return QPA/GAP's stored representative of quotient coordinates.
    ///
    /// # Errors
    ///
    /// Rejects a row with the wrong quotient dimension.
    pub fn lift_coordinates(
        &self,
        coordinates: &[F::Element],
    ) -> Result<Vec<F::Element>, StructureConstantsError> {
        if coordinates.len() != self.algebra.dimension() {
            return Err(StructureConstantsError::ElementLength {
                actual: coordinates.len(),
                expected: self.algebra.dimension(),
            });
        }
        Ok(multiply_row_matrix(coordinates, &self.lifts))
    }

    /// Translate QPA's `LiftingCompleteSetOfOrthogonalIdempotents` for this
    /// natural quotient map.
    ///
    /// QPA first applies `LiftIdempotent` to the first quotient idempotent and
    /// then calls `LiftTwoOrthogonalIdempotents` successively against the sum
    /// already lifted. Both Anderson--Fuller/Lam nilpotent-series formulas are
    /// followed literally.
    ///
    /// # Errors
    ///
    /// Rejects an empty entered set, non-idempotents/nonorthogonal elements,
    /// or a nonnilpotent kernel series, and propagates exact coordinates.
    pub fn lift_complete_set_of_orthogonal_idempotents(
        &self,
        idempotents: &[Vec<F::Element>],
    ) -> Result<Vec<Vec<F::Element>>, StructureConstantsError> {
        if idempotents.is_empty() {
            return Err(StructureConstantsError::EmptyIdempotentSet);
        }
        let quotient_identity = self
            .algebra
            .find_multiplicative_identity()?
            .ok_or(StructureConstantsError::MissingMultiplicativeIdentity)?;
        let mut sum = vec![self.algebra.field().zero(); self.algebra.dimension()];
        for (position, idempotent) in idempotents.iter().enumerate() {
            self.validate_quotient_idempotent(idempotent, position + 1)?;
            for previous in &idempotents[..position] {
                let previous_then_current = self.algebra.coordinate_product(previous, idempotent);
                let current_then_previous = self.algebra.coordinate_product(idempotent, previous);
                if !is_zero_row(self.algebra.field(), &previous_then_current)
                    || !is_zero_row(self.algebra.field(), &current_then_previous)
                {
                    return Err(StructureConstantsError::EnteredIdempotentsNotOrthogonal);
                }
            }
            sum = add_rows(self.algebra.field(), &sum, idempotent);
        }
        if sum != quotient_identity.coordinates() {
            return Err(StructureConstantsError::EnteredIdempotentsNotComplete);
        }

        let first = self.lift_idempotent(&idempotents[0])?;
        let mut lifted = vec![first];
        for idempotent in &idempotents[1..] {
            let accumulated = lifted.iter().fold(
                vec![self.original.field().zero(); self.original.dimension()],
                |sum, idempotent| add_rows(self.original.field(), &sum, idempotent),
            );
            let pair = self.lift_two_orthogonal_idempotents(&accumulated, idempotent)?;
            lifted.push(pair[1].clone());
        }
        Ok(lifted)
    }

    /// Translate QPA's `LiftingIdempotent(f, e)` for this natural quotient
    /// map.
    ///
    /// This is the older Anderson--Fuller implementation in `algebra.gi`: it
    /// stops when the spaces `A h^i` and `A h^(i+1)` have equal dimension. A
    /// zero stabilized power produces the lifted idempotent; `None` preserves
    /// QPA's `false` result when the algorithm cannot construct a lift.
    ///
    /// # Errors
    ///
    /// Rejects an incorrectly sized or nonidempotent quotient element and
    /// propagates exact matrix and identity failures.
    pub fn lifting_idempotent(
        &self,
        idempotent: &[F::Element],
    ) -> Result<Option<Vec<F::Element>>, StructureConstantsError> {
        self.validate_quotient_idempotent(idempotent, 1)?;
        let representative = self.lift_coordinates(idempotent)?;
        let square = self
            .original
            .coordinate_product(&representative, &representative);
        let defect = subtract_rows(self.original.field(), &representative, &square);
        if is_zero_row(self.original.field(), &defect) {
            return Ok(Some(representative));
        }
        let identity = self
            .original
            .find_multiplicative_identity()?
            .ok_or(StructureConstantsError::MissingMultiplicativeIdentity)?;
        let mut index = 1usize;
        let mut power = defect.clone();
        loop {
            let next = self.original.coordinate_product(&power, &defect);
            let dimension = self.right_principal_space_dimension(&power)?;
            let next_dimension = self.right_principal_space_dimension(&next)?;
            if dimension == next_dimension {
                break;
            }
            index = index
                .checked_add(1)
                .ok_or(StructureConstantsError::NilpotentSeriesOverflow)?;
            power = next;
        }
        if !is_zero_row(self.original.field(), &power) {
            return Ok(None);
        }
        Ok(Some(self.idempotent_lifting_polynomial(
            &representative,
            index,
            identity.coordinates(),
        )?))
    }

    /// Translate QPA's newer `LiftIdempotent(f, v)` operation for this natural
    /// quotient map.
    ///
    /// Unlike [`Self::lifting_idempotent`], the source operation assumes that
    /// the quotient kernel is nilpotent and determines the nilpotence index by
    /// multiplying the defect until zero. It is also the primitive used by
    /// QPA's complete orthogonal-set lifting method.
    ///
    /// # Errors
    ///
    /// Rejects an incorrectly sized or nonidempotent quotient element. A
    /// nonnilpotent kernel is reported instead of allowing QPA's repeat loop to
    /// run forever.
    pub fn lift_idempotent(
        &self,
        idempotent: &[F::Element],
    ) -> Result<Vec<F::Element>, StructureConstantsError> {
        self.validate_quotient_idempotent(idempotent, 1)?;
        let representative = self.lift_coordinates(idempotent)?;
        let square = self
            .original
            .coordinate_product(&representative, &representative);
        let defect = subtract_rows(self.original.field(), &representative, &square);
        if is_zero_row(self.original.field(), &defect) {
            return Ok(representative);
        }
        let mut nilindex = 1usize;
        let mut power = defect.clone();
        while !is_zero_row(self.original.field(), &power) {
            nilindex = nilindex
                .checked_add(1)
                .ok_or(StructureConstantsError::NilpotentSeriesOverflow)?;
            power = self.original.coordinate_product(&defect, &power);
            if nilindex > self.original.dimension().saturating_add(1) {
                return Err(StructureConstantsError::QuotientKernelNotNilpotent);
            }
        }
        let identity = self
            .original
            .find_multiplicative_identity()?
            .ok_or(StructureConstantsError::MissingMultiplicativeIdentity)?;
        self.idempotent_lifting_polynomial(&representative, nilindex, identity.coordinates())
    }

    fn validate_quotient_idempotent(
        &self,
        idempotent: &[F::Element],
        position: usize,
    ) -> Result<(), StructureConstantsError> {
        if idempotent.len() != self.algebra.dimension() {
            return Err(StructureConstantsError::ElementLength {
                actual: idempotent.len(),
                expected: self.algebra.dimension(),
            });
        }
        if self.algebra.coordinate_product(idempotent, idempotent) != idempotent {
            return Err(StructureConstantsError::EnteredElementNotIdempotent { position });
        }
        Ok(())
    }

    fn right_principal_space_dimension(
        &self,
        right: &[F::Element],
    ) -> Result<usize, StructureConstantsError> {
        let products = self
            .original
            .basis_elements()
            .iter()
            .map(|basis| {
                self.original
                    .multiply_coordinates(basis.coordinates(), right)
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(span_rows(self.original.field(), self.original.dimension(), products)?.nrows())
    }

    fn idempotent_lifting_polynomial(
        &self,
        representative: &[F::Element],
        nilindex: usize,
        identity: &[F::Element],
    ) -> Result<Vec<F::Element>, StructureConstantsError> {
        let mut polynomial_value = vec![self.original.field().zero(); self.original.dimension()];
        for exponent in 0..nilindex {
            let binomial = binomial_integer(nilindex, exponent + 1);
            let coefficient = self.original.field().nth(if exponent.is_multiple_of(2) {
                binomial
            } else {
                -binomial
            });
            let term = algebra_power(
                &self.original,
                representative,
                u64::try_from(exponent)
                    .map_err(|_| StructureConstantsError::NilpotentSeriesOverflow)?,
                identity,
            );
            for (entry, term_entry) in polynomial_value.iter_mut().zip(term) {
                self.original
                    .field()
                    .add_mul_assign(entry, &coefficient, &term_entry);
            }
        }
        let representative_power = algebra_power(
            &self.original,
            representative,
            u64::try_from(nilindex)
                .map_err(|_| StructureConstantsError::NilpotentSeriesOverflow)?,
            identity,
        );
        let polynomial_power = algebra_power(
            &self.original,
            &polynomial_value,
            u64::try_from(nilindex)
                .map_err(|_| StructureConstantsError::NilpotentSeriesOverflow)?,
            identity,
        );
        Ok(self
            .original
            .coordinate_product(&representative_power, &polynomial_power))
    }

    /// Translate QPA's `LiftTwoOrthogonalIdempotents(f, v, w)` for this
    /// natural quotient map.
    ///
    /// The first element is an already lifted idempotent in the original
    /// algebra, while the second is an orthogonal quotient idempotent. The
    /// returned pair retains the first element and applies QPA's finite
    /// geometric-series correction to a lift of the second.
    ///
    /// # Errors
    ///
    /// Rejects incorrectly sized, nonidempotent, or nonorthogonal inputs and
    /// a nonnilpotent kernel series, and propagates exact coordinate failures.
    pub fn lift_two_orthogonal_idempotents(
        &self,
        first: &[F::Element],
        second_quotient: &[F::Element],
    ) -> Result<[Vec<F::Element>; 2], StructureConstantsError> {
        if first.len() != self.original.dimension() {
            return Err(StructureConstantsError::ElementLength {
                actual: first.len(),
                expected: self.original.dimension(),
            });
        }
        if self.original.coordinate_product(first, first) != first {
            return Err(StructureConstantsError::EnteredElementNotIdempotent { position: 1 });
        }
        self.validate_quotient_idempotent(second_quotient, 2)?;
        let first_quotient = self.project_coordinates(first)?;
        if !is_zero_row(
            self.algebra.field(),
            &self
                .algebra
                .coordinate_product(&first_quotient, second_quotient),
        ) || !is_zero_row(
            self.algebra.field(),
            &self
                .algebra
                .coordinate_product(second_quotient, &first_quotient),
        ) {
            return Err(StructureConstantsError::EnteredIdempotentsNotOrthogonal);
        }
        let second = self.lift_idempotent(second_quotient)?;
        let product = self.original.coordinate_product(&second, first);
        let mut power = product.clone();
        let identity = self
            .original
            .find_multiplicative_identity()?
            .ok_or(StructureConstantsError::MissingMultiplicativeIdentity)?;
        let mut series = identity.coordinates().to_vec();
        let mut steps = 0usize;
        while !is_zero_row(self.original.field(), &power) {
            series = add_rows(self.original.field(), &series, &power);
            power = self.original.coordinate_product(&product, &power);
            steps += 1;
            if steps > self.original.dimension() {
                return Err(StructureConstantsError::QuotientKernelNotNilpotent);
            }
        }
        let one_minus_first = subtract_rows(self.original.field(), identity.coordinates(), first);
        let one_minus_product =
            subtract_rows(self.original.field(), identity.coordinates(), &product);
        let result = self.original.coordinate_product(&one_minus_first, &series);
        let result = self.original.coordinate_product(&result, &second);
        let corrected = self
            .original
            .coordinate_product(&result, &one_minus_product);
        Ok([first.to_vec(), corrected])
    }
}

fn binomial_integer(n: usize, k: usize) -> Integer {
    let k = k.min(n - k);
    let mut result = Integer::from(1);
    for step in 1..=k {
        result *= Integer::from(n - k + step);
        result /= Integer::from(step);
    }
    result
}

fn is_zero_row<F: Field>(field: &F, row: &[F::Element]) -> bool {
    row.iter().all(|entry| field.is_zero(entry))
}

pub(crate) fn span_rows<F: Field>(
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

pub(crate) fn is_commutative<F: Field>(algebra: &StructureConstantsAlgebra<F>) -> bool {
    for left in 0..algebra.dimension() {
        for right in left + 1..algebra.dimension() {
            if algebra.basis_product(left, right) != algebra.basis_product(right, left) {
                return false;
            }
        }
    }
    true
}

fn multiply_row_matrix<F: Field>(row: &[F::Element], matrix: &DenseMatrix<F>) -> Vec<F::Element> {
    debug_assert_eq!(row.len(), matrix.nrows());
    let field = matrix.field();
    let mut result = vec![field.zero(); matrix.ncols()];
    for (inner, coefficient) in row.iter().enumerate() {
        if field.is_zero(coefficient) {
            continue;
        }
        for column in 0..matrix.ncols() {
            field.add_mul_assign(&mut result[column], coefficient, &matrix[(inner, column)]);
        }
    }
    result
}

impl<F: Field> StructureConstantsSubalgebra<F> {
    pub(crate) fn from_basis(
        ambient: StructureConstantsAlgebra<F>,
        basis: DenseMatrix<F>,
    ) -> Result<Self, StructureConstantsError> {
        if basis.ncols() != ambient.dimension() {
            return Err(StructureConstantsError::SubalgebraBasisWidth {
                actual: basis.ncols(),
                expected: ambient.dimension(),
            });
        }
        if basis.rank() != basis.nrows() {
            return Err(StructureConstantsError::DependentSubalgebraBasis);
        }
        let dimension = basis.nrows();
        let tensor_length = dimension
            .checked_mul(dimension)
            .and_then(|square| square.checked_mul(dimension))
            .ok_or(StructureConstantsError::TensorShapeOverflow { dimension })?;
        let mut products = Vec::with_capacity(tensor_length);
        for left in 0..dimension {
            for right in 0..dimension {
                let product = ambient.coordinate_product(basis.row(left), basis.row(right));
                let coordinates = basis
                    .solution_row(&product)?
                    .ok_or(StructureConstantsError::BasisNotClosedUnderMultiplication)?;
                products.extend(coordinates);
            }
        }
        let algebra = StructureConstantsAlgebra::new(ambient.field().clone(), dimension, products)?;
        Ok(Self {
            ambient,
            basis,
            algebra,
        })
    }

    /// Return the ambient algebra.
    #[must_use]
    pub const fn ambient(&self) -> &StructureConstantsAlgebra<F> {
        &self.ambient
    }

    /// Return basis rows in ambient standard coordinates.
    #[must_use]
    pub const fn basis(&self) -> &DenseMatrix<F> {
        &self.basis
    }

    /// Return the induced structure-constant algebra in this basis.
    #[must_use]
    pub const fn algebra(&self) -> &StructureConstantsAlgebra<F> {
        &self.algebra
    }

    /// Lift a coordinate row from the subalgebra basis to the ambient basis.
    ///
    /// # Errors
    ///
    /// Rejects a row with the wrong subalgebra dimension.
    pub fn lift_coordinates(
        &self,
        coordinates: &[F::Element],
    ) -> Result<Vec<F::Element>, StructureConstantsError> {
        if coordinates.len() != self.algebra.dimension() {
            return Err(StructureConstantsError::ElementLength {
                actual: coordinates.len(),
                expected: self.algebra.dimension(),
            });
        }
        let row = DenseMatrix::from_flat(
            1,
            coordinates.len(),
            coordinates.to_vec(),
            self.ambient.field().clone(),
        )?;
        Ok(row.multiply(&self.basis)?.row(0).to_vec())
    }
}

impl<F: Field> PartialEq for StructureConstantsAlgebra<F> {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }
}

impl<F: Field> Eq for StructureConstantsAlgebra<F> {}

impl<F: Field> OppositeStructureConstantsAlgebra<F> {
    /// Return whether QPA constructed an actual opposite-algebra wrapper.
    #[must_use]
    pub const fn is_opposite_algebra(&self) -> bool {
        matches!(self, Self::Opposite(_))
    }

    /// Return whether QPA's commutative fast path returned the original object.
    #[must_use]
    pub const fn is_identical_to_underlying(&self) -> bool {
        matches!(self, Self::Original(_))
    }

    /// Return QPA's retained `UnderlyingAlgebra` object.
    #[must_use]
    pub fn underlying_algebra(&self) -> &StructureConstantsAlgebra<F> {
        match self {
            Self::Original(algebra) => algebra,
            Self::Opposite(data) => &data.underlying,
        }
    }

    /// Return the structure-constant algebra implementing this object's
    /// multiplication.
    #[must_use]
    pub fn as_structure_constants_algebra(&self) -> &StructureConstantsAlgebra<F> {
        match self {
            Self::Original(algebra) => algebra,
            Self::Opposite(data) => &data.reversed,
        }
    }

    /// Return the exact coefficient field.
    #[must_use]
    pub fn field(&self) -> &F {
        self.underlying_algebra().field()
    }

    /// Return QPA's delegated opposite-algebra dimension.
    #[must_use]
    pub fn dimension(&self) -> usize {
        self.underlying_algebra().dimension()
    }

    /// Finite structure-constant algebras remain finite dimensional under
    /// QPA's opposite construction.
    #[must_use]
    pub const fn is_finite_dimensional(&self) -> bool {
        true
    }

    /// Wrap an element of the identical underlying algebra.
    ///
    /// This is QPA's `ObjByExtRep` family check specialized to the typed
    /// structure-constant representation.
    ///
    /// # Errors
    ///
    /// Rejects an element from a different underlying parent object.
    pub fn element_from_underlying(
        &self,
        underlying: StructureConstantsAlgebraElement<F>,
    ) -> Result<OppositeStructureConstantsAlgebraElement<F>, StructureConstantsError> {
        if underlying.algebra() != self.underlying_algebra() {
            return Err(StructureConstantsError::DifferentParents);
        }
        Ok(OppositeStructureConstantsAlgebraElement {
            algebra: self.clone(),
            underlying,
        })
    }

    /// Construct an opposite element from underlying-basis coordinates.
    ///
    /// # Errors
    ///
    /// Rejects a coordinate row of the wrong dimension.
    pub fn element(
        &self,
        coordinates: Vec<F::Element>,
    ) -> Result<OppositeStructureConstantsAlgebraElement<F>, StructureConstantsError> {
        self.element_from_underlying(self.underlying_algebra().element(coordinates)?)
    }

    /// Translate QPA's generic opposite-algebra `A.name` lookup.
    ///
    /// The source delegates the record-name lookup to `UnderlyingAlgebra(A)`
    /// and wraps the returned element in the opposite family.
    ///
    /// # Errors
    ///
    /// Propagates checked named lookup and opposite-family wrapping.
    pub fn element_by_name(
        &self,
        name: &str,
    ) -> Result<Option<OppositeStructureConstantsAlgebraElement<F>>, StructureConstantsError> {
        self.underlying_algebra()
            .element_by_name(name)?
            .map(|element| self.element_from_underlying(element))
            .transpose()
    }

    /// Return QPA's wrapped `BasisVectors(Basis(A))` in underlying basis order.
    #[must_use]
    pub fn basis_elements(&self) -> Vec<OppositeStructureConstantsAlgebraElement<F>> {
        self.underlying_algebra()
            .basis_elements()
            .into_iter()
            .map(|underlying| OppositeStructureConstantsAlgebraElement {
                algebra: self.clone(),
                underlying,
            })
            .collect()
    }

    /// Return QPA's wrapped canonical basis vectors.
    ///
    /// Structure-constant algebras use their ordered standard basis as the
    /// canonical basis, so QPA's generic opposite-basis wrapper retains the
    /// same list as [`Self::basis_elements`].
    #[must_use]
    pub fn canonical_basis_elements(&self) -> Vec<OppositeStructureConstantsAlgebraElement<F>> {
        self.basis_elements()
    }

    /// Return coefficients in QPA's wrapped underlying basis.
    ///
    /// # Errors
    ///
    /// Rejects an element from a different opposite-algebra object.
    pub fn coefficients(
        &self,
        element: &OppositeStructureConstantsAlgebraElement<F>,
    ) -> Result<Vec<F::Element>, StructureConstantsError> {
        if !self.contains(element) {
            return Err(StructureConstantsError::DifferentParents);
        }
        Ok(element.coordinates().to_vec())
    }

    /// Return the additive zero opposite element.
    #[must_use]
    pub fn zero(&self) -> OppositeStructureConstantsAlgebraElement<F> {
        OppositeStructureConstantsAlgebraElement {
            algebra: self.clone(),
            underlying: self.underlying_algebra().zero(),
        }
    }

    /// Return the wrapped multiplicative identity, or QPA's `fail` as `None`.
    ///
    /// # Errors
    ///
    /// Propagates the underlying exact identity computation.
    pub fn one(
        &self,
    ) -> Result<Option<OppositeStructureConstantsAlgebraElement<F>>, StructureConstantsError> {
        self.underlying_algebra()
            .find_multiplicative_identity()?
            .map(|underlying| self.element_from_underlying(underlying))
            .transpose()
    }

    /// Test QPA family membership in this exact opposite-algebra object.
    #[must_use]
    pub fn contains(&self, element: &OppositeStructureConstantsAlgebraElement<F>) -> bool {
        self == element.algebra()
    }
}

impl<F: Field> PartialEq for OppositeStructureConstantsAlgebra<F> {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Original(left), Self::Original(right)) => left == right,
            (Self::Opposite(left), Self::Opposite(right)) => Arc::ptr_eq(left, right),
            (Self::Original(_), Self::Opposite(_)) | (Self::Opposite(_), Self::Original(_)) => {
                false
            }
        }
    }
}

impl<F: Field> Eq for OppositeStructureConstantsAlgebra<F> {}

impl<F: Field> fmt::Display for OppositeStructureConstantsAlgebra<F> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "OppositeAlgebra( <structure-constants algebra of dimension {}> )",
            self.dimension()
        )
    }
}

/// An element of a [`StructureConstantsAlgebra`].
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct StructureConstantsAlgebraElement<F: Field> {
    algebra: StructureConstantsAlgebra<F>,
    coordinates: Vec<F::Element>,
}

impl<F: Field> StructureConstantsAlgebraElement<F> {
    /// Return the parent algebra.
    #[must_use]
    pub const fn algebra(&self) -> &StructureConstantsAlgebra<F> {
        &self.algebra
    }

    /// Return standard-basis coordinates.
    #[must_use]
    pub fn coordinates(&self) -> &[F::Element] {
        &self.coordinates
    }

    /// Multiply elements using the stored structure constants.
    ///
    /// # Errors
    ///
    /// Rejects elements from distinct parent objects.
    pub fn multiply(&self, right: &Self) -> Result<Self, StructureConstantsError> {
        if self.algebra != right.algebra {
            return Err(StructureConstantsError::DifferentParents);
        }
        let field = self.algebra.field();
        let mut coordinates = vec![field.zero(); self.algebra.dimension()];
        for (left_index, left_coefficient) in self.coordinates.iter().enumerate() {
            if field.is_zero(left_coefficient) {
                continue;
            }
            for (right_index, right_coefficient) in right.coordinates.iter().enumerate() {
                if field.is_zero(right_coefficient) {
                    continue;
                }
                let scalar = field.mul(left_coefficient, right_coefficient);
                for (coordinate, product_coefficient) in coordinates
                    .iter_mut()
                    .zip(self.algebra.basis_product(left_index, right_index))
                {
                    field.add_mul_assign(coordinate, &scalar, product_coefficient);
                }
            }
        }
        self.algebra.element(coordinates)
    }
}

impl<F: Field> OppositeStructureConstantsAlgebraElement<F> {
    /// Return the parent opposite algebra.
    #[must_use]
    pub const fn algebra(&self) -> &OppositeStructureConstantsAlgebra<F> {
        &self.algebra
    }

    /// Return whether this is one of QPA's wrapped opposite elements.
    #[must_use]
    pub const fn is_opposite_algebra_element(&self) -> bool {
        self.algebra.is_opposite_algebra()
    }

    /// Return the exact element wrapped from QPA's underlying algebra.
    #[must_use]
    pub const fn underlying_element(&self) -> &StructureConstantsAlgebraElement<F> {
        &self.underlying
    }

    /// Return coefficients in the wrapped underlying basis.
    #[must_use]
    pub fn coordinates(&self) -> &[F::Element] {
        self.underlying.coordinates()
    }

    fn check_same_parent(&self, other: &Self) -> Result<(), StructureConstantsError> {
        if self.algebra != other.algebra {
            return Err(StructureConstantsError::DifferentParents);
        }
        Ok(())
    }

    /// Add two elements by wrapping QPA's underlying sum.
    ///
    /// # Errors
    ///
    /// Rejects elements from distinct opposite-algebra objects.
    pub fn add(&self, right: &Self) -> Result<Self, StructureConstantsError> {
        self.check_same_parent(right)?;
        let coordinates = self
            .coordinates()
            .iter()
            .zip(right.coordinates())
            .map(|(left, right)| self.algebra.field().add(left, right))
            .collect();
        self.algebra.element(coordinates)
    }

    /// Multiply using QPA's `op(a) * op(b) = op(b * a)` rule.
    ///
    /// # Errors
    ///
    /// Rejects elements from distinct opposite-algebra objects.
    pub fn multiply(&self, right: &Self) -> Result<Self, StructureConstantsError> {
        self.check_same_parent(right)?;
        let underlying = if self.algebra.is_opposite_algebra() {
            right.underlying.multiply(&self.underlying)?
        } else {
            self.underlying.multiply(&right.underlying)?
        };
        self.algebra.element_from_underlying(underlying)
    }

    /// Multiply by a coefficient-field scalar on either side.
    ///
    /// # Errors
    ///
    /// Propagates checked construction of the scaled element.
    pub fn scale(&self, scalar: &F::Element) -> Result<Self, StructureConstantsError> {
        self.algebra.element(
            self.coordinates()
                .iter()
                .map(|coefficient| self.algebra.field().mul(coefficient, scalar))
                .collect(),
        )
    }

    /// Return QPA's same-mutability additive inverse.
    ///
    /// # Errors
    ///
    /// Propagates checked construction of the negated element.
    pub fn negate(&self) -> Result<Self, StructureConstantsError> {
        self.algebra.element(
            self.coordinates()
                .iter()
                .map(|coefficient| self.algebra.field().neg(coefficient))
                .collect(),
        )
    }

    /// Return QPA's `ZeroOp` result in the identical opposite parent.
    #[must_use]
    pub fn zero_same_parent(&self) -> Self {
        self.algebra.zero()
    }

    /// Return QPA's `ExtRepOfObj` payload for an opposite element.
    #[must_use]
    pub fn external_representation(&self) -> &[F::Element] {
        self.underlying.coordinates()
    }
}

impl<F: QpaCoefficientField> OppositeStructureConstantsAlgebraElement<F> {
    /// Compare by QPA's delegated underlying-element order.
    ///
    /// # Errors
    ///
    /// Rejects elements from distinct opposite-algebra objects and propagates
    /// exact coefficient-order construction.
    pub fn compare(&self, other: &Self) -> Result<Ordering, StructureConstantsError> {
        self.check_same_parent(other)?;
        for (left, right) in self.coordinates().iter().zip(other.coordinates()) {
            match self.algebra.field().qpa_compare(left, right)? {
                Ordering::Equal => {}
                ordering => return Ok(ordering),
            }
        }
        Ok(Ordering::Equal)
    }
}

impl<F: Field> PartialEq for OppositeStructureConstantsAlgebraElement<F> {
    fn eq(&self, other: &Self) -> bool {
        self.algebra == other.algebra && self.underlying == other.underlying
    }
}

impl<F: Field> Eq for OppositeStructureConstantsAlgebraElement<F> {}

impl<F> fmt::Display for OppositeStructureConstantsAlgebraElement<F>
where
    F: QpaCoefficientField,
    F::Element: fmt::Display,
{
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let field = self.algebra.field();
        let nonzero = self
            .coordinates()
            .iter()
            .enumerate()
            .filter(|(_, coefficient)| !field.is_zero(coefficient))
            .collect::<Vec<_>>();
        if nonzero.is_empty() {
            return if self.coordinates().is_empty() {
                write!(formatter, "op( 0 )")
            } else {
                write!(formatter, "op( 0*v.1 )")
            };
        }

        let characteristic = field.qpa_characteristic();
        let finite_elements = match characteristic {
            QpaFieldCharacteristic::Zero => None,
            QpaFieldCharacteristic::Prime(_) => field.qpa_elements().map_err(|_| fmt::Error)?,
        };
        write!(formatter, "op( ")?;
        for (term, (index, coefficient)) in nonzero.into_iter().enumerate() {
            if term != 0 {
                write!(formatter, "+")?;
            }
            if field.is_one(coefficient) {
                write!(formatter, "v.{}", index + 1)?;
                continue;
            }
            write!(formatter, "(")?;
            match characteristic {
                QpaFieldCharacteristic::Zero => write!(formatter, "{coefficient}")?,
                QpaFieldCharacteristic::Prime(prime) => {
                    let position = finite_elements
                        .as_ref()
                        .and_then(|elements| {
                            elements.iter().position(|element| element == coefficient)
                        })
                        .ok_or(fmt::Error)?;
                    let exponent = position.checked_sub(1).ok_or(fmt::Error)?;
                    if exponent == 1 {
                        write!(formatter, "Z({prime})")?;
                    } else {
                        write!(formatter, "Z({prime})^{exponent}")?;
                    }
                }
            }
            write!(formatter, ")*v.{}", index + 1)?;
        }
        write!(formatter, " )")
    }
}

/// Failure to construct or operate on a structure-constant algebra.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum StructureConstantsError {
    /// The requested multiplication tensor cube overflowed `usize`.
    #[error("the structure-constant tensor shape for dimension {dimension} overflows")]
    TensorShapeOverflow {
        /// Requested algebra dimension.
        dimension: usize,
    },
    /// The supplied flattened multiplication tensor has the wrong length.
    #[error(
        "the structure-constant tensor has {actual} entries, expected {expected} for dimension {dimension}"
    )]
    TensorLength {
        /// Supplied coefficient count.
        actual: usize,
        /// Required coefficient count.
        expected: usize,
        /// Requested algebra dimension.
        dimension: usize,
    },
    /// An element coordinate row has the wrong length.
    #[error("the algebra element has {actual} coordinates, expected {expected}")]
    ElementLength {
        /// Supplied coordinate count.
        actual: usize,
        /// Required coordinate count.
        expected: usize,
    },
    /// The block matrix used by QPA's identity method overflowed `usize`.
    #[error("the multiplicative-identity matrix shape for dimension {dimension} overflows")]
    IdentityMatrixShapeOverflow {
        /// Algebra dimension.
        dimension: usize,
    },
    /// GAP's `BlockMatrix` call has no block-shape information in dimension zero.
    #[error("need <rpb>, <cpb>, <zero> if <blocks> is empty")]
    EmptyIdentityBlockMatrix,
    /// The entered elements have different parent objects.
    #[error("the algebra elements belong to different structure-constant algebras")]
    DifferentParents,
    /// A subspace equation shape overflowed addressable storage.
    #[error("the algebra subspace matrix shape for dimension {dimension} overflows")]
    SubspaceMatrixShapeOverflow {
        /// Algebra dimension.
        dimension: usize,
    },
    /// A subalgebra basis row has the wrong ambient width.
    #[error("the subalgebra basis has width {actual}, expected {expected}")]
    SubalgebraBasisWidth {
        /// Supplied row width.
        actual: usize,
        /// Ambient dimension.
        expected: usize,
    },
    /// GAP's supplied `"basis"` list must be linearly independent.
    #[error("the supplied subalgebra basis is linearly dependent")]
    DependentSubalgebraBasis,
    /// A purported subalgebra basis was not closed under multiplication.
    #[error("the supplied basis is not closed under algebra multiplication")]
    BasisNotClosedUnderMultiplication,
    /// QPA's faithful radical representation requires an algebra identity.
    #[error("the algebra does not have a multiplicative identity")]
    MissingMultiplicativeIdentity,
    /// The radical algorithm requested a characteristic-polynomial coefficient outside its degree.
    #[error(
        "radical coefficient index is invalid for dimension {dimension} and prime power {prime_power}"
    )]
    RadicalCoefficientIndex {
        /// Matrix dimension.
        dimension: usize,
        /// Current characteristic power.
        prime_power: u64,
    },
    /// An ideal basis row has the wrong original width.
    #[error("the ideal basis has width {actual}, expected {expected}")]
    IdealBasisWidth {
        /// Supplied row width.
        actual: usize,
        /// Original algebra dimension.
        expected: usize,
    },
    /// GAP's entered ideal basis must be linearly independent.
    #[error("the supplied ideal basis is linearly dependent")]
    DependentIdealBasis,
    /// Extending the ideal basis by standard vectors failed to span the original algebra.
    #[error("the natural quotient basis does not span the original algebra")]
    IncompleteNaturalQuotientBasis,
    /// `CentralIdempotentsOfAlgebra` requires a finite coefficient field.
    #[error("the entered algebra is not an algebra over a finite field")]
    FiniteFieldRequired,
    /// GAP's `IsPosInt` filter excludes zero before `PowerSubalgebraOp` runs.
    #[error("the power-subalgebra degree must be positive")]
    PowerSubalgebraDegreeNotPositive,
    /// QPA requires the requested power degree to divide the field degree.
    #[error("degree {degree} is not a factor of field degree {field_degree}")]
    PowerSubalgebraDegreeNotDivisor {
        /// Entered power-subalgebra degree.
        degree: usize,
        /// Degree of the coefficient field over its prime field.
        field_degree: usize,
    },
    /// The deterministic traversal of GAP's splitting-element coefficient set overflowed.
    #[error("the splitting-element search space overflows addressable storage")]
    SplittingSearchOverflow,
    /// No element in GAP's prescribed splitting coefficient set generated the algebra.
    #[error("no splitting element was found in GAP's prescribed coefficient set")]
    MissingSplittingElement,
    /// The Frobenius-fixed element was expected to decompose a nonsimple split ideal.
    #[error("the Frobenius-fixed element did not split the commutative ideal")]
    FrobeniusElementDidNotSplit,
    /// An ideal product unexpectedly fell outside the retained ideal basis.
    #[error("an ideal element has no coordinates in its retained basis")]
    IdealBasisCoordinateFailure,
    /// Initial prime-field support does not yet factor nonsplit extension polynomials.
    #[error(
        "the minimal polynomial has degree {degree} but only {linear_factors} linear factors over the prime field"
    )]
    NonsplitMinimalPolynomial {
        /// Polynomial degree.
        degree: usize,
        /// Number of roots in the entered prime field.
        linear_factors: usize,
    },
    /// QPA's complete-set lifting indexes the first entered idempotent.
    #[error("the entered complete idempotent set is empty")]
    EmptyIdempotentSet,
    /// One entered quotient element was not idempotent.
    #[error("entered element {position} is not an idempotent")]
    EnteredElementNotIdempotent {
        /// One-based list position.
        position: usize,
    },
    /// The entered quotient idempotents were not pairwise orthogonal.
    #[error("the entered idempotents are not orthogonal")]
    EnteredIdempotentsNotOrthogonal,
    /// The entered quotient idempotents did not sum to one.
    #[error("the entered idempotents are not complete")]
    EnteredIdempotentsNotComplete,
    /// A nilpotence-index counter overflowed addressable storage.
    #[error("the idempotent-lifting nilpotence index overflows usize")]
    NilpotentSeriesOverflow,
    /// QPA's lifting formula requires a nilpotent quotient kernel.
    #[error("the natural quotient kernel is not nilpotent")]
    QuotientKernelNotNilpotent,
    /// Exact matrix construction or elimination failed.
    #[error(transparent)]
    Matrix(#[from] MatrixError),
    /// Coefficient-field enumeration failed.
    #[error(transparent)]
    Coefficient(#[from] CoefficientDomainError),
}

#[cfg(test)]
mod tests {
    use std::cmp::Ordering;

    use numerica::domains::{
        Field, Ring,
        finite_field::{FiniteFieldCore, Z2, Zp},
        rational::Q,
    };
    use serde_json::{Value, json};

    use super::StructureConstantsAlgebra;
    use crate::{DenseMatrix, QpaCoefficientField};

    #[test]
    fn finds_gap_left_identity_and_rejects_a_nonunital_algebra() {
        // e1^2=e1, e1e2=e2, e2e1=e2, e2^2=0.
        let q = |value: i64| Q.nth(value.into());
        let dual_numbers = StructureConstantsAlgebra::new(
            Q,
            2,
            [
                1, 0, // e1 e1
                0, 1, // e1 e2
                0, 1, // e2 e1
                0, 0, // e2 e2
            ]
            .into_iter()
            .map(q)
            .collect(),
        )
        .unwrap();
        assert_eq!(
            dual_numbers
                .find_multiplicative_identity()
                .unwrap()
                .unwrap()
                .coordinates(),
            &[q(1), q(0)]
        );

        let zero = StructureConstantsAlgebra::new(Z2, 1, vec![Z2.zero()]).unwrap();
        assert!(zero.find_multiplicative_identity().unwrap().is_none());

        let empty = StructureConstantsAlgebra::new(Q, 0, Vec::new()).unwrap();
        assert!(matches!(
            empty.find_multiplicative_identity(),
            Err(super::StructureConstantsError::EmptyIdentityBlockMatrix)
        ));
    }

    #[test]
    fn pinned_qpa_identity_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/algebra_structure_constants.schema-v1.json"
        ))
        .unwrap();

        let rational = StructureConstantsAlgebra::new(
            Q,
            2,
            [1, 0, 0, 1, 0, 0, 0, 0]
                .into_iter()
                .map(|value| Q.nth(value.into()))
                .collect(),
        )
        .unwrap();
        let rational_identity = rational.find_multiplicative_identity().unwrap().unwrap();
        assert_eq!(
            json!(
                rational_identity
                    .coordinates()
                    .iter()
                    .map(|coefficient| coefficient.to_string().parse::<i64>().unwrap())
                    .collect::<Vec<_>>()
            ),
            expected["rational"]["identity"]
        );

        let binary = StructureConstantsAlgebra::new(
            Z2,
            2,
            vec![
                Z2.one(),
                Z2.zero(),
                Z2.zero(),
                Z2.one(),
                Z2.zero(),
                Z2.zero(),
                Z2.zero(),
                Z2.zero(),
            ],
        )
        .unwrap();
        let binary_identity = binary.find_multiplicative_identity().unwrap().unwrap();
        assert_eq!(
            json!(
                binary_identity
                    .coordinates()
                    .iter()
                    .map(|coefficient| Z2
                        .to_integer(coefficient)
                        .to_string()
                        .parse::<i64>()
                        .unwrap())
                    .collect::<Vec<_>>()
            ),
            expected["f2"]["identity"]
        );
    }

    #[test]
    fn multiplies_in_standard_basis_coordinates() {
        let algebra = StructureConstantsAlgebra::new(
            Z2,
            2,
            vec![
                Z2.one(),
                Z2.zero(),
                Z2.zero(),
                Z2.one(),
                Z2.zero(),
                Z2.one(),
                Z2.zero(),
                Z2.zero(),
            ],
        )
        .unwrap();
        let basis = algebra.basis_elements();
        assert_eq!(basis[1].multiply(&basis[0]).unwrap(), basis[1]);
        assert_eq!(basis[1].multiply(&basis[1]).unwrap(), algebra.zero());
    }

    fn product_of_fields<F: Field>(field: F, dimension: usize) -> StructureConstantsAlgebra<F> {
        let mut products = vec![field.zero(); dimension * dimension * dimension];
        for coordinate in 0..dimension {
            products[(coordinate * dimension + coordinate) * dimension + coordinate] = field.one();
        }
        StructureConstantsAlgebra::new(field, dimension, products).unwrap()
    }

    fn upper_triangular_two_by_two<F: Field>(field: F) -> StructureConstantsAlgebra<F> {
        let dimension = 3;
        let mut products = vec![field.zero(); dimension * dimension * dimension];
        let mut set = |left: usize, right: usize, output: usize| {
            products[(left * dimension + right) * dimension + output] = field.one();
        };
        set(0, 0, 0);
        set(0, 1, 1);
        set(1, 2, 1);
        set(2, 2, 2);
        StructureConstantsAlgebra::new(field, dimension, products).unwrap()
    }

    fn dual_numbers<F: Field>(field: F) -> StructureConstantsAlgebra<F> {
        let dimension = 2;
        let mut products = vec![field.zero(); dimension * dimension * dimension];
        let mut set = |left: usize, right: usize, output: usize| {
            products[(left * dimension + right) * dimension + output] = field.one();
        };
        set(0, 0, 0);
        set(0, 1, 1);
        set(1, 0, 1);
        StructureConstantsAlgebra::new(field, dimension, products).unwrap()
    }

    fn shifted_dual_numbers<F: Field>(field: F) -> StructureConstantsAlgebra<F> {
        let dimension = 2;
        let mut products = vec![field.zero(); dimension * dimension * dimension];
        products[0] = field.one(); // b1 b1 = b1 + b2
        products[1] = field.one();
        products[3] = field.one(); // b1 b2 = b2
        products[5] = field.one(); // b2 b1 = b2
        StructureConstantsAlgebra::new(field, dimension, products).unwrap()
    }

    fn shifted_upper_triangular_two_by_two<F: Field>(field: F) -> StructureConstantsAlgebra<F> {
        let dimension = 3;
        let mut products = vec![field.zero(); dimension * dimension * dimension];
        products[0] = field.one(); // b1 b1 = b1
        products[5] = field.add(&field.one(), &field.one()); // b1 b2 = 2 b3
        products[8] = field.one(); // b1 b3 = b3
        products[13] = field.one(); // b2 b2 = b2
        products[23] = field.one(); // b3 b2 = b3
        StructureConstantsAlgebra::new(field, dimension, products).unwrap()
    }

    fn dense_rows_json<F: Field>(
        _field: &F,
        matrix: &DenseMatrix<F>,
        coefficient_json: &impl Fn(&F::Element) -> Value,
    ) -> Value {
        json!(
            (0..matrix.nrows())
                .map(|row| matrix
                    .row(row)
                    .iter()
                    .map(coefficient_json)
                    .collect::<Vec<_>>())
                .collect::<Vec<_>>()
        )
    }

    fn row_json<F: Field>(
        _field: &F,
        row: &[F::Element],
        coefficient_json: &impl Fn(&F::Element) -> Value,
    ) -> Value {
        json!(row.iter().map(coefficient_json).collect::<Vec<_>>())
    }

    fn algebra_properties_json<F: QpaCoefficientField>(
        algebra: &StructureConstantsAlgebra<F>,
        coefficient_json: &impl Fn(&F::Element) -> Value,
    ) -> Value {
        let sparse_basis = |basis: &DenseMatrix<F>| {
            (0..basis.nrows())
                .map(|row| {
                    basis
                        .row(row)
                        .iter()
                        .enumerate()
                        .filter(|(_, coefficient)| !algebra.field().is_zero(coefficient))
                        .map(|(coordinate, coefficient)| {
                            json!([coordinate + 1, coefficient_json(coefficient)])
                        })
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        };
        let series = algebra
            .radical_series()
            .unwrap()
            .iter()
            .map(sparse_basis)
            .collect::<Vec<_>>();
        let blocks = algebra
            .blocks()
            .unwrap()
            .iter()
            .map(|block| sparse_basis(block.basis()))
            .collect::<Vec<_>>();
        json!({
            "radical_series": series,
            "radical_square_zero": algebra.is_radical_square_zero().unwrap(),
            "basic": algebra.is_basic_algebra().unwrap().unwrap(),
            "elementary": algebra.is_elementary_algebra().unwrap().unwrap(),
            "blocks": blocks
        })
    }

    #[test]
    fn centre_radical_and_natural_quotient_match_gap_basis_order() {
        let algebra = upper_triangular_two_by_two(Z2);
        let centre = algebra.centre().unwrap();
        assert_eq!(centre.basis().nrows(), 1);
        assert_eq!(centre.basis().row(0), &[Z2.one(), Z2.zero(), Z2.one()]);
        let radical = algebra.radical().unwrap();
        assert_eq!(radical.basis().nrows(), 1);
        assert_eq!(radical.basis().row(0), &[Z2.zero(), Z2.one(), Z2.zero()]);
        let quotient = algebra
            .quotient_by_ideal_basis(radical.basis().clone())
            .unwrap();
        assert_eq!(quotient.algebra().dimension(), 2);
        assert_eq!(
            quotient.lifts(),
            &DenseMatrix::from_rows(
                vec![
                    vec![Z2.one(), Z2.zero(), Z2.zero()],
                    vec![Z2.zero(), Z2.zero(), Z2.one()]
                ],
                Z2
            )
            .unwrap()
        );

        let rational = upper_triangular_two_by_two(Q);
        assert_eq!(
            rational.radical().unwrap().basis().row(0),
            &[Q.nth(0.into()), Q.nth(1.into()), Q.nth(0.into())]
        );
    }

    fn check_product_idempotents<F: QpaCoefficientField>(algebra: &StructureConstantsAlgebra<F>) {
        let field = algebra.field();
        let zero = field.zero();
        let one = field.one();
        assert_eq!(
            algebra.central_idempotents_split_prime().unwrap(),
            vec![
                vec![zero.clone(), zero.clone(), one.clone()],
                vec![zero.clone(), one.clone(), zero.clone()],
                vec![one, zero.clone(), zero],
            ]
        );
    }

    #[test]
    fn central_idempotents_follow_small_and_large_gap_splitting_branches() {
        check_product_idempotents(&product_of_fields(Z2, 3));
        check_product_idempotents(&product_of_fields(Zp::new(5), 3));
        check_product_idempotents(&product_of_fields(Zp::new(19), 3));

        let triangular = upper_triangular_two_by_two(Zp::new(5));
        let field = triangular.field();
        assert_eq!(
            triangular.central_idempotents_split_prime().unwrap(),
            vec![vec![field.one(), field.zero(), field.one()]]
        );
    }

    #[test]
    fn pinned_qpa_radical_properties_and_blocks_match() {
        let fixture: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_module_preprojective_structure.schema-v1.json"
        ))
        .unwrap();
        assert_eq!(
            algebra_properties_json(&upper_triangular_two_by_two(Z2), &|coefficient| {
                Z2.to_integer(coefficient)
                    .to_string()
                    .parse::<i64>()
                    .unwrap()
                    .into()
            }),
            fixture["algebra_properties"]["upper_triangular_f2"]
        );
        let field = Zp::new(5);
        assert_eq!(
            algebra_properties_json(&dual_numbers(field.clone()), &|coefficient| {
                field
                    .to_integer(coefficient)
                    .to_string()
                    .parse::<i64>()
                    .unwrap()
                    .into()
            }),
            fixture["algebra_properties"]["dual_numbers_f5"]
        );
        assert_eq!(
            algebra_properties_json(&product_of_fields(Z2, 2), &|coefficient| {
                Z2.to_integer(coefficient)
                    .to_string()
                    .parse::<i64>()
                    .unwrap()
                    .into()
            }),
            fixture["algebra_properties"]["product_fields_f2"]
        );
        assert_eq!(dual_numbers(Q).is_basic_algebra().unwrap(), None);
        assert_eq!(dual_numbers(Q).is_elementary_algebra().unwrap(), None);
    }

    #[test]
    fn pinned_qpa_power_subalgebra_matches() {
        let fixture: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/algebra_structure_constants.schema-v1.json"
        ))
        .unwrap();
        let f5 = Zp::new(5);
        let f5_json = |coefficient: &_| -> Value {
            f5.to_integer(coefficient)
                .to_string()
                .parse::<i64>()
                .unwrap()
                .into()
        };
        let dual_power = dual_numbers(f5.clone()).power_subalgebra(1).unwrap();
        assert_eq!(
            dense_rows_json(&f5, dual_power.basis(), &f5_json),
            fixture["power_subalgebra"]["dual_numbers_f5"]
        );
        let triangular_power = upper_triangular_two_by_two(Z2).power_subalgebra(1).unwrap();
        assert_eq!(
            dense_rows_json(&Z2, triangular_power.basis(), &|coefficient| {
                Z2.to_integer(coefficient)
                    .to_string()
                    .parse::<i64>()
                    .unwrap()
                    .into()
            }),
            fixture["power_subalgebra"]["upper_triangular_f2"]
        );
        let product_power = product_of_fields(Z2, 2).power_subalgebra(1).unwrap();
        assert_eq!(
            dense_rows_json(&Z2, product_power.basis(), &|coefficient| {
                Z2.to_integer(coefficient)
                    .to_string()
                    .parse::<i64>()
                    .unwrap()
                    .into()
            }),
            fixture["power_subalgebra"]["product_f2"]
        );

        assert!(matches!(
            dual_numbers(f5.clone()).power_subalgebra(2),
            Err(
                super::StructureConstantsError::PowerSubalgebraDegreeNotDivisor {
                    degree: 2,
                    field_degree: 1
                }
            )
        ));
        assert!(matches!(
            dual_numbers(Q).power_subalgebra(1),
            Err(super::StructureConstantsError::FiniteFieldRequired)
        ));
    }

    #[test]
    fn pinned_qpa_idempotent_lifting_matches() {
        let fixture: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/algebra_structure_constants.schema-v1.json"
        ))
        .unwrap();
        let f5 = Zp::new(5);
        let f5_json = |coefficient: &_| -> Value {
            f5.to_integer(coefficient)
                .to_string()
                .parse::<i64>()
                .unwrap()
                .into()
        };
        let shifted = shifted_dual_numbers(f5.clone());
        let ideal = DenseMatrix::from_rows(
            vec![vec![shifted.field().zero(), shifted.field().one()]],
            f5.clone(),
        )
        .unwrap();
        let quotient = shifted.quotient_by_ideal_basis(ideal).unwrap();
        let quotient_idempotent = vec![f5.one()];
        let expected = &fixture["idempotent_lifting"]["shifted_dual_numbers_f5"];
        assert_eq!(
            json!(quotient_idempotent.iter().map(&f5_json).collect::<Vec<_>>()),
            expected["quotient_idempotent"]
        );
        let lifting = quotient
            .lifting_idempotent(&quotient_idempotent)
            .unwrap()
            .unwrap();
        assert_eq!(
            json!(lifting.iter().map(&f5_json).collect::<Vec<_>>()),
            expected["lifting"]
        );
        let lift = quotient.lift_idempotent(&quotient_idempotent).unwrap();
        assert_eq!(
            json!(lift.iter().map(&f5_json).collect::<Vec<_>>()),
            expected["lift"]
        );
        let complete = quotient
            .lift_complete_set_of_orthogonal_idempotents(&[quotient_idempotent])
            .unwrap();
        assert_eq!(
            json!(
                complete
                    .iter()
                    .map(|row| row.iter().map(&f5_json).collect::<Vec<_>>())
                    .collect::<Vec<_>>()
            ),
            expected["complete"]
        );

        let shifted_triangular = shifted_upper_triangular_two_by_two(f5.clone());
        let ideal = DenseMatrix::from_rows(
            vec![vec![
                shifted_triangular.field().zero(),
                shifted_triangular.field().zero(),
                shifted_triangular.field().one(),
            ]],
            f5.clone(),
        )
        .unwrap();
        let quotient = shifted_triangular.quotient_by_ideal_basis(ideal).unwrap();
        let idempotents = vec![vec![f5.zero(), f5.one()], vec![f5.one(), f5.zero()]];
        let expected = &fixture["idempotent_lifting"]["shifted_upper_triangular_f5"];
        assert_eq!(
            json!(
                idempotents
                    .iter()
                    .map(|row| row.iter().map(&f5_json).collect::<Vec<_>>())
                    .collect::<Vec<_>>()
            ),
            expected["idempotents"]
        );
        let first_lift = vec![f5.zero(), f5.one(), f5.zero()];
        let pair = quotient
            .lift_two_orthogonal_idempotents(&first_lift, &idempotents[1])
            .unwrap();
        assert_eq!(
            json!(
                pair.iter()
                    .map(|row| row.iter().map(&f5_json).collect::<Vec<_>>())
                    .collect::<Vec<_>>()
            ),
            expected["pair"]
        );
        let complete = quotient
            .lift_complete_set_of_orthogonal_idempotents(&idempotents)
            .unwrap();
        assert_eq!(
            json!(
                complete
                    .iter()
                    .map(|row| row.iter().map(&f5_json).collect::<Vec<_>>())
                    .collect::<Vec<_>>()
            ),
            expected["complete"]
        );
    }

    #[test]
    fn pinned_qpa_generic_opposite_algebra_matches() {
        let fixture: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/algebra_opposite.schema-v1.json"
        ))
        .unwrap();
        let field = Zp::new(5);
        let coefficient_json = |coefficient: &_| -> Value {
            field
                .to_integer(coefficient)
                .to_string()
                .parse::<i64>()
                .unwrap()
                .into()
        };
        let algebra = upper_triangular_two_by_two(field.clone());
        let underlying_basis = algebra.basis_elements();
        let opposite = algebra.opposite().unwrap();
        let opposite2 = algebra.opposite().unwrap();
        let basis = opposite.basis_elements();
        let canonical_basis = opposite.canonical_basis_elements();
        let [x, y, z] = [&basis[0], &basis[1], &basis[2]];
        let recreated = opposite
            .element_from_underlying(underlying_basis[1].clone())
            .unwrap();
        algebra
            .bind_named_element("y", &underlying_basis[1])
            .unwrap();
        let named = opposite.element_by_name("y").unwrap().unwrap();
        let products = [
            x.multiply(y).unwrap(),
            y.multiply(x).unwrap(),
            y.multiply(z).unwrap(),
            z.multiply(y).unwrap(),
        ];
        let upper = json!({
            "is_opposite": opposite.is_opposite_algebra(),
            "is_algebra": true,
            "underlying_with_one": false,
            "opposite_with_one": false,
            "underlying_identical": opposite.underlying_algebra() == &algebra,
            "fresh": opposite != opposite2,
            "dimension": opposite.dimension(),
            "finite_dimensional": opposite.is_finite_dimensional(),
            "basis": basis.iter().map(|element| row_json(&field, element.coordinates(), &coefficient_json)).collect::<Vec<_>>(),
            "canonical_basis": canonical_basis.iter().map(|element| row_json(&field, element.coordinates(), &coefficient_json)).collect::<Vec<_>>(),
            "coefficients": row_json(&field, &opposite.coefficients(y).unwrap(), &coefficient_json),
            "member": opposite.contains(y),
            "products": products.iter().map(|element| row_json(&field, element.coordinates(), &coefficient_json)).collect::<Vec<_>>(),
            "sum": row_json(&field, y.add(z).unwrap().coordinates(), &coefficient_json),
            "left_scalar": row_json(&field, y.scale(&field.nth(2.into())).unwrap().coordinates(), &coefficient_json),
            "right_scalar": row_json(&field, y.scale(&field.nth(2.into())).unwrap().coordinates(), &coefficient_json),
            "negative": row_json(&field, y.negate().unwrap().coordinates(), &coefficient_json),
            "zero": row_json(&field, y.zero_same_parent().coordinates(), &coefficient_json),
            "equal": recreated == *y,
            "less": [x.compare(y).unwrap() == Ordering::Less, y.compare(x).unwrap() == Ordering::Less],
            "underlying_less": [underlying_basis[0].coordinates() < underlying_basis[1].coordinates(), underlying_basis[1].coordinates() < underlying_basis[0].coordinates()],
            "external_equal": y.external_representation() == underlying_basis[1].coordinates(),
            "recreated": row_json(&field, recreated.coordinates(), &coefficient_json),
            "element_print": y.to_string(),
            "sum_print": y.add(z).unwrap().to_string(),
            "scalar_print": y.scale(&field.nth(2.into())).unwrap().to_string(),
            "negative_print": y.negate().unwrap().to_string(),
            "zero_print": y.zero_same_parent().to_string()
        });
        assert_eq!(upper, fixture["upper_triangular_f5"]);
        assert_eq!(named, *y);
        assert!(opposite.element_by_name("missing").unwrap().is_none());
        assert_eq!(
            opposite.to_string(),
            "OppositeAlgebra( <structure-constants algebra of dimension 3> )"
        );

        let rational = upper_triangular_two_by_two(Q).opposite().unwrap();
        let rational_basis = rational.basis_elements();
        let [r1, r2, _] = [&rational_basis[0], &rational_basis[1], &rational_basis[2]];
        assert_eq!(
            json!({
                "sum": r1.add(r2).unwrap().to_string(),
                "double": r2.scale(&Q.nth(2.into())).unwrap().to_string(),
                "half": r2
                    .scale(&Q.div(&Q.one(), &Q.nth(2.into())))
                    .unwrap()
                    .to_string(),
                "negative": r2.negate().unwrap().to_string(),
                "mixed": r1
                    .add(&r2.scale(&Q.nth((-2).into())).unwrap())
                    .unwrap()
                    .to_string(),
                "zero": r1.zero_same_parent().to_string()
            }),
            fixture["rational_print"]
        );

        let commutative = dual_numbers(field).opposite().unwrap();
        assert_eq!(
            json!({
                "identical": commutative.is_identical_to_underlying(),
                "has_is_opposite": commutative.is_opposite_algebra()
            }),
            fixture["dual_numbers_f5"]
        );
    }
}
