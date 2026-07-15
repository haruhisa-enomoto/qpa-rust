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
    DenseMatrix, Generator, MatrixError, Path, PathEndpoint, QuiverAlgebra, TensorProductAlgebra,
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
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        Ring,
        finite_field::Zp,
        rational::{Q, Rational},
    };

    use crate::{
        ArrowSpec, DenseMatrix, PathAlgebra, Quiver, QuiverAlgebra, QuotientPathAlgebra,
        VertexIndex,
    };

    use super::{
        NamedArrowMatrix, PathAlgebraModuleError, PathAlgebraModuleMatrix,
        right_module_over_path_algebra,
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
}
