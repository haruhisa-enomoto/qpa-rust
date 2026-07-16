//! QPA's representation-type decision methods.
//!
//! This module translates `lib/pathalgreptype.gi` one installed overload at a
//! time. Methods whose source path reaches modules remain deferred until that
//! exact substrate exists.

use numerica::domains::Field;
use qpa_core::{
    AlgebraDimension, DenseMatrix, MatrixError, PathAlgebra, QpaCoefficientField, Quiver,
    QuiverAlgebra, QuiverError, QuotientError, QuotientPathAlgebra, StructureConstantsAlgebra,
    StructureConstantsError,
};
use thiserror::Error;

pub use qpa_core::DYNKIN_UNORIENTED_CYCLE_MESSAGE as DYNKIN_ORIENTED_CYCLE_MESSAGE;

use crate::{
    module::{
        EnvelopingModuleError, d_tr, indec_injective_modules, indec_projective_modules,
        simple_modules, tr_d,
    },
    pathalg::{
        PathalgGbnpError, quotient_dimension, quotient_is_admissible, quotient_is_special_biserial,
    },
};

/// QPA's final semisimple message from the path-algebra finite-type method.
pub const FINITE_TYPE_SEMISIMPLE_MESSAGE: &str =
    "Finite type!\nSemisimple algebra (Q has no arrows).\n";

/// QPA's final Dynkin-union message from the path-algebra finite-type method.
pub const FINITE_TYPE_DYNKIN_MESSAGE: &str =
    "Finite type!\nQuiver is a (union of) Dynkin quiver(s).\n";

/// QPA's final non-Dynkin message from the path-algebra finite-type method.
pub const INFINITE_TYPE_NON_DYNKIN_MESSAGE: &str =
    "Infinite type!\nQuiver is not a (union of) Dynkin quiver(s).\n";

/// QPA's diagnostic when the separated quiver proves infinite type modulo
/// the square of the radical.
pub const INFINITE_TYPE_SEPARATED_QUIVER_MESSAGE: &str =
    "Have checked representation type of the algebra modulo the square of the radical.\n";

/// QPA's finite-type diagnostic for the acyclic special-biserial branch.
pub const FINITE_TYPE_SPECIAL_BISERIAL_MESSAGE: &str =
    "Finite type!\nSpecial biserial algebra with no unoriented cycles in Q.\n";

/// QPA's final diagnostic when every implemented criterion is inconclusive.
pub const FINITE_TYPE_INCONCLUSIVE_MESSAGE: &str = "Can not determine the representation type.\n";

/// QPA's unconditional diagnostic at the start of `BongartzTest` after its
/// admissibility check.
pub const BONGARTZ_ENTERING_MESSAGE: &str = "Entering Bongartz' test for infinite type.\n";

/// QPA's `false`/`fail` result from `BongartzTest`, with its printed message.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BongartzTestDecision {
    value: Option<bool>,
    messages: Vec<Box<str>>,
}

impl BongartzTestDecision {
    /// Return `Some(false)` when the dimension bound proves infinite type, or
    /// `None` for QPA's inconclusive `fail`.
    #[must_use]
    pub const fn value(&self) -> Option<bool> {
        self.value
    }

    /// Return QPA's ordered printed side effects.
    #[must_use]
    pub fn messages(&self) -> &[Box<str>] {
        &self.messages
    }
}

/// The source branch deciding finite or infinite path-algebra type.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum PathAlgebraFiniteTypeReason {
    /// QPA returned before component construction because the quiver has no arrows.
    Semisimple,
    /// Every connected component was recognized as Dynkin.
    DynkinComponents,
    /// `ForAll` stopped at the first component not recognized as Dynkin.
    NonDynkinComponent,
}

/// QPA's Boolean result together with its ordered printed side effects.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PathAlgebraFiniteTypeDecision {
    finite: bool,
    reason: PathAlgebraFiniteTypeReason,
    messages: Vec<Box<str>>,
}

/// The source branch deciding a quotient path algebra's representation type.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum QuotientAlgebraFiniteTypeReason {
    /// QPA returned before component construction because the quiver has no arrows.
    Semisimple,
    /// Every connected component of the original quiver was Dynkin.
    DynkinComponents,
    /// A separated-quiver component was not Dynkin.
    SeparatedQuiverNonDynkin,
    /// The separated quiver was Dynkin and the algebra had radical square zero.
    RadicalSquareZero,
    /// Jans' distributivity criterion rejected finite representation type.
    NonDistributive,
    /// The algebra was special biserial with an unoriented-acyclic quiver.
    AcyclicSpecialBiserial,
    /// Bongartz' bounded test proved infinite type.
    BongartzInfinite,
    /// Every implemented source criterion, including Bongartz, was inconclusive.
    Inconclusive,
}

/// QPA's `true`/`false`/`fail` quotient decision and ordered printed effects.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct QuotientAlgebraFiniteTypeDecision {
    value: Option<bool>,
    reason: QuotientAlgebraFiniteTypeReason,
    messages: Vec<Box<str>>,
}

struct QuotientStructureModel<F: Field> {
    algebra: StructureConstantsAlgebra<F>,
    primitive_idempotents: Vec<Vec<F::Element>>,
}

impl QuotientAlgebraFiniteTypeDecision {
    /// Return `Some(true)`, `Some(false)`, or `None` for QPA's `fail`.
    #[must_use]
    pub const fn value(&self) -> Option<bool> {
        self.value
    }

    /// Return the source branch that produced the result.
    #[must_use]
    pub const fn reason(&self) -> QuotientAlgebraFiniteTypeReason {
        self.reason
    }

    /// Return nested Dynkin/Bongartz output and the final branch diagnostic in
    /// QPA call order.
    #[must_use]
    pub fn messages(&self) -> &[Box<str>] {
        &self.messages
    }
}

impl PathAlgebraFiniteTypeDecision {
    /// Return QPA's cached `IsFiniteTypeAlgebra` Boolean value.
    #[must_use]
    pub const fn is_finite(&self) -> bool {
        self.finite
    }

    /// Return the source branch that produced the Boolean.
    #[must_use]
    pub const fn reason(&self) -> PathAlgebraFiniteTypeReason {
        self.reason
    }

    /// Return every message printed by nested `IsDynkinQuiver` calls and the
    /// final method branch, in source order.
    #[must_use]
    pub fn messages(&self) -> &[Box<str>] {
        &self.messages
    }
}

/// Translate QPA's `IsFiniteTypeAlgebra` overload for path algebras.
///
/// The source first recognizes arrow-free semisimple algebras. Otherwise it
/// builds QPA-ordered connected components and calls `IsDynkinQuiver` through
/// GAP's short-circuiting `ForAll`.
///
/// # Errors
///
/// Propagates the source-translated full-subquiver construction failures used
/// by `ConnectedComponentsOfQuiver`.
pub fn path_algebra_finite_type<F: Field>(
    algebra: &PathAlgebra<F>,
) -> Result<PathAlgebraFiniteTypeDecision, QuiverError> {
    let quiver = algebra.quiver();
    if quiver.number_of_arrows() == 0 {
        return Ok(PathAlgebraFiniteTypeDecision {
            finite: true,
            reason: PathAlgebraFiniteTypeReason::Semisimple,
            messages: vec![FINITE_TYPE_SEMISIMPLE_MESSAGE.into()],
        });
    }

    let mut messages = Vec::new();
    for component in quiver.connected_components()? {
        let dynkin = component.dynkin_decision();
        messages.extend(dynkin.messages().iter().cloned());
        if dynkin.is_dynkin() {
            continue;
        }
        messages.push(INFINITE_TYPE_NON_DYNKIN_MESSAGE.into());
        return Ok(PathAlgebraFiniteTypeDecision {
            finite: false,
            reason: PathAlgebraFiniteTypeReason::NonDynkinComponent,
            messages,
        });
    }

    messages.push(FINITE_TYPE_DYNKIN_MESSAGE.into());
    Ok(PathAlgebraFiniteTypeDecision {
        finite: true,
        reason: PathAlgebraFiniteTypeReason::DynkinComponents,
        messages,
    })
}

/// Translate QPA's `IsDistributiveAlgebra` fallback for path algebras.
///
/// The installed method delegates exactly to `IsTreeQuiver` on the defining
/// quiver; unlike the quotient overload, it performs no finite-dimensionality
/// or radical computation.
#[must_use]
pub fn path_algebra_is_distributive<F: Field>(algebra: &PathAlgebra<F>) -> bool {
    algebra.quiver().is_tree()
}

/// Translate QPA's `IsDistributiveAlgebra` overload for an admissible quotient.
///
/// The implementation follows the source's primitive-vertex corner rings,
/// radical-series layer test, and left-then-right uniserial tests for every
/// off-diagonal corner. Quotient multiplication is first written in QPA's
/// canonical basis and every subsequent subspace uses GAP's semi-echelon row
/// order.
///
/// # Errors
///
/// Reproduces the finite-dimensional and admissible preconditions and
/// propagates exact quotient, matrix, structure-constant, and field failures.
pub fn quotient_is_distributive<F: QpaCoefficientField>(
    quotient: &QuotientPathAlgebra<F>,
) -> Result<bool, RepresentationTypeError> {
    if !matches!(quotient_dimension(quotient)?, AlgebraDimension::Finite(_)) {
        return Err(RepresentationTypeError::DistributiveAlgebraInfinite);
    }
    if !quotient_is_admissible(quotient)? {
        return Err(RepresentationTypeError::DistributiveRequiresAdmissibleQuotient);
    }
    let model = quotient_structure_constants(quotient)?;
    structure_algebra_is_distributive(&model.algebra, &model.primitive_idempotents)
}

/// Translate QPA's `IsFiniteTypeAlgebra` overload for quotient path algebras.
///
/// The source order is preserved: admissibility, semisimplicity, original
/// Dynkin components, separated-quiver components, radical-square-zero,
/// distributivity, acyclic special biseriality, and finally Bongartz with
/// bound 100. Nested printed diagnostics are retained in call order.
///
/// # Errors
///
/// QPA accepts only admissible quotients. Exact quiver, quotient,
/// structure-constant, distributivity, and AR-functor failures propagate.
pub fn quotient_algebra_finite_type<F: QpaCoefficientField>(
    quotient: &QuotientPathAlgebra<F>,
) -> Result<QuotientAlgebraFiniteTypeDecision, RepresentationTypeError> {
    if !quotient_is_admissible(quotient)? {
        return Err(RepresentationTypeError::FiniteTypeRequiresAdmissibleQuotient);
    }
    let quiver = quotient.original_path_algebra().quiver();
    if quiver.number_of_arrows() == 0 {
        return Ok(QuotientAlgebraFiniteTypeDecision {
            value: Some(true),
            reason: QuotientAlgebraFiniteTypeReason::Semisimple,
            messages: vec![FINITE_TYPE_SEMISIMPLE_MESSAGE.into()],
        });
    }

    let mut messages = Vec::new();
    if quiver_components_are_dynkin(quiver, &mut messages)? {
        messages.push(FINITE_TYPE_DYNKIN_MESSAGE.into());
        return Ok(QuotientAlgebraFiniteTypeDecision {
            value: Some(true),
            reason: QuotientAlgebraFiniteTypeReason::DynkinComponents,
            messages,
        });
    }

    let separated = quiver.separated()?;
    if !quiver_components_are_dynkin(&separated, &mut messages)? {
        messages.push(INFINITE_TYPE_SEPARATED_QUIVER_MESSAGE.into());
        return Ok(QuotientAlgebraFiniteTypeDecision {
            value: Some(false),
            reason: QuotientAlgebraFiniteTypeReason::SeparatedQuiverNonDynkin,
            messages,
        });
    }

    let model = quotient_structure_constants(quotient)?;
    if model.algebra.is_radical_square_zero()? {
        return Ok(QuotientAlgebraFiniteTypeDecision {
            value: Some(true),
            reason: QuotientAlgebraFiniteTypeReason::RadicalSquareZero,
            messages,
        });
    }
    if !structure_algebra_is_distributive(&model.algebra, &model.primitive_idempotents)? {
        return Ok(QuotientAlgebraFiniteTypeDecision {
            value: Some(false),
            reason: QuotientAlgebraFiniteTypeReason::NonDistributive,
            messages,
        });
    }
    if quotient_is_special_biserial(quotient)? && quiver.is_unoriented_acyclic() {
        messages.push(FINITE_TYPE_SPECIAL_BISERIAL_MESSAGE.into());
        return Ok(QuotientAlgebraFiniteTypeDecision {
            value: Some(true),
            reason: QuotientAlgebraFiniteTypeReason::AcyclicSpecialBiserial,
            messages,
        });
    }

    let bongartz = bongartz_test(&QuiverAlgebra::Quotient(quotient.clone()), 100)?;
    messages.extend(bongartz.messages().iter().cloned());
    if bongartz.value() == Some(false) {
        return Ok(QuotientAlgebraFiniteTypeDecision {
            value: Some(false),
            reason: QuotientAlgebraFiniteTypeReason::BongartzInfinite,
            messages,
        });
    }
    messages.push(FINITE_TYPE_INCONCLUSIVE_MESSAGE.into());
    Ok(QuotientAlgebraFiniteTypeDecision {
        value: None,
        reason: QuotientAlgebraFiniteTypeReason::Inconclusive,
        messages,
    })
}

fn quiver_components_are_dynkin(
    quiver: &Quiver,
    messages: &mut Vec<Box<str>>,
) -> Result<bool, QuiverError> {
    for component in quiver.connected_components()? {
        let decision = component.dynkin_decision();
        messages.extend(decision.messages().iter().cloned());
        if !decision.is_dynkin() {
            return Ok(false);
        }
    }
    Ok(true)
}

fn quotient_structure_constants<F: QpaCoefficientField>(
    quotient: &QuotientPathAlgebra<F>,
) -> Result<QuotientStructureModel<F>, RepresentationTypeError> {
    let basis = quotient
        .canonical_basis()?
        .ok_or(RepresentationTypeError::MissingFiniteQuotientBasis)?;
    let dimension = basis.len();
    let mut products = Vec::new();
    for left in &basis {
        for right in &basis {
            let product = left.multiply(right)?;
            products.extend(quotient.canonical_basis_coefficients(&product)?);
        }
    }
    let algebra = StructureConstantsAlgebra::new(
        quotient.original_path_algebra().field().clone(),
        dimension,
        products,
    )?;
    let primitive_idempotents = quotient
        .original_path_algebra()
        .quiver()
        .vertices()
        .map(|vertex| {
            let element = quotient.path(vertex.path())?;
            Ok(quotient.canonical_basis_coefficients(&element)?)
        })
        .collect::<Result<Vec<_>, RepresentationTypeError>>()?;
    Ok(QuotientStructureModel {
        algebra,
        primitive_idempotents,
    })
}

fn structure_algebra_is_distributive<F: QpaCoefficientField>(
    algebra: &StructureConstantsAlgebra<F>,
    primitive_idempotents: &[Vec<F::Element>],
) -> Result<bool, RepresentationTypeError> {
    let local_rings = primitive_idempotents
        .iter()
        .map(|idempotent| {
            let basis = corner_basis(algebra, idempotent, idempotent)?;
            let local = induced_structure_algebra(algebra, &basis)?;
            Ok((basis, local))
        })
        .collect::<Result<Vec<_>, RepresentationTypeError>>()?;

    for (_, local) in &local_rings {
        let series = local.radical_series()?;
        if series
            .windows(2)
            .any(|pair| pair[0].nrows().saturating_sub(pair[1].nrows()) != 1)
        {
            return Ok(false);
        }
    }

    let radical_bases = local_rings
        .iter()
        .map(|(ambient_basis, local)| Ok(local.radical()?.basis().multiply(ambient_basis)?))
        .collect::<Result<Vec<_>, RepresentationTypeError>>()?;
    for (left_index, left_idempotent) in primitive_idempotents.iter().enumerate() {
        for (right_index, right_idempotent) in primitive_idempotents.iter().enumerate() {
            if left_index == right_index {
                continue;
            }
            let corner = corner_basis(algebra, left_idempotent, right_idempotent)?;
            if corner_is_uniserial(algebra, &corner, &radical_bases[left_index], true)? {
                continue;
            }
            if !corner_is_uniserial(algebra, &corner, &radical_bases[right_index], false)? {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

fn corner_basis<F: Field>(
    algebra: &StructureConstantsAlgebra<F>,
    left: &[F::Element],
    right: &[F::Element],
) -> Result<DenseMatrix<F>, RepresentationTypeError> {
    let mut rows = Vec::with_capacity(algebra.dimension());
    for basis_index in 0..algebra.dimension() {
        let mut basis = vec![algebra.field().zero(); algebra.dimension()];
        basis[basis_index] = algebra.field().one();
        let product = algebra.multiply_coordinates(left, &basis)?;
        rows.push(algebra.multiply_coordinates(&product, right)?);
    }
    row_space_basis(algebra.field(), algebra.dimension(), rows)
}

fn induced_structure_algebra<F: Field>(
    ambient: &StructureConstantsAlgebra<F>,
    basis: &DenseMatrix<F>,
) -> Result<StructureConstantsAlgebra<F>, RepresentationTypeError> {
    let mut products = Vec::new();
    for left in 0..basis.nrows() {
        for right in 0..basis.nrows() {
            let product = ambient.multiply_coordinates(basis.row(left), basis.row(right))?;
            products.extend(
                basis
                    .solution_row(&product)?
                    .ok_or(RepresentationTypeError::DistributiveSubspaceNotClosed)?,
            );
        }
    }
    Ok(StructureConstantsAlgebra::new(
        ambient.field().clone(),
        basis.nrows(),
        products,
    )?)
}

fn corner_is_uniserial<F: Field>(
    algebra: &StructureConstantsAlgebra<F>,
    corner: &DenseMatrix<F>,
    radical: &DenseMatrix<F>,
    radical_acts_on_left: bool,
) -> Result<bool, RepresentationTypeError> {
    let mut current = corner.clone();
    while current.nrows() != 0 {
        let mut products = Vec::new();
        for row in 0..current.nrows() {
            for radical_row in 0..radical.nrows() {
                let product = if radical_acts_on_left {
                    algebra.multiply_coordinates(radical.row(radical_row), current.row(row))?
                } else {
                    algebra.multiply_coordinates(current.row(row), radical.row(radical_row))?
                };
                if product
                    .iter()
                    .any(|coefficient| !algebra.field().is_zero(coefficient))
                {
                    products.push(product);
                }
            }
        }
        let next = row_space_basis(algebra.field(), algebra.dimension(), products)?;
        if next.nrows() == 0 {
            if current.nrows() != 1 {
                return Ok(false);
            }
        } else if current.nrows().saturating_sub(next.nrows()) != 1 {
            return Ok(false);
        }
        current = next;
    }
    Ok(true)
}

fn row_space_basis<F: Field>(
    field: &F,
    width: usize,
    rows: Vec<Vec<F::Element>>,
) -> Result<DenseMatrix<F>, RepresentationTypeError> {
    let row_count = rows.len();
    Ok(DenseMatrix::from_flat(
        row_count,
        width,
        rows.into_iter().flatten().collect(),
        field.clone(),
    )?
    .semi_echelon_transformation()
    .vectors)
}

/// Translate QPA's `BongartzTest(A, bound)` method.
///
/// The method initializes simples, projectives, and injectives in source
/// order, sets `maxdim = Maximum(2 * Dimension(A), 29)`, and then updates the
/// four families `DTr(S)`, `TrD(S)`, `TrD(P)`, and `DTr(I)` one step at a time.
/// The first module exceeding `maxdim` returns `Some(false)`; exhausting the
/// integer bound (including a nonpositive bound) returns `None`, QPA's `fail`.
///
/// # Errors
///
/// QPA accepts only an admissible quotient of a path algebra. This function
/// propagates the configured GBNP admissibility/dimension path, projective and
/// injective construction, and exact AR-functor computations.
pub fn bongartz_test<F: Field>(
    algebra: &QuiverAlgebra<F>,
    bound: i64,
) -> Result<BongartzTestDecision, RepresentationTypeError> {
    let QuiverAlgebra::Quotient(quotient) = algebra else {
        return Err(RepresentationTypeError::BongartzRequiresAdmissibleQuotient);
    };
    if !quotient_is_admissible(quotient)? {
        return Err(RepresentationTypeError::BongartzRequiresAdmissibleQuotient);
    }

    let simples = simple_modules(algebra)?;
    let projectives = indec_projective_modules(algebra)?;
    let injectives = indec_injective_modules(algebra)?;
    let AlgebraDimension::Finite(dimension) = quotient_dimension(quotient)? else {
        return Err(RepresentationTypeError::AdmissibleBongartzAlgebraInfinite);
    };
    let dimension = usize::try_from(dimension)
        .map_err(|_| RepresentationTypeError::BongartzDimensionOverflow)?;
    let max_dimension = dimension
        .checked_mul(2)
        .ok_or(RepresentationTypeError::BongartzDimensionOverflow)?
        .max(29);

    let mut d_tr_simples = simples.clone();
    let mut tr_d_simples = simples;
    let mut tr_d_projectives = projectives;
    let mut d_tr_injectives = injectives;
    let mut completed = 0_i64;
    while completed < bound {
        for index in 0..d_tr_simples.len() {
            d_tr_simples[index] = d_tr(&d_tr_simples[index])?;
            if d_tr_simples[index].dimension() > max_dimension {
                return Ok(BongartzTestDecision {
                    value: Some(false),
                    messages: vec![BONGARTZ_ENTERING_MESSAGE.into()],
                });
            }
            tr_d_simples[index] = tr_d(&tr_d_simples[index])?;
            if tr_d_simples[index].dimension() > max_dimension {
                return Ok(BongartzTestDecision {
                    value: Some(false),
                    messages: vec![BONGARTZ_ENTERING_MESSAGE.into()],
                });
            }
            tr_d_projectives[index] = tr_d(&tr_d_projectives[index])?;
            if tr_d_projectives[index].dimension() > max_dimension {
                return Ok(BongartzTestDecision {
                    value: Some(false),
                    messages: vec![BONGARTZ_ENTERING_MESSAGE.into()],
                });
            }
            d_tr_injectives[index] = d_tr(&d_tr_injectives[index])?;
            if d_tr_injectives[index].dimension() > max_dimension {
                return Ok(BongartzTestDecision {
                    value: Some(false),
                    messages: vec![BONGARTZ_ENTERING_MESSAGE.into()],
                });
            }
        }
        completed += 1;
    }
    Ok(BongartzTestDecision {
        value: None,
        messages: vec![BONGARTZ_ENTERING_MESSAGE.into()],
    })
}

/// Failure in QPA's representation-type decision methods.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum RepresentationTypeError {
    /// The quotient finite-type overload requires QPA's admissibility filter.
    #[error("the entered algebra is not an admissible quotient of a path algebra")]
    FiniteTypeRequiresAdmissibleQuotient,
    /// QPA's quotient distributivity method errors on an infinite algebra.
    #[error("the entered algebra is not finite dimensional")]
    DistributiveAlgebraInfinite,
    /// QPA tries another method when its quotient is not admissible.
    #[error("the quotient is not an admissible quotient of a path algebra")]
    DistributiveRequiresAdmissibleQuotient,
    /// An admissible quotient must expose its finite canonical basis.
    #[error("the admissible quotient has no finite canonical basis")]
    MissingFiniteQuotientBasis,
    /// A corner-space basis must remain closed under multiplication.
    #[error("a distributivity-test corner subspace is not closed under multiplication")]
    DistributiveSubspaceNotClosed,
    /// `BongartzTest` dispatches only on admissible path-algebra quotients.
    #[error("the entered algebra is not an admissible quotient of a path algebra")]
    BongartzRequiresAdmissibleQuotient,
    /// Admissibility should imply finite dimension in this method domain.
    #[error("an admissible Bongartz-test algebra unexpectedly has infinite dimension")]
    AdmissibleBongartzAlgebraInfinite,
    /// Twice the algebra dimension must fit Rust's module-dimension type.
    #[error("the Bongartz-test dimension bound overflows usize")]
    BongartzDimensionOverflow,
    /// Configured QPA/GBNP admissibility or dimension computation failed.
    #[error(transparent)]
    Backend(#[from] PathalgGbnpError),
    /// Exact matrix/subspace construction failed.
    #[error(transparent)]
    Matrix(#[from] MatrixError),
    /// Exact quotient construction, normalization, or coordinates failed.
    #[error(transparent)]
    Quotient(#[from] QuotientError),
    /// Exact structure-constant conversion or radical computation failed.
    #[error(transparent)]
    StructureConstants(#[from] StructureConstantsError),
    /// QPA component or separated-quiver construction failed.
    #[error(transparent)]
    Quiver(#[from] QuiverError),
    /// Projective, injective, or AR-functor construction failed.
    #[error(transparent)]
    Module(#[from] EnvelopingModuleError),
}

#[cfg(test)]
mod tests {
    use numerica::domains::{Field, finite_field::Z2, rational::Q};
    use qpa_core::{
        ArrowSpec, NamedArrowSpec, PathAlgebra, Quiver, QuiverAlgebra, QuotientPathAlgebra,
        VertexIndex,
    };
    use serde_json::{Value, json};

    use super::{
        BONGARTZ_ENTERING_MESSAGE, DYNKIN_ORIENTED_CYCLE_MESSAGE, FINITE_TYPE_DYNKIN_MESSAGE,
        FINITE_TYPE_INCONCLUSIVE_MESSAGE, FINITE_TYPE_SEMISIMPLE_MESSAGE,
        FINITE_TYPE_SPECIAL_BISERIAL_MESSAGE, INFINITE_TYPE_NON_DYNKIN_MESSAGE,
        INFINITE_TYPE_SEPARATED_QUIVER_MESSAGE, PathAlgebraFiniteTypeReason,
        QuotientAlgebraFiniteTypeDecision, QuotientAlgebraFiniteTypeReason,
        RepresentationTypeError, bongartz_test, path_algebra_finite_type,
        path_algebra_is_distributive, quotient_algebra_finite_type, quotient_is_distributive,
    };
    use crate::quotient_by_relations;

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn algebra<F: Field>(
        field: F,
        number_of_vertices: usize,
        arrows: &[(usize, usize, &str)],
    ) -> PathAlgebra<F> {
        let arrows = arrows
            .iter()
            .map(|(source, target, name)| ArrowSpec::named(vertex(*source), vertex(*target), *name))
            .collect();
        PathAlgebra::new(
            field,
            Quiver::from_indexed(number_of_vertices, arrows).unwrap(),
        )
        .unwrap()
    }

    fn messages(decision: &super::PathAlgebraFiniteTypeDecision) -> Vec<&str> {
        decision.messages().iter().map(Box::as_ref).collect()
    }

    #[test]
    fn source_branches_and_for_all_short_circuit_are_explicit() {
        let semisimple = PathAlgebra::new(Q, Quiver::from_indexed(3, Vec::new()).unwrap()).unwrap();
        let decision = path_algebra_finite_type(&semisimple).unwrap();
        assert!(decision.is_finite());
        assert_eq!(decision.reason(), PathAlgebraFiniteTypeReason::Semisimple);
        assert_eq!(messages(&decision), [FINITE_TYPE_SEMISIMPLE_MESSAGE]);

        let a3 = PathAlgebra::new(
            Z2,
            Quiver::from_indexed(
                3,
                vec![
                    ArrowSpec::named(vertex(1), vertex(2), "a"),
                    ArrowSpec::named(vertex(2), vertex(3), "b"),
                ],
            )
            .unwrap(),
        )
        .unwrap();
        let decision = path_algebra_finite_type(&a3).unwrap();
        assert!(decision.is_finite());
        assert_eq!(
            decision.reason(),
            PathAlgebraFiniteTypeReason::DynkinComponents
        );
        assert_eq!(messages(&decision), ["A_3\n", FINITE_TYPE_DYNKIN_MESSAGE]);

        let first_dynkin_then_cycle = PathAlgebra::new(
            Q,
            Quiver::from_indexed(
                4,
                vec![
                    ArrowSpec::named(vertex(1), vertex(2), "a"),
                    ArrowSpec::named(vertex(3), vertex(4), "b"),
                    ArrowSpec::named(vertex(4), vertex(3), "c"),
                ],
            )
            .unwrap(),
        )
        .unwrap();
        let decision = path_algebra_finite_type(&first_dynkin_then_cycle).unwrap();
        assert!(!decision.is_finite());
        assert_eq!(
            decision.reason(),
            PathAlgebraFiniteTypeReason::NonDynkinComponent
        );
        assert_eq!(
            messages(&decision),
            [
                "A_2\n",
                DYNKIN_ORIENTED_CYCLE_MESSAGE,
                INFINITE_TYPE_NON_DYNKIN_MESSAGE
            ]
        );
    }

    #[test]
    fn pinned_qpa_component_messages_match() {
        let disjoint_dynkin = algebra(Q, 6, &[(1, 2, "a"), (3, 4, "b"), (3, 5, "c"), (3, 6, "d")]);
        let decision = path_algebra_finite_type(&disjoint_dynkin).unwrap();
        assert!(decision.is_finite());
        assert_eq!(
            decision.reason(),
            PathAlgebraFiniteTypeReason::DynkinComponents
        );
        assert_eq!(
            messages(&decision),
            ["A_2\n", "D_4\n", FINITE_TYPE_DYNKIN_MESSAGE]
        );

        let nondynkin_tree = algebra(Z2, 5, &[(1, 2, "a"), (1, 3, "b"), (1, 4, "c"), (1, 5, "d")]);
        let decision = path_algebra_finite_type(&nondynkin_tree).unwrap();
        assert!(!decision.is_finite());
        assert_eq!(messages(&decision), [INFINITE_TYPE_NON_DYNKIN_MESSAGE]);
    }

    #[test]
    fn pinned_qpa_bongartz_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/reptype_bongartz.schema-v1.json"
        ))
        .unwrap();

        let loop_quiver =
            Quiver::from_named(vec!["v".into()], vec![NamedArrowSpec::named("v", "v", "x")])
                .unwrap();
        let loop_parent = PathAlgebra::new(Q, loop_quiver).unwrap();
        let x = loop_parent
            .element(loop_parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let dual_numbers = QuiverAlgebra::Quotient(
            quotient_by_relations(loop_parent.clone(), &[x.multiply(&x).unwrap()]).unwrap(),
        );

        let case = |bound| {
            let decision = bongartz_test(&dual_numbers, bound).unwrap();
            assert_eq!(
                decision.messages(),
                [Box::<str>::from(BONGARTZ_ENTERING_MESSAGE)]
            );
            json!({"success": true, "value": decision.value()})
        };
        assert_eq!(case(0), expected["bound_zero"]);
        assert_eq!(case(1), expected["bound_one"]);
        assert_eq!(case(-1), expected["negative_bound"]);

        let path_input = bongartz_test(&QuiverAlgebra::Path(loop_parent), 0);
        assert_eq!(
            path_input,
            Err(RepresentationTypeError::BongartzRequiresAdmissibleQuotient)
        );
        assert_eq!(
            json!({"success": path_input.is_ok(), "value": Value::Null}),
            expected["path_input"]
        );

        let wild_quiver = Quiver::from_named(
            vec!["w".into()],
            vec![
                NamedArrowSpec::named("w", "w", "x"),
                NamedArrowSpec::named("w", "w", "y"),
            ],
        )
        .unwrap();
        let wild_parent = PathAlgebra::new(Z2, wild_quiver).unwrap();
        let arrows = wild_parent
            .quiver()
            .arrows()
            .map(|arrow| wild_parent.element(arrow.path()).unwrap())
            .collect::<Vec<_>>();
        let relations = vec![
            arrows[0].multiply(&arrows[0]).unwrap(),
            arrows[0].multiply(&arrows[1]).unwrap(),
            arrows[1].multiply(&arrows[0]).unwrap(),
            arrows[1].multiply(&arrows[1]).unwrap(),
        ];
        let wild = QuiverAlgebra::Quotient(quotient_by_relations(wild_parent, &relations).unwrap());
        let wild_decision = bongartz_test(&wild, 20).unwrap();
        assert_eq!(wild_decision.value(), Some(false));
        assert_eq!(
            json!({"success": true, "value": wild_decision.value()}),
            expected["wild"]
        );
    }

    fn quotient_messages(decision: &QuotientAlgebraFiniteTypeDecision) -> Vec<&str> {
        decision.messages().iter().map(Box::as_ref).collect()
    }

    struct QuotientReptypeCases {
        semi: QuotientPathAlgebra<Q>,
        dynkin: QuotientPathAlgebra<Z2>,
        dual_numbers: QuotientPathAlgebra<Q>,
        truncated_cube: QuotientPathAlgebra<Q>,
        nonadmissible: QuotientPathAlgebra<Q>,
        wild: QuotientPathAlgebra<Z2>,
        diamond: QuotientPathAlgebra<Q>,
        star: QuotientPathAlgebra<Q>,
    }

    fn quotient_reptype_cases() -> QuotientReptypeCases {
        let semi_parent = algebra(Q, 1, &[]);
        let semi = quotient_by_relations(semi_parent, &[]).unwrap();
        let dynkin_parent = algebra(Z2, 3, &[(1, 2, "a"), (2, 3, "b")]);
        let dynkin_arrows = dynkin_parent
            .quiver()
            .arrows()
            .map(|arrow| dynkin_parent.element(arrow.path()).unwrap())
            .collect::<Vec<_>>();
        let dynkin = quotient_by_relations(
            dynkin_parent,
            &[dynkin_arrows[0].multiply(&dynkin_arrows[1]).unwrap()],
        )
        .unwrap();
        let loop_parent = algebra(Q, 1, &[(1, 1, "x")]);
        let x = loop_parent
            .element(loop_parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let dual_numbers =
            quotient_by_relations(loop_parent.clone(), &[x.pow(2).unwrap()]).unwrap();
        let truncated_cube =
            quotient_by_relations(loop_parent.clone(), &[x.pow(3).unwrap()]).unwrap();
        let nonadmissible = quotient_by_relations(loop_parent, std::slice::from_ref(&x)).unwrap();
        let wild_parent = algebra(Z2, 1, &[(1, 1, "x"), (1, 1, "y")]);
        let wild_arrows = wild_parent
            .quiver()
            .arrows()
            .map(|arrow| wild_parent.element(arrow.path()).unwrap())
            .collect::<Vec<_>>();
        let wild_relations = [
            wild_arrows[0].multiply(&wild_arrows[0]).unwrap(),
            wild_arrows[0].multiply(&wild_arrows[1]).unwrap(),
            wild_arrows[1].multiply(&wild_arrows[0]).unwrap(),
            wild_arrows[1].multiply(&wild_arrows[1]).unwrap(),
        ];
        let wild = quotient_by_relations(wild_parent, &wild_relations).unwrap();
        let diamond_parent = algebra(
            Q,
            7,
            &[
                (1, 2, "a"),
                (1, 3, "b"),
                (2, 4, "c"),
                (3, 4, "d"),
                (5, 6, "u"),
                (6, 7, "z"),
            ],
        );
        let diamond_arrows = diamond_parent
            .quiver()
            .arrows()
            .map(|arrow| diamond_parent.element(arrow.path()).unwrap())
            .collect::<Vec<_>>();
        let diamond_relation = diamond_arrows[4].multiply(&diamond_arrows[5]).unwrap();
        let diamond = quotient_by_relations(diamond_parent, &[diamond_relation]).unwrap();
        let star_parent = algebra(Q, 5, &[(1, 3, "a"), (2, 3, "b"), (3, 4, "d"), (3, 5, "e")]);
        let star_arrows = star_parent
            .quiver()
            .arrows()
            .map(|arrow| star_parent.element(arrow.path()).unwrap())
            .collect::<Vec<_>>();
        let star_relations = [
            star_arrows[0].multiply(&star_arrows[3]).unwrap(),
            star_arrows[1].multiply(&star_arrows[2]).unwrap(),
        ];
        let star = quotient_by_relations(star_parent, &star_relations).unwrap();

        QuotientReptypeCases {
            semi,
            dynkin,
            dual_numbers,
            truncated_cube,
            nonadmissible,
            wild,
            diamond,
            star,
        }
    }

    struct QuotientReptypeDecisions {
        semisimple: QuotientAlgebraFiniteTypeDecision,
        dynkin: QuotientAlgebraFiniteTypeDecision,
        radical_square_zero: QuotientAlgebraFiniteTypeDecision,
        separated_nondynkin: QuotientAlgebraFiniteTypeDecision,
        nondistributive: QuotientAlgebraFiniteTypeDecision,
        special_biserial: QuotientAlgebraFiniteTypeDecision,
        inconclusive: QuotientAlgebraFiniteTypeDecision,
        invalid: Result<QuotientAlgebraFiniteTypeDecision, RepresentationTypeError>,
    }

    fn quotient_reptype_decisions(cases: &QuotientReptypeCases) -> QuotientReptypeDecisions {
        QuotientReptypeDecisions {
            semisimple: quotient_algebra_finite_type(&cases.semi).unwrap(),
            dynkin: quotient_algebra_finite_type(&cases.dynkin).unwrap(),
            radical_square_zero: quotient_algebra_finite_type(&cases.dual_numbers).unwrap(),
            separated_nondynkin: quotient_algebra_finite_type(&cases.wild).unwrap(),
            nondistributive: quotient_algebra_finite_type(&cases.diamond).unwrap(),
            special_biserial: quotient_algebra_finite_type(&cases.star).unwrap(),
            inconclusive: quotient_algebra_finite_type(&cases.truncated_cube).unwrap(),
            invalid: quotient_algebra_finite_type(&cases.nonadmissible),
        }
    }

    fn assert_quotient_reptype_messages(decisions: &QuotientReptypeDecisions) {
        assert_eq!(
            decisions.semisimple.reason(),
            QuotientAlgebraFiniteTypeReason::Semisimple
        );
        assert_eq!(
            quotient_messages(&decisions.semisimple),
            [FINITE_TYPE_SEMISIMPLE_MESSAGE]
        );
        assert_eq!(
            decisions.dynkin.reason(),
            QuotientAlgebraFiniteTypeReason::DynkinComponents
        );
        assert_eq!(
            quotient_messages(&decisions.dynkin),
            ["A_3\n", FINITE_TYPE_DYNKIN_MESSAGE]
        );
        assert_eq!(
            decisions.radical_square_zero.reason(),
            QuotientAlgebraFiniteTypeReason::RadicalSquareZero
        );
        assert_eq!(
            quotient_messages(&decisions.radical_square_zero),
            [DYNKIN_ORIENTED_CYCLE_MESSAGE, "A_2\n"]
        );
        assert_eq!(
            decisions.separated_nondynkin.reason(),
            QuotientAlgebraFiniteTypeReason::SeparatedQuiverNonDynkin
        );
        assert_eq!(
            quotient_messages(&decisions.separated_nondynkin),
            [
                DYNKIN_ORIENTED_CYCLE_MESSAGE,
                DYNKIN_ORIENTED_CYCLE_MESSAGE,
                INFINITE_TYPE_SEPARATED_QUIVER_MESSAGE
            ]
        );
        assert_eq!(
            decisions.nondistributive.reason(),
            QuotientAlgebraFiniteTypeReason::NonDistributive
        );
        assert_eq!(
            quotient_messages(&decisions.nondistributive),
            [
                DYNKIN_ORIENTED_CYCLE_MESSAGE,
                "A_3\n",
                "A_3\n",
                "A_1\n",
                "A_2\n",
                "A_2\n",
                "A_1\n",
                "A_1\n",
                "A_1\n"
            ]
        );
        assert_eq!(
            decisions.special_biserial.reason(),
            QuotientAlgebraFiniteTypeReason::AcyclicSpecialBiserial
        );
        assert_eq!(
            quotient_messages(&decisions.special_biserial),
            [
                "A_3\n",
                "A_3\n",
                "A_1\n",
                "A_1\n",
                "A_1\n",
                "A_1\n",
                FINITE_TYPE_SPECIAL_BISERIAL_MESSAGE
            ]
        );
        assert_eq!(
            decisions.inconclusive.reason(),
            QuotientAlgebraFiniteTypeReason::Inconclusive
        );
        assert_eq!(
            quotient_messages(&decisions.inconclusive),
            [
                DYNKIN_ORIENTED_CYCLE_MESSAGE,
                "A_2\n",
                BONGARTZ_ENTERING_MESSAGE,
                FINITE_TYPE_INCONCLUSIVE_MESSAGE
            ]
        );
        assert_eq!(
            decisions.invalid,
            Err(RepresentationTypeError::FiniteTypeRequiresAdmissibleQuotient),
        );
    }

    fn quotient_reptype_values(decisions: &QuotientReptypeDecisions) -> Value {
        json!({
            "semisimple": decisions.semisimple.value(),
            "dynkin": decisions.dynkin.value(),
            "radical_square_zero": decisions.radical_square_zero.value(),
            "separated_nondynkin": decisions.separated_nondynkin.value(),
            "nondistributive": decisions.nondistributive.value(),
            "special_biserial": decisions.special_biserial.value(),
            "inconclusive": decisions.inconclusive.value(),
            "nonadmissible": {"success": decisions.invalid.is_ok(), "value": Value::Null}
        })
    }

    #[test]
    fn pinned_qpa_quotient_representation_type_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/reptype_quotient.schema-v1.json"
        ))
        .unwrap();
        let cases = quotient_reptype_cases();
        let decisions = quotient_reptype_decisions(&cases);
        assert_quotient_reptype_messages(&decisions);
        assert_eq!(quotient_reptype_values(&decisions), expected["values"]);
        assert_eq!(
            json!({
                "diamond": quotient_is_distributive(&cases.diamond).unwrap(),
                "star": quotient_is_distributive(&cases.star).unwrap(),
                "truncated_cube": quotient_is_distributive(&cases.truncated_cube).unwrap(),
                "wild_gf2": quotient_is_distributive(&cases.wild).unwrap(),
                "path_tree": path_algebra_is_distributive(cases.star.original_path_algebra()),
                "path_disconnected": path_algebra_is_distributive(
                    cases.diamond.original_path_algebra()
                )
            }),
            expected["distributive"]
        );
        assert_eq!(expected["is_quotient"], json!([true, true]));
    }
}
