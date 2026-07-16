//! QPA's representation-type decision methods.
//!
//! This module translates `lib/pathalgreptype.gi` one installed overload at a
//! time. Methods whose source path reaches modules remain deferred until that
//! exact substrate exists.

use numerica::domains::Field;
use qpa_core::{AlgebraDimension, PathAlgebra, QuiverAlgebra, QuiverError};
use thiserror::Error;

pub use qpa_core::DYNKIN_UNORIENTED_CYCLE_MESSAGE as DYNKIN_ORIENTED_CYCLE_MESSAGE;

use crate::{
    module::{
        EnvelopingModuleError, d_tr, indec_injective_modules, indec_projective_modules,
        simple_modules, tr_d,
    },
    pathalg::{PathalgGbnpError, quotient_dimension, quotient_is_admissible},
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
    /// Projective, injective, or AR-functor construction failed.
    #[error(transparent)]
    Module(#[from] EnvelopingModuleError),
}

#[cfg(test)]
mod tests {
    use numerica::domains::{Field, finite_field::Z2, rational::Q};
    use qpa_core::{ArrowSpec, NamedArrowSpec, PathAlgebra, Quiver, QuiverAlgebra, VertexIndex};
    use serde_json::{Value, json};

    use super::{
        BONGARTZ_ENTERING_MESSAGE, DYNKIN_ORIENTED_CYCLE_MESSAGE, FINITE_TYPE_DYNKIN_MESSAGE,
        FINITE_TYPE_SEMISIMPLE_MESSAGE, INFINITE_TYPE_NON_DYNKIN_MESSAGE,
        PathAlgebraFiniteTypeReason, RepresentationTypeError, bongartz_test,
        path_algebra_finite_type,
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
}
