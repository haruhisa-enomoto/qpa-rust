//! QPA's representation-type decision methods.
//!
//! This module translates `lib/pathalgreptype.gi` one installed overload at a
//! time. Methods whose source path reaches modules remain deferred until that
//! exact substrate exists.

use numerica::domains::Field;
use qpa_core::{PathAlgebra, QuiverError};

/// QPA's final semisimple message from the path-algebra finite-type method.
pub const FINITE_TYPE_SEMISIMPLE_MESSAGE: &str =
    "Finite type!\nSemisimple algebra (Q has no arrows).\n";

/// QPA's final Dynkin-union message from the path-algebra finite-type method.
pub const FINITE_TYPE_DYNKIN_MESSAGE: &str =
    "Finite type!\nQuiver is a (union of) Dynkin quiver(s).\n";

/// QPA's final non-Dynkin message from the path-algebra finite-type method.
pub const INFINITE_TYPE_NON_DYNKIN_MESSAGE: &str =
    "Infinite type!\nQuiver is not a (union of) Dynkin quiver(s).\n";

/// QPA's message when `IsDynkinQuiver` rejects a connected component by its
/// arrow count.
pub const DYNKIN_ORIENTED_CYCLE_MESSAGE: &str = "Quiver contains an (un)oriented cycle.\n";

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
        if let Some(dynkin) = component.dynkin_type() {
            messages.push(format!("{}_{}\n", dynkin.family, dynkin.rank).into_boxed_str());
            continue;
        }
        if component.number_of_arrows() != component.number_of_vertices().saturating_sub(1) {
            messages.push(DYNKIN_ORIENTED_CYCLE_MESSAGE.into());
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

#[cfg(test)]
mod tests {
    use numerica::domains::{Field, finite_field::Z2, rational::Q};
    use qpa_core::{ArrowSpec, PathAlgebra, Quiver, VertexIndex};

    use super::{
        DYNKIN_ORIENTED_CYCLE_MESSAGE, FINITE_TYPE_DYNKIN_MESSAGE, FINITE_TYPE_SEMISIMPLE_MESSAGE,
        INFINITE_TYPE_NON_DYNKIN_MESSAGE, PathAlgebraFiniteTypeReason, path_algebra_finite_type,
    };

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
}
