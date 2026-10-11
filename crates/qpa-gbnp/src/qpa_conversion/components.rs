//! Structural preconditioning for checked path-algebra completion.
//!
//! A path with endpoints in one strongly connected component never leaves
//! that component. Relations internal to components can therefore be completed
//! first, independently of the relations crossing the condensation DAG.
//! Replacing each crossing relation by its normal form modulo that internal
//! ideal preserves the full ideal. The usual `SGrobner` final reduction still
//! supplies the ordered, monic reduced basis; no presentation-specific data or
//! new monomial order is used.

use numerica::domains::Field;
use qpa_core::{PathAlgebra, PathAlgebraElement, PathEndpoint, Quiver};

use super::{QpaGbnpError, path_to_cohen};
use crate::{NpPolynomial, SgrobnerError, sgrobner, strong_normal_form_np};

/// None means there is no internal/crossing split to exploit.
pub(super) fn precondition<F: Field>(
    uniform: &[PathAlgebraElement<F>],
    algebra: &PathAlgebra<F>,
) -> Result<Option<Vec<NpPolynomial<F>>>, QpaGbnpError> {
    if uniform.is_empty() || algebra.quiver().number_of_vertices() < 2 {
        return Ok(None);
    }
    let components = component_ids(algebra.quiver());
    let mut internal = Vec::new();
    let mut crossing = Vec::new();
    for relation in uniform {
        let Some(term) = relation.terms().first() else {
            continue;
        };
        let (PathEndpoint::Vertex(source), PathEndpoint::Vertex(target)) =
            (term.path().source(), term.path().target())
        else {
            // Preserve the existing converter's error boundary.
            return Ok(None);
        };
        if components[source.index().get() - 1] == components[target.index().get() - 1] {
            internal.push(relation.clone());
        } else {
            crossing.push(relation.clone());
        }
    }
    if internal.is_empty() || crossing.is_empty() {
        return Ok(None);
    }
    let input = path_to_cohen(&internal)?
        .iter()
        .map(super::CohenPolynomial::clean)
        .collect::<Result<Vec<_>, _>>()?;
    let basis = sgrobner(&input)?;
    let mut remainders = Vec::new();
    // Stay in NP coordinates: use SGrobner's order even if the parent carries
    // a different path ordering. No extra completion or order conversion.
    for relation in path_to_cohen(&crossing)? {
        let remainder =
            strong_normal_form_np(&relation.clean()?, &basis).map_err(SgrobnerError::from)?;
        if !remainder.is_zero() {
            remainders.push(remainder);
        }
    }
    let mut completed = basis;
    completed.extend(remainders);
    Ok(Some(completed))
}

/// Iterative Kosaraju traversal: O(vertices + arrows), without recursion or
/// an all-pairs reachability matrix. Component IDs have no algebraic meaning.
fn component_ids(quiver: &Quiver) -> Vec<usize> {
    let count = quiver.number_of_vertices();
    let mut edges = vec![Vec::new(); count];
    let mut reversed = vec![Vec::new(); count];
    for arrow in quiver.arrows() {
        let source = arrow.source().index().get() - 1;
        let target = arrow.target().index().get() - 1;
        edges[source].push(target);
        reversed[target].push(source);
    }
    let mut visited = vec![false; count];
    let mut order = Vec::with_capacity(count);
    let mut stack = Vec::new();
    for start in 0..count {
        if visited[start] {
            continue;
        }
        visited[start] = true;
        stack.push((start, 0));
        while let Some((vertex, next)) = stack.last_mut() {
            let current = *vertex;
            if *next < edges[current].len() {
                let target = edges[current][*next];
                *next += 1;
                if !visited[target] {
                    visited[target] = true;
                    stack.push((target, 0));
                }
            } else {
                order.push(current);
                stack.pop();
            }
        }
    }
    let mut components = vec![usize::MAX; count];
    let mut pending = Vec::new();
    for start in order.into_iter().rev() {
        if components[start] != usize::MAX {
            continue;
        }
        components[start] = start;
        pending.push(start);
        while let Some(vertex) = pending.pop() {
            for &source in &reversed[vertex] {
                if components[source] == usize::MAX {
                    components[source] = start;
                    pending.push(source);
                }
            }
        }
    }
    components
}

#[cfg(test)]
mod tests;
