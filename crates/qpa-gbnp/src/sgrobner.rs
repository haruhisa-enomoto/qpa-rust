//! GBNP's strong noncommutative Gröbner-basis driver.
//!
//! This module follows `GBNP.AllObs`, `GBNP.ObsTall`, the tree strong-normal
//! forms, `GBNP.CentralT`, `GBNP.SGrobnerLoops`, and public `SGrobner` in the
//! pinned GBNP sources.

use std::collections::BTreeSet;

use numerica::domains::Field;
use thiserror::Error;

use crate::{
    GbnpOptions, NpMonomial, NpPolynomial, Obstruction, ObstructionError, OccurTree,
    OccurTreeError, PreprocessingError, ReductionError, TodoHeap, TodoHeapError, TreeDirection,
    calculate_prefix_generators, left_obstructions_tree, mk_monic_np, occur,
    preprocessing::{
        reduce_cancellation, reduce_polynomials_in_place, reduce_polynomials_with_options,
        strong_normal_form_except,
    },
    reduction::reduce_indexed_term,
    right_obstructions_tree, spolynomial,
};

/// Translate `GBNP.StrongNormalForm2Tall` using two occurrence trees.
///
/// When `strong` is false, only the leading term is reduced, matching the
/// source record switch. `SGrobner` always passes `true`.
///
/// # Errors
///
/// Returns an error for mismatched tree/list indices or reduction failures.
pub fn strong_normal_form_trees<F: Field>(
    polynomial: &NpPolynomial<F>,
    primary: &[NpPolynomial<F>],
    secondary: &[NpPolynomial<F>],
    primary_tree: &OccurTree,
    secondary_tree: &OccurTree,
    strong: bool,
    options: &GbnpOptions,
) -> Result<NpPolynomial<F>, SgrobnerError> {
    strong_normal_form_trees_skipping_secondary(
        polynomial,
        primary,
        secondary,
        primary_tree,
        secondary_tree,
        strong,
        None,
        options,
    )
}

#[allow(clippy::too_many_arguments)]
fn strong_normal_form_trees_skipping_secondary<F: Field>(
    polynomial: &NpPolynomial<F>,
    primary: &[NpPolynomial<F>],
    secondary: &[NpPolynomial<F>],
    primary_tree: &OccurTree,
    secondary_tree: &OccurTree,
    strong: bool,
    skipped_secondary: Option<usize>,
    options: &GbnpOptions,
) -> Result<NpPolynomial<F>, SgrobnerError> {
    if polynomial.is_zero() {
        return Ok(polynomial.clone());
    }
    let mut result = polynomial.clone();
    let mut term_index = 0;
    loop {
        let Some(term) = result.terms().get(term_index) else {
            break;
        };
        let mut monomial = term.monomial().clone();
        let mut primary_match = primary_tree.first_occurrence(&monomial, TreeDirection::Forward);
        let mut secondary_occurrence =
            secondary_match(secondary_tree, &monomial, skipped_secondary);
        while primary_match.is_some() || secondary_occurrence.is_some() {
            let (reducer, position) = if let Some((index, position)) = primary_match {
                (
                    primary
                        .get(index)
                        .ok_or(SgrobnerError::TreeIndexOutOfRange {
                            index,
                            list_length: primary.len(),
                        })?,
                    position,
                )
            } else {
                let (index, position) =
                    secondary_occurrence.ok_or(SgrobnerError::MissingTreeReductionMatch)?;
                (
                    secondary
                        .get(index)
                        .ok_or(SgrobnerError::TreeIndexOutOfRange {
                            index,
                            list_length: secondary.len(),
                        })?,
                    position,
                )
            };
            result = reduce_indexed_term(&result, term_index, reducer, position + 1)?;
            if result.is_zero() || term_index >= result.terms().len() {
                return Ok(reduce_cancellation(&result, options)?);
            }
            monomial = result.terms()[term_index].monomial().clone();
            primary_match = primary_tree.first_occurrence(&monomial, TreeDirection::Forward);
            secondary_occurrence = secondary_match(secondary_tree, &monomial, skipped_secondary);
        }
        term_index += 1;
        if !strong || term_index >= result.terms().len() {
            break;
        }
    }
    Ok(reduce_cancellation(&result, options)?)
}

fn secondary_match(
    tree: &OccurTree,
    monomial: &NpMonomial,
    skipped: Option<usize>,
) -> Option<(usize, usize)> {
    if let Some(skipped) = skipped {
        tree.all_occurrences(monomial, TreeDirection::Forward)
            .into_iter()
            .find(|(index, _)| *index != skipped)
    } else {
        tree.first_occurrence(monomial, TreeDirection::Forward)
    }
}

/// Translate `GBNP.AllObs` for a preprocessed monic basis.
///
/// # Errors
///
/// Returns an error for invalid basis/tree invariants or NP arithmetic.
pub fn all_obstructions<F: Field>(
    basis: &[NpPolynomial<F>],
    options: &GbnpOptions,
) -> Result<Vec<NpPolynomial<F>>, SgrobnerError> {
    if is_trivial_basis(basis) {
        return Ok(Vec::new());
    }
    let leading = leading_monomials(basis)?;
    let prefix_generators = calculate_prefix_generators(&leading).max(options.prefix_generators());
    let basis_left = OccurTree::new(&leading, prefix_generators, TreeDirection::Forward)?;
    let mut partial_left = OccurTree::new(&[], prefix_generators, TreeDirection::Forward)?;
    let mut partial_right = OccurTree::new(&[], prefix_generators, TreeDirection::Reverse)?;
    let mut todo = Vec::new();
    let mut todo_tree = OccurTree::new(&[], prefix_generators, TreeDirection::Forward)?;

    for (basis_index, leading_monomial) in leading.iter().enumerate().take(basis.len()) {
        partial_left.add(leading_monomial, basis_index, TreeDirection::Forward)?;
        partial_right.add(leading_monomial, basis_index, TreeDirection::Reverse)?;
        add_obstructions_to_list(
            basis_index,
            basis,
            &mut todo,
            &basis_left,
            &partial_left,
            &partial_right,
            &mut todo_tree,
            options,
        )?;
    }

    let empty_tree = OccurTree::new(&[], prefix_generators, TreeDirection::Forward)?;
    for polynomial in &mut todo {
        *polynomial = strong_normal_form_trees(
            polynomial,
            basis,
            &[],
            &basis_left,
            &empty_tree,
            true,
            options,
        )?;
    }
    sort_polynomials(&mut todo);
    Ok(todo)
}

/// Translate public GBNP `SGrobner` with default options.
///
/// Like the source routine, this may not terminate for arbitrary input.
///
/// # Errors
///
/// Returns an error if any source-translated Gröbner phase fails.
pub fn sgrobner<F: Field>(
    input: &[NpPolynomial<F>],
) -> Result<Vec<NpPolynomial<F>>, SgrobnerError> {
    Ok(sgrobner_with_options(input, &GbnpOptions::default())?.into_basis())
}

/// Translate public GBNP `SGrobner` with typed options and a result record.
///
/// Set [`GbnpOptions::set_max_iterations`] to obtain the source's interruptible
/// Gröbner-pair result.
///
/// # Errors
///
/// Returns an error if any source-translated Gröbner phase fails.
pub fn sgrobner_with_options<F: Field>(
    input: &[NpPolynomial<F>],
    options: &GbnpOptions,
) -> Result<SgrobnerResult<F>, SgrobnerError> {
    let basis = reduce_polynomials_with_options(input, options)?;
    let mut loop_options = options.clone();
    let prefix_generators =
        calculate_prefix_generators(&leading_monomials(&basis)?).max(options.prefix_generators());
    loop_options.set_prefix_generators(prefix_generators);
    let todo = all_obstructions(&basis, &loop_options)?;
    let mut result = sgrobner_loops(basis, todo, &loop_options)?;

    let tree = reduce_polynomials_in_place(&mut result.basis, &loop_options)?;
    for index in 0..result.basis.len() {
        result.basis[index] = strong_normal_form_except(&result.basis, index, &tree)?;
    }
    Ok(result)
}

/// Translate `GBNP.SGrobnerLoops` for an existing Gröbner pair.
///
/// # Errors
///
/// Returns an error for invalid pair/tree invariants or exact NP failures.
#[allow(clippy::too_many_lines)]
pub fn sgrobner_loops<F: Field>(
    mut basis: Vec<NpPolynomial<F>>,
    todo: Vec<NpPolynomial<F>>,
    options: &GbnpOptions,
) -> Result<SgrobnerResult<F>, SgrobnerError> {
    if is_trivial_basis(&basis) {
        return Ok(SgrobnerResult {
            basis,
            todo: Vec::new(),
            completed: true,
            iterations: 0,
        });
    }

    let leading = leading_monomials(&basis)?;
    let prefix_generators = calculate_prefix_generators(&leading).max(options.prefix_generators());
    let mut basis_left = OccurTree::new(&leading, prefix_generators, TreeDirection::Forward)?;
    let mut basis_right = OccurTree::new(&leading, prefix_generators, TreeDirection::Reverse)?;
    let mut symbol_sets = leading.iter().map(symbol_set).collect::<Vec<_>>();
    let mut heap = TodoHeap::new(todo, prefix_generators)?;
    let mut iterations = 0;

    while !heap.is_empty() {
        iterations += 1;

        let spoly = heap
            .minimum()
            .cloned()
            .ok_or(SgrobnerError::EmptyHeapMinimum)?;
        let spoly_tip = spoly
            .leading_monomial()
            .cloned()
            .ok_or(SgrobnerError::ZeroBasisPolynomial)?;
        let spoly_symbols = symbol_set(&spoly_tip);
        basis.push(spoly);
        symbol_sets.push(spoly_symbols.clone());
        heap.remove(0)?;
        basis_left.append(&spoly_tip, TreeDirection::Forward)?;
        basis_right.append(&spoly_tip, TreeDirection::Reverse)?;

        let new_index = basis.len() - 1;
        add_obstructions_to_heap(
            new_index,
            &basis,
            &mut heap,
            &basis_left,
            &basis_right,
            options,
        )?;

        let mut index = 0;
        let mut loop_bound = basis.len();
        while index + 1 < loop_bound {
            let old_tip = basis[index]
                .leading_monomial()
                .cloned()
                .ok_or(SgrobnerError::ZeroBasisPolynomial)?;
            let reducible =
                spoly_symbols.is_subset(&symbol_sets[index]) && occur(&spoly_tip, &old_tip) > 0;
            if !reducible {
                index += 1;
                continue;
            }

            let old = basis.remove(index);
            symbol_sets.remove(index);
            basis_left.remove(&old_tip, index, TreeDirection::Forward)?;
            basis_right.remove(&old_tip, index, TreeDirection::Reverse)?;
            let reduced = strong_normal_form_trees(
                &old,
                &basis,
                heap.polynomials(),
                &basis_left,
                heap.occurrence_tree(),
                true,
                options,
            )?;
            if !reduced.is_zero() {
                let monic = mk_monic_np(&reduced);
                let tip = monic
                    .leading_monomial()
                    .cloned()
                    .ok_or(SgrobnerError::ZeroBasisPolynomial)?;
                basis.push(monic);
                symbol_sets.push(symbol_set(&tip));
                basis_left.append(&tip, TreeDirection::Forward)?;
                basis_right.append(&tip, TreeDirection::Reverse)?;
                let replacement_index = basis.len() - 1;
                add_central_obstructions_to_heap(
                    replacement_index,
                    &basis,
                    &mut heap,
                    &basis_left,
                    options,
                )?;
                add_obstructions_to_heap(
                    replacement_index,
                    &basis,
                    &mut heap,
                    &basis_left,
                    &basis_right,
                    options,
                )?;
            }
            loop_bound -= 1;
        }

        let mut todo_index = 0;
        while todo_index < heap.len() {
            let polynomial =
                heap.get(todo_index)
                    .cloned()
                    .ok_or(SgrobnerError::TreeIndexOutOfRange {
                        index: todo_index,
                        list_length: heap.len(),
                    })?;
            let tip = polynomial
                .leading_monomial()
                .ok_or(SgrobnerError::ZeroBasisPolynomial)?;
            let reducible = basis_left
                .first_occurrence(tip, TreeDirection::Forward)
                .is_some()
                || secondary_match(heap.occurrence_tree(), tip, Some(todo_index)).is_some();
            if !reducible {
                todo_index += 1;
                continue;
            }
            let reduced = strong_normal_form_trees_skipping_secondary(
                &polynomial,
                &basis,
                heap.polynomials(),
                &basis_left,
                heap.occurrence_tree(),
                true,
                Some(todo_index),
                options,
            )?;
            if reduced.is_zero() {
                heap.remove(todo_index)?;
            } else {
                heap.replace(todo_index, mk_monic_np(&reduced))?;
                todo_index += 1;
            }
        }

        if options
            .max_iterations()
            .is_some_and(|limit| iterations >= limit)
        {
            return Ok(SgrobnerResult {
                basis,
                todo: heap.into_polynomials(),
                completed: false,
                iterations,
            });
        }
    }

    Ok(SgrobnerResult {
        basis,
        todo: Vec::new(),
        completed: true,
        iterations,
    })
}

fn add_obstructions_to_heap<F: Field>(
    basis_index: usize,
    basis: &[NpPolynomial<F>],
    heap: &mut TodoHeap<F>,
    basis_left: &OccurTree,
    basis_right: &OccurTree,
    options: &GbnpOptions,
) -> Result<(), SgrobnerError> {
    let leading = leading_monomials(basis)?;
    let obstructions = collect_tree_obstructions(basis_index, &leading, basis_left, basis_right)?;
    for obstruction in obstructions {
        let spoly = spolynomial(&obstruction, basis)?;
        if spoly.is_zero() {
            continue;
        }
        let reduced = strong_normal_form_trees(
            &spoly,
            basis,
            heap.polynomials(),
            basis_left,
            heap.occurrence_tree(),
            true,
            options,
        )?;
        if !reduced.is_zero() {
            heap.add(mk_monic_np(&reduced))?;
        }
    }
    Ok(())
}

fn add_central_obstructions_to_heap<F: Field>(
    basis_index: usize,
    basis: &[NpPolynomial<F>],
    heap: &mut TodoHeap<F>,
    basis_left: &OccurTree,
    options: &GbnpOptions,
) -> Result<(), SgrobnerError> {
    let leading = leading_monomials(basis)?;
    let word = &leading[basis_index];
    for (other_index, other) in leading[..basis_index].iter().enumerate() {
        let occurrence = occur(word, other);
        if occurrence <= 1 || occurrence + word.degree() > other.degree() {
            continue;
        }
        let start = occurrence - 1;
        let obstruction = Obstruction::new(
            NpMonomial::from(&other.generators()[..start]),
            basis_index,
            NpMonomial::from(&other.generators()[start + word.degree()..]),
            NpMonomial::from(&[][..]),
            other_index,
            NpMonomial::from(&[][..]),
        );
        let spoly = spolynomial(&obstruction, basis)?;
        if spoly.is_zero() {
            continue;
        }
        let reduced = strong_normal_form_trees(
            &spoly,
            basis,
            heap.polynomials(),
            basis_left,
            heap.occurrence_tree(),
            true,
            options,
        )?;
        if !reduced.is_zero() {
            heap.add(mk_monic_np(&reduced))?;
        }
    }
    Ok(())
}

fn symbol_set(monomial: &NpMonomial) -> BTreeSet<i32> {
    monomial.generators().iter().copied().collect()
}

#[allow(clippy::too_many_arguments)]
fn add_obstructions_to_list<F: Field>(
    basis_index: usize,
    basis: &[NpPolynomial<F>],
    todo: &mut Vec<NpPolynomial<F>>,
    basis_left: &OccurTree,
    obstruction_left: &OccurTree,
    obstruction_right: &OccurTree,
    todo_tree: &mut OccurTree,
    options: &GbnpOptions,
) -> Result<(), SgrobnerError> {
    let leading = leading_monomials(basis)?;
    let obstructions =
        collect_tree_obstructions(basis_index, &leading, obstruction_left, obstruction_right)?;
    for obstruction in obstructions {
        let spoly = spolynomial(&obstruction, basis)?;
        if spoly.is_zero() {
            continue;
        }
        let reduced =
            strong_normal_form_trees(&spoly, basis, todo, basis_left, todo_tree, true, options)?;
        if !reduced.is_zero() {
            let monic = mk_monic_np(&reduced);
            todo_tree.append(
                monic
                    .leading_monomial()
                    .ok_or(SgrobnerError::ZeroBasisPolynomial)?,
                TreeDirection::Forward,
            )?;
            todo.push(monic);
        }
    }
    sort_polynomials_and_rebuild_tree(todo, todo_tree)?;
    Ok(())
}

fn collect_tree_obstructions(
    basis_index: usize,
    leading: &[NpMonomial],
    left_tree: &OccurTree,
    right_tree: &OccurTree,
) -> Result<Vec<Obstruction>, SgrobnerError> {
    let (mut left, left_self) =
        left_obstructions_tree(basis_index, leading, left_tree)?.into_parts();
    let (mut right, right_self) =
        right_obstructions_tree(basis_index, leading, right_tree)?.into_parts();
    if let Some(index) = left_self {
        left.remove(index);
    } else if let Some(index) = right_self {
        right.remove(index);
    }
    left.extend(right);
    Ok(left)
}

fn sort_polynomials<F: Field>(polynomials: &mut [NpPolynomial<F>]) {
    polynomials.sort_by(|left, right| left.leading_monomial().cmp(&right.leading_monomial()));
}

fn sort_polynomials_and_rebuild_tree<F: Field>(
    polynomials: &mut [NpPolynomial<F>],
    tree: &mut OccurTree,
) -> Result<(), SgrobnerError> {
    sort_polynomials(polynomials);
    let leading = leading_monomials(polynomials)?;
    *tree = OccurTree::new(&leading, tree.prefix_generators(), TreeDirection::Forward)?;
    Ok(())
}

fn leading_monomials<F: Field>(
    polynomials: &[NpPolynomial<F>],
) -> Result<Vec<NpMonomial>, SgrobnerError> {
    polynomials
        .iter()
        .map(|polynomial| {
            polynomial
                .leading_monomial()
                .cloned()
                .ok_or(SgrobnerError::ZeroBasisPolynomial)
        })
        .collect()
}

fn is_trivial_basis<F: Field>(basis: &[NpPolynomial<F>]) -> bool {
    basis.is_empty()
        || (basis.len() == 1
            && basis[0].terms().len() == 1
            && basis[0].terms()[0].monomial().generators().is_empty()
            && basis[0].field().is_one(basis[0].terms()[0].coefficient()))
}

/// Result record returned by an iteration-limited `SGrobner` computation.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SgrobnerResult<F: Field> {
    basis: Vec<NpPolynomial<F>>,
    todo: Vec<NpPolynomial<F>>,
    completed: bool,
    iterations: usize,
}

impl<F: Field> SgrobnerResult<F> {
    /// Current Gröbner basis list.
    #[must_use]
    pub fn basis(&self) -> &[NpPolynomial<F>] {
        &self.basis
    }

    /// Current todo list when an iteration limit interrupted the computation.
    #[must_use]
    pub fn todo(&self) -> &[NpPolynomial<F>] {
        &self.todo
    }

    /// Whether the todo heap was exhausted.
    #[must_use]
    pub const fn completed(&self) -> bool {
        self.completed
    }

    /// Number of completed main-loop iterations.
    #[must_use]
    pub const fn iterations(&self) -> usize {
        self.iterations
    }

    /// Consume the result and return its basis.
    #[must_use]
    pub fn into_basis(self) -> Vec<NpPolynomial<F>> {
        self.basis
    }
}

/// Failure in GBNP's strong Gröbner machinery.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum SgrobnerError {
    /// A tree returned an index outside its synchronized polynomial list.
    #[error("occurrence-tree index {index} is outside a list of length {list_length}")]
    TreeIndexOutOfRange {
        /// Invalid zero-based index.
        index: usize,
        /// Synchronized list length.
        list_length: usize,
    },
    /// A nonzero tree match unexpectedly disappeared.
    #[error("a secondary occurrence-tree match was unexpectedly missing")]
    MissingTreeReductionMatch,
    /// A source routine requiring nonzero polynomials received zero.
    #[error("GBNP Gröbner machinery received a zero basis polynomial")]
    ZeroBasisPolynomial,
    /// The loop observed an empty heap after checking it was nonempty.
    #[error("GBNP todo heap unexpectedly had no minimum")]
    EmptyHeapMinimum,
    /// NP reduction failed.
    #[error(transparent)]
    Reduction(#[from] ReductionError),
    /// Basis preprocessing failed.
    #[error(transparent)]
    Preprocessing(#[from] PreprocessingError),
    /// Obstruction construction or S-polynomial arithmetic failed.
    #[error(transparent)]
    Obstruction(#[from] ObstructionError),
    /// Occurrence-tree construction or mutation failed.
    #[error(transparent)]
    OccurTree(#[from] OccurTreeError),
    /// Todo-heap construction or mutation failed.
    #[error(transparent)]
    TodoHeap(#[from] TodoHeapError),
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        Field, Ring,
        finite_field::{FiniteFieldCore, Z2, Zp},
        integer::Integer,
        rational::{Q, RationalField},
    };
    use serde_json::{Value, json};

    use crate::{GbnpOptions, NpMonomial, NpPolynomial, clean_np, reduce_polynomials};

    use super::{all_obstructions, sgrobner, sgrobner_with_options};

    fn polynomial(
        monomials: &[&[i32]],
        coefficients: &[i64],
    ) -> NpPolynomial<numerica::domains::rational::RationalField> {
        clean_np(
            Q,
            monomials
                .iter()
                .map(|monomial| NpMonomial::from(*monomial))
                .collect(),
            coefficients
                .iter()
                .map(|&coefficient| Q.nth(coefficient.into()))
                .collect(),
        )
        .unwrap()
    }

    fn small_integer(value: &Integer) -> i64 {
        value
            .to_string()
            .parse()
            .expect("the SGrobner fixture uses only small integers")
    }

    fn polynomial_json<F: Field>(
        polynomial: &NpPolynomial<F>,
        coefficient_json: impl Fn(&F::Element) -> Value,
    ) -> Value {
        Value::Array(
            polynomial
                .terms()
                .iter()
                .map(|term| {
                    json!({
                        "monomial": term.monomial().generators(),
                        "coefficient": coefficient_json(term.coefficient())
                    })
                })
                .collect(),
        )
    }

    fn rational_list_json(polynomials: &[NpPolynomial<RationalField>]) -> Value {
        Value::Array(
            polynomials
                .iter()
                .map(|polynomial| {
                    polynomial_json(polynomial, |coefficient| {
                        json!([
                            small_integer(coefficient.numerator_ref()),
                            small_integer(coefficient.denominator_ref())
                        ])
                    })
                })
                .collect(),
        )
    }

    fn zp_list_json(polynomials: &[NpPolynomial<Zp>], field: &Zp) -> Value {
        Value::Array(
            polynomials
                .iter()
                .map(|polynomial| {
                    polynomial_json(polynomial, |coefficient| {
                        small_integer(&field.to_integer(coefficient)).into()
                    })
                })
                .collect(),
        )
    }

    fn z2_list_json(polynomials: &[NpPolynomial<Z2>]) -> Value {
        Value::Array(
            polynomials
                .iter()
                .map(|polynomial| {
                    polynomial_json(polynomial, |coefficient| {
                        small_integer(&Z2.to_integer(coefficient)).into()
                    })
                })
                .collect(),
        )
    }

    #[test]
    fn all_obs_builds_reduced_spolynomial_todo() {
        let basis = [
            polynomial(&[&[1, 2, 1], &[1]], &[1, -1]),
            polynomial(&[&[2, 1, 2], &[2]], &[1, 2]),
        ];
        let todo = all_obstructions(&basis, &GbnpOptions::default()).unwrap();
        assert!(todo.iter().all(|polynomial| !polynomial.is_zero()));
    }

    #[test]
    fn sgrobner_matches_the_documented_two_relation_example() {
        let input = [
            polynomial(&[&[1, 1, 2], &[]], &[1, -1]),
            polynomial(&[&[1, 2, 2], &[]], &[1, -1]),
        ];
        let basis = sgrobner(&input).unwrap();
        assert_eq!(basis.len(), 2);
        assert_eq!(basis[0].leading_monomial().unwrap().generators(), [2]);
        assert_eq!(basis[1].leading_monomial().unwrap().generators(), [1, 1, 1]);

        let mut options = GbnpOptions::default();
        options.set_max_iterations(Some(1));
        let interrupted = sgrobner_with_options(&input, &options).unwrap();
        assert!(!interrupted.completed());
        assert_eq!(interrupted.iterations(), 1);
        assert_eq!(interrupted.basis().len(), 1);
        assert_eq!(interrupted.todo().len(), 1);
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn pinned_gbnp_sgrobner_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/gbnp_sgrobner.schema-v1.json"
        ))
        .unwrap();
        assert_eq!(
            expected["oracle"],
            json!({
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            })
        );

        let rational_input = [
            polynomial(&[&[1, 1, 2], &[]], &[1, -1]),
            polynomial(&[&[1, 2, 2], &[]], &[1, -1]),
        ];
        let rational_phase1 = reduce_polynomials(&rational_input).unwrap();
        let rational_initial_todo =
            all_obstructions(&rational_phase1, &GbnpOptions::default()).unwrap();
        let rational_full = sgrobner(&rational_input).unwrap();
        let mut one_options = GbnpOptions::default();
        one_options.set_max_iterations(Some(1));
        let rational_one = sgrobner_with_options(&rational_input, &one_options).unwrap();
        let rational_empty = sgrobner(&Vec::<NpPolynomial<RationalField>>::new()).unwrap();
        let rational_identity = sgrobner(&[NpPolynomial::one(Q)]).unwrap();

        let field5 = Zp::new(5);
        let f5_input = [
            clean_np(
                field5.clone(),
                vec![NpMonomial::from(&[1, 1, 2][..]), NpMonomial::from(&[][..])],
                vec![field5.one(), field5.to_element(4)],
            )
            .unwrap(),
            clean_np(
                field5.clone(),
                vec![NpMonomial::from(&[1, 2, 2][..]), NpMonomial::from(&[][..])],
                vec![field5.one(), field5.to_element(4)],
            )
            .unwrap(),
        ];
        let f5_full = sgrobner(&f5_input).unwrap();

        let one2 = Z2.one();
        let f2_input = [
            clean_np(
                Z2,
                vec![NpMonomial::from(&[1, 1, 2][..]), NpMonomial::from(&[][..])],
                vec![one2, one2],
            )
            .unwrap(),
            clean_np(
                Z2,
                vec![NpMonomial::from(&[1, 2, 2][..]), NpMonomial::from(&[][..])],
                vec![one2, one2],
            )
            .unwrap(),
        ];
        let f2_full = sgrobner(&f2_input).unwrap();

        assert_eq!(
            json!({
                "schema_version": 1,
                "kind": "gbnp_sgrobner",
                "oracle": expected["oracle"].clone(),
                "rational": {
                    "phase1": rational_list_json(&rational_phase1),
                    "initial_todo": rational_list_json(&rational_initial_todo),
                    "full": rational_list_json(&rational_full),
                    "one_iteration": {
                        "basis": rational_list_json(rational_one.basis()),
                        "todo": rational_list_json(rational_one.todo()),
                        "iterations": rational_one.iterations(),
                        "completed": rational_one.completed()
                    },
                    "empty": rational_list_json(&rational_empty),
                    "identity": rational_list_json(&rational_identity)
                },
                "f5": {"full": zp_list_json(&f5_full, &field5)},
                "f2": {"full": z2_list_json(&f2_full)}
            }),
            expected
        );
    }
}
