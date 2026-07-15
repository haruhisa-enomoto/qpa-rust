//! GBNP's ternary todo heap coupled to an occurrence tree.
//!
//! This is a direct typed translation of `lib/kary-heap.gi`. Heap comparison
//! uses only the NP leading monomial and every swap updates the occurrence
//! tree's parallel array mapping.

use numerica::domains::Field;
use thiserror::Error;

use crate::{NpPolynomial, OccurTree, OccurTreeError, TreeDirection, calculate_prefix_generators};

/// GBNP `THeapOT`: a ternary min-heap with a synchronized occurrence tree.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TodoHeap<F: Field> {
    list: Vec<NpPolynomial<F>>,
    tree: OccurTree,
}

impl<F: Field> TodoHeap<F> {
    /// Translate the `THeapOT` constructor.
    ///
    /// # Errors
    ///
    /// Returns an error for zero polynomials, mixed fields, or invalid prefix
    /// generator symbols.
    pub fn new(
        mut todo: Vec<NpPolynomial<F>>,
        prefix_generators: usize,
    ) -> Result<Self, TodoHeapError> {
        validate_polynomials(&todo)?;
        let strictly_sorted = todo
            .windows(2)
            .all(|pair| pair[0].leading_monomial() < pair[1].leading_monomial());
        if !strictly_sorted {
            todo.sort_by(|left, right| left.leading_monomial().cmp(&right.leading_monomial()));
        }
        let tips = todo
            .iter()
            .filter_map(|polynomial| polynomial.leading_monomial().cloned())
            .collect::<Vec<_>>();
        let prefix_generators = prefix_generators.max(calculate_prefix_generators(&tips));
        let tree = OccurTree::new(&tips, prefix_generators, TreeDirection::Forward)?;
        Ok(Self { list: todo, tree })
    }

    /// Number of pending polynomials.
    #[must_use]
    pub fn len(&self) -> usize {
        self.list.len()
    }

    /// Whether the todo heap is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.list.is_empty()
    }

    /// Heap-ordered polynomial storage, matching the source `heap!.list`.
    #[must_use]
    pub fn polynomials(&self) -> &[NpPolynomial<F>] {
        &self.list
    }

    /// Occurrence tree synchronized with `polynomials`.
    #[must_use]
    pub const fn occurrence_tree(&self) -> &OccurTree {
        &self.tree
    }

    /// Return a polynomial at a zero-based heap position.
    #[must_use]
    pub fn get(&self, position: usize) -> Option<&NpPolynomial<F>> {
        self.list.get(position)
    }

    /// Translate `HeapMin`.
    #[must_use]
    pub fn minimum(&self) -> Option<&NpPolynomial<F>> {
        self.list.first()
    }

    /// Translate heap `Add`.
    ///
    /// # Errors
    ///
    /// Returns an error for zero, a different field, or an invalid tree word.
    pub fn add(&mut self, polynomial: NpPolynomial<F>) -> Result<(), TodoHeapError> {
        self.validate_new(&polynomial)?;
        let tip = polynomial
            .leading_monomial()
            .ok_or(TodoHeapError::ZeroPolynomial)?
            .clone();
        self.list.push(polynomial);
        self.tree.append(&tip, TreeDirection::Forward)?;
        self.balance_up(self.list.len() - 1)?;
        Ok(())
    }

    /// Translate heap `Remove` at a zero-based position.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid position or tree invariant failure.
    pub fn remove(&mut self, position: usize) -> Result<NpPolynomial<F>, TodoHeapError> {
        if position >= self.len() {
            return Err(TodoHeapError::HeapIndexOutOfRange {
                index: position,
                heap_length: self.len(),
            });
        }
        let last = self.len() - 1;
        self.swap(position, last)?;
        let tip = self.list[last]
            .leading_monomial()
            .ok_or(TodoHeapError::ZeroPolynomial)?
            .clone();
        self.tree.remove(&tip, last, TreeDirection::Forward)?;
        let removed = self.list.pop().ok_or(TodoHeapError::HeapIndexOutOfRange {
            index: position,
            heap_length: 0,
        })?;
        if position < self.len() {
            let parent = parent_index(position);
            if position > 0 && leading_is_less(&self.list[position], &self.list[parent]) {
                self.balance_up(position)?;
            } else {
                self.balance_down(position)?;
            }
        }
        Ok(removed)
    }

    /// Translate heap `Replace` at a zero-based position.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid position, zero/different-field input,
    /// or tree invariant failure.
    pub fn replace(
        &mut self,
        position: usize,
        polynomial: NpPolynomial<F>,
    ) -> Result<(), TodoHeapError> {
        if position >= self.len() {
            return Err(TodoHeapError::HeapIndexOutOfRange {
                index: position,
                heap_length: self.len(),
            });
        }
        self.validate_new(&polynomial)?;
        let old_tip = self.list[position]
            .leading_monomial()
            .ok_or(TodoHeapError::ZeroPolynomial)?
            .clone();
        let new_tip = polynomial
            .leading_monomial()
            .ok_or(TodoHeapError::ZeroPolynomial)?
            .clone();
        self.tree
            .remove(&old_tip, position, TreeDirection::Forward)?;
        self.tree.add(&new_tip, position, TreeDirection::Forward)?;
        self.list[position] = polynomial;
        let parent = parent_index(position);
        if position > 0 && leading_is_less(&self.list[position], &self.list[parent]) {
            self.balance_up(position)?;
        } else {
            self.balance_down(position)?;
        }
        Ok(())
    }

    /// Consume the heap and return its internal heap-order list.
    #[must_use]
    pub fn into_polynomials(self) -> Vec<NpPolynomial<F>> {
        self.list
    }

    fn validate_new(&self, polynomial: &NpPolynomial<F>) -> Result<(), TodoHeapError> {
        if polynomial.is_zero() {
            return Err(TodoHeapError::ZeroPolynomial);
        }
        if let Some(first) = self.list.first()
            && first.field() != polynomial.field()
        {
            return Err(TodoHeapError::DifferentFields);
        }
        Ok(())
    }

    fn balance_down(&mut self, mut index: usize) -> Result<(), TodoHeapError> {
        loop {
            let mut minimum_child = None;
            for child in children_indices(index) {
                if child < self.len()
                    && minimum_child.map_or_else(
                        || leading_is_less(&self.list[child], &self.list[index]),
                        |current| leading_is_less(&self.list[child], &self.list[current]),
                    )
                {
                    minimum_child = Some(child);
                }
            }
            let Some(child) = minimum_child else {
                return Ok(());
            };
            self.swap(index, child)?;
            index = child;
        }
    }

    fn balance_up(&mut self, mut index: usize) -> Result<(), TodoHeapError> {
        while index > 0 {
            let parent = parent_index(index);
            if !leading_is_less(&self.list[index], &self.list[parent]) {
                break;
            }
            self.swap(index, parent)?;
            index = parent;
        }
        Ok(())
    }

    fn swap(&mut self, first: usize, second: usize) -> Result<(), TodoHeapError> {
        if first != second {
            self.list.swap(first, second);
            self.tree.swap_array_indices(first, second)?;
        }
        Ok(())
    }
}

fn validate_polynomials<F: Field>(polynomials: &[NpPolynomial<F>]) -> Result<(), TodoHeapError> {
    if polynomials.iter().any(NpPolynomial::is_zero) {
        return Err(TodoHeapError::ZeroPolynomial);
    }
    if let Some(first) = polynomials.first()
        && polynomials
            .iter()
            .skip(1)
            .any(|polynomial| polynomial.field() != first.field())
    {
        return Err(TodoHeapError::DifferentFields);
    }
    Ok(())
}

fn leading_is_less<F: Field>(left: &NpPolynomial<F>, right: &NpPolynomial<F>) -> bool {
    left.leading_monomial() < right.leading_monomial()
}

const fn children_indices(index: usize) -> [usize; 3] {
    [3 * index + 1, 3 * index + 2, 3 * index + 3]
}

const fn parent_index(index: usize) -> usize {
    index.saturating_sub(1) / 3
}

/// Failure in GBNP's source-translated todo heap.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum TodoHeapError {
    /// Heap entries must be nonzero.
    #[error("GBNP todo heaps cannot contain the zero polynomial")]
    ZeroPolynomial,
    /// All heap entries must share a coefficient field.
    #[error("GBNP todo heap received different coefficient fields")]
    DifferentFields,
    /// A heap position lies outside the list.
    #[error("heap index {index} is outside a heap of length {heap_length}")]
    HeapIndexOutOfRange {
        /// Invalid zero-based index.
        index: usize,
        /// Current heap length.
        heap_length: usize,
    },
    /// Occurrence-tree synchronization failed.
    #[error(transparent)]
    OccurTree(#[from] OccurTreeError),
}

#[cfg(test)]
mod tests {
    use numerica::domains::{Ring, rational::Q};

    use crate::{NpMonomial, NpPolynomial, TreeDirection, clean_np};

    use super::TodoHeap;

    fn polynomial(tip: &[i32]) -> NpPolynomial<numerica::domains::rational::RationalField> {
        clean_np(Q, vec![NpMonomial::from(tip)], vec![Q.one()]).unwrap()
    }

    #[test]
    fn ternary_heap_operations_keep_the_occurrence_tree_synchronized() {
        let mut heap = TodoHeap::new(
            vec![polynomial(&[2, 1]), polynomial(&[1]), polynomial(&[1, 2])],
            0,
        )
        .unwrap();
        assert_eq!(
            heap.minimum()
                .unwrap()
                .leading_monomial()
                .unwrap()
                .generators(),
            [1]
        );
        heap.add(polynomial(&[])).unwrap();
        assert!(
            heap.minimum()
                .unwrap()
                .leading_monomial()
                .unwrap()
                .generators()
                .is_empty()
        );
        heap.remove(0).unwrap();
        heap.replace(1, polynomial(&[1, 1, 1])).unwrap();
        for (index, polynomial) in heap.polynomials().iter().enumerate() {
            assert!(
                heap.occurrence_tree()
                    .all_occurrences(
                        polynomial.leading_monomial().unwrap(),
                        TreeDirection::Forward
                    )
                    .iter()
                    .any(|(matched, _)| *matched == index)
            );
        }
    }
}
