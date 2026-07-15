//! GBNP's mutable left/right occurrence tree.
//!
//! This translates the trie and parallel-index bookkeeping in
//! `lib/occurtree.gi`. GAP uses sparse lists offset by the number of negative
//! prefix generators; Rust uses a typed ordered map while preserving lookup
//! and array-index selection semantics.

use std::collections::BTreeMap;

use thiserror::Error;

use crate::NpMonomial;

/// Direction in which a GBNP occurrence tree consumes a word.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TreeDirection {
    /// Consume the monomial from left to right.
    Forward,
    /// Consume the monomial from right to left.
    Reverse,
}

/// Mutable occurrence trie with GBNP's stable tree-number/array-index maps.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct OccurTree {
    root: Node,
    tree_to_array: Vec<Option<usize>>,
    array_to_tree: Vec<usize>,
    free_tree_numbers: Vec<usize>,
    prefix_generators: usize,
}

impl OccurTree {
    /// Translate `GBNP.CreateOccurTreePTSLRimpl`.
    ///
    /// # Errors
    ///
    /// Returns an error if a negative symbol lies outside the declared prefix
    /// generator range.
    pub fn new(
        monomials: &[NpMonomial],
        prefix_generators: usize,
        direction: TreeDirection,
    ) -> Result<Self, OccurTreeError> {
        let mut tree = Self {
            root: Node::default(),
            tree_to_array: Vec::new(),
            array_to_tree: Vec::new(),
            free_tree_numbers: Vec::new(),
            prefix_generators,
        };
        for (array_index, monomial) in monomials.iter().enumerate() {
            tree.add(monomial, array_index, direction)?;
        }
        Ok(tree)
    }

    /// Number of array entries represented by the tree.
    #[must_use]
    pub fn len(&self) -> usize {
        self.array_to_tree.len()
    }

    /// Whether the represented array is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.array_to_tree.is_empty()
    }

    /// Number of negative prefix/module generators admitted by the trie.
    #[must_use]
    pub const fn prefix_generators(&self) -> usize {
        self.prefix_generators
    }

    /// Translate `GBNP.AddMonToTreePTSLR` with a zero-based insertion index.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid insertion index or prefix symbol.
    pub fn add(
        &mut self,
        monomial: &NpMonomial,
        array_index: usize,
        direction: TreeDirection,
    ) -> Result<(), OccurTreeError> {
        if array_index > self.array_to_tree.len() {
            return Err(OccurTreeError::ArrayIndexOutOfRange {
                index: array_index,
                array_length: self.array_to_tree.len(),
            });
        }
        self.validate_monomial(monomial)?;
        let tree_number = self.allocate_tree_number();
        self.array_to_tree.insert(array_index, tree_number);
        self.rebuild_tree_to_array();

        let mut node = &mut self.root;
        for symbol in directed_symbols(monomial, direction) {
            node = node.children.entry(symbol).or_default();
        }
        if node.terminal.is_none() {
            node.terminal = Some(tree_number);
        }
        Ok(())
    }

    /// Append a monomial, corresponding to source index `-1`.
    ///
    /// # Errors
    ///
    /// Returns an error if a prefix symbol is outside the declared range.
    pub fn append(
        &mut self,
        monomial: &NpMonomial,
        direction: TreeDirection,
    ) -> Result<(), OccurTreeError> {
        self.add(monomial, self.len(), direction)
    }

    /// Translate `GBNP.RemoveMonFromTreePTSLR` with a zero-based array index.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid index, prefix symbol, or missing path.
    pub fn remove(
        &mut self,
        monomial: &NpMonomial,
        array_index: usize,
        direction: TreeDirection,
    ) -> Result<(), OccurTreeError> {
        let Some(&tree_number) = self.array_to_tree.get(array_index) else {
            return Err(OccurTreeError::ArrayIndexOutOfRange {
                index: array_index,
                array_length: self.array_to_tree.len(),
            });
        };
        self.validate_monomial(monomial)?;
        let symbols = directed_symbols(monomial, direction).collect::<Vec<_>>();
        if !remove_terminal(&mut self.root, &symbols, tree_number) {
            return Err(OccurTreeError::MissingMonomialPath);
        }
        self.array_to_tree.remove(array_index);
        self.tree_to_array[tree_number] = None;
        self.free_tree_numbers.push(tree_number);
        self.rebuild_tree_to_array();
        Ok(())
    }

    /// Translate `GBNP.LookUpOccurTreePTSLRPos`.
    ///
    /// Positions and returned array indices are zero-based.
    #[must_use]
    pub fn lookup_at(
        &self,
        monomial: &NpMonomial,
        direction: TreeDirection,
        start_position: usize,
    ) -> Option<usize> {
        let symbols = directed_symbols(monomial, direction).collect::<Vec<_>>();
        let mut node = &self.root;
        for &symbol in symbols.get(start_position..)? {
            if let Some(array_index) = self.terminal_array_index(node) {
                return Some(array_index);
            }
            node = node.children.get(&symbol)?;
        }
        self.terminal_array_index(node)
    }

    /// Translate `GBNP.LookUpOccurTreeAllPTSLRPos`.
    #[must_use]
    pub fn lookup_all_at(
        &self,
        monomial: &NpMonomial,
        direction: TreeDirection,
        start_position: usize,
    ) -> Vec<usize> {
        let symbols = directed_symbols(monomial, direction).collect::<Vec<_>>();
        let Some(suffix) = symbols.get(start_position..) else {
            return Vec::new();
        };
        let mut answer = Vec::new();
        let mut node = &self.root;
        for &symbol in suffix {
            if let Some(array_index) = self.terminal_array_index(node) {
                answer.push(array_index);
            }
            let Some(next) = node.children.get(&symbol) else {
                return answer;
            };
            node = next;
        }
        if let Some(array_index) = self.terminal_array_index(node) {
            answer.push(array_index);
        }
        answer
    }

    /// Translate `GBNP.OccurInLstPTSLR`.
    #[must_use]
    pub fn first_occurrence(
        &self,
        monomial: &NpMonomial,
        direction: TreeDirection,
    ) -> Option<(usize, usize)> {
        for start_position in 0..monomial.degree() {
            if let Some(array_index) = self.lookup_at(monomial, direction, start_position) {
                return Some((array_index, start_position));
            }
        }
        None
    }

    /// Translate `GBNP.LookUpOccurTreeAllLstPTSLR`.
    ///
    /// Results are sorted by array index and retain the first position at
    /// which each represented monomial was found.
    #[must_use]
    pub fn all_occurrences(
        &self,
        monomial: &NpMonomial,
        direction: TreeDirection,
    ) -> Vec<(usize, usize)> {
        let mut first_positions = vec![None; self.len()];
        for start_position in 0..monomial.degree().max(1) {
            for array_index in self.lookup_all_at(monomial, direction, start_position) {
                if first_positions[array_index].is_none() {
                    first_positions[array_index] = Some(start_position);
                }
            }
        }
        first_positions
            .into_iter()
            .enumerate()
            .filter_map(|(array_index, position)| position.map(|position| (array_index, position)))
            .collect()
    }

    /// Translate `GBNP.LookUpOccurTreeForObsPTSLR`.
    ///
    /// `excluded_index` implements the source's special exclusion of a full
    /// match at the first position. Results retain sparse-trie traversal order.
    #[must_use]
    pub fn obstruction_matches(
        &self,
        monomial: &NpMonomial,
        excluded_index: Option<usize>,
        direction: TreeDirection,
    ) -> Vec<(usize, usize)> {
        let mut answer = Vec::new();
        for start_position in 0..monomial.degree() {
            for array_index in self.obstruction_matches_at(monomial, direction, start_position) {
                if !(start_position == 0 && excluded_index == Some(array_index)) {
                    answer.push((array_index, start_position));
                }
            }
        }
        answer
    }

    /// Swap two represented array indices as `GBNP.THeapOTSwap` does.
    ///
    /// # Errors
    ///
    /// Returns an error if either index is outside the represented array.
    pub fn swap_array_indices(
        &mut self,
        first: usize,
        second: usize,
    ) -> Result<(), OccurTreeError> {
        let length = self.len();
        for index in [first, second] {
            if index >= length {
                return Err(OccurTreeError::ArrayIndexOutOfRange {
                    index,
                    array_length: length,
                });
            }
        }
        self.array_to_tree.swap(first, second);
        let first_tree = self.array_to_tree[first];
        let second_tree = self.array_to_tree[second];
        self.tree_to_array[first_tree] = Some(first);
        self.tree_to_array[second_tree] = Some(second);
        Ok(())
    }

    fn obstruction_matches_at(
        &self,
        monomial: &NpMonomial,
        direction: TreeDirection,
        start_position: usize,
    ) -> Vec<usize> {
        let symbols = directed_symbols(monomial, direction).collect::<Vec<_>>();
        let Some(suffix) = symbols.get(start_position..) else {
            return Vec::new();
        };
        let mut node = &self.root;
        for &symbol in suffix {
            let Some(next) = node.children.get(&symbol) else {
                return Vec::new();
            };
            node = next;
        }
        let mut answer = Vec::new();
        self.flatten_terminal_indices(node, &mut answer);
        answer
    }

    fn flatten_terminal_indices(&self, node: &Node, answer: &mut Vec<usize>) {
        for child in node.children.range(..0).map(|(_, child)| child) {
            self.flatten_terminal_indices(child, answer);
        }
        if let Some(array_index) = self.terminal_array_index(node) {
            answer.push(array_index);
        }
        for child in node.children.range(0..).map(|(_, child)| child) {
            self.flatten_terminal_indices(child, answer);
        }
    }

    fn terminal_array_index(&self, node: &Node) -> Option<usize> {
        node.terminal
            .and_then(|number| self.tree_to_array.get(number).copied().flatten())
    }

    fn validate_monomial(&self, monomial: &NpMonomial) -> Result<(), OccurTreeError> {
        if monomial.generators().contains(&0) {
            return Err(OccurTreeError::ZeroGenerator);
        }
        if let Some(&symbol) = monomial
            .generators()
            .iter()
            .find(|&&symbol| symbol < 0 && symbol.unsigned_abs() as usize > self.prefix_generators)
        {
            Err(OccurTreeError::PrefixGeneratorOutOfRange {
                symbol,
                prefix_generators: self.prefix_generators,
            })
        } else {
            Ok(())
        }
    }

    fn allocate_tree_number(&mut self) -> usize {
        if let Some(number) = self.free_tree_numbers.pop() {
            number
        } else {
            let number = self.tree_to_array.len();
            self.tree_to_array.push(None);
            number
        }
    }

    fn rebuild_tree_to_array(&mut self) {
        self.tree_to_array.fill(None);
        for (array_index, &tree_number) in self.array_to_tree.iter().enumerate() {
            self.tree_to_array[tree_number] = Some(array_index);
        }
    }
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
struct Node {
    terminal: Option<usize>,
    children: BTreeMap<i32, Self>,
}

fn directed_symbols(
    monomial: &NpMonomial,
    direction: TreeDirection,
) -> impl Iterator<Item = i32> + '_ {
    let generators = monomial.generators();
    (0..generators.len()).map(move |position| match direction {
        TreeDirection::Forward => generators[position],
        TreeDirection::Reverse => generators[generators.len() - position - 1],
    })
}

fn remove_terminal(node: &mut Node, symbols: &[i32], tree_number: usize) -> bool {
    if let Some((&symbol, rest)) = symbols.split_first() {
        let Some(child) = node.children.get_mut(&symbol) else {
            return false;
        };
        if !remove_terminal(child, rest, tree_number) {
            return false;
        }
        if child.terminal.is_none() && child.children.is_empty() {
            node.children.remove(&symbol);
        }
        true
    } else if node.terminal == Some(tree_number) {
        node.terminal = None;
        true
    } else {
        // The source still remaps its parallel arrays for a duplicate path
        // whose terminal points at another stable tree number.
        true
    }
}

/// Failure in GBNP's source-translated occurrence tree.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum OccurTreeError {
    /// A mutable operation used an index outside the represented array.
    #[error("array index {index} is outside an occurrence array of length {array_length}")]
    ArrayIndexOutOfRange {
        /// Invalid zero-based array index.
        index: usize,
        /// Current represented array length.
        array_length: usize,
    },
    /// A negative generator lies outside the prefix-generator offset.
    #[error(
        "prefix generator {symbol} is outside the declared range of {prefix_generators} generators"
    )]
    PrefixGeneratorOutOfRange {
        /// Invalid negative generator number.
        symbol: i32,
        /// Declared number of prefix generators.
        prefix_generators: usize,
    },
    /// Zero collides with the terminal slot in GBNP's sparse trie layout.
    #[error("generator zero is not valid in a GBNP NP monomial")]
    ZeroGenerator,
    /// The requested word did not follow an existing trie path.
    #[error("the monomial is not represented by this occurrence tree path")]
    MissingMonomialPath,
}

#[cfg(test)]
mod tests {
    use crate::NpMonomial;

    use super::{OccurTree, TreeDirection};

    fn words(values: &[&[i32]]) -> Vec<NpMonomial> {
        values
            .iter()
            .map(|value| NpMonomial::from(*value))
            .collect()
    }

    #[test]
    fn lookup_prefers_shorter_prefix_at_the_leftmost_position() {
        let patterns = words(&[&[1, 2], &[1], &[2, 1]]);
        let tree = OccurTree::new(&patterns, 0, TreeDirection::Forward).unwrap();
        let word = NpMonomial::from(&[2, 1, 2][..]);
        assert_eq!(
            tree.first_occurrence(&word, TreeDirection::Forward),
            Some((2, 0))
        );
        assert_eq!(
            tree.all_occurrences(&word, TreeDirection::Forward),
            [(0, 1), (1, 1), (2, 0)]
        );
    }

    #[test]
    fn insertion_and_removal_remap_parallel_array_indices() {
        let patterns = words(&[&[1], &[2]]);
        let mut tree = OccurTree::new(&patterns, 0, TreeDirection::Forward).unwrap();
        let twelve = NpMonomial::from(&[1, 2][..]);
        tree.add(&twelve, 1, TreeDirection::Forward).unwrap();
        assert_eq!(
            tree.all_occurrences(&twelve, TreeDirection::Forward),
            [(0, 0), (1, 0), (2, 1)]
        );
        tree.remove(&twelve, 1, TreeDirection::Forward).unwrap();
        assert_eq!(
            tree.all_occurrences(&twelve, TreeDirection::Forward),
            [(0, 0), (1, 1)]
        );
    }

    #[test]
    fn reverse_tree_and_negative_prefix_generators_are_supported() {
        let patterns = words(&[&[-2, 1], &[1, 2]]);
        let tree = OccurTree::new(&patterns, 2, TreeDirection::Reverse).unwrap();
        assert_eq!(
            tree.first_occurrence(&NpMonomial::from(&[3, -2, 1][..]), TreeDirection::Reverse),
            Some((0, 0))
        );
    }

    #[test]
    fn obstruction_lookup_returns_extensions_after_consuming_the_suffix() {
        let patterns = words(&[&[1, 2, 1], &[2, 1, 2]]);
        let tree = OccurTree::new(&patterns, 0, TreeDirection::Forward).unwrap();
        assert_eq!(
            tree.obstruction_matches(&patterns[1], Some(1), TreeDirection::Forward),
            [(0, 1), (1, 2)]
        );
    }
}
