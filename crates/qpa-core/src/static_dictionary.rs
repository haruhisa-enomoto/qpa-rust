//! QPA's static dictionary for forbidden quiver paths.
//!
//! This is a direct Rust model of `lib/groebnerstadict.gi`. QPA builds an
//! Aho--Corasick-style trie from Groebner tips, completes it with quiver
//! transitions, detects cycles in the nonaccepting automaton, and uses the
//! resulting DAG to count or enumerate nontip paths.

use std::collections::VecDeque;

use numerica::domains::integer::Integer;
use thiserror::Error;

use crate::{ArrowIndex, Path, PathEndpoint, Quiver, QuiverError, VertexIndex};

#[derive(Clone, Debug)]
struct AutomatonState {
    target: VertexIndex,
    trie_transitions: Vec<Option<usize>>,
    transitions: Vec<Option<usize>>,
    failure: Option<usize>,
    own_patterns: Vec<usize>,
    outputs: Vec<usize>,
}

impl AutomatonState {
    fn new(target: VertexIndex, number_of_arrows: usize) -> Self {
        Self {
            target,
            trie_transitions: vec![None; number_of_arrows],
            transitions: vec![None; number_of_arrows],
            failure: None,
            own_patterns: Vec::new(),
            outputs: Vec::new(),
        }
    }

    fn is_accepting(&self) -> bool {
        !self.outputs.is_empty()
    }
}

/// One QPA static-dictionary search match.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DictionaryMatch {
    /// One-based arrow position where the match ends; vertex-prefix matches use
    /// zero, matching QPA's `PrefixSearch`.
    pub end_position: usize,
    /// One-based indices into [`QuiverStaticDictionary::patterns`].
    pub pattern_indices: Vec<usize>,
}

/// A QPA quiver static dictionary built from forbidden path patterns.
#[derive(Clone, Debug)]
pub struct QuiverStaticDictionary {
    quiver: Quiver,
    patterns: Vec<Path>,
    states: Vec<AutomatonState>,
    vertex_states: Vec<usize>,
}

impl QuiverStaticDictionary {
    /// Translate QPA's `QuiverStaticDictionary(quiver, patterns)`.
    ///
    /// # Errors
    ///
    /// Returns an error for a pattern from another quiver or QPA's zero path.
    #[allow(clippy::too_many_lines)]
    pub fn new(quiver: &Quiver, patterns: &[Path]) -> Result<Self, StaticDictionaryError> {
        let number_of_arrows = quiver.number_of_arrows();
        let mut states = Vec::with_capacity(quiver.number_of_vertices() + patterns.len());
        let mut vertex_states = Vec::with_capacity(quiver.number_of_vertices());
        for vertex in quiver.vertices() {
            vertex_states.push(states.len());
            states.push(AutomatonState::new(vertex.index(), number_of_arrows));
        }

        for (pattern_index, pattern) in patterns.iter().enumerate() {
            if !pattern.quiver().is_identical(quiver) {
                return Err(StaticDictionaryError::DifferentQuivers);
            }
            let source = endpoint_index(pattern.source())?;
            let walk = pattern.walk().ok_or(StaticDictionaryError::ZeroPattern)?;
            let mut current = vertex_states[source.zero_based()];
            for arrow_index in walk {
                let position = arrow_index.zero_based();
                let next = if let Some(next) = states[current].trie_transitions[position] {
                    next
                } else {
                    let arrow = quiver
                        .arrow(*arrow_index)
                        .ok_or(StaticDictionaryError::InvalidPattern)?;
                    let next = states.len();
                    states.push(AutomatonState::new(
                        arrow.target().index(),
                        number_of_arrows,
                    ));
                    states[current].trie_transitions[position] = Some(next);
                    next
                };
                current = next;
            }
            states[current].own_patterns.push(pattern_index + 1);
        }

        let mut breadth_first = Vec::new();
        let mut queue = VecDeque::new();
        for &vertex_state in &vertex_states {
            breadth_first.push(vertex_state);
            let children = states[vertex_state]
                .trie_transitions
                .iter()
                .enumerate()
                .filter_map(|(arrow, child)| child.map(|child| (arrow, child)))
                .collect::<Vec<_>>();
            for (arrow, child) in children {
                let target = quiver
                    .arrow(ArrowIndex::from_zero_based(arrow))
                    .ok_or(StaticDictionaryError::InvalidPattern)?
                    .target()
                    .index();
                states[child].failure = Some(vertex_states[target.zero_based()]);
                queue.push_back(child);
            }
        }

        while let Some(current) = queue.pop_front() {
            breadth_first.push(current);
            let failure_outputs = states[current]
                .failure
                .map(|failure| states[failure].outputs.clone())
                .unwrap_or_default();
            let own_patterns = states[current].own_patterns.clone();
            states[current].outputs.clone_from(&own_patterns);
            states[current].outputs.extend(failure_outputs);

            let children = states[current]
                .trie_transitions
                .iter()
                .enumerate()
                .filter_map(|(arrow, child)| child.map(|child| (arrow, child)))
                .collect::<Vec<_>>();
            for (arrow, child) in children {
                let failure = failure_transition(
                    &states,
                    &vertex_states,
                    quiver,
                    states[current].failure,
                    ArrowIndex::from_zero_based(arrow),
                )?;
                states[child].failure = Some(failure);
                queue.push_back(child);
            }
        }

        for &vertex_state in &vertex_states {
            let own_patterns = states[vertex_state].own_patterns.clone();
            states[vertex_state].outputs.clone_from(&own_patterns);
        }
        for &state_index in &breadth_first {
            let trie_transitions = states[state_index].trie_transitions.clone();
            states[state_index]
                .transitions
                .clone_from(&trie_transitions);
            let vertex = quiver
                .vertex(states[state_index].target)
                .ok_or(StaticDictionaryError::InvalidPattern)?;
            for arrow in quiver.outgoing_arrows(&vertex)? {
                let position = arrow.index().zero_based();
                if states[state_index].transitions[position].is_some() {
                    continue;
                }
                let fallback = states[state_index].failure.map_or(
                    vertex_states[arrow.target().index().zero_based()],
                    |failure| {
                        states[failure].transitions[position]
                            .unwrap_or(vertex_states[arrow.target().index().zero_based()])
                    },
                );
                states[state_index].transitions[position] = Some(fallback);
            }
        }

        Ok(Self {
            quiver: quiver.clone(),
            patterns: patterns.to_vec(),
            states,
            vertex_states,
        })
    }

    /// Return QPA's immutable pattern list.
    #[must_use]
    pub fn patterns(&self) -> &[Path] {
        &self.patterns
    }

    /// Translate QPA's `Search(dict, path)`.
    ///
    /// # Errors
    ///
    /// Returns an error for a path from another quiver or the zero path.
    pub fn search(&self, path: &Path) -> Result<Vec<DictionaryMatch>, StaticDictionaryError> {
        self.validate_path(path)?;
        let source = endpoint_index(path.source())?;
        let walk = path.walk().ok_or(StaticDictionaryError::ZeroPattern)?;
        let mut current = self.vertex_states[source.zero_based()];
        let mut matches = Vec::new();
        for (position, arrow) in walk.iter().enumerate() {
            current = self.states[current].transitions[arrow.zero_based()]
                .ok_or(StaticDictionaryError::InvalidPattern)?;
            if !self.states[current].outputs.is_empty() {
                matches.push(DictionaryMatch {
                    end_position: position + 1,
                    pattern_indices: self.states[current].outputs.clone(),
                });
            }
        }
        Ok(matches)
    }

    /// Translate QPA's `PrefixSearch(dict, path)`.
    ///
    /// # Errors
    ///
    /// Returns an error for a path from another quiver or the zero path.
    pub fn prefix_search(
        &self,
        path: &Path,
    ) -> Result<Vec<DictionaryMatch>, StaticDictionaryError> {
        self.validate_path(path)?;
        let source = endpoint_index(path.source())?;
        let walk = path.walk().ok_or(StaticDictionaryError::ZeroPattern)?;
        let mut current = self.vertex_states[source.zero_based()];
        let mut matches = Vec::new();
        if !self.states[current].own_patterns.is_empty() {
            matches.push(DictionaryMatch {
                end_position: 0,
                pattern_indices: self.states[current].own_patterns.clone(),
            });
        }
        for (position, arrow) in walk.iter().enumerate() {
            let Some(next) = self.states[current].trie_transitions[arrow.zero_based()] else {
                break;
            };
            current = next;
            if !self.states[current].own_patterns.is_empty() {
                matches.push(DictionaryMatch {
                    end_position: position + 1,
                    pattern_indices: self.states[current].own_patterns.clone(),
                });
            }
        }
        Ok(matches)
    }

    /// Translate QPA's `IsFiniteDifference` cycle test.
    #[must_use]
    pub fn is_finite_difference(&self) -> bool {
        self.topological_order().is_some()
    }

    /// Translate QPA's `DifferenceSize`; `None` represents `infinity`.
    #[must_use]
    pub fn difference_size(&self) -> Option<Integer> {
        let order = self.topological_order()?;
        let mut counts = vec![Integer::from(0); self.states.len()];
        for &vertex in &self.vertex_states {
            counts[vertex] += 1_u8;
        }
        let mut total = Integer::from(0);
        for state in order {
            if self.states[state].is_accepting() {
                continue;
            }
            let contribution = counts[state].clone();
            total += &contribution;
            for target in self.nonaccepting_transitions(state) {
                counts[target] += &contribution;
            }
        }
        Some(total)
    }

    /// Translate QPA's `DifferenceWords`; `None` represents `fail` for an
    /// infinite difference language.
    #[must_use]
    pub fn difference_words(&self) -> Option<Vec<Path>> {
        let order = self.topological_order()?;
        let mut state_paths = vec![Vec::new(); self.states.len()];
        for (position, &state) in self.vertex_states.iter().enumerate() {
            if let Some(vertex) = self.quiver.vertex(VertexIndex::from_zero_based(position)) {
                state_paths[state].push(vertex.path());
            }
        }
        let mut paths = Vec::new();
        for state in order {
            if self.states[state].is_accepting() {
                continue;
            }
            paths.extend(state_paths[state].iter().cloned());
            let current_paths = state_paths[state].clone();
            for (arrow_index, target) in self.indexed_nonaccepting_transitions(state) {
                let Some(arrow) = self.quiver.arrow(arrow_index) else {
                    continue;
                };
                for path in &current_paths {
                    if let Ok(extension) = path.multiply(&arrow.path()) {
                        state_paths[target].push(extension);
                    }
                }
            }
        }
        paths.sort_by(|left, right| {
            self.quiver
                .compare_paths(left, right)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        Some(paths)
    }

    fn validate_path(&self, path: &Path) -> Result<(), StaticDictionaryError> {
        if path.quiver().is_identical(&self.quiver) {
            Ok(())
        } else {
            Err(StaticDictionaryError::DifferentQuivers)
        }
    }

    fn nonaccepting_transitions(&self, state: usize) -> impl Iterator<Item = usize> + '_ {
        self.indexed_nonaccepting_transitions(state)
            .map(|(_, target)| target)
    }

    fn indexed_nonaccepting_transitions(
        &self,
        state: usize,
    ) -> impl Iterator<Item = (ArrowIndex, usize)> + '_ {
        let accepting = self.states[state].is_accepting();
        self.states[state]
            .transitions
            .iter()
            .enumerate()
            .filter_map(move |(arrow, target)| {
                (!accepting)
                    .then(|| target.map(|target| (ArrowIndex::from_zero_based(arrow), target)))
                    .flatten()
            })
    }

    fn topological_order(&self) -> Option<Vec<usize>> {
        let mut indegree = vec![0_usize; self.states.len()];
        for state in 0..self.states.len() {
            for target in self.nonaccepting_transitions(state) {
                indegree[target] += 1;
            }
        }
        let mut queue = indegree
            .iter()
            .enumerate()
            .filter_map(|(state, &degree)| (degree == 0).then_some(state))
            .collect::<VecDeque<_>>();
        let mut order = Vec::with_capacity(self.states.len());
        while let Some(state) = queue.pop_front() {
            order.push(state);
            for target in self.nonaccepting_transitions(state) {
                indegree[target] -= 1;
                if indegree[target] == 0 {
                    queue.push_back(target);
                }
            }
        }
        (order.len() == self.states.len()).then_some(order)
    }
}

fn failure_transition(
    states: &[AutomatonState],
    vertex_states: &[usize],
    quiver: &Quiver,
    mut failure: Option<usize>,
    arrow_index: ArrowIndex,
) -> Result<usize, StaticDictionaryError> {
    while let Some(state) = failure {
        if let Some(target) = states[state].trie_transitions[arrow_index.zero_based()] {
            return Ok(target);
        }
        failure = states[state].failure;
    }
    let target = quiver
        .arrow(arrow_index)
        .ok_or(StaticDictionaryError::InvalidPattern)?
        .target()
        .index();
    Ok(vertex_states[target.zero_based()])
}

fn endpoint_index(endpoint: PathEndpoint) -> Result<VertexIndex, StaticDictionaryError> {
    match endpoint {
        PathEndpoint::Vertex(vertex) => Ok(vertex.index()),
        PathEndpoint::Zero(_) => Err(StaticDictionaryError::ZeroPattern),
    }
}

/// Failure to construct or query QPA's quiver static dictionary.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum StaticDictionaryError {
    /// Patterns and searched paths must belong to the dictionary's quiver.
    #[error("a dictionary path belongs to a different quiver")]
    DifferentQuivers,
    /// QPA's zero path has no source or walk for dictionary construction.
    #[error("the zero path cannot be a static-dictionary pattern or query")]
    ZeroPattern,
    /// A pattern contains an invalid quiver transition.
    #[error("a static-dictionary pattern is not a valid quiver path")]
    InvalidPattern,
    /// A QPA quiver adjacency query failed.
    #[error(transparent)]
    Quiver(#[from] QuiverError),
}

#[cfg(test)]
mod tests {
    use crate::{ArrowIndex, ArrowSpec, Quiver, VertexIndex};

    use super::{DictionaryMatch, QuiverStaticDictionary};

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn arrow(index: usize) -> ArrowIndex {
        ArrowIndex::new(index).unwrap()
    }

    fn two_loops() -> Quiver {
        Quiver::from_indexed(
            1,
            vec![
                ArrowSpec::named(vertex(1), vertex(1), "x"),
                ArrowSpec::named(vertex(1), vertex(1), "y"),
            ],
        )
        .unwrap()
    }

    #[test]
    fn finite_difference_enumerates_qpa_nontips() {
        let quiver = two_loops();
        let x = quiver.arrow(arrow(1)).unwrap().path();
        let y = quiver.arrow(arrow(2)).unwrap().path();
        let patterns = [
            x.multiply(&x).unwrap(),
            y.multiply(&x).unwrap(),
            y.multiply(&y).unwrap(),
        ];
        let dictionary = QuiverStaticDictionary::new(&quiver, &patterns).unwrap();
        assert!(dictionary.is_finite_difference());
        assert_eq!(dictionary.difference_size().unwrap(), 4);
        assert_eq!(
            dictionary
                .difference_words()
                .unwrap()
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>(),
            ["v1", "x", "y", "x*y"]
        );
    }

    #[test]
    fn unrestricted_loop_keeps_difference_infinite() {
        let quiver = two_loops();
        let x = quiver.arrow(arrow(1)).unwrap().path();
        let dictionary = QuiverStaticDictionary::new(&quiver, &[x.multiply(&x).unwrap()]).unwrap();
        assert!(!dictionary.is_finite_difference());
        assert!(dictionary.difference_size().is_none());
        assert!(dictionary.difference_words().is_none());
    }

    #[test]
    fn search_reports_one_based_pattern_and_end_positions() {
        let quiver = two_loops();
        let x = quiver.arrow(arrow(1)).unwrap().path();
        let y = quiver.arrow(arrow(2)).unwrap().path();
        let xx = x.multiply(&x).unwrap();
        let yxx = y.multiply(&xx).unwrap();
        let dictionary = QuiverStaticDictionary::new(&quiver, &[xx]).unwrap();
        assert_eq!(
            dictionary.search(&yxx).unwrap(),
            vec![DictionaryMatch {
                end_position: 3,
                pattern_indices: vec![1]
            }]
        );
    }
}
