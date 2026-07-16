//! QPA finite posets.
//!
//! Construction translates `lib/qpaposet.gi`: QPA first forms a GAP domain,
//! adds the generating pairs, takes reflexive and transitive closures, checks
//! antisymmetry, and obtains the minimal relations from GAP's Hasse-diagram
//! routine. The closure loop below is the modified Floyd-Warshall method in
//! GAP's `lib/relation.gi` used by that pipeline.

use std::{fmt, sync::Arc};

use thiserror::Error;

/// A finite QPA poset with a source-compatible closed binary relation.
#[derive(Clone, Debug)]
pub struct Poset<Element>(Arc<PosetData<Element>>);

#[derive(Debug)]
struct PosetData<Element> {
    underlying_set: Box<[Element]>,
    relation: Box<[Box<[bool]>]>,
    minimal_relations: Box<[(Element, Element)]>,
}

impl<Element> Poset<Element>
where
    Element: Clone + Ord,
{
    /// Translate QPA's `Poset(points, relations)` constructor.
    ///
    /// Each nonempty relation row `[x, y, z, ...]` generates `x <= y`,
    /// `x <= z`, and so on. As in QPA, the underlying set is sorted and
    /// duplicate-free, while Hasse pairs retain first-occurrence source order.
    ///
    /// # Errors
    ///
    /// Returns an error for an empty underlying set, a relation element outside
    /// that set, or a generated preorder that is not antisymmetric.
    pub fn new(points: Vec<Element>, relations: &[Vec<Element>]) -> Result<Self, PosetError> {
        if points.is_empty() {
            return Err(PosetError::EmptyUnderlyingSet);
        }

        let mut source_order = Vec::new();
        for point in &points {
            if !source_order.contains(point) {
                source_order.push(point.clone());
            }
        }
        let mut underlying_set = points;
        underlying_set.sort();
        underlying_set.dedup();

        for element in relations.iter().flatten() {
            if underlying_set.binary_search(element).is_err() {
                return Err(PosetError::RelationElementOutsideUnderlyingSet);
            }
        }

        let size = underlying_set.len();
        let mut relation = vec![vec![false; size]; size];
        for row in relations {
            let Some(first) = row.first() else {
                continue;
            };
            let source = underlying_set
                .binary_search(first)
                .map_err(|_| PosetError::RelationElementOutsideUnderlyingSet)?;
            for target in &row[1..] {
                let target = underlying_set
                    .binary_search(target)
                    .map_err(|_| PosetError::RelationElementOutsideUnderlyingSet)?;
                relation[source][target] = true;
            }
        }

        // GAP's ReflexiveClosureBinaryRelation.
        for (index, successors) in relation.iter_mut().enumerate() {
            successors[index] = true;
        }

        // GAP relation.gi:544-552, its modified Floyd-Warshall closure.
        for intermediate in 0..size {
            for source in 0..size {
                if relation[source][intermediate] {
                    for target in 0..size {
                        relation[source][target] |= relation[intermediate][target];
                    }
                }
            }
        }

        for left in 0..size {
            for right in left + 1..size {
                if relation[left][right] && relation[right][left] {
                    return Err(PosetError::NotPartialOrder);
                }
            }
        }

        // GAP's HDBRListCovers iterates AsList(original Domain), hence the
        // first-occurrence source order, while ImagesElm supplies sorted targets.
        let mut minimal_relations = Vec::new();
        for source_element in source_order {
            let source = underlying_set
                .binary_search(&source_element)
                .map_err(|_| PosetError::ElementOutsideUnderlyingSet)?;
            let greater = (0..size)
                .filter(|&target| target != source && relation[source][target])
                .collect::<Vec<_>>();
            for &target in &greater {
                let has_intermediate = greater
                    .iter()
                    .copied()
                    .any(|candidate| candidate != target && relation[candidate][target]);
                if !has_intermediate {
                    minimal_relations
                        .push((source_element.clone(), underlying_set[target].clone()));
                }
            }
        }

        Ok(Self(Arc::new(PosetData {
            underlying_set: underlying_set.into_boxed_slice(),
            relation: relation
                .into_iter()
                .map(Vec::into_boxed_slice)
                .collect::<Vec<_>>()
                .into_boxed_slice(),
            minimal_relations: minimal_relations.into_boxed_slice(),
        })))
    }

    /// Return QPA's `Size` attribute.
    #[must_use]
    pub fn size(&self) -> usize {
        self.0.underlying_set.len()
    }

    /// Return QPA's sorted, duplicate-free `UnderlyingSet`.
    #[must_use]
    pub fn underlying_set(&self) -> &[Element] {
        &self.0.underlying_set
    }

    /// Return the minimal relations stored from QPA's Hasse diagram.
    #[must_use]
    pub fn minimal_relations(&self) -> &[(Element, Element)] {
        &self.0.minimal_relations
    }

    /// Apply QPA's reflexive `PartialOrderOfPoset` function.
    ///
    /// # Errors
    ///
    /// Returns an error when either argument is outside the underlying set.
    pub fn less_than_or_equal(&self, left: &Element, right: &Element) -> Result<bool, PosetError> {
        let left = self
            .0
            .underlying_set
            .binary_search(left)
            .map_err(|_| PosetError::ElementOutsideUnderlyingSet)?;
        let right = self
            .0
            .underlying_set
            .binary_search(right)
            .map_err(|_| PosetError::ElementOutsideUnderlyingSet)?;
        Ok(self.0.relation[left][right])
    }

    /// Return the callable represented by QPA's `PartialOrderOfPoset`.
    pub fn partial_order_of_poset(
        &self,
    ) -> impl Fn(&Element, &Element) -> Result<bool, PosetError> + '_ {
        move |left, right| self.less_than_or_equal(left, right)
    }

    /// Return QPA's `Display(P)` text, including its final newline.
    ///
    /// QPA prints the underlying GAP list with `PrintObj`, but the entries of
    /// each stored `DirectProductElement` with their compact view.  Rust's
    /// `Debug` and `Display` conventions reproduce that distinction for the
    /// scalar and string point types supported by this port.
    #[must_use]
    pub fn display_string(&self) -> String
    where
        Element: fmt::Debug + fmt::Display,
    {
        let underlying = self
            .underlying_set()
            .iter()
            .map(|element| format!("{element:?}"))
            .collect::<Vec<_>>()
            .join(", ");
        let minimal = self
            .minimal_relations()
            .iter()
            .map(|(left, right)| format!("DirectProductElement( [ {left}, {right} ] )"))
            .collect::<Vec<_>>()
            .join(", ");
        format!(
            "A poset with underlying set:\n[ {underlying} ]\nand partial order:\n[ {minimal} ]\n"
        )
    }

    /// Return whether two handles denote the identical QPA poset object.
    #[must_use]
    pub fn is_identical(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }
}

impl<Element> fmt::Display for Poset<Element> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "<A poset on {} points>",
            self.0.underlying_set.len()
        )
    }
}

/// Failure to construct or query a QPA poset.
#[derive(Clone, Copy, Debug, Eq, Error, PartialEq)]
pub enum PosetError {
    /// QPA rejects an empty first argument.
    #[error("the underlying set of a poset cannot be empty")]
    EmptyUnderlyingSet,
    /// At least one relation entry is absent from the underlying set.
    #[error("a relation element is outside the underlying set")]
    RelationElementOutsideUnderlyingSet,
    /// The reflexive/transitive closure is not antisymmetric.
    #[error("the entered relations do not define a partial order")]
    NotPartialOrder,
    /// A query argument is absent from the underlying set.
    #[error("the entered element is outside the underlying set")]
    ElementOutsideUnderlyingSet,
}

#[cfg(test)]
mod tests {
    use serde_json::{Value, json};

    use super::{Poset, PosetError};

    fn matrix<Element>(poset: &Poset<Element>) -> Value
    where
        Element: Clone + Ord,
    {
        json!(
            poset
                .underlying_set()
                .iter()
                .map(|left| {
                    poset
                        .underlying_set()
                        .iter()
                        .map(|right| poset.less_than_or_equal(left, right).unwrap())
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        )
    }

    #[test]
    fn pinned_qpa_poset_fixture_matches() {
        let expected: Value =
            serde_json::from_str(include_str!("../../../tests/fixtures/poset.schema-v1.json"))
                .unwrap();

        let example = Poset::new(
            vec!["d", "b", "a", "c"],
            &[vec!["a", "b", "c"], vec!["b", "d"], vec!["c", "d"]],
        )
        .unwrap();
        assert_eq!(example.size(), expected["example"]["size"]);
        assert_eq!(
            json!(example.underlying_set()),
            expected["example"]["underlying"]
        );
        assert_eq!(
            json!(example.minimal_relations()),
            expected["example"]["minimal"]
        );
        assert_eq!(matrix(&example), expected["example"]["matrix"]);

        let duplicates = Poset::new(vec![3, 1, 2, 1], &[vec![1, 2], vec![2, 3]]).unwrap();
        assert_eq!(duplicates.size(), expected["duplicates"]["size"]);
        assert_eq!(
            json!(
                duplicates
                    .underlying_set()
                    .iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>()
            ),
            expected["duplicates"]["underlying"]
        );
        assert_eq!(
            json!(
                duplicates
                    .minimal_relations()
                    .iter()
                    .map(|(left, right)| [left.to_string(), right.to_string()])
                    .collect::<Vec<_>>()
            ),
            expected["duplicates"]["minimal"]
        );
        assert_eq!(matrix(&duplicates), expected["duplicates"]["matrix"]);

        let self_relation = Poset::new(vec![2, 1], &[vec![1, 1]]).unwrap();
        assert_eq!(
            json!(
                self_relation
                    .underlying_set()
                    .iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>()
            ),
            expected["self_relation"]["underlying"]
        );
        assert!(self_relation.minimal_relations().is_empty());
        assert_eq!(matrix(&self_relation), expected["self_relation"]["matrix"]);
    }

    #[test]
    fn qpa_poset_errors_and_identity_are_explicit() {
        assert!(matches!(
            Poset::<i32>::new(vec![], &[]),
            Err(PosetError::EmptyUnderlyingSet)
        ));
        assert!(matches!(
            Poset::new(vec![1, 2], &[vec![1, 3]]),
            Err(PosetError::RelationElementOutsideUnderlyingSet)
        ));
        assert!(matches!(
            Poset::new(vec![1, 2], &[vec![1, 2], vec![2, 1]]),
            Err(PosetError::NotPartialOrder)
        ));

        let poset = Poset::new(vec![1, 2], &[vec![1, 2]]).unwrap();
        let clone = poset.clone();
        let equal_data = Poset::new(vec![1, 2], &[vec![1, 2]]).unwrap();
        assert!(poset.is_identical(&clone));
        assert!(!poset.is_identical(&equal_data));
        assert_eq!(poset.to_string(), "<A poset on 2 points>");
        assert_eq!(
            poset.less_than_or_equal(&1, &3),
            Err(PosetError::ElementOutsideUnderlyingSet)
        );
    }

    #[test]
    fn pinned_qpa_poset_output_fixture_matches() {
        let example = Poset::new(
            vec!["d", "b", "a", "c"],
            &[vec!["a", "b", "c"], vec!["b", "d"], vec!["c", "d"]],
        )
        .unwrap();
        assert_eq!(
            format!("{example}\n{}", example.display_string()),
            include_str!("../../../tests/fixtures/poset_output.txt")
        );
    }
}
