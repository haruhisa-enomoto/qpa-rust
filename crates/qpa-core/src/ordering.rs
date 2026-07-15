//! QPA quiver orderings.
//!
//! The comparison routines translate the constructors and
//! `SetComparisonFunction` closures in QPA 1.37
//! `lib/groebnerordering.gi`. They compare QPA path external
//! representations; no generic Rust ordering is substituted.

use std::{cmp::Ordering, sync::Arc};

use numerica::domains::integer::Integer;
use thiserror::Error;

use crate::{Generator, Path, Quiver};

/// A source-translated QPA ordering tied to a quiver's generator positions.
#[derive(Clone, Debug)]
pub struct QuiverOrdering {
    quiver: Quiver,
    kind: Arc<OrderingKind>,
    lexicographic_table: Option<Arc<LexicographicTable>>,
}

#[derive(Clone, Debug)]
enum OrderingKind {
    Length {
        next: Option<QuiverOrdering>,
    },
    LeftLexicographic {
        table: Arc<LexicographicTable>,
    },
    RightLexicographic {
        table: Arc<LexicographicTable>,
    },
    Reverse {
        next: QuiverOrdering,
    },
    LeftVector {
        table: Arc<LexicographicTable>,
        next: Option<QuiverOrdering>,
    },
    RightVector {
        table: Arc<LexicographicTable>,
        next: Option<QuiverOrdering>,
    },
    Weight {
        weights: Arc<[Integer]>,
        next: Option<QuiverOrdering>,
    },
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct LexicographicTable {
    /// Generator positions in increasing lexicographic order.
    generators: Box<[usize]>,
    /// `ranks[generator_position]` is its one-based lexicographic rank.
    ranks: Box<[usize]>,
}

impl QuiverOrdering {
    /// Construct QPA's default
    /// `LengthOrdering(Q, LeftLexicographicOrdering(Q, generators))`.
    ///
    /// # Errors
    ///
    /// Construction is infallible for a valid quiver; the result type remains
    /// explicit so it composes with the other checked constructors.
    pub fn default_for(quiver: &Quiver) -> Result<Self, OrderingError> {
        let lexicographic = Self::left_lexicographic(quiver, &[])?;
        Self::length(quiver, Some(lexicographic))
    }

    /// Translate QPA's `LengthOrdering(Q)` or `LengthOrdering(Q, next)`.
    ///
    /// # Errors
    ///
    /// Returns an error if `next` belongs to a different quiver.
    pub fn length(quiver: &Quiver, next: Option<Self>) -> Result<Self, OrderingError> {
        validate_next(quiver, next.as_ref())?;
        let lexicographic_table = next
            .as_ref()
            .and_then(|ordering| ordering.lexicographic_table.clone());
        Ok(Self {
            quiver: quiver.clone(),
            kind: Arc::new(OrderingKind::Length { next }),
            lexicographic_table,
        })
    }

    /// Translate QPA's `LeftLexicographicOrdering(Q, preferred)`.
    ///
    /// Preferred vertices are moved first among vertices and preferred arrows
    /// first among arrows; duplicates are ignored, matching the GAP source.
    ///
    /// # Errors
    ///
    /// Returns an error if a preferred generator belongs to another quiver.
    pub fn left_lexicographic(
        quiver: &Quiver,
        preferred: &[Generator],
    ) -> Result<Self, OrderingError> {
        let table = Arc::new(build_lexicographic_table(
            quiver,
            preferred,
            &default_generator_positions(quiver),
        )?);
        Ok(Self {
            quiver: quiver.clone(),
            kind: Arc::new(OrderingKind::LeftLexicographic {
                table: table.clone(),
            }),
            lexicographic_table: Some(table),
        })
    }

    /// Translate QPA's `RightLexicographicOrdering(Q, preferred)`.
    ///
    /// # Errors
    ///
    /// Returns an error if a preferred generator belongs to another quiver.
    pub fn right_lexicographic(
        quiver: &Quiver,
        preferred: &[Generator],
    ) -> Result<Self, OrderingError> {
        let table = Arc::new(build_lexicographic_table(
            quiver,
            preferred,
            &default_generator_positions(quiver),
        )?);
        Ok(Self {
            quiver: quiver.clone(),
            kind: Arc::new(OrderingKind::RightLexicographic {
                table: table.clone(),
            }),
            lexicographic_table: Some(table),
        })
    }

    /// Translate QPA's `ReverseOrdering(Q, next)`.
    ///
    /// QPA reverses comparisons only when both inputs have positive path
    /// length; zero paths and vertices are delegated unchanged.
    ///
    /// # Errors
    ///
    /// Returns an error if `next` belongs to a different quiver.
    pub fn reverse(quiver: &Quiver, next: Self) -> Result<Self, OrderingError> {
        validate_next(quiver, Some(&next))?;
        Ok(Self {
            quiver: quiver.clone(),
            lexicographic_table: next.lexicographic_table.clone(),
            kind: Arc::new(OrderingKind::Reverse { next }),
        })
    }

    /// Translate QPA's `LeftVectorOrdering(Q, preferred, next)`.
    ///
    /// # Errors
    ///
    /// Returns an error for foreign preferred generators, a foreign next
    /// ordering, or QPA's "Must use same lexicographic table" condition.
    pub fn left_vector(
        quiver: &Quiver,
        preferred: &[Generator],
        next: Option<Self>,
    ) -> Result<Self, OrderingError> {
        Self::vector(quiver, preferred, next, false)
    }

    /// Translate QPA's `RightVectorOrdering(Q, preferred, next)`.
    ///
    /// # Errors
    ///
    /// Returns an error for foreign preferred generators, a foreign next
    /// ordering, or QPA's "Must use same lexicographic table" condition.
    pub fn right_vector(
        quiver: &Quiver,
        preferred: &[Generator],
        next: Option<Self>,
    ) -> Result<Self, OrderingError> {
        Self::vector(quiver, preferred, next, true)
    }

    fn vector(
        quiver: &Quiver,
        preferred: &[Generator],
        next: Option<Self>,
        right: bool,
    ) -> Result<Self, OrderingError> {
        validate_next(quiver, next.as_ref())?;
        let fallback = next
            .as_ref()
            .and_then(|ordering| ordering.lexicographic_table.as_ref())
            .map_or_else(
                || default_generator_positions(quiver),
                |table| table.generators.to_vec(),
            );
        let table = Arc::new(build_lexicographic_table(quiver, preferred, &fallback)?);
        if let Some(next_table) = next
            .as_ref()
            .and_then(|ordering| ordering.lexicographic_table.as_ref())
            && **next_table != *table
        {
            return Err(OrderingError::ConflictingLexicographicTables);
        }

        let kind = if right {
            OrderingKind::RightVector {
                table: table.clone(),
                next,
            }
        } else {
            OrderingKind::LeftVector {
                table: table.clone(),
                next,
            }
        };
        Ok(Self {
            quiver: quiver.clone(),
            kind: Arc::new(kind),
            lexicographic_table: Some(table),
        })
    }

    /// Translate QPA's `WeightOrdering(Q, weights, next)` for exact integer
    /// weights.
    ///
    /// Missing arrow weights default to zero and excess weights are ignored by
    /// comparisons, reproducing how QPA fills and subsequently indexes its
    /// weight list.
    ///
    /// # Errors
    ///
    /// Returns an error if `next` belongs to a different quiver.
    pub fn weight(
        quiver: &Quiver,
        weights: &[Integer],
        next: Option<Self>,
    ) -> Result<Self, OrderingError> {
        validate_next(quiver, next.as_ref())?;
        let mut stored = vec![Integer::from(0); quiver.number_of_arrows()];
        for (target, source) in stored.iter_mut().zip(weights) {
            *target = source.clone();
        }
        let lexicographic_table = next
            .as_ref()
            .and_then(|ordering| ordering.lexicographic_table.clone());
        Ok(Self {
            quiver: quiver.clone(),
            kind: Arc::new(OrderingKind::Weight {
                weights: stored.into(),
                next,
            }),
            lexicographic_table,
        })
    }

    /// Return QPA's ordering name.
    #[must_use]
    pub fn name(&self) -> &'static str {
        match &*self.kind {
            OrderingKind::Length { .. } => "length",
            OrderingKind::LeftLexicographic { .. } => "left lexicographic",
            OrderingKind::RightLexicographic { .. } => "right lexicographic",
            OrderingKind::Reverse { .. } => "reverse",
            OrderingKind::LeftVector { .. } => "left vector",
            OrderingKind::RightVector { .. } => "right vector",
            OrderingKind::Weight { .. } => "weight",
        }
    }

    /// Return the containing quiver.
    #[must_use]
    pub const fn quiver(&self) -> &Quiver {
        &self.quiver
    }

    /// Compare two paths using QPA's source closure.
    ///
    /// # Errors
    ///
    /// Returns an error if either path belongs to another quiver.
    pub fn compare(&self, left: &Path, right: &Path) -> Result<Ordering, OrderingError> {
        if !left.quiver().is_identical(&self.quiver) || !right.quiver().is_identical(&self.quiver) {
            return Err(OrderingError::DifferentQuivers);
        }
        Ok(compare_kind(&self.kind, left, right))
    }

    /// Return QPA's `LessThanByOrdering` result.
    ///
    /// # Errors
    ///
    /// Returns an error if either path belongs to another quiver.
    pub fn less_than(&self, left: &Path, right: &Path) -> Result<bool, OrderingError> {
        self.compare(left, right).map(Ordering::is_lt)
    }
}

fn compare_kind(kind: &OrderingKind, left: &Path, right: &Path) -> Ordering {
    match kind {
        OrderingKind::Length { next } => left
            .length()
            .unwrap_or(0)
            .cmp(&right.length().unwrap_or(0))
            .then_with(|| compare_optional(next.as_ref(), left, right)),
        OrderingKind::LeftLexicographic { table } => {
            compare_lexicographic(table, left, right, false)
        }
        OrderingKind::RightLexicographic { table } => {
            compare_lexicographic(table, left, right, true)
        }
        OrderingKind::Reverse { next } => {
            let comparison = compare_kind(&next.kind, left, right);
            if left.is_zero()
                || right.is_zero()
                || (left.is_vertex() && !right.is_vertex())
                || (!left.is_vertex() && right.is_vertex())
            {
                comparison
            } else {
                comparison.reverse()
            }
        }
        OrderingKind::LeftVector { table, next } => {
            compare_vector(table, next.as_ref(), left, right, false)
        }
        OrderingKind::RightVector { table, next } => {
            compare_vector(table, next.as_ref(), left, right, true)
        }
        OrderingKind::Weight { weights, next } => {
            compare_weight(weights, next.as_ref(), left, right)
        }
    }
}

fn compare_lexicographic(
    table: &LexicographicTable,
    left: &Path,
    right: &Path,
    reverse_words: bool,
) -> Ordering {
    if left.is_zero() {
        return if right.is_zero() {
            Ordering::Equal
        } else {
            Ordering::Less
        };
    }
    if right.is_zero() {
        return Ordering::Greater;
    }

    let left_representation = left.external_representation();
    let right_representation = right.external_representation();
    if reverse_words {
        compare_ranked_iterators(
            table,
            left_representation.iter().rev(),
            right_representation.iter().rev(),
        )
    } else {
        compare_ranked_iterators(
            table,
            left_representation.iter(),
            right_representation.iter(),
        )
    }
}

fn compare_ranked_iterators<'a>(
    table: &LexicographicTable,
    mut left: impl ExactSizeIterator<Item = &'a usize>,
    mut right: impl ExactSizeIterator<Item = &'a usize>,
) -> Ordering {
    let left_length = left.len();
    let right_length = right.len();
    for (left, right) in left.by_ref().zip(right.by_ref()) {
        let comparison = table.ranks[*left].cmp(&table.ranks[*right]);
        if comparison != Ordering::Equal {
            return comparison;
        }
    }
    left_length.cmp(&right_length)
}

fn compare_vector(
    table: &LexicographicTable,
    next: Option<&QuiverOrdering>,
    left: &Path,
    right: &Path,
    reverse_vector: bool,
) -> Ordering {
    let left_short = left.is_zero() || left.is_vertex();
    let right_short = right.is_zero() || right.is_vertex();
    match (left_short, right_short) {
        (true, true) => return compare_optional(next, left, right),
        (true, false) => return Ordering::Less,
        (false, true) => return Ordering::Greater,
        (false, false) => {}
    }

    let vertices = left.quiver().number_of_vertices();
    let arrows = left.quiver().number_of_arrows();
    let counts = |path: &Path| {
        let mut result = vec![0_usize; arrows];
        for generator in path.external_representation() {
            let position = table.ranks[generator] - vertices - 1;
            result[position] += 1;
        }
        result
    };
    let left_counts = counts(left);
    let right_counts = counts(right);
    let comparison = if reverse_vector {
        left_counts.iter().rev().cmp(right_counts.iter().rev())
    } else {
        left_counts.iter().cmp(right_counts.iter())
    };
    comparison.then_with(|| compare_optional(next, left, right))
}

fn compare_weight(
    weights: &[Integer],
    next: Option<&QuiverOrdering>,
    left: &Path,
    right: &Path,
) -> Ordering {
    let left_short = left.is_zero() || left.is_vertex();
    let right_short = right.is_zero() || right.is_vertex();
    match (left_short, right_short) {
        (true, true) => return compare_optional(next, left, right),
        (true, false) => return Ordering::Less,
        (false, true) => return Ordering::Greater,
        (false, false) => {}
    }

    let vertex_count = left.quiver().number_of_vertices();
    let total = |path: &Path| {
        path.external_representation()
            .into_iter()
            .fold(Integer::from(0), |sum, generator| {
                sum + weights[generator - vertex_count - 1].clone()
            })
    };
    total(left)
        .cmp(&total(right))
        .then_with(|| compare_optional(next, left, right))
}

fn compare_optional(next: Option<&QuiverOrdering>, left: &Path, right: &Path) -> Ordering {
    next.map_or(Ordering::Equal, |ordering| {
        compare_kind(&ordering.kind, left, right)
    })
}

fn build_lexicographic_table(
    quiver: &Quiver,
    preferred: &[Generator],
    fallback: &[usize],
) -> Result<LexicographicTable, OrderingError> {
    for generator in preferred {
        if !generator.path().quiver().is_identical(quiver) {
            return Err(OrderingError::DifferentQuivers);
        }
    }
    let vertex_count = quiver.number_of_vertices();
    let generator_count = vertex_count + quiver.number_of_arrows();
    let mut generators = Vec::with_capacity(generator_count);

    append_preferred(&mut generators, preferred, |generator| {
        matches!(generator, Generator::Vertex(_))
    });
    append_positions(
        &mut generators,
        fallback
            .iter()
            .copied()
            .filter(|position| *position <= vertex_count),
    );
    append_preferred(&mut generators, preferred, |generator| {
        matches!(generator, Generator::Arrow(_))
    });
    append_positions(
        &mut generators,
        fallback
            .iter()
            .copied()
            .filter(|position| *position > vertex_count),
    );

    let mut ranks = vec![0; generator_count + 1];
    for (rank, &position) in generators.iter().enumerate() {
        ranks[position] = rank + 1;
    }
    Ok(LexicographicTable {
        generators: generators.into_boxed_slice(),
        ranks: ranks.into_boxed_slice(),
    })
}

fn append_preferred(
    positions: &mut Vec<usize>,
    preferred: &[Generator],
    predicate: impl Fn(&Generator) -> bool,
) {
    append_positions(
        positions,
        preferred
            .iter()
            .filter(|generator| predicate(generator))
            .map(generator_position),
    );
}

fn append_positions(positions: &mut Vec<usize>, source: impl IntoIterator<Item = usize>) {
    for position in source {
        if !positions.contains(&position) {
            positions.push(position);
        }
    }
}

fn generator_position(generator: &Generator) -> usize {
    match generator {
        Generator::Vertex(vertex) => vertex.index().get(),
        Generator::Arrow(arrow) => arrow.generator_index(),
    }
}

fn default_generator_positions(quiver: &Quiver) -> Vec<usize> {
    (1..=quiver.number_of_vertices() + quiver.number_of_arrows()).collect()
}

fn validate_next(quiver: &Quiver, next: Option<&QuiverOrdering>) -> Result<(), OrderingError> {
    if next.is_none_or(|ordering| ordering.quiver.is_identical(quiver)) {
        Ok(())
    } else {
        Err(OrderingError::DifferentQuivers)
    }
}

/// Failure to construct or apply a QPA quiver ordering.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum OrderingError {
    /// Ordering data and paths must use one identical quiver family.
    #[error("the ordering and paths belong to different quivers")]
    DifferentQuivers,
    /// QPA requires stacked vector orderings to share one lexicographic table.
    #[error("conflicting lexicographic tables")]
    ConflictingLexicographicTables,
}

#[cfg(test)]
mod tests {
    use numerica::domains::integer::Integer;
    use serde_json::{Value, json};

    use crate::{ArrowIndex, ArrowSpec, Generator, Quiver, VertexIndex};

    use super::QuiverOrdering;

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn arrow(index: usize) -> ArrowIndex {
        ArrowIndex::new(index).unwrap()
    }

    fn quiver() -> Quiver {
        Quiver::from_indexed(
            2,
            vec![
                ArrowSpec::named(vertex(1), vertex(1), "a"),
                ArrowSpec::named(vertex(1), vertex(2), "b"),
            ],
        )
        .unwrap()
    }

    fn comparison_matrix(ordering: &QuiverOrdering, paths: &[crate::Path]) -> Value {
        json!(
            paths
                .iter()
                .map(|left| {
                    paths
                        .iter()
                        .map(|right| match ordering.compare(left, right).unwrap() {
                            std::cmp::Ordering::Less => -1,
                            std::cmp::Ordering::Equal => 0,
                            std::cmp::Ordering::Greater => 1,
                        })
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        )
    }

    #[test]
    fn default_ordering_matches_quiver_comparison() {
        let quiver = quiver();
        let ordering = QuiverOrdering::default_for(&quiver).unwrap();
        let a = quiver.arrow(arrow(1)).unwrap().path();
        let aa = a.multiply(&a).unwrap();
        assert!(ordering.less_than(&quiver.zero(), &a).unwrap());
        assert!(ordering.less_than(&a, &aa).unwrap());
        assert_eq!(
            ordering.compare(&a, &aa).unwrap(),
            quiver.compare_paths(&a, &aa).unwrap()
        );
    }

    #[test]
    fn lexicographic_preference_and_right_direction_are_source_faithful() {
        let quiver = quiver();
        let a = quiver.arrow(arrow(1)).unwrap();
        let b = quiver.arrow(arrow(2)).unwrap();
        let preferred = [Generator::Arrow(b.clone()), Generator::Arrow(a.clone())];
        let left = QuiverOrdering::left_lexicographic(&quiver, &preferred).unwrap();
        assert!(left.less_than(&b.path(), &a.path()).unwrap());

        let aa = a.path().multiply(&a.path()).unwrap();
        let ba = b.path().multiply(&quiver.zero()).unwrap();
        assert!(ba.is_zero());
        let right = QuiverOrdering::right_lexicographic(&quiver, &preferred).unwrap();
        assert!(right.less_than(&b.path(), &aa).unwrap());
    }

    #[test]
    fn reverse_preserves_zero_and_vertex_comparisons() {
        let quiver = quiver();
        let base = QuiverOrdering::default_for(&quiver).unwrap();
        let reverse = QuiverOrdering::reverse(&quiver, base).unwrap();
        let a = quiver.arrow(arrow(1)).unwrap().path();
        let aa = a.multiply(&a).unwrap();
        assert!(reverse.less_than(&quiver.zero(), &a).unwrap());
        assert!(reverse.less_than(&aa, &a).unwrap());
    }

    #[test]
    fn vector_and_weight_orderings_count_qpa_arrow_positions() {
        let quiver = quiver();
        let a = quiver.arrow(arrow(1)).unwrap();
        let b = quiver.arrow(arrow(2)).unwrap();
        let aa = a.path().multiply(&a.path()).unwrap();

        let vector = QuiverOrdering::left_vector(&quiver, &[], None).unwrap();
        assert!(vector.less_than(&b.path(), &aa).unwrap());

        let weight =
            QuiverOrdering::weight(&quiver, &[Integer::from(1), Integer::from(5)], None).unwrap();
        assert!(weight.less_than(&aa, &b.path()).unwrap());
    }

    #[test]
    fn pinned_qpa_ordering_matrices_match() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/orderings.schema-v1.json"
        ))
        .unwrap();
        let quiver = quiver();
        let a = quiver.arrow(arrow(1)).unwrap();
        let b = quiver.arrow(arrow(2)).unwrap();
        let aa = a.path().multiply(&a.path()).unwrap();
        let ab = a.path().multiply(&b.path()).unwrap();
        let aab = aa.multiply(&b.path()).unwrap();
        let paths = vec![
            quiver.zero(),
            quiver.vertex(vertex(1)).unwrap().path(),
            quiver.vertex(vertex(2)).unwrap().path(),
            a.path(),
            b.path(),
            aa,
            ab,
            aab,
        ];
        assert_eq!(
            json!(paths.iter().map(ToString::to_string).collect::<Vec<_>>()),
            expected["paths"]
        );

        let preferred = [Generator::Arrow(b), Generator::Arrow(a)];
        let default = QuiverOrdering::default_for(&quiver).unwrap();
        let orderings = [
            ("default", default.clone()),
            (
                "left_lex_ba",
                QuiverOrdering::left_lexicographic(&quiver, &preferred).unwrap(),
            ),
            (
                "right_lex_ba",
                QuiverOrdering::right_lexicographic(&quiver, &preferred).unwrap(),
            ),
            (
                "reverse_default",
                QuiverOrdering::reverse(&quiver, default).unwrap(),
            ),
            (
                "left_vector_ba",
                QuiverOrdering::left_vector(&quiver, &preferred, None).unwrap(),
            ),
            (
                "right_vector_ba",
                QuiverOrdering::right_vector(&quiver, &preferred, None).unwrap(),
            ),
            (
                "weight_1_5",
                QuiverOrdering::weight(&quiver, &[Integer::from(1), Integer::from(5)], None)
                    .unwrap(),
            ),
        ];
        for (name, ordering) in orderings {
            assert_eq!(
                comparison_matrix(&ordering, &paths),
                expected["matrices"][name],
                "QPA comparison matrix differs for {name}"
            );
        }
    }
}
