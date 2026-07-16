//! QPA quiver orderings.
//!
//! The comparison routines translate the constructors and
//! `SetComparisonFunction` closures in QPA 1.37
//! `lib/groebnerordering.gi`. They compare QPA path external
//! representations; no generic Rust ordering is substituted.

use std::{cmp::Ordering, fmt, sync::Arc};

use numerica::domains::rational::Rational;
use thiserror::Error;

use crate::{Generator, Path, Quiver};

/// A source-translated QPA ordering tied to a quiver's generator positions.
#[derive(Clone, Debug)]
pub struct QuiverOrdering {
    quiver: Quiver,
    kind: Arc<OrderingKind>,
    lexicographic_table: Option<Arc<LexicographicTable>>,
    identity: Arc<()>,
}

/// The QPA category installed on a quiver-ordering object.
///
/// The two lexicographic variants also satisfy `IsLexicographicOrdering`, and
/// the two vector variants also satisfy `IsVectorOrdering`.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum QuiverOrderingKind {
    /// `IsLengthOrdering`.
    Length,
    /// `IsLeftLexicographicOrdering` and `IsLexicographicOrdering`.
    LeftLexicographic,
    /// `IsRightLexicographicOrdering` and `IsLexicographicOrdering`.
    RightLexicographic,
    /// `IsReverseOrdering`.
    Reverse,
    /// `IsLeftVectorOrdering` and `IsVectorOrdering`.
    LeftVector,
    /// `IsRightVectorOrdering` and `IsVectorOrdering`.
    RightVector,
    /// `IsWeightOrdering`.
    Weight,
    /// `IsBlockOrdering`.
    Block,
    /// `IsWreathOrdering`.
    Wreath,
}

impl QuiverOrderingKind {
    /// Return QPA's `IsLexicographicOrdering` category result.
    #[must_use]
    pub const fn is_lexicographic(self) -> bool {
        matches!(self, Self::LeftLexicographic | Self::RightLexicographic)
    }

    /// Return QPA's `IsVectorOrdering` category result.
    #[must_use]
    pub const fn is_vector(self) -> bool {
        matches!(self, Self::LeftVector | Self::RightVector)
    }
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
        weights: Arc<[Rational]>,
        next: Option<QuiverOrdering>,
    },
    Block {
        blocks: Arc<[usize]>,
        orders: Arc<[QuiverOrdering]>,
        next: Option<QuiverOrdering>,
    },
    Wreath {
        blocks: Arc<[usize]>,
        orders: Arc<[QuiverOrdering]>,
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

    pub(crate) fn default_for_with_identity(
        quiver: &Quiver,
        identity: Arc<()>,
    ) -> Result<Self, OrderingError> {
        let mut ordering = Self::default_for(quiver)?;
        ordering.identity = identity;
        Ok(ordering)
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
            identity: Arc::new(()),
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
            identity: Arc::new(()),
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
            identity: Arc::new(()),
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
            identity: Arc::new(()),
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
            identity: Arc::new(()),
        })
    }

    /// Translate QPA's `WeightOrdering(Q, weights, next)` for exact rational
    /// weights admitted by Numerica.
    ///
    /// Missing arrow weights default to zero and excess weights are ignored by
    /// comparisons, reproducing how QPA fills and subsequently indexes its
    /// weight list.
    ///
    /// # Errors
    ///
    /// Returns an error if `next` belongs to a different quiver.
    pub fn weight<Weight>(
        quiver: &Quiver,
        weights: &[Weight],
        next: Option<Self>,
    ) -> Result<Self, OrderingError>
    where
        Weight: Clone + Into<Rational>,
    {
        validate_next(quiver, next.as_ref())?;
        let mut stored = vec![Rational::from(0); quiver.number_of_arrows()];
        for (target, source) in stored.iter_mut().zip(weights) {
            *target = source.clone().into();
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
            identity: Arc::new(()),
        })
    }

    /// Translate QPA's `BlockOrdering(Q, blocks, next)`.
    ///
    /// Each pair is one ordering and its explicitly assigned generators.
    /// Unassigned generators belong to the last block, matching QPA.
    ///
    /// # Errors
    ///
    /// Returns an error for an empty block list, foreign data, or conflicting
    /// lexicographic tables among component/fallback orderings.
    pub fn block(
        quiver: &Quiver,
        blocks: Vec<(Self, Vec<Generator>)>,
        next: Option<Self>,
    ) -> Result<Self, OrderingError> {
        Self::block_or_wreath(quiver, blocks, next, false)
    }

    /// Translate QPA's `WreathOrdering(Q, blocks, next)`.
    ///
    /// # Errors
    ///
    /// Returns the same checked-construction errors as [`Self::block`].
    pub fn wreath(
        quiver: &Quiver,
        blocks: Vec<(Self, Vec<Generator>)>,
        next: Option<Self>,
    ) -> Result<Self, OrderingError> {
        Self::block_or_wreath(quiver, blocks, next, true)
    }

    fn block_or_wreath(
        quiver: &Quiver,
        specifications: Vec<(Self, Vec<Generator>)>,
        next: Option<Self>,
        wreath: bool,
    ) -> Result<Self, OrderingError> {
        if specifications.is_empty() {
            return Err(OrderingError::EmptyBlockList);
        }
        validate_next(quiver, next.as_ref())?;
        let number_of_orders = specifications.len();
        let mut assignments = vec![number_of_orders - 1; generator_count(quiver) + 1];
        let mut orders = Vec::with_capacity(number_of_orders);
        let mut lexicographic_table = None;

        for (block_index, (ordering, generators)) in specifications.into_iter().enumerate() {
            validate_next(quiver, Some(&ordering))?;
            merge_lexicographic_table(&mut lexicographic_table, &ordering)?;
            for generator in generators {
                if !generator.path().quiver().is_identical(quiver) {
                    return Err(OrderingError::DifferentQuivers);
                }
                assignments[generator_position(&generator)] = block_index;
            }
            orders.push(ordering);
        }
        if let Some(next) = &next {
            merge_lexicographic_table(&mut lexicographic_table, next)?;
        }

        let kind = if wreath {
            OrderingKind::Wreath {
                blocks: assignments.into(),
                orders: orders.into(),
                next,
            }
        } else {
            OrderingKind::Block {
                blocks: assignments.into(),
                orders: orders.into(),
                next,
            }
        };
        Ok(Self {
            quiver: quiver.clone(),
            kind: Arc::new(kind),
            lexicographic_table,
            identity: Arc::new(()),
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
            OrderingKind::Block { .. } => "block",
            OrderingKind::Wreath { .. } => "wreath",
        }
    }

    /// Return the specific QPA ordering category.
    #[must_use]
    pub fn kind(&self) -> QuiverOrderingKind {
        match &*self.kind {
            OrderingKind::Length { .. } => QuiverOrderingKind::Length,
            OrderingKind::LeftLexicographic { .. } => QuiverOrderingKind::LeftLexicographic,
            OrderingKind::RightLexicographic { .. } => QuiverOrderingKind::RightLexicographic,
            OrderingKind::Reverse { .. } => QuiverOrderingKind::Reverse,
            OrderingKind::LeftVector { .. } => QuiverOrderingKind::LeftVector,
            OrderingKind::RightVector { .. } => QuiverOrderingKind::RightVector,
            OrderingKind::Weight { .. } => QuiverOrderingKind::Weight,
            OrderingKind::Block { .. } => QuiverOrderingKind::Block,
            OrderingKind::Wreath { .. } => QuiverOrderingKind::Wreath,
        }
    }

    /// Return the containing quiver.
    #[must_use]
    pub const fn quiver(&self) -> &Quiver {
        &self.quiver
    }

    /// Return whether two handles denote the identical QPA ordering object.
    #[must_use]
    pub fn is_identical(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.identity, &other.identity)
    }

    /// Compare paths from QPA's fresh `OrderedBy(Q, O)` family.
    pub(crate) fn compare_on_equally_generated_quiver(
        &self,
        quiver: &Quiver,
        left: &Path,
        right: &Path,
    ) -> Result<Ordering, OrderingError> {
        if !left.quiver().is_identical(quiver) || !right.quiver().is_identical(quiver) {
            return Err(OrderingError::DifferentQuivers);
        }
        self.ensure_equal_generator_count(quiver)?;
        Ok(compare_kind(&self.kind, left, right))
    }

    /// Rebind an ordering to QPA's position-preserving `OrderedBy` copy.
    pub(crate) fn rebase_to(&self, quiver: &Quiver) -> Result<Self, OrderingError> {
        self.ensure_equal_generator_count(quiver)?;
        let kind = match &*self.kind {
            OrderingKind::Length { next } => OrderingKind::Length {
                next: next
                    .as_ref()
                    .map(|ordering| ordering.rebase_to(quiver))
                    .transpose()?,
            },
            OrderingKind::LeftLexicographic { table } => OrderingKind::LeftLexicographic {
                table: table.clone(),
            },
            OrderingKind::RightLexicographic { table } => OrderingKind::RightLexicographic {
                table: table.clone(),
            },
            OrderingKind::Reverse { next } => OrderingKind::Reverse {
                next: next.rebase_to(quiver)?,
            },
            OrderingKind::LeftVector { table, next } => OrderingKind::LeftVector {
                table: table.clone(),
                next: next
                    .as_ref()
                    .map(|ordering| ordering.rebase_to(quiver))
                    .transpose()?,
            },
            OrderingKind::RightVector { table, next } => OrderingKind::RightVector {
                table: table.clone(),
                next: next
                    .as_ref()
                    .map(|ordering| ordering.rebase_to(quiver))
                    .transpose()?,
            },
            OrderingKind::Weight { weights, next } => OrderingKind::Weight {
                weights: weights.clone(),
                next: next
                    .as_ref()
                    .map(|ordering| ordering.rebase_to(quiver))
                    .transpose()?,
            },
            OrderingKind::Block {
                blocks,
                orders,
                next,
            } => OrderingKind::Block {
                blocks: blocks.clone(),
                orders: orders
                    .iter()
                    .map(|ordering| ordering.rebase_to(quiver))
                    .collect::<Result<Vec<_>, _>>()?
                    .into(),
                next: next
                    .as_ref()
                    .map(|ordering| ordering.rebase_to(quiver))
                    .transpose()?,
            },
            OrderingKind::Wreath {
                blocks,
                orders,
                next,
            } => OrderingKind::Wreath {
                blocks: blocks.clone(),
                orders: orders
                    .iter()
                    .map(|ordering| ordering.rebase_to(quiver))
                    .collect::<Result<Vec<_>, _>>()?
                    .into(),
                next: next
                    .as_ref()
                    .map(|ordering| ordering.rebase_to(quiver))
                    .transpose()?,
            },
        };
        Ok(Self {
            quiver: quiver.clone(),
            kind: Arc::new(kind),
            lexicographic_table: self.lexicographic_table.clone(),
            identity: self.identity.clone(),
        })
    }

    fn ensure_equal_generator_count(&self, quiver: &Quiver) -> Result<(), OrderingError> {
        let ordering_count = generator_count(&self.quiver);
        let quiver_count = generator_count(quiver);
        if ordering_count == quiver_count {
            Ok(())
        } else {
            Err(OrderingError::OrderedQuiverGeneratorCountMismatch {
                ordering_count,
                quiver_count,
            })
        }
    }

    /// Return QPA's optional `NextOrdering` attribute.
    #[must_use]
    pub fn next_ordering(&self) -> Option<&Self> {
        match &*self.kind {
            OrderingKind::Length { next }
            | OrderingKind::LeftVector { next, .. }
            | OrderingKind::RightVector { next, .. }
            | OrderingKind::Weight { next, .. }
            | OrderingKind::Block { next, .. }
            | OrderingKind::Wreath { next, .. } => next.as_ref(),
            OrderingKind::Reverse { next } => Some(next),
            OrderingKind::LeftLexicographic { .. } | OrderingKind::RightLexicographic { .. } => {
                None
            }
        }
    }

    /// Return QPA's optional `LexicographicTable` in increasing order.
    #[must_use]
    pub fn lexicographic_table(&self) -> Option<Vec<Generator>> {
        self.lexicographic_table.as_ref().map(|table| {
            table
                .generators
                .iter()
                .map(|&position| generator_at_position(&self.quiver, position))
                .collect()
        })
    }

    /// Return whether two orderings carry the identical QPA
    /// `LexicographicTable` attribute object.
    #[must_use]
    pub fn shares_lexicographic_table_with(&self, other: &Self) -> bool {
        match (&self.lexicographic_table, &other.lexicographic_table) {
            (Some(left), Some(right)) => Arc::ptr_eq(left, right),
            _ => false,
        }
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

    /// Return the exact value produced by QPA's `ComparisonFunction`.
    ///
    /// Most constructors return `-1`, `0`, or `1`; weight orderings return the
    /// exact rational weight difference, including through composed block,
    /// wreath, reverse, and fallback orderings.
    ///
    /// # Errors
    ///
    /// Returns an error if either path belongs to another quiver.
    pub fn comparison_value(&self, left: &Path, right: &Path) -> Result<Rational, OrderingError> {
        if !left.quiver().is_identical(&self.quiver) || !right.quiver().is_identical(&self.quiver) {
            return Err(OrderingError::DifferentQuivers);
        }
        Ok(compare_word_value(
            &self.kind,
            &OrderingWord::from_path(left),
            &OrderingWord::from_path(right),
            self.quiver.number_of_vertices(),
        ))
    }

    /// Return QPA's `LessThanByOrdering` result.
    ///
    /// # Errors
    ///
    /// Returns an error if either path belongs to another quiver.
    pub fn less_than(&self, left: &Path, right: &Path) -> Result<bool, OrderingError> {
        self.compare(left, right).map(Ordering::is_lt)
    }

    /// Translate QPA's `CompareArrowLists`, including its formal,
    /// not-necessarily-composable `PathFromArrowList` behavior.
    ///
    /// # Errors
    ///
    /// Returns an error if either list contains a foreign generator.
    pub fn compare_arrow_lists(
        &self,
        left: &[Generator],
        right: &[Generator],
    ) -> Result<Ordering, OrderingError> {
        validate_subset(&self.quiver, left)?;
        validate_subset(&self.quiver, right)?;
        Ok(compare_word_kind(
            &self.kind,
            &OrderingWord::from_generators(left.iter().map(generator_position).collect::<Vec<_>>()),
            &OrderingWord::from_generators(
                right.iter().map(generator_position).collect::<Vec<_>>(),
            ),
            self.quiver.number_of_vertices(),
        ))
    }

    /// Return QPA's exact `CompareArrowLists` value.
    ///
    /// # Errors
    ///
    /// Returns an error if either list contains a foreign generator.
    pub fn compare_arrow_lists_value(
        &self,
        left: &[Generator],
        right: &[Generator],
    ) -> Result<Rational, OrderingError> {
        validate_subset(&self.quiver, left)?;
        validate_subset(&self.quiver, right)?;
        Ok(compare_word_value(
            &self.kind,
            &OrderingWord::from_generators(left.iter().map(generator_position).collect::<Vec<_>>()),
            &OrderingWord::from_generators(
                right.iter().map(generator_position).collect::<Vec<_>>(),
            ),
            self.quiver.number_of_vertices(),
        ))
    }

    /// Translate QPA's `IsWellOrdering` property.
    ///
    /// # Errors
    ///
    /// Returns a typed error for an observable failure in the pinned QPA
    /// implementation, notably its immutable adjacency-matrix mutation or a
    /// wreath ordering's missing subset method.
    pub fn is_well_ordering(&self) -> Result<bool, OrderingError> {
        match &*self.kind {
            OrderingKind::LeftLexicographic { table }
            | OrderingKind::RightLexicographic { table } => {
                lexicographic_is_well_ordering(&self.quiver, table)
            }
            OrderingKind::Reverse { next } => next.is_well_reversed_ordering(),
            _ => self.is_well_subset(&all_generators(&self.quiver)),
        }
    }

    /// Translate QPA's `IsWellReversedOrdering` property.
    ///
    /// # Errors
    ///
    /// Returns a typed error when the corresponding pinned QPA method fails.
    pub fn is_well_reversed_ordering(&self) -> Result<bool, OrderingError> {
        match &*self.kind {
            OrderingKind::Length { .. }
            | OrderingKind::LeftVector { .. }
            | OrderingKind::RightVector { .. } => Ok(self.quiver.is_acyclic()),
            OrderingKind::Reverse { next } => next.is_well_ordering(),
            OrderingKind::Weight { weights, next } => {
                if self.quiver.is_acyclic() {
                    return Ok(true);
                }
                let nonzero = self
                    .quiver
                    .arrows()
                    .enumerate()
                    .filter(|(index, _)| weights[*index] != 0)
                    .map(|(_, arrow)| Generator::Arrow(arrow))
                    .collect::<Vec<_>>();
                if !nonzero.is_empty() && !check_quiver_subset_for_cycles(&self.quiver, &nonzero)? {
                    return Ok(false);
                }
                next.as_ref().map_or(Ok(true), |ordering| {
                    ordering.is_well_reversed_subset(&nonzero)
                })
            }
            _ => self.is_well_reversed_subset(&all_generators(&self.quiver)),
        }
    }

    /// Translate QPA's `IsTotalOrdering` property.
    ///
    /// # Errors
    ///
    /// Returns a typed error when a transitive QPA property method fails.
    pub fn is_total_ordering(&self) -> Result<bool, OrderingError> {
        match &*self.kind {
            OrderingKind::Length { next }
            | OrderingKind::LeftVector { next, .. }
            | OrderingKind::RightVector { next, .. } => {
                if self.quiver.number_of_vertices() == 0 {
                    Ok(true)
                } else {
                    next.as_ref()
                        .map_or(Ok(false), QuiverOrdering::is_total_ordering)
                }
            }
            OrderingKind::Reverse { next } => next.is_total_ordering(),
            OrderingKind::Weight { next, .. } => {
                if self.quiver.number_of_vertices() == 0 {
                    Ok(true)
                } else {
                    next.as_ref()
                        .map_or(Ok(true), QuiverOrdering::is_total_ordering)
                }
            }
            OrderingKind::Block {
                blocks,
                orders,
                next,
            } => block_is_total_ordering(&self.quiver, blocks, orders, next.as_ref()),
            _ => self.is_total_subset(&all_generators(&self.quiver)),
        }
    }

    /// Translate QPA's `IsAdmissibleOrdering`: totality is tested first.
    ///
    /// # Errors
    ///
    /// Returns a typed error from either transitive property computation.
    pub fn is_admissible_ordering(&self) -> Result<bool, OrderingError> {
        if self.is_total_ordering()? {
            self.is_well_ordering()
        } else {
            Ok(false)
        }
    }

    /// Translate QPA's `IsWellSubset` operation.
    ///
    /// # Errors
    ///
    /// Returns an error for foreign generators or an observable pinned QPA
    /// method failure.
    pub fn is_well_subset(&self, subset: &[Generator]) -> Result<bool, OrderingError> {
        validate_subset(&self.quiver, subset)?;
        match &*self.kind {
            OrderingKind::Length { .. }
            | OrderingKind::LeftVector { .. }
            | OrderingKind::RightVector { .. }
            | OrderingKind::Weight { .. } => Ok(true),
            OrderingKind::LeftLexicographic { table }
            | OrderingKind::RightLexicographic { table } => Ok(lexicographic_subset_is_well(
                &self.quiver,
                table,
                subset,
                false,
            )),
            OrderingKind::Reverse { next } => next.is_well_reversed_subset(subset),
            OrderingKind::Block { blocks, orders, .. } => {
                block_subset_property(blocks, orders, subset, |ordering, generators| {
                    ordering.is_well_subset(generators)
                })
            }
            OrderingKind::Wreath { .. } => Err(OrderingError::MissingWreathSubsetMethod),
        }
    }

    /// Translate QPA's `IsWellReversedSubset` operation.
    ///
    /// # Errors
    ///
    /// Returns an error for foreign generators or an observable pinned QPA
    /// method failure.
    pub fn is_well_reversed_subset(&self, subset: &[Generator]) -> Result<bool, OrderingError> {
        validate_subset(&self.quiver, subset)?;
        match &*self.kind {
            OrderingKind::Length { .. }
            | OrderingKind::LeftVector { .. }
            | OrderingKind::RightVector { .. } => {
                check_quiver_subset_for_cycles(&self.quiver, subset)
            }
            OrderingKind::LeftLexicographic { table }
            | OrderingKind::RightLexicographic { table } => Ok(lexicographic_subset_is_well(
                &self.quiver,
                table,
                subset,
                true,
            )),
            OrderingKind::Reverse { next } => next.is_well_subset(subset),
            OrderingKind::Weight { weights, next } => {
                if self.quiver.is_acyclic() {
                    return Ok(true);
                }
                let vertices = self.quiver.number_of_vertices();
                let (nonzero, zero): (Vec<_>, Vec<_>) = subset
                    .iter()
                    .filter(|generator| matches!(generator, Generator::Arrow(_)))
                    .cloned()
                    .partition(|generator| {
                        weights[generator_position(generator) - vertices - 1] != 0
                    });
                if !nonzero.is_empty() && !check_quiver_subset_for_cycles(&self.quiver, &nonzero)? {
                    return Ok(false);
                }
                if !zero.is_empty() && !check_quiver_subset_for_cycles(&self.quiver, &zero)? {
                    return next
                        .as_ref()
                        .map_or(Ok(true), |ordering| ordering.is_well_reversed_subset(&zero));
                }
                Ok(true)
            }
            OrderingKind::Block { blocks, orders, .. } => {
                block_subset_property(blocks, orders, subset, |ordering, generators| {
                    ordering.is_well_reversed_subset(generators)
                })
            }
            OrderingKind::Wreath { .. } => Err(OrderingError::MissingWreathSubsetMethod),
        }
    }

    /// Translate QPA's `IsTotalSubset` operation.
    ///
    /// # Errors
    ///
    /// Returns an error for foreign generators or one of the pinned QPA
    /// failure branches documented by [`OrderingError`].
    pub fn is_total_subset(&self, subset: &[Generator]) -> Result<bool, OrderingError> {
        validate_subset(&self.quiver, subset)?;
        match &*self.kind {
            OrderingKind::Length { next }
            | OrderingKind::LeftVector { next, .. }
            | OrderingKind::RightVector { next, .. } => {
                if subset.is_empty() {
                    Ok(true)
                } else {
                    next.as_ref()
                        .map_or(Ok(false), |ordering| ordering.is_total_subset(subset))
                }
            }
            OrderingKind::LeftLexicographic { .. } | OrderingKind::RightLexicographic { .. } => {
                Ok(true)
            }
            OrderingKind::Reverse { next } => next.is_total_subset(subset),
            OrderingKind::Weight { weights, next } => {
                let singleton_nonzero_arrow = subset.len() == 1
                    && matches!(&subset[0], Generator::Arrow(_))
                    && weights
                        [generator_position(&subset[0]) - self.quiver.number_of_vertices() - 1]
                        != 0;
                if subset.is_empty() || singleton_nonzero_arrow {
                    Ok(true)
                } else if next.is_some() {
                    Err(OrderingError::PinnedWeightTotalSubsetNextCall)
                } else {
                    Ok(true)
                }
            }
            OrderingKind::Block {
                blocks,
                orders,
                next,
            } => block_is_total_subset(&self.quiver, blocks, orders, next.as_ref(), subset),
            OrderingKind::Wreath { .. } => Err(OrderingError::MissingWreathSubsetMethod),
        }
    }
}

impl fmt::Display for QuiverOrdering {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("<")?;
        let mut ordering = Some(self);
        while let Some(current) = ordering {
            formatter.write_str(current.name())?;
            formatter.write_str(" ")?;
            ordering = current.next_ordering();
        }
        formatter.write_str("ordering>")
    }
}

fn compare_kind(kind: &OrderingKind, left: &Path, right: &Path) -> Ordering {
    compare_word_kind(
        kind,
        &OrderingWord::from_path(left),
        &OrderingWord::from_path(right),
        left.quiver().number_of_vertices(),
    )
}

#[derive(Clone, Debug)]
struct OrderingWord {
    generators: Vec<usize>,
    zero: bool,
}

impl OrderingWord {
    fn from_path(path: &Path) -> Self {
        Self {
            generators: if path.is_zero() {
                Vec::new()
            } else {
                path.external_representation()
            },
            zero: path.is_zero(),
        }
    }

    fn from_generators(generators: Vec<usize>) -> Self {
        Self {
            zero: generators.is_empty(),
            generators,
        }
    }

    fn is_vertex(&self, vertices: usize) -> bool {
        !self.zero && self.generators.len() == 1 && self.generators[0] <= vertices
    }

    fn path_length(&self, vertices: usize) -> usize {
        if self.zero || self.is_vertex(vertices) {
            0
        } else {
            self.generators.len()
        }
    }
}

fn compare_word_kind(
    kind: &OrderingKind,
    left: &OrderingWord,
    right: &OrderingWord,
    vertices: usize,
) -> Ordering {
    compare_word_value(kind, left, right, vertices).cmp(&Rational::from(0))
}

fn compare_word_value(
    kind: &OrderingKind,
    left: &OrderingWord,
    right: &OrderingWord,
    vertices: usize,
) -> Rational {
    match kind {
        OrderingKind::Length { next } => {
            let comparison = left.path_length(vertices).cmp(&right.path_length(vertices));
            if comparison.is_eq() {
                compare_optional_words_value(next.as_ref(), left, right, vertices)
            } else {
                ordering_value(comparison)
            }
        }
        OrderingKind::LeftLexicographic { table } => {
            ordering_value(compare_lexicographic(table, left, right, false))
        }
        OrderingKind::RightLexicographic { table } => {
            ordering_value(compare_lexicographic(table, left, right, true))
        }
        OrderingKind::Reverse { next } => {
            let comparison = compare_word_value(&next.kind, left, right, vertices);
            if left.zero
                || right.zero
                || (left.is_vertex(vertices) && !right.is_vertex(vertices))
                || (!left.is_vertex(vertices) && right.is_vertex(vertices))
            {
                comparison
            } else {
                -comparison
            }
        }
        OrderingKind::LeftVector { table, next } => {
            compare_vector_value(table, next.as_ref(), left, right, false, vertices)
        }
        OrderingKind::RightVector { table, next } => {
            compare_vector_value(table, next.as_ref(), left, right, true, vertices)
        }
        OrderingKind::Weight { weights, next } => {
            compare_weight_value(weights, next.as_ref(), left, right, vertices)
        }
        OrderingKind::Block {
            blocks,
            orders,
            next,
        } => compare_block_value(blocks, orders, next.as_ref(), left, right, vertices),
        OrderingKind::Wreath {
            blocks,
            orders,
            next,
        } => compare_wreath_value(blocks, orders, next.as_ref(), left, right, vertices),
    }
}

fn compare_block_value(
    blocks: &[usize],
    orders: &[QuiverOrdering],
    next: Option<&QuiverOrdering>,
    left: &OrderingWord,
    right: &OrderingWord,
    vertices: usize,
) -> Rational {
    let mut left_blocks = vec![Vec::new(); orders.len()];
    let mut right_blocks = vec![Vec::new(); orders.len()];
    for &generator in &left.generators {
        left_blocks[blocks[generator]].push(generator);
    }
    for &generator in &right.generators {
        right_blocks[blocks[generator]].push(generator);
    }
    for index in (0..orders.len()).rev() {
        let comparison = compare_word_value(
            &orders[index].kind,
            &OrderingWord::from_generators(left_blocks[index].clone()),
            &OrderingWord::from_generators(right_blocks[index].clone()),
            vertices,
        );
        if comparison != 0 {
            return comparison;
        }
    }
    compare_optional_words_value(next, left, right, vertices)
}

fn compare_wreath_value(
    blocks: &[usize],
    orders: &[QuiverOrdering],
    next: Option<&QuiverOrdering>,
    left: &OrderingWord,
    right: &OrderingWord,
    vertices: usize,
) -> Rational {
    let common = left
        .generators
        .iter()
        .zip(&right.generators)
        .take_while(|(left, right)| left == right)
        .count();
    if common == left.generators.len() && common == right.generators.len() {
        return Rational::from(0);
    }

    let (mut left_groups, left_lower) = wreath_groups(&left.generators[common..], blocks);
    let (mut right_groups, right_lower) = wreath_groups(&right.generators[common..], blocks);
    while let (Some((left_block, left_word)), Some((right_block, right_word))) =
        (left_groups.last(), right_groups.last())
    {
        let (ordering_index, left_comparison, right_comparison, pop_left, pop_right) =
            match left_block.cmp(right_block) {
                Ordering::Equal => (
                    *left_block,
                    left_word.clone(),
                    right_word.clone(),
                    true,
                    true,
                ),
                Ordering::Greater => (*left_block, left_word.clone(), Vec::new(), true, false),
                Ordering::Less => (*right_block, Vec::new(), right_word.clone(), false, true),
            };
        let comparison = compare_word_value(
            &orders[ordering_index].kind,
            &OrderingWord::from_generators(left_comparison),
            &OrderingWord::from_generators(right_comparison),
            vertices,
        );
        if comparison != 0 {
            return comparison;
        }
        if pop_left {
            left_groups.pop();
        }
        if pop_right {
            right_groups.pop();
        }
    }
    while let Some((block, word)) = left_groups.pop() {
        let comparison = compare_word_value(
            &orders[block].kind,
            &OrderingWord::from_generators(word),
            &OrderingWord::from_generators(Vec::new()),
            vertices,
        );
        if comparison != 0 {
            return comparison;
        }
    }
    while let Some((block, word)) = right_groups.pop() {
        let comparison = compare_word_value(
            &orders[block].kind,
            &OrderingWord::from_generators(Vec::new()),
            &OrderingWord::from_generators(word),
            vertices,
        );
        if comparison != 0 {
            return comparison;
        }
    }

    let lower_comparison = compare_wreath_value(
        blocks,
        orders,
        next,
        &OrderingWord::from_generators(left_lower),
        &OrderingWord::from_generators(right_lower),
        vertices,
    );
    if lower_comparison == 0 {
        compare_optional_words_value(next, left, right, vertices)
    } else {
        lower_comparison
    }
}

fn wreath_groups(generators: &[usize], blocks: &[usize]) -> (Vec<(usize, Vec<usize>)>, Vec<usize>) {
    let mut groups: Vec<(usize, Vec<usize>)> = Vec::new();
    let mut lower = Vec::new();
    let mut maximum_block = None;
    for &generator in generators {
        let block = blocks[generator];
        match maximum_block {
            None => {
                groups.push((block, vec![generator]));
                maximum_block = Some(block);
            }
            Some(maximum) if block > maximum => {
                groups.push((block, vec![generator]));
                maximum_block = Some(block);
            }
            Some(maximum) if block == maximum => {
                groups
                    .last_mut()
                    .expect("a maximum block has a corresponding wreath group")
                    .1
                    .push(generator);
            }
            Some(_) => lower.push(generator),
        }
    }
    (groups, lower)
}

fn compare_lexicographic(
    table: &LexicographicTable,
    left: &OrderingWord,
    right: &OrderingWord,
    reverse_words: bool,
) -> Ordering {
    if left.zero {
        return if right.zero {
            Ordering::Equal
        } else {
            Ordering::Less
        };
    }
    if right.zero {
        return Ordering::Greater;
    }

    if reverse_words {
        compare_ranked_iterators(
            table,
            left.generators.iter().rev(),
            right.generators.iter().rev(),
        )
    } else {
        compare_ranked_iterators(table, left.generators.iter(), right.generators.iter())
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

fn compare_vector_value(
    table: &LexicographicTable,
    next: Option<&QuiverOrdering>,
    left: &OrderingWord,
    right: &OrderingWord,
    reverse_vector: bool,
    vertices: usize,
) -> Rational {
    let left_short = left.zero || left.is_vertex(vertices);
    let right_short = right.zero || right.is_vertex(vertices);
    match (left_short, right_short) {
        (true, true) => return compare_optional_words_value(next, left, right, vertices),
        (true, false) => return Rational::from(-1),
        (false, true) => return Rational::from(1),
        (false, false) => {}
    }

    let arrows = table.generators.len() - vertices;
    let counts = |word: &OrderingWord| {
        let mut result = vec![0_usize; arrows];
        for &generator in &word.generators {
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
    if comparison.is_eq() {
        compare_optional_words_value(next, left, right, vertices)
    } else {
        ordering_value(comparison)
    }
}

fn compare_weight_value(
    weights: &[Rational],
    next: Option<&QuiverOrdering>,
    left: &OrderingWord,
    right: &OrderingWord,
    vertices: usize,
) -> Rational {
    let left_short = left.zero || left.is_vertex(vertices);
    let right_short = right.zero || right.is_vertex(vertices);
    match (left_short, right_short) {
        (true, true) => return compare_optional_words_value(next, left, right, vertices),
        (true, false) => return Rational::from(-1),
        (false, true) => return Rational::from(1),
        (false, false) => {}
    }

    let total = |word: &OrderingWord| {
        word.generators
            .iter()
            .copied()
            .fold(Rational::from(0), |sum, generator| {
                sum + weights[generator - vertices - 1].clone()
            })
    };
    let difference = total(left) - total(right);
    if difference == 0 {
        compare_optional_words_value(next, left, right, vertices)
    } else {
        difference
    }
}

fn compare_optional_words_value(
    next: Option<&QuiverOrdering>,
    left: &OrderingWord,
    right: &OrderingWord,
    vertices: usize,
) -> Rational {
    next.map_or_else(
        || Rational::from(0),
        |ordering| compare_word_value(&ordering.kind, left, right, vertices),
    )
}

fn ordering_value(ordering: Ordering) -> Rational {
    Rational::from(match ordering {
        Ordering::Less => -1,
        Ordering::Equal => 0,
        Ordering::Greater => 1,
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

fn generator_at_position(quiver: &Quiver, position: usize) -> Generator {
    if position <= quiver.number_of_vertices() {
        Generator::Vertex(
            quiver
                .vertices()
                .nth(position - 1)
                .expect("lexicographic vertex position belongs to the quiver"),
        )
    } else {
        Generator::Arrow(
            quiver
                .arrows()
                .nth(position - quiver.number_of_vertices() - 1)
                .expect("lexicographic arrow position belongs to the quiver"),
        )
    }
}

fn default_generator_positions(quiver: &Quiver) -> Vec<usize> {
    (1..=generator_count(quiver)).collect()
}

fn generator_count(quiver: &Quiver) -> usize {
    quiver.number_of_vertices() + quiver.number_of_arrows()
}

fn merge_lexicographic_table(
    table: &mut Option<Arc<LexicographicTable>>,
    ordering: &QuiverOrdering,
) -> Result<(), OrderingError> {
    let Some(candidate) = &ordering.lexicographic_table else {
        return Ok(());
    };
    if table
        .as_ref()
        .is_some_and(|existing| **existing != **candidate)
    {
        Err(OrderingError::ConflictingLexicographicTables)
    } else {
        if table.is_none() {
            *table = Some(candidate.clone());
        }
        Ok(())
    }
}

fn validate_next(quiver: &Quiver, next: Option<&QuiverOrdering>) -> Result<(), OrderingError> {
    if next.is_none_or(|ordering| ordering.quiver.is_identical(quiver)) {
        Ok(())
    } else {
        Err(OrderingError::DifferentQuivers)
    }
}

fn all_generators(quiver: &Quiver) -> Vec<Generator> {
    quiver
        .vertices()
        .map(Generator::Vertex)
        .chain(quiver.arrows().map(Generator::Arrow))
        .collect()
}

fn validate_subset(quiver: &Quiver, subset: &[Generator]) -> Result<(), OrderingError> {
    if subset
        .iter()
        .all(|generator| generator.path().quiver().is_identical(quiver))
    {
        Ok(())
    } else {
        Err(OrderingError::DifferentQuivers)
    }
}

fn lexicographic_is_well_ordering(
    quiver: &Quiver,
    table: &LexicographicTable,
) -> Result<bool, OrderingError> {
    if quiver.is_acyclic() {
        return Ok(true);
    }

    // The second of QPA's two identically-filtered installed methods is the
    // active one. It chooses the first arrow in the lexicographic table.
    let vertices = quiver.number_of_vertices();
    let generator = table.generators[vertices];
    let arrow = quiver
        .arrows()
        .nth(generator - vertices - 1)
        .expect("a cyclic finite quiver has an arrow");
    if arrow.source() != arrow.target() || arrow.source().in_degree() > 1 {
        return Ok(false);
    }

    // groebnerordering.gi:327-330 attempts to mutate the immutable adjacency
    // matrix returned by QPA under the pinned GAP revision.
    Err(OrderingError::PinnedImmutableAdjacencyMutation)
}

#[derive(Clone, Copy, Debug)]
struct QuiverSubsetCycleNode {
    label: usize,
    active: bool,
    cycle_root: Option<usize>,
}

/// Translate QPA's global `CheckQuiverSubsetForCycles` DFS.
///
/// # Errors
///
/// Returns an error if the subset contains a generator from another quiver.
pub fn check_quiver_subset_for_cycles(
    quiver: &Quiver,
    subset: &[Generator],
) -> Result<bool, OrderingError> {
    validate_subset(quiver, subset)?;
    Ok(check_quiver_subset_for_cycles_unchecked(quiver, subset))
}

/// Translate QPA's public `CheckBlocksForCycles(O, block)` helper.
///
/// `blocks` is in generator order and contains QPA's one-based block number,
/// or zero for a generator outside the tested subset.
///
/// # Errors
///
/// On a cyclic quiver, returns a shape error unless there is at least one entry
/// for every quiver generator, or the typed pinned-QPA failure reached by the
/// block-cycle DFS. QPA ignores excess entries.
pub fn check_blocks_for_cycles(
    ordering: &QuiverOrdering,
    blocks: &[usize],
) -> Result<bool, OrderingError> {
    if ordering.quiver().is_acyclic() {
        return Ok(true);
    }
    let expected = generator_count(ordering.quiver());
    if blocks.len() < expected {
        return Err(OrderingError::BlockAssignmentLength {
            actual: blocks.len(),
            expected,
        });
    }
    let mut indexed = Vec::with_capacity(expected + 1);
    indexed.push(0);
    indexed.extend_from_slice(blocks);
    check_blocks_for_cycles_internal(ordering.quiver(), &indexed)
}

fn check_quiver_subset_for_cycles_unchecked(quiver: &Quiver, subset: &[Generator]) -> bool {
    if quiver.is_acyclic() {
        return true;
    }

    let mut included = vec![false; generator_count(quiver) + 1];
    let mut starts = vec![false; quiver.number_of_vertices()];
    for generator in subset {
        if let Generator::Arrow(arrow) = generator {
            included[arrow.generator_index()] = true;
            starts[arrow.source().index().zero_based()] = true;
        }
    }
    if !included.iter().any(|included| *included) {
        return true;
    }

    let mut states = vec![None; quiver.number_of_vertices()];
    let mut next_label = 1;
    for start in 0..starts.len() {
        if starts[start]
            && states[start].is_none()
            && !visit_quiver_subset_cycles(
                quiver,
                start,
                0,
                &included,
                &mut states,
                &mut next_label,
            )
        {
            return false;
        }
    }
    true
}

fn visit_quiver_subset_cycles(
    quiver: &Quiver,
    vertex: usize,
    last_subset_source_label: usize,
    included: &[bool],
    states: &mut [Option<QuiverSubsetCycleNode>],
    next_label: &mut usize,
) -> bool {
    states[vertex] = Some(QuiverSubsetCycleNode {
        label: *next_label,
        active: true,
        cycle_root: None,
    });
    *next_label += 1;

    let current_vertex = quiver
        .vertices()
        .nth(vertex)
        .expect("the DFS vertex index belongs to the quiver");
    for arrow in quiver
        .outgoing_arrows(&current_vertex)
        .expect("the DFS vertex belongs to the quiver")
    {
        let target = arrow.target().index().zero_based();
        if included[arrow.generator_index()] {
            if states[target].is_none() {
                let source_label = states[vertex].expect("current vertex is visited").label;
                if !visit_quiver_subset_cycles(
                    quiver,
                    target,
                    source_label,
                    included,
                    states,
                    next_label,
                ) {
                    return false;
                }
            } else {
                let target_state = states[target].expect("target is visited");
                let root_is_active = target_state
                    .cycle_root
                    .is_some_and(|root| states[root].is_some_and(|state| state.active));
                if target_state.active || root_is_active {
                    return false;
                }
            }
        } else if states[target].is_none() {
            if !visit_quiver_subset_cycles(
                quiver,
                target,
                last_subset_source_label,
                included,
                states,
                next_label,
            ) {
                return false;
            }
            let target_root = states[target].expect("target is visited").cycle_root;
            if let Some(root) = target_root
                && states[root].is_some_and(|state| state.active)
            {
                let replace = states[vertex]
                    .expect("current vertex is visited")
                    .cycle_root
                    .is_none_or(|current| {
                        states[current].expect("cycle root is visited").label
                            > states[root].expect("cycle root is visited").label
                    });
                if replace {
                    states[vertex]
                        .as_mut()
                        .expect("current vertex is visited")
                        .cycle_root = Some(root);
                }
            }
        } else {
            let target_state = states[target].expect("target is visited");
            let root_is_active = target_state
                .cycle_root
                .is_some_and(|root| states[root].is_some_and(|state| state.active));
            if target_state.active || root_is_active {
                let root = target_state.cycle_root.unwrap_or(target);
                if states[root].expect("cycle root is visited").label <= last_subset_source_label {
                    return false;
                }
                let replace = states[vertex]
                    .expect("current vertex is visited")
                    .cycle_root
                    .is_none_or(|current| {
                        states[current].expect("cycle root is visited").label
                            > states[root].expect("cycle root is visited").label
                    });
                if replace {
                    states[vertex]
                        .as_mut()
                        .expect("current vertex is visited")
                        .cycle_root = Some(root);
                }
            }
        }
    }
    states[vertex]
        .as_mut()
        .expect("current vertex is visited")
        .active = false;
    true
}

#[derive(Clone, Copy, Debug)]
struct LexicographicCycleNode {
    label: usize,
    status: i8,
    cycle_root: Option<usize>,
}

fn lexicographic_subset_is_well(
    quiver: &Quiver,
    table: &LexicographicTable,
    subset: &[Generator],
    reversed: bool,
) -> bool {
    if quiver.is_acyclic() {
        return true;
    }

    let mut included = vec![false; generator_count(quiver) + 1];
    let mut starts = vec![false; quiver.number_of_vertices()];
    let mut selected: Option<(usize, usize)> = None;
    for generator in subset {
        if let Generator::Arrow(arrow) = generator {
            let position = arrow.generator_index();
            included[position] = true;
            starts[arrow.source().index().zero_based()] = true;
            let rank = table.ranks[position];
            let replace = selected.is_none_or(|(_, selected_rank)| {
                if reversed {
                    rank < selected_rank
                } else {
                    rank > selected_rank
                }
            });
            if replace {
                selected = Some((arrow.index().zero_based(), rank));
            }
        }
    }
    let Some((selected_arrow, _)) = selected else {
        return true;
    };
    let selected_generator = quiver.number_of_vertices() + selected_arrow + 1;
    included[selected_generator] = false;

    let mut states = vec![None; quiver.number_of_vertices()];
    let mut next_label = 1;
    let selected_source = quiver
        .arrows()
        .nth(selected_arrow)
        .expect("the selected arrow belongs to the quiver")
        .source()
        .index()
        .zero_based();
    if !visit_lexicographic_subset(
        quiver,
        selected_source,
        0,
        selected_arrow,
        &included,
        &mut states,
        &mut next_label,
    ) {
        return false;
    }
    for start in 0..starts.len() {
        if starts[start]
            && states[start].is_none()
            && !visit_lexicographic_subset(
                quiver,
                start,
                0,
                selected_arrow,
                &included,
                &mut states,
                &mut next_label,
            )
        {
            return false;
        }
    }
    true
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn visit_lexicographic_subset(
    quiver: &Quiver,
    vertex: usize,
    last_parent_source_label: i64,
    selected_arrow: usize,
    included: &[bool],
    states: &mut [Option<LexicographicCycleNode>],
    next_label: &mut usize,
) -> bool {
    states[vertex] = Some(LexicographicCycleNode {
        label: *next_label,
        status: -1,
        cycle_root: None,
    });
    *next_label += 1;

    let current_vertex = quiver
        .vertices()
        .nth(vertex)
        .expect("the DFS vertex index belongs to the quiver");
    for arrow in quiver
        .outgoing_arrows(&current_vertex)
        .expect("the DFS vertex belongs to the quiver")
    {
        let target = arrow.target().index().zero_based();
        let arrow_index = arrow.index().zero_based();
        if included[arrow.generator_index()] {
            if states[target].is_none() {
                let source_label = states[vertex].expect("current vertex is visited").label;
                if !visit_lexicographic_subset(
                    quiver,
                    target,
                    i64::try_from(source_label).expect("DFS label fits in i64"),
                    selected_arrow,
                    included,
                    states,
                    next_label,
                ) {
                    return false;
                }
            } else {
                let target_state = states[target].expect("target is visited");
                let root_is_active = target_state.cycle_root.is_some_and(|root| {
                    states[root].is_some_and(|state| state.status.is_negative())
                });
                if (target_state.label == 1 && target_state.status == 2)
                    || target_state.status.is_negative()
                    || root_is_active
                {
                    return false;
                }
            }
        } else if states[target].is_none() {
            let next_parent = if arrow_index == selected_arrow {
                -1
            } else {
                last_parent_source_label
            };
            if !visit_lexicographic_subset(
                quiver,
                target,
                next_parent,
                selected_arrow,
                included,
                states,
                next_label,
            ) {
                return false;
            }
            let target_root = states[target].expect("target is visited").cycle_root;
            if let Some(root) = target_root
                && states[root].is_some_and(|state| state.status.is_negative())
            {
                let replace = states[vertex]
                    .expect("current vertex is visited")
                    .cycle_root
                    .is_none_or(|current| {
                        states[current].expect("cycle root is visited").label
                            > states[root].expect("cycle root is visited").label
                    });
                if replace {
                    states[vertex]
                        .as_mut()
                        .expect("current vertex is visited")
                        .cycle_root = Some(root);
                }
            }
        } else {
            let target_state = states[target].expect("target is visited");
            let root_is_active = target_state
                .cycle_root
                .is_some_and(|root| states[root].is_some_and(|state| state.status.is_negative()));
            if target_state.status.is_negative() || root_is_active {
                let root = target_state.cycle_root.unwrap_or(target);
                let root_label = states[root].expect("cycle root is visited").label;
                if i64::try_from(root_label).expect("DFS label fits in i64")
                    <= last_parent_source_label
                {
                    return false;
                }
                let replace = states[vertex]
                    .expect("current vertex is visited")
                    .cycle_root
                    .is_none_or(|current| {
                        states[current].expect("cycle root is visited").label > root_label
                    });
                if replace {
                    states[vertex]
                        .as_mut()
                        .expect("current vertex is visited")
                        .cycle_root = Some(root);
                }
                if last_parent_source_label < 0 || arrow_index == selected_arrow {
                    states[root].as_mut().expect("cycle root is visited").status = -2;
                }
            }
        }
    }
    let node = states[vertex].as_mut().expect("current vertex is visited");
    node.status = -node.status;
    true
}

fn block_subset_property(
    blocks: &[usize],
    orders: &[QuiverOrdering],
    subset: &[Generator],
    property: impl Fn(&QuiverOrdering, &[Generator]) -> Result<bool, OrderingError>,
) -> Result<bool, OrderingError> {
    let mut grouped = vec![Vec::new(); orders.len()];
    for generator in subset {
        grouped[blocks[generator_position(generator)]].push(generator.clone());
    }
    for (ordering, generators) in orders.iter().zip(&grouped) {
        if !property(ordering, generators)? {
            return Ok(false);
        }
    }
    Ok(true)
}

fn block_is_total_ordering(
    quiver: &Quiver,
    blocks: &[usize],
    orders: &[QuiverOrdering],
    next: Option<&QuiverOrdering>,
) -> Result<bool, OrderingError> {
    let generators = all_generators(quiver);
    let mut grouped = vec![Vec::new(); orders.len()];
    for generator in &generators {
        grouped[blocks[generator_position(generator)]].push(generator.clone());
    }
    for (ordering, generators) in orders.iter().zip(&grouped) {
        if !ordering.is_total_subset(generators)? {
            return next.map_or(Ok(false), QuiverOrdering::is_total_ordering);
        }
    }

    if generator_count(quiver) > 0 {
        let source_blocks = blocks.iter().map(|block| block + 1).collect::<Vec<_>>();
        if !check_blocks_for_cycles_internal(quiver, &source_blocks)? {
            return next.map_or(Ok(false), QuiverOrdering::is_total_ordering);
        }
    }
    Ok(true)
}

fn block_is_total_subset(
    quiver: &Quiver,
    blocks: &[usize],
    orders: &[QuiverOrdering],
    next: Option<&QuiverOrdering>,
    subset: &[Generator],
) -> Result<bool, OrderingError> {
    let mut grouped = vec![Vec::new(); orders.len()];
    let mut source_blocks = vec![0; blocks.len()];
    for generator in subset {
        if let Generator::Arrow(arrow) = generator {
            let source = arrow.source().index().get();
            source_blocks[source] = blocks[source] + 1;
        }
        let position = generator_position(generator);
        source_blocks[position] = blocks[position] + 1;
        grouped[blocks[position]].push(generator.clone());
    }
    for (ordering, generators) in orders.iter().zip(&grouped) {
        if !ordering.is_total_subset(generators)? {
            return next.map_or(Ok(false), |ordering| ordering.is_total_subset(subset));
        }
    }

    if generator_count(quiver) > 0 && !check_blocks_for_cycles_internal(quiver, &source_blocks)? {
        // This is intentionally `IsTotalOrdering`, not `IsTotalSubset`, in
        // QPA's source branch at groebnerordering.gi:1795.
        return next.map_or(Ok(false), QuiverOrdering::is_total_ordering);
    }
    Ok(true)
}

#[derive(Clone, Copy, Debug)]
struct BlockCycleNode {
    label: usize,
    active: bool,
    cycle_root: Option<usize>,
    cycle_block: usize,
}

fn check_blocks_for_cycles_internal(
    quiver: &Quiver,
    blocks: &[usize],
) -> Result<bool, OrderingError> {
    if quiver.is_acyclic() {
        return Ok(true);
    }

    let mut states = vec![None; quiver.number_of_vertices()];
    let mut next_label = 1;
    for vertex in 0..quiver.number_of_vertices() {
        if blocks[vertex + 1] > 0
            && states[vertex].is_none()
            && !visit_block_cycles(
                quiver,
                vertex,
                0,
                0,
                0,
                blocks,
                &mut states,
                &mut next_label,
            )?
        {
            return Ok(false);
        }
    }
    Ok(true)
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn visit_block_cycles(
    quiver: &Quiver,
    vertex: usize,
    transition_source_label: usize,
    current_block: usize,
    preceding_zero_source_label: usize,
    blocks: &[usize],
    states: &mut [Option<BlockCycleNode>],
    next_label: &mut usize,
) -> Result<bool, OrderingError> {
    states[vertex] = Some(BlockCycleNode {
        label: *next_label,
        active: true,
        cycle_root: None,
        cycle_block: 0,
    });
    *next_label += 1;

    let current_vertex = quiver
        .vertices()
        .nth(vertex)
        .expect("the DFS vertex index belongs to the quiver");
    for arrow in quiver
        .outgoing_arrows(&current_vertex)
        .expect("the DFS vertex belongs to the quiver")
    {
        let mut arrow_block = blocks[arrow.generator_index()];
        let (next_transition, next_zero_source) = if arrow_block > 0 {
            let transition = if current_block > 0 && arrow_block != current_block {
                preceding_zero_source_label
            } else {
                transition_source_label
            };
            (transition, *next_label)
        } else {
            arrow_block = current_block;
            (transition_source_label, preceding_zero_source_label)
        };
        let target = arrow.target().index().zero_based();
        if states[target].is_none() {
            if !visit_block_cycles(
                quiver,
                target,
                next_transition,
                arrow_block,
                next_zero_source,
                blocks,
                states,
                next_label,
            )? {
                return Ok(false);
            }
            let target_state = states[target].expect("target is visited");
            if let Some(root) = target_state.cycle_root
                && states[root].is_some_and(|state| state.active)
            {
                let replace = states[vertex]
                    .expect("current vertex is visited")
                    .cycle_root
                    .is_none_or(|current| {
                        states[current].expect("cycle root is visited").label
                            > states[root].expect("cycle root is visited").label
                    });
                if replace {
                    states[vertex]
                        .as_mut()
                        .expect("current vertex is visited")
                        .cycle_root = Some(root);
                }
                if target_state.cycle_block != 0 {
                    states[vertex]
                        .as_mut()
                        .expect("current vertex is visited")
                        .cycle_block = target_state.cycle_block;
                }
            }
        } else {
            let target_state = states[target].expect("target is visited");
            let root_is_active = target_state
                .cycle_root
                .is_some_and(|root| states[root].is_some_and(|state| state.active));
            if target_state.active || root_is_active {
                // QPA assigns `m := [k-2]` in this branch when the root is
                // bound, then indexes `n[m]`; GAP rejects the list position.
                if target_state.cycle_root.is_some() {
                    return Err(OrderingError::PinnedBlockCycleListIndex);
                }
                let root = target;
                if states[root].expect("cycle root is visited").label < next_transition
                    || (target_state.cycle_block > 0
                        && arrow_block > 0
                        && target_state.cycle_block != arrow_block)
                {
                    return Ok(false);
                }
                let replace = states[vertex]
                    .expect("current vertex is visited")
                    .cycle_root
                    .is_none_or(|current| {
                        states[current].expect("cycle root is visited").label
                            > states[root].expect("cycle root is visited").label
                    });
                if replace {
                    states[vertex]
                        .as_mut()
                        .expect("current vertex is visited")
                        .cycle_root = Some(root);
                }
            }
        }
    }
    states[vertex]
        .as_mut()
        .expect("current vertex is visited")
        .active = false;
    Ok(true)
}

/// Failure to construct or apply a QPA quiver ordering.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum OrderingError {
    /// Ordering data and paths must use one identical quiver family.
    #[error("the ordering and paths belong to different quivers")]
    DifferentQuivers,
    /// An `OrderedBy` copy has a generator count incompatible with its
    /// retained ordering. QPA accepts construction but cannot safely apply
    /// that ordering's position table to such paths.
    #[error("ordered quiver has {quiver_count} generators but its ordering has {ordering_count}")]
    OrderedQuiverGeneratorCountMismatch {
        /// Number of generator positions in the retained ordering.
        ordering_count: usize,
        /// Number of generator positions in the fresh quiver.
        quiver_count: usize,
    },
    /// QPA requires stacked vector orderings to share one lexicographic table.
    #[error("conflicting lexicographic tables")]
    ConflictingLexicographicTables,
    /// QPA rejects an empty list of block-ordering specifications.
    #[error("block and wreath ordering subset lists cannot be empty")]
    EmptyBlockList,
    /// `CheckBlocksForCycles` requires at least one block number per generator
    /// once its acyclic-quiver shortcut does not apply.
    #[error("block assignment has length {actual}, expected at least {expected}")]
    BlockAssignmentLength {
        /// Supplied number of block entries.
        actual: usize,
        /// Number of quiver generators.
        expected: usize,
    },
    /// QPA installs no subset-property methods for wreath orderings.
    #[error("QPA has no applicable ordering-subset property method for wreath orderings")]
    MissingWreathSubsetMethod,
    /// The active QPA lexicographic well-order method mutates an immutable GAP
    /// adjacency matrix on this branch.
    #[error("pinned QPA attempts to mutate an immutable adjacency matrix")]
    PinnedImmutableAdjacencyMutation,
    /// QPA's weight-order fallback calls `NextOrdering(O, L)` with two args.
    #[error("pinned QPA calls NextOrdering with two arguments in IsTotalSubset")]
    PinnedWeightTotalSubsetNextCall,
    /// QPA's block-cycle DFS uses a singleton list as a list position.
    #[error("pinned QPA uses a singleton list as an index in CheckBlocksForCycles")]
    PinnedBlockCycleListIndex,
}

#[cfg(test)]
mod tests {
    use numerica::domains::{integer::Integer, rational::Rational};
    use serde_json::{Value, json};

    use crate::{ArrowIndex, ArrowSpec, Generator, Path, PathEndpoint, Quiver, VertexIndex};

    use super::{
        OrderingError, QuiverOrdering, QuiverOrderingKind, all_generators, check_blocks_for_cycles,
        check_quiver_subset_for_cycles,
    };

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

    fn ordering_properties(ordering: &QuiverOrdering) -> Value {
        json!([
            ordering.is_well_ordering().unwrap(),
            ordering.is_well_reversed_ordering().unwrap(),
            ordering.is_total_ordering().unwrap(),
            ordering.is_admissible_ordering().unwrap(),
        ])
    }

    fn subset_properties(ordering: &QuiverOrdering, subset: &[Generator]) -> Value {
        json!([
            ordering.is_well_subset(subset).unwrap(),
            ordering.is_well_reversed_subset(subset).unwrap(),
            ordering.is_total_subset(subset).unwrap(),
        ])
    }

    fn endpoint_name(endpoint: PathEndpoint) -> String {
        match endpoint {
            PathEndpoint::Vertex(vertex) => vertex.name().to_owned(),
            PathEndpoint::Zero(path) => path.to_string(),
        }
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
    fn ordering_attributes_and_formal_arrow_lists_follow_qpa() {
        let quiver = quiver();
        let ordering = QuiverOrdering::default_for(&quiver).unwrap();
        assert_eq!(ordering.to_string(), "<length left lexicographic ordering>");
        assert_eq!(
            ordering.next_ordering().map(QuiverOrdering::name),
            Some("left lexicographic")
        );
        assert_eq!(
            ordering
                .lexicographic_table()
                .unwrap()
                .iter()
                .map(|generator| generator.path().to_string())
                .collect::<Vec<_>>(),
            ["v1", "v2", "a", "b"]
        );

        let a = Generator::Arrow(quiver.arrow(arrow(1)).unwrap());
        let b = Generator::Arrow(quiver.arrow(arrow(2)).unwrap());
        assert_eq!(
            ordering.compare_arrow_lists(&[b.clone(), a.clone()], &[a, b]),
            Ok(std::cmp::Ordering::Greater)
        );
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn pinned_qpa_ordering_surface_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/ordering_surface.schema-v1.json"
        ))
        .unwrap();
        let quiver = quiver();
        let a = quiver.arrow(arrow(1)).unwrap();
        let b = quiver.arrow(arrow(2)).unwrap();
        let default = QuiverOrdering::default_for(&quiver).unwrap();
        let length = QuiverOrdering::length(&quiver, None).unwrap();
        let left_lexicographic = QuiverOrdering::left_lexicographic(
            &quiver,
            &[Generator::Arrow(b.clone()), Generator::Arrow(a.clone())],
        )
        .unwrap();
        let right_lexicographic = QuiverOrdering::right_lexicographic(
            &quiver,
            &[Generator::Arrow(b.clone()), Generator::Arrow(a.clone())],
        )
        .unwrap();
        let reverse = QuiverOrdering::reverse(&quiver, default.clone()).unwrap();
        let left_vector = QuiverOrdering::left_vector(
            &quiver,
            &[Generator::Arrow(b.clone()), Generator::Arrow(a.clone())],
            None,
        )
        .unwrap();
        let right_vector = QuiverOrdering::right_vector(
            &quiver,
            &[Generator::Arrow(b.clone()), Generator::Arrow(a.clone())],
            None,
        )
        .unwrap();
        let weight =
            QuiverOrdering::weight(&quiver, &[Integer::from(1), Integer::from(5)], None).unwrap();
        let block = QuiverOrdering::block(
            &quiver,
            vec![(weight.clone(), all_generators(&quiver))],
            None,
        )
        .unwrap();
        let wreath = QuiverOrdering::wreath(
            &quiver,
            vec![(weight.clone(), all_generators(&quiver))],
            None,
        )
        .unwrap();
        let orderings = [
            length,
            left_lexicographic,
            right_lexicographic,
            reverse.clone(),
            left_vector,
            right_vector,
            weight.clone(),
            block.clone(),
            wreath.clone(),
        ];
        let categories = orderings
            .iter()
            .map(|ordering| {
                let kind = ordering.kind();
                json!([
                    true,
                    kind.is_lexicographic(),
                    kind == QuiverOrderingKind::Length,
                    kind == QuiverOrderingKind::Reverse,
                    kind.is_vector(),
                    kind == QuiverOrderingKind::Weight,
                    kind == QuiverOrderingKind::Block,
                    kind == QuiverOrderingKind::Wreath,
                    kind == QuiverOrderingKind::LeftVector,
                    kind == QuiverOrderingKind::RightVector,
                    kind == QuiverOrderingKind::LeftLexicographic,
                    kind == QuiverOrderingKind::RightLexicographic,
                ])
            })
            .collect::<Vec<_>>();
        assert_eq!(json!(categories), expected["categories"]);
        assert_eq!(
            json!(
                orderings
                    .iter()
                    .map(QuiverOrdering::name)
                    .collect::<Vec<_>>()
            ),
            expected["names"]
        );
        assert_eq!(
            json!(
                orderings
                    .iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>()
            ),
            expected["views"]
        );
        assert_eq!(
            json!({
                "default_view": default.to_string(),
                "default_name": default.name(),
                "default_has_next": default.next_ordering().is_some(),
                "default_next_identity": default.next_ordering().is_some_and(|next| next.is_identical(next)),
                "default_lex_table": default.lexicographic_table().unwrap().iter().map(|generator| generator.path().to_string()).collect::<Vec<_>>(),
                "reverse_next_identity": reverse.next_ordering().is_some_and(|next| next.is_identical(&default)),
                "reverse_lex_identity": reverse.shares_lexicographic_table_with(&default),
            }),
            expected["attributes"]
        );
        assert_eq!(
            json!({
                "weight_a_b": weight.comparison_value(&a.path(), &b.path()).unwrap().to_string(),
                "weight_aa_b": weight.compare_arrow_lists_value(
                    &[Generator::Arrow(a.clone()), Generator::Arrow(a.clone())],
                    &[Generator::Arrow(b.clone())],
                ).unwrap().to_string(),
                "block_a_b": block.comparison_value(&a.path(), &b.path()).unwrap().to_string(),
                "wreath_a_b": wreath.comparison_value(&a.path(), &b.path()).unwrap().to_string(),
            }),
            expected["comparison_values"]
        );

        let empty = Path::from_generator_list(&quiver, &[]).unwrap();
        let singleton = Path::from_generator_list(
            &quiver,
            &[Generator::Vertex(quiver.vertex(vertex(1)).unwrap())],
        )
        .unwrap();
        let formal_vertices = Path::from_generator_list(
            &quiver,
            &[
                Generator::Vertex(quiver.vertex(vertex(1)).unwrap()),
                Generator::Vertex(quiver.vertex(vertex(2)).unwrap()),
            ],
        )
        .unwrap();
        let formal_noncomposable = Path::from_generator_list(
            &quiver,
            &[Generator::Arrow(b.clone()), Generator::Arrow(a.clone())],
        )
        .unwrap();
        let path_record = |path: &Path| {
            json!({
                "view": path.to_string(),
                "extrep": path.external_representation(),
                "walk": path.generator_walk().unwrap().iter().map(|generator| generator.path().to_string()).collect::<Vec<_>>(),
                "source": endpoint_name(path.source()),
                "target": endpoint_name(path.target()),
                "length": path.length().unwrap(),
            })
        };
        assert_eq!(
            json!({
                "empty_is_zero": empty.is_zero(),
                "singleton_vertex_identity": singleton == quiver.vertex(vertex(1)).unwrap().path(),
                "formal_vertices": path_record(&formal_vertices),
                "formal_noncomposable": path_record(&formal_noncomposable),
            }),
            expected["path_from_arrow_list"]
        );

        let cycle_quiver = Quiver::from_indexed(
            2,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "x"),
                ArrowSpec::named(vertex(2), vertex(1), "y"),
            ],
        )
        .unwrap();
        let cycle_ordering = QuiverOrdering::length(&cycle_quiver, None).unwrap();
        assert_eq!(
            json!([
                check_blocks_for_cycles(&cycle_ordering, &[1, 1, 1, 1]).unwrap(),
                check_blocks_for_cycles(&cycle_ordering, &[1, 1, 1, 2]).unwrap(),
                check_blocks_for_cycles(&cycle_ordering, &[0, 0, 1, 2]).unwrap(),
            ]),
            expected["block_cycle_checks"]
        );
    }

    #[test]
    fn pinned_qpa_ordering_dispatch_failures_are_typed() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/ordering_dispatch.schema-v1.json"
        ))
        .unwrap();
        assert_eq!(expected["success"]["empty_block"], false);
        assert_eq!(expected["success"]["empty_wreath"], false);
        assert_eq!(expected["success"]["foreign_length"], true);

        let quiver = quiver();
        let a = quiver.arrow(arrow(1)).unwrap();
        let b = quiver.arrow(arrow(2)).unwrap();
        let default = QuiverOrdering::default_for(&quiver).unwrap();
        assert_eq!(
            QuiverOrdering::block(&quiver, Vec::new(), None).unwrap_err(),
            OrderingError::EmptyBlockList
        );
        assert_eq!(
            QuiverOrdering::wreath(&quiver, Vec::new(), None).unwrap_err(),
            OrderingError::EmptyBlockList
        );
        assert_eq!(
            QuiverOrdering::left_vector(
                &quiver,
                &[Generator::Arrow(b.clone()), Generator::Arrow(a.clone())],
                Some(default.clone()),
            )
            .unwrap_err(),
            OrderingError::ConflictingLexicographicTables
        );
        assert_eq!(
            check_blocks_for_cycles(&default, &[1]),
            Err(OrderingError::BlockAssignmentLength {
                actual: 1,
                expected: 4,
            })
        );

        let foreign =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "z")]).unwrap();
        let foreign_default = QuiverOrdering::default_for(&foreign).unwrap();
        assert_eq!(
            QuiverOrdering::length(&quiver, Some(foreign_default.clone())).unwrap_err(),
            OrderingError::DifferentQuivers
        );
        assert_eq!(
            QuiverOrdering::reverse(&quiver, foreign_default.clone()).unwrap_err(),
            OrderingError::DifferentQuivers
        );
        assert_eq!(
            QuiverOrdering::weight(
                &quiver,
                &[Integer::from(1), Integer::from(2)],
                Some(foreign_default),
            )
            .unwrap_err(),
            OrderingError::DifferentQuivers
        );
        assert_eq!(
            Path::from_generator_list(
                &quiver,
                &[
                    Generator::Vertex(quiver.vertex(vertex(1)).unwrap()),
                    Generator::Vertex(foreign.vertex(vertex(1)).unwrap()),
                ],
            ),
            Err(crate::PathError::DifferentQuivers)
        );
    }

    #[test]
    #[allow(clippy::too_many_lines)]
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

        let block_specifications = vec![
            (
                QuiverOrdering::default_for(&quiver).unwrap(),
                vec![Generator::Arrow(a.clone())],
            ),
            (
                QuiverOrdering::default_for(&quiver).unwrap(),
                vec![Generator::Arrow(b.clone())],
            ),
        ];
        let block = QuiverOrdering::block(&quiver, block_specifications.clone(), None).unwrap();
        let wreath = QuiverOrdering::wreath(&quiver, block_specifications, None).unwrap();
        let formal_left = [Generator::Arrow(b.clone()), Generator::Arrow(a.clone())];
        let formal_right = [Generator::Arrow(a.clone()), Generator::Arrow(b.clone())];
        let comparison_sign = |ordering: &QuiverOrdering| match ordering
            .compare_arrow_lists(&formal_left, &formal_right)
            .unwrap()
        {
            std::cmp::Ordering::Less => -1,
            std::cmp::Ordering::Equal => 0,
            std::cmp::Ordering::Greater => 1,
        };
        assert_eq!(
            json!({
                "default_ba_ab": comparison_sign(&QuiverOrdering::default_for(&quiver).unwrap()),
                "block_ba_ab": comparison_sign(&block),
                "wreath_ba_ab": comparison_sign(&wreath),
            }),
            expected["formal_arrow_lists"]
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
            (
                "weight_half_negative_third",
                QuiverOrdering::weight(
                    &quiver,
                    &[
                        Rational::from((1_i64, 2_i64)),
                        Rational::from((-1_i64, 3_i64)),
                    ],
                    None,
                )
                .unwrap(),
            ),
            ("block_a_b", block),
            ("wreath_a_b", wreath),
        ];
        for (name, ordering) in orderings {
            assert_eq!(
                comparison_matrix(&ordering, &paths),
                expected["matrices"][name],
                "QPA comparison matrix differs for {name}"
            );
        }

        let two_loop_quiver = Quiver::from_indexed(
            1,
            vec![
                ArrowSpec::named(vertex(1), vertex(1), "x"),
                ArrowSpec::named(vertex(1), vertex(1), "y"),
            ],
        )
        .unwrap();
        let x = two_loop_quiver.arrow(arrow(1)).unwrap();
        let y = two_loop_quiver.arrow(arrow(2)).unwrap();
        let xy = x.path().multiply(&y.path()).unwrap();
        let yx = y.path().multiply(&x.path()).unwrap();
        let xyx = xy.multiply(&x.path()).unwrap();
        let yxy = yx.multiply(&y.path()).unwrap();
        let xyxy = xyx.multiply(&y.path()).unwrap();
        let yxyx = yxy.multiply(&x.path()).unwrap();
        let paths = vec![
            two_loop_quiver.zero(),
            two_loop_quiver.vertex(vertex(1)).unwrap().path(),
            x.path(),
            y.path(),
            xy,
            yx,
            xyx,
            yxy,
            xyxy,
            yxyx,
        ];
        assert_eq!(
            json!(paths.iter().map(ToString::to_string).collect::<Vec<_>>()),
            expected["block_wreath"]["paths"]
        );
        let specifications = vec![
            (
                QuiverOrdering::default_for(&two_loop_quiver).unwrap(),
                vec![Generator::Arrow(x)],
            ),
            (
                QuiverOrdering::default_for(&two_loop_quiver).unwrap(),
                vec![Generator::Arrow(y)],
            ),
        ];
        for (name, ordering) in [
            (
                "block",
                QuiverOrdering::block(&two_loop_quiver, specifications.clone(), None).unwrap(),
            ),
            (
                "wreath",
                QuiverOrdering::wreath(&two_loop_quiver, specifications, None).unwrap(),
            ),
        ] {
            assert_eq!(
                comparison_matrix(&ordering, &paths),
                expected["block_wreath"][name],
                "QPA comparison matrix differs for {name}"
            );
        }
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn pinned_qpa_ordering_properties_match() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/ordering_properties.schema-v1.json"
        ))
        .unwrap();

        let acyclic =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let acyclic_default = QuiverOrdering::default_for(&acyclic).unwrap();
        let acyclic_length = QuiverOrdering::length(&acyclic, None).unwrap();
        let acyclic_lex = QuiverOrdering::left_lexicographic(&acyclic, &[]).unwrap();
        let acyclic_reverse = QuiverOrdering::reverse(&acyclic, acyclic_length.clone()).unwrap();
        let acyclic_vector = QuiverOrdering::left_vector(&acyclic, &[], None).unwrap();
        let acyclic_weight = QuiverOrdering::weight(&acyclic, &[Integer::from(0)], None).unwrap();
        let acyclic_block = QuiverOrdering::block(
            &acyclic,
            vec![(acyclic_default.clone(), all_generators(&acyclic))],
            None,
        )
        .unwrap();
        for (name, ordering) in [
            ("length", acyclic_length),
            ("default", acyclic_default),
            ("lex", acyclic_lex),
            ("reverse_length", acyclic_reverse),
            ("vector", acyclic_vector),
            ("weight", acyclic_weight),
            ("block", acyclic_block),
        ] {
            assert_eq!(
                ordering_properties(&ordering),
                expected["acyclic"][name],
                "QPA ordering properties differ for acyclic {name}"
            );
        }

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let loop_default = QuiverOrdering::default_for(&loop_quiver).unwrap();
        let loop_length = QuiverOrdering::length(&loop_quiver, None).unwrap();
        let loop_vector = QuiverOrdering::right_vector(&loop_quiver, &[], None).unwrap();
        let loop_weight_zero =
            QuiverOrdering::weight(&loop_quiver, &[Integer::from(0)], None).unwrap();
        let loop_weight_one =
            QuiverOrdering::weight(&loop_quiver, &[Integer::from(1)], None).unwrap();
        let loop_block = QuiverOrdering::block(
            &loop_quiver,
            vec![(loop_default.clone(), all_generators(&loop_quiver))],
            None,
        )
        .unwrap();
        for (name, ordering) in [
            ("length", loop_length),
            ("default", loop_default),
            ("vector", loop_vector),
            ("weight_zero", loop_weight_zero),
            ("weight_one", loop_weight_one),
            ("block", loop_block),
        ] {
            assert_eq!(
                ordering_properties(&ordering),
                expected["loop"][name],
                "QPA ordering properties differ for loop {name}"
            );
        }

        let subset_quiver = quiver();
        let a = subset_quiver.arrow(arrow(1)).unwrap();
        let b = subset_quiver.arrow(arrow(2)).unwrap();
        let v = subset_quiver.vertex(vertex(1)).unwrap();
        let subset_default = QuiverOrdering::default_for(&subset_quiver).unwrap();
        let subset_length = QuiverOrdering::length(&subset_quiver, None).unwrap();
        let subset_lex = QuiverOrdering::left_lexicographic(
            &subset_quiver,
            &[Generator::Arrow(b.clone()), Generator::Arrow(a.clone())],
        )
        .unwrap();
        let subset_weight_zero =
            QuiverOrdering::weight(&subset_quiver, &[Integer::from(0), Integer::from(0)], None)
                .unwrap();
        let subset_weight_mixed =
            QuiverOrdering::weight(&subset_quiver, &[Integer::from(1), Integer::from(0)], None)
                .unwrap();
        let subset_block = QuiverOrdering::block(
            &subset_quiver,
            vec![
                (subset_default.clone(), vec![Generator::Arrow(a.clone())]),
                (subset_default, vec![Generator::Arrow(b.clone())]),
            ],
            None,
        )
        .unwrap();
        let subsets = [
            vec![],
            vec![Generator::Vertex(v)],
            vec![Generator::Arrow(a.clone())],
            vec![Generator::Arrow(b.clone())],
            vec![Generator::Arrow(a), Generator::Arrow(b)],
            all_generators(&subset_quiver),
        ];
        assert_eq!(
            json!(
                subsets
                    .iter()
                    .map(|subset| {
                        check_quiver_subset_for_cycles(&subset_quiver, subset).unwrap()
                    })
                    .collect::<Vec<_>>()
            ),
            expected["subsets"]["cycle_checks"]
        );
        for (name, ordering) in [
            ("length", subset_length),
            ("lex_ba", subset_lex),
            ("weight_zero", subset_weight_zero),
            ("weight_mixed", subset_weight_mixed),
            ("block", subset_block),
        ] {
            assert_eq!(
                json!(
                    subsets
                        .iter()
                        .map(|subset| subset_properties(&ordering, subset))
                        .collect::<Vec<_>>()
                ),
                expected["subsets"][name],
                "QPA subset properties differ for {name}"
            );
        }

        let universe = [
            ArrowSpec::named(vertex(1), vertex(1), "a"),
            ArrowSpec::named(vertex(1), vertex(2), "b"),
            ArrowSpec::named(vertex(2), vertex(1), "c"),
            ArrowSpec::named(vertex(2), vertex(2), "d"),
        ];
        let mut exhaustive = Vec::new();
        for quiver_mask in 0_u8..16 {
            let specifications = universe
                .iter()
                .enumerate()
                .filter(|(index, _)| quiver_mask & (1 << index) != 0)
                .map(|(_, specification)| specification.clone())
                .collect();
            let quiver = Quiver::from_indexed(2, specifications).unwrap();
            let arrows = quiver.arrows().collect::<Vec<_>>();
            let lexicographic = QuiverOrdering::left_lexicographic(&quiver, &[]).unwrap();
            let preferred = arrows
                .iter()
                .rev()
                .cloned()
                .map(Generator::Arrow)
                .collect::<Vec<_>>();
            let reverse_lexicographic =
                QuiverOrdering::left_lexicographic(&quiver, &preferred).unwrap();
            let weights = (1..=arrows.len())
                .map(|index| Integer::from(index % 2))
                .collect::<Vec<_>>();
            let weight = QuiverOrdering::weight(&quiver, &weights, None).unwrap();
            for subset_mask in 0..(1_usize << arrows.len()) {
                let subset = arrows
                    .iter()
                    .enumerate()
                    .filter(|(index, _)| subset_mask & (1 << index) != 0)
                    .map(|(_, arrow)| Generator::Arrow(arrow.clone()))
                    .collect::<Vec<_>>();
                exhaustive.push(json!([
                    quiver_mask,
                    subset_mask,
                    check_quiver_subset_for_cycles(&quiver, &subset).unwrap(),
                    lexicographic.is_well_subset(&subset).unwrap(),
                    lexicographic.is_well_reversed_subset(&subset).unwrap(),
                    reverse_lexicographic.is_well_subset(&subset).unwrap(),
                    reverse_lexicographic
                        .is_well_reversed_subset(&subset)
                        .unwrap(),
                    weight.is_well_reversed_subset(&subset).unwrap(),
                ]));
            }
        }
        assert_eq!(json!(exhaustive), expected["exhaustive"]);
    }

    #[test]
    fn pinned_qpa_ordering_property_failures_are_typed() {
        let quiver = quiver();
        let lexicographic = QuiverOrdering::left_lexicographic(&quiver, &[]).unwrap();
        assert_eq!(
            lexicographic.is_well_ordering(),
            Err(OrderingError::PinnedImmutableAdjacencyMutation)
        );

        let default = QuiverOrdering::default_for(&quiver).unwrap();
        let a = quiver.arrow(arrow(1)).unwrap();
        let wreath = QuiverOrdering::wreath(
            &quiver,
            vec![(default.clone(), vec![Generator::Arrow(a.clone())])],
            None,
        )
        .unwrap();
        assert_eq!(
            wreath.is_total_ordering(),
            Err(OrderingError::MissingWreathSubsetMethod)
        );

        let weight = QuiverOrdering::weight(
            &quiver,
            &[Integer::from(0), Integer::from(0)],
            Some(default),
        )
        .unwrap();
        assert_eq!(
            weight.is_total_subset(&[Generator::Arrow(a)]),
            Err(OrderingError::PinnedWeightTotalSubsetNextCall)
        );

        let two_loops = Quiver::from_indexed(
            1,
            vec![
                ArrowSpec::named(vertex(1), vertex(1), "x"),
                ArrowSpec::named(vertex(1), vertex(1), "y"),
            ],
        )
        .unwrap();
        let x = two_loops.arrow(arrow(1)).unwrap();
        let y = two_loops.arrow(arrow(2)).unwrap();
        let two_loop_default = QuiverOrdering::default_for(&two_loops).unwrap();
        let block = QuiverOrdering::block(
            &two_loops,
            vec![
                (two_loop_default.clone(), vec![Generator::Arrow(x)]),
                (two_loop_default, vec![Generator::Arrow(y)]),
            ],
            None,
        )
        .unwrap();
        assert_eq!(
            block.is_total_ordering(),
            Err(OrderingError::PinnedBlockCycleListIndex)
        );
    }
}
