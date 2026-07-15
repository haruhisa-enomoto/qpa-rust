//! QPA's noncommutative Groebner-basis routines and basis objects.
//!
//! The algorithms here translate `lib/gbhighlevel.gi` and `lib/groebner.gi`
//! from QPA 1.37. They intentionally operate on QPA paths and QPA's monomial
//! ordering rather than substituting a generic polynomial implementation.

use numerica::domains::Field;
use thiserror::Error;

use crate::{
    ArrowIndex, Path, PathAlgebra, PathAlgebraElement, PathAlgebraError, PathEndpoint,
    QuiverStaticDictionary, StaticDictionaryError,
};

/// Whether QPA treats a basis as two-sided or right-sided.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum GroebnerBasisKind {
    /// Ordinary two-sided Gröbner basis.
    TwoSided,
    /// Right Gröbner basis, whose tips reduce only at a word prefix.
    Right,
}

/// QPA's default Gröbner-basis object representation.
#[derive(Clone, Debug)]
pub struct GroebnerBasis<F: Field> {
    algebra: PathAlgebra<F>,
    relations: Vec<PathAlgebraElement<F>>,
    kind: GroebnerBasisKind,
    complete: bool,
    tip_reduced: Option<bool>,
    completely_reduced: Option<bool>,
}

impl<F: Field> GroebnerBasis<F> {
    /// Translate `GroebnerBasis(I, relations)` for a complete two-sided basis.
    ///
    /// QPA sorts the supplied collection, marks completeness, and initially
    /// leaves its tip-reduced and completely-reduced cached properties unset.
    ///
    /// # Errors
    ///
    /// Returns an error for relations from another path algebra or a failed
    /// QPA comparison while sorting.
    pub fn from_complete_relations(
        algebra: &PathAlgebra<F>,
        relations: &[PathAlgebraElement<F>],
    ) -> Result<Self, GroebnerError> {
        ensure_elements(relations, algebra)?;
        let mut relations = relations.to_vec();
        sort_elements(&mut relations)?;
        Ok(Self {
            algebra: algebra.clone(),
            relations,
            kind: GroebnerBasisKind::TwoSided,
            complete: true,
            tip_reduced: None,
            completely_reduced: None,
        })
    }

    /// Run QPA's `HighLevelGroebnerBasis` implementation.
    ///
    /// # Errors
    ///
    /// Returns an error for relations from another algebra, relations outside
    /// the arrow ideal, or a failure in source-translated arithmetic.
    pub fn high_level(
        algebra: &PathAlgebra<F>,
        relations: &[PathAlgebraElement<F>],
    ) -> Result<Self, GroebnerError> {
        let relations = high_level_groebner_basis(relations, algebra)?;
        Self::from_complete_relations(algebra, &relations)
    }

    /// Return the containing path algebra.
    #[must_use]
    pub const fn algebra(&self) -> &PathAlgebra<F> {
        &self.algebra
    }

    /// Return QPA's final reduced relation list.
    #[must_use]
    pub fn relations(&self) -> &[PathAlgebraElement<F>] {
        &self.relations
    }

    /// Iterate over relations in QPA's stored enumerator order.
    pub fn iter(&self) -> std::slice::Iter<'_, PathAlgebraElement<F>> {
        self.relations.iter()
    }

    /// Return whether this is QPA's two-sided or right basis category.
    #[must_use]
    pub const fn kind(&self) -> GroebnerBasisKind {
        self.kind
    }

    /// Return QPA's known `IsCompleteGroebnerBasis` property.
    #[must_use]
    pub const fn is_complete(&self) -> bool {
        self.complete
    }

    /// Return the cached `IsTipReducedGroebnerBasis` property.
    ///
    /// `None` corresponds to GAP's property not yet being known.
    #[must_use]
    pub const fn tip_reduced_property(&self) -> Option<bool> {
        self.tip_reduced
    }

    /// Return the cached `IsCompletelyReducedGroebnerBasis` property.
    ///
    /// `None` corresponds to GAP's property not yet being known.
    #[must_use]
    pub const fn completely_reduced_property(&self) -> Option<bool> {
        self.completely_reduced
    }

    /// Translate QPA `TipReduceGroebnerBasis` for a complete two-sided basis.
    ///
    /// Rust returns a transformed clone instead of mutating the original GAP
    /// object. The pending-set and append-only pattern-slot behavior follows
    /// QPA's dynamic-suffix-tree method.
    ///
    /// # Errors
    ///
    /// Returns an error for an incomplete/right basis or arithmetic failure.
    pub fn tip_reduced_basis(&self) -> Result<Self, GroebnerError> {
        if self.tip_reduced == Some(true) {
            return Ok(self.clone());
        }
        if !self.complete {
            return Err(GroebnerError::IncompleteBasis);
        }
        if self.kind != GroebnerBasisKind::TwoSided {
            return Err(GroebnerError::UnsupportedRightTipReductionConstruction);
        }

        let mut pending = Vec::new();
        for relation in &self.relations {
            add_set(&mut pending, relation.clone())?;
        }
        let mut slots: Vec<Option<PathAlgebraElement<F>>> = Vec::new();

        while !pending.is_empty() {
            let mut relation = pending.remove(0);
            if relation.is_zero() {
                continue;
            }
            let relation_walk = tip_walk(&relation)?;
            if let Some((slot, start, end)) = first_tip_slot_occurrence(&slots, &relation_walk)? {
                let reducer = slots[slot]
                    .as_ref()
                    .ok_or(GroebnerError::InvalidBasisState)?;
                relation = reduce_relation_at(&relation, reducer, start, end)?;
                if !relation.is_zero() {
                    add_set(&mut pending, relation)?;
                }
                continue;
            }

            for slot in &mut slots {
                let Some(existing) = slot.as_ref() else {
                    continue;
                };
                let existing_walk = tip_walk(existing)?;
                if leftmost_occurrence(&existing_walk, &relation_walk).is_some() {
                    let displaced = slot.take().ok_or(GroebnerError::InvalidBasisState)?;
                    add_set(&mut pending, displaced)?;
                }
            }
            slots.push(Some(relation));
        }

        if slots.iter().any(Option::is_none) {
            return Err(GroebnerError::PinnedSparseTipReduction);
        }
        let mut relations = slots.into_iter().flatten().collect::<Vec<_>>();
        sort_elements(&mut relations)?;
        Ok(Self {
            algebra: self.algebra.clone(),
            relations,
            kind: self.kind,
            complete: self.complete,
            tip_reduced: Some(true),
            completely_reduced: self.completely_reduced,
        })
    }

    /// Translate QPA `TipReduce(GB, relation)`.
    ///
    /// A two-sided basis reduces any occurrence in the current tip; a right
    /// basis reduces only a prefix occurrence.
    ///
    /// # Errors
    ///
    /// Requires a complete basis whose tip-reduced property is known true.
    pub fn tip_reduce(
        &self,
        mut relation: PathAlgebraElement<F>,
    ) -> Result<PathAlgebraElement<F>, GroebnerError> {
        ensure_algebra(&relation, &self.algebra)?;
        self.require_tip_reduced()?;
        if relation.is_zero() {
            return Ok(relation);
        }
        let dictionary = self.static_dictionary_from_relations()?;

        loop {
            let tip = relation.leading_monomial();
            let matches = match self.kind {
                GroebnerBasisKind::TwoSided => dictionary.search(&tip)?,
                GroebnerBasisKind::Right => dictionary.prefix_search(&tip)?,
            };
            let Some(first_match) = matches.first() else {
                break;
            };
            let pattern_index = first_match
                .pattern_indices
                .first()
                .copied()
                .ok_or(GroebnerError::InvalidBasisState)?
                .checked_sub(1)
                .ok_or(GroebnerError::InvalidBasisState)?;
            let reducer = self
                .relations
                .get(pattern_index)
                .ok_or(GroebnerError::InvalidBasisState)?;
            let pattern_length = reducer
                .leading_monomial()
                .walk()
                .ok_or(GroebnerError::InvalidTipWalk)?
                .len();
            let end = first_match.end_position;
            let start = end
                .checked_sub(pattern_length)
                .ok_or(GroebnerError::InvalidBasisState)?;
            relation = reduce_relation_at(&relation, reducer, start, end)?;
            if relation.is_zero() {
                break;
            }
        }
        Ok(relation)
    }

    /// Translate QPA `CompletelyReduce(GB, relation)`.
    ///
    /// # Errors
    ///
    /// Requires a complete tip-reduced basis and a relation in its algebra.
    pub fn completely_reduce(
        &self,
        relation: PathAlgebraElement<F>,
    ) -> Result<PathAlgebraElement<F>, GroebnerError> {
        ensure_algebra(&relation, &self.algebra)?;
        self.require_tip_reduced()?;
        let mut remaining = relation;
        let mut reduced = self.algebra.zero();
        while !remaining.is_zero() {
            remaining = self.tip_reduce(remaining)?;
            if remaining.is_zero() {
                break;
            }
            let tip = remaining.leading_term();
            remaining = remaining.sub(&tip)?;
            reduced = reduced.add(&tip)?;
        }
        Ok(reduced)
    }

    /// Translate QPA `CompletelyReduceGroebnerBasis`.
    ///
    /// # Errors
    ///
    /// Returns an error for an incomplete basis or arithmetic failure.
    pub fn completely_reduced_basis(&self) -> Result<Self, GroebnerError> {
        if self.completely_reduced == Some(true) {
            return Ok(self.clone());
        }
        let mut basis = self.tip_reduced_basis()?;
        for index in 0..basis.relations.len() {
            let relation = basis.relations[index].clone();
            let tip = relation.leading_term();
            let tail = relation.sub(&tip)?;
            let reduced_tail = basis.completely_reduce(tail)?;
            basis.relations[index] = tip.add(&reduced_tail)?;
        }
        basis.completely_reduced = Some(true);
        basis.tip_reduced = Some(true);
        Ok(basis)
    }

    /// Translate `IsPrefixOfTipInTipIdeal`.
    ///
    /// # Errors
    ///
    /// Returns an error for a foreign element or dictionary failure.
    pub fn is_prefix_of_tip_in_tip_ideal(
        &self,
        relation: &PathAlgebraElement<F>,
    ) -> Result<bool, GroebnerError> {
        ensure_algebra(relation, &self.algebra)?;
        if relation.is_zero() {
            return Ok(true);
        }
        let reduced = self.tip_reduced_basis()?;
        let tip = relation.leading_monomial();
        let tip_length = tip.walk().ok_or(GroebnerError::InvalidTipWalk)?.len();
        Ok(reduced
            .static_dictionary_from_relations()?
            .search(&tip)?
            .first()
            .is_some_and(|matched| matched.end_position < tip_length))
    }

    /// Translate QPA `RightGroebnerBasis` from a finite-nontip complete basis.
    ///
    /// # Errors
    ///
    /// Returns an error unless this is a complete two-sided basis admitting
    /// finitely many nontips, or if source-translated arithmetic fails.
    pub fn right_basis(&self) -> Result<Self, GroebnerError> {
        if !self.complete {
            return Err(GroebnerError::IncompleteBasis);
        }
        if self.kind != GroebnerBasisKind::TwoSided {
            return Err(GroebnerError::NotTwoSidedBasis);
        }
        let reduced = self.completely_reduced_basis()?;
        let nontips = reduced.nontips()?.ok_or(GroebnerError::InfiniteNontips)?;
        let mut relations = Vec::new();
        for path in nontips {
            let prefix = reduced.algebra.element(path)?;
            for generator in &reduced.relations {
                let candidate = prefix.multiply(generator)?;
                if !reduced.is_prefix_of_tip_in_tip_ideal(&candidate)? {
                    add_set(&mut relations, candidate)?;
                }
            }
        }
        sort_elements(&mut relations)?;
        Ok(Self {
            algebra: self.algebra.clone(),
            relations,
            kind: GroebnerBasisKind::Right,
            complete: true,
            tip_reduced: Some(true),
            completely_reduced: None,
        })
    }

    /// Compute QPA's `RemainderOfDivision` by this basis.
    ///
    /// # Errors
    ///
    /// Propagates relation-family and arithmetic failures.
    pub fn remainder(
        &self,
        element: PathAlgebraElement<F>,
    ) -> Result<PathAlgebraElement<F>, GroebnerError> {
        remainder_of_division(element, &self.relations, &self.algebra)
    }

    /// Build QPA's `QuiverStaticDictionary` from leading monomials.
    ///
    /// # Errors
    ///
    /// Propagates static-dictionary validation failures.
    pub fn static_dictionary(&self) -> Result<QuiverStaticDictionary, GroebnerError> {
        let reduced = self.tip_reduced_basis()?;
        reduced.static_dictionary_from_relations()
    }

    fn static_dictionary_from_relations(&self) -> Result<QuiverStaticDictionary, GroebnerError> {
        let patterns = self
            .relations
            .iter()
            .map(PathAlgebraElement::leading_monomial)
            .collect::<Vec<_>>();
        Ok(QuiverStaticDictionary::new(
            self.algebra.quiver(),
            &patterns,
        )?)
    }

    fn require_tip_reduced(&self) -> Result<(), GroebnerError> {
        if !self.complete {
            Err(GroebnerError::IncompleteBasis)
        } else if self.tip_reduced != Some(true) {
            Err(GroebnerError::BasisNotTipReduced)
        } else {
            Ok(())
        }
    }

    /// Translate QPA's `AdmitsFinitelyManyNontips`.
    ///
    /// # Errors
    ///
    /// Propagates static-dictionary construction failures.
    pub fn admits_finitely_many_nontips(&self) -> Result<bool, GroebnerError> {
        Ok(self.static_dictionary()?.is_finite_difference())
    }

    /// Translate QPA's `NontipSize`; `None` represents `infinity`.
    ///
    /// # Errors
    ///
    /// Propagates static-dictionary construction failures.
    pub fn nontip_size(
        &self,
    ) -> Result<Option<numerica::domains::integer::Integer>, GroebnerError> {
        Ok(self.static_dictionary()?.difference_size())
    }

    /// Translate QPA's `Nontips`; `None` represents `fail` for an infinite
    /// nontip language.
    ///
    /// # Errors
    ///
    /// Propagates static-dictionary construction failures.
    pub fn nontips(&self) -> Result<Option<Vec<Path>>, GroebnerError> {
        Ok(self.static_dictionary()?.difference_words())
    }
}

impl<'a, F: Field> IntoIterator for &'a GroebnerBasis<F> {
    type Item = &'a PathAlgebraElement<F>;
    type IntoIter = std::slice::Iter<'a, PathAlgebraElement<F>>;

    fn into_iter(self) -> Self::IntoIter {
        self.relations.iter()
    }
}

impl<F: Field> AsRef<[PathAlgebraElement<F>]> for GroebnerBasis<F> {
    fn as_ref(&self) -> &[PathAlgebraElement<F>] {
        &self.relations
    }
}

/// Translate QPA's `HighLevelGroebnerBasis` overlap-completion algorithm.
///
/// This routine can be nonterminating when the ideal has no finite Groebner
/// basis for the chosen ordering, exactly as its QPA source can be.
///
/// # Errors
///
/// Returns an error for relations from another algebra, relations outside the
/// arrow ideal, or a failure in source-translated arithmetic.
pub fn high_level_groebner_basis<F: Field>(
    relations: &[PathAlgebraElement<F>],
    algebra: &PathAlgebra<F>,
) -> Result<Vec<PathAlgebraElement<F>>, GroebnerError> {
    if !qpa_in_arrow_ideal(relations, algebra)? {
        return Err(GroebnerError::OutsideArrowIdeal);
    }

    let mut pending = make_uniform_list(relations, algebra)?;
    let mut basis = Vec::new();

    while !pending.is_empty() {
        for relation in pending {
            let leading_coefficient = relation.leading_coefficient();
            basis.push(relation.divide_scalar(&leading_coefficient)?);
        }
        basis = reduced_list_qpa(&basis, algebra)?;

        let number_of_relations = basis.len();
        pending = Vec::new();
        for left_index in 0..number_of_relations {
            let left_walk = tip_walk(&basis[left_index])?;
            let left_length = left_walk.len();

            for right_index in 0..number_of_relations {
                let right_walk = tip_walk(&basis[right_index])?;
                let right_length = right_walk.len();
                let first_overlap = left_length.saturating_sub(right_length);

                for split in first_overlap..left_length {
                    let overlap_length = left_length - split;
                    if left_walk[split..] != right_walk[..overlap_length] {
                        continue;
                    }

                    let left_context = context_from_walk(algebra, &left_walk[..split])?;
                    let right_context = context_from_walk(algebra, &right_walk[overlap_length..])?;
                    let overlap = basis[left_index]
                        .multiply(&right_context)?
                        .sub(&left_context.multiply(&basis[right_index])?)?;
                    let remainder = remainder_of_division(overlap, &basis, algebra)?;
                    if !remainder.is_zero() {
                        add_set(&mut pending, remainder)?;
                    }
                }
            }
        }
    }

    let basis = tip_reduced_list(&basis, algebra)?;
    reduced_list_qpa(&basis, algebra)
}

/// Translate QPA's noncommutative `RemainderOfDivision`.
///
/// The first divisor whose tip occurs in the current tip is used, and the
/// leftmost occurrence is selected.
///
/// # Errors
///
/// Returns an error for elements from another algebra, zero divisors in the
/// divisor list, or source-translated arithmetic failures.
pub fn remainder_of_division<F: Field>(
    mut dividend: PathAlgebraElement<F>,
    divisors: &[PathAlgebraElement<F>],
    algebra: &PathAlgebra<F>,
) -> Result<PathAlgebraElement<F>, GroebnerError> {
    ensure_algebra(&dividend, algebra)?;
    for divisor in divisors {
        ensure_algebra(divisor, algebra)?;
        if divisor.is_zero() {
            return Err(GroebnerError::ZeroDivisor);
        }
    }

    let mut remainder = algebra.zero();
    while !dividend.is_zero() {
        let leading_term = dividend.leading_term();
        let dividend_walk = tip_walk(&leading_term)?;
        let mut reduction_applied = false;

        for divisor in divisors {
            let divisor_walk = tip_walk(divisor)?;
            let Some((start, end)) = leftmost_occurrence(&dividend_walk, &divisor_walk) else {
                continue;
            };

            let left = context_from_walk(algebra, &dividend_walk[..start])?;
            let right = context_from_walk(algebra, &dividend_walk[end + 1..])?;
            let coefficient = algebra.field().div(
                &leading_term.leading_coefficient(),
                &divisor.leading_coefficient(),
            );
            let reduction = left
                .multiply(divisor)?
                .multiply(&right)?
                .scale_left(&coefficient)?;
            dividend = dividend.sub(&reduction)?;
            reduction_applied = true;
            break;
        }

        if !reduction_applied {
            remainder = remainder.add(&leading_term)?;
            dividend = dividend.sub(&leading_term)?;
        }
    }

    Ok(remainder)
}

/// Translate QPA's `ReducedListQPA`.
///
/// # Errors
///
/// Returns an error for elements from another algebra or arithmetic failures.
pub fn reduced_list_qpa<F: Field>(
    elements: &[PathAlgebraElement<F>],
    algebra: &PathAlgebra<F>,
) -> Result<Vec<PathAlgebraElement<F>>, GroebnerError> {
    ensure_elements(elements, algebra)?;
    let mut result = elements
        .iter()
        .filter(|element| !element.is_zero())
        .cloned()
        .collect::<Vec<_>>();

    let mut index = result.len();
    while index > 0 {
        let current = index - 1;
        let others = result
            .iter()
            .enumerate()
            .filter(|(candidate, _)| *candidate != current)
            .map(|(_, element)| element.clone())
            .collect::<Vec<_>>();
        let reduced = remainder_of_division(result[current].clone(), &others, algebra)?;
        if reduced.is_zero() {
            result.remove(current);
        } else {
            result[current] = reduced;
        }
        index -= 1;
    }

    Ok(result)
}

/// Translate QPA's `TipReducedList`.
///
/// # Errors
///
/// Returns an error for elements from another algebra or invalid zero tips.
pub fn tip_reduced_list<F: Field>(
    elements: &[PathAlgebraElement<F>],
    algebra: &PathAlgebra<F>,
) -> Result<Vec<PathAlgebraElement<F>>, GroebnerError> {
    ensure_elements(elements, algebra)?;
    let mut result = Vec::new();
    for element in elements.iter().filter(|element| !element.is_zero()) {
        add_set(&mut result, element.clone())?;
    }

    let mut index = result.len();
    while index > 0 {
        let current = index - 1;
        let current_tip = tip_walk(&result[current])?;
        let mut contained = false;
        for prior in &result[..current] {
            if leftmost_occurrence(&current_tip, &tip_walk(prior)?).is_some() {
                contained = true;
                break;
            }
        }
        if contained {
            result.remove(current);
        }
        index -= 1;
    }
    Ok(result)
}

/// Translate QPA's `LeftmostOccurrence` using zero-based inclusive indices.
#[must_use]
pub fn leftmost_occurrence<T: PartialEq>(sequence: &[T], pattern: &[T]) -> Option<(usize, usize)> {
    if pattern.is_empty() || pattern.len() > sequence.len() {
        return None;
    }
    sequence
        .windows(pattern.len())
        .position(|window| window == pattern)
        .map(|start| (start, start + pattern.len() - 1))
}

/// Translate `QPA_InArrowIdeal` for a list of relations.
///
/// # Errors
///
/// Returns an error when a relation belongs to another algebra.
pub fn qpa_in_arrow_ideal<F: Field>(
    elements: &[PathAlgebraElement<F>],
    algebra: &PathAlgebra<F>,
) -> Result<bool, GroebnerError> {
    ensure_elements(elements, algebra)?;
    Ok(elements
        .iter()
        .flat_map(PathAlgebraElement::terms)
        .all(|term| !term.path().is_vertex()))
}

/// Translate QPA's `MakeUniform` for one path-algebra element.
///
/// # Errors
///
/// Returns an error for an element from another algebra or arithmetic failure.
pub fn make_uniform<F: Field>(
    element: &PathAlgebraElement<F>,
    algebra: &PathAlgebra<F>,
) -> Result<Vec<PathAlgebraElement<F>>, GroebnerError> {
    ensure_algebra(element, algebra)?;
    let mut sources = Vec::new();
    let mut targets = Vec::new();
    for term in element.terms() {
        if let PathEndpoint::Vertex(source) = term.path().source()
            && !sources.contains(&source)
        {
            sources.push(source);
        }
        if let PathEndpoint::Vertex(target) = term.path().target()
            && !targets.contains(&target)
        {
            targets.push(target);
        }
    }

    let mut uniform = Vec::new();
    for source in &sources {
        let source_idempotent = algebra.element(source.path())?;
        for target in &targets {
            let target_idempotent = algebra.element(target.path())?;
            let component = source_idempotent
                .multiply(element)?
                .multiply(&target_idempotent)?;
            if !component.is_zero() {
                add_set(&mut uniform, component)?;
            }
        }
    }
    Ok(uniform)
}

fn make_uniform_list<F: Field>(
    elements: &[PathAlgebraElement<F>],
    algebra: &PathAlgebra<F>,
) -> Result<Vec<PathAlgebraElement<F>>, GroebnerError> {
    let mut result = Vec::new();
    for element in elements {
        for component in make_uniform(element, algebra)? {
            add_set(&mut result, component)?;
        }
    }
    Ok(result)
}

fn context_from_walk<F: Field>(
    algebra: &PathAlgebra<F>,
    walk: &[ArrowIndex],
) -> Result<PathAlgebraElement<F>, GroebnerError> {
    if walk.is_empty() {
        return Ok(algebra.one());
    }
    let arrows = walk
        .iter()
        .map(|index| {
            algebra
                .quiver()
                .arrow(*index)
                .ok_or(GroebnerError::InvalidTipWalk)
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(algebra.element(Path::from_arrow_list(algebra.quiver(), &arrows)?)?)
}

fn tip_walk<F: Field>(element: &PathAlgebraElement<F>) -> Result<Vec<ArrowIndex>, GroebnerError> {
    if element.is_zero() {
        return Err(GroebnerError::ZeroDivisor);
    }
    element
        .leading_monomial()
        .walk()
        .map(<[ArrowIndex]>::to_vec)
        .ok_or(GroebnerError::InvalidTipWalk)
}

fn first_tip_slot_occurrence<F: Field>(
    slots: &[Option<PathAlgebraElement<F>>],
    sequence: &[ArrowIndex],
) -> Result<Option<(usize, usize, usize)>, GroebnerError> {
    for start in 0..sequence.len() {
        for (slot_index, slot) in slots.iter().enumerate() {
            let Some(relation) = slot else {
                continue;
            };
            let pattern = tip_walk(relation)?;
            let end = start + pattern.len();
            if !pattern.is_empty() && end <= sequence.len() && sequence[start..end] == pattern {
                return Ok(Some((slot_index, start, end)));
            }
        }
    }
    Ok(None)
}

fn reduce_relation_at<F: Field>(
    relation: &PathAlgebraElement<F>,
    reducer: &PathAlgebraElement<F>,
    start: usize,
    end: usize,
) -> Result<PathAlgebraElement<F>, GroebnerError> {
    ensure_algebra(reducer, relation.algebra())?;
    let tip = tip_walk(relation)?;
    if start > end || end > tip.len() {
        return Err(GroebnerError::InvalidBasisState);
    }
    let left = context_from_walk(relation.algebra(), &tip[..start])?;
    let right = context_from_walk(relation.algebra(), &tip[end..])?;
    let coefficient = relation.algebra().field().div(
        &relation.leading_coefficient(),
        &reducer.leading_coefficient(),
    );
    let reduction = left
        .multiply(reducer)?
        .multiply(&right)?
        .scale_left(&coefficient)?;
    Ok(relation.sub(&reduction)?)
}

fn ensure_algebra<F: Field>(
    element: &PathAlgebraElement<F>,
    algebra: &PathAlgebra<F>,
) -> Result<(), GroebnerError> {
    if element.algebra().is_identical(algebra) {
        Ok(())
    } else {
        Err(GroebnerError::DifferentAlgebras)
    }
}

fn ensure_elements<F: Field>(
    elements: &[PathAlgebraElement<F>],
    algebra: &PathAlgebra<F>,
) -> Result<(), GroebnerError> {
    for element in elements {
        ensure_algebra(element, algebra)?;
    }
    Ok(())
}

fn add_set<F: Field>(
    elements: &mut Vec<PathAlgebraElement<F>>,
    element: PathAlgebraElement<F>,
) -> Result<(), GroebnerError> {
    if elements.contains(&element) {
        return Ok(());
    }
    elements.push(element);
    let mut comparison_error = None;
    elements.sort_by(|left, right| match left.compare(right) {
        Ok(ordering) => ordering,
        Err(error) => {
            comparison_error = Some(error);
            std::cmp::Ordering::Equal
        }
    });
    comparison_error.map_or(Ok(()), |error| Err(error.into()))
}

fn sort_elements<F: Field>(elements: &mut [PathAlgebraElement<F>]) -> Result<(), GroebnerError> {
    let mut comparison_error = None;
    elements.sort_by(|left, right| match left.compare(right) {
        Ok(ordering) => ordering,
        Err(error) => {
            comparison_error = Some(error);
            std::cmp::Ordering::Equal
        }
    });
    comparison_error.map_or(Ok(()), |error| Err(error.into()))
}

/// Failure in QPA's source-translated Groebner-basis routines.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum GroebnerError {
    /// A relation belongs to another path-algebra family.
    #[error("a relation belongs to a different path algebra")]
    DifferentAlgebras,
    /// `HighLevelGroebnerBasis` accepts only relations in the arrow ideal.
    #[error("relations do not belong to the arrow ideal of the path algebra")]
    OutsideArrowIdeal,
    /// A zero relation cannot serve as a division rule.
    #[error("a Groebner divisor relation is zero")]
    ZeroDivisor,
    /// A tip was not a valid path walk.
    #[error("a relation tip does not have a valid arrow walk")]
    InvalidTipWalk,
    /// A method requiring a complete Gröbner basis received a partial basis.
    #[error("the Groebner basis is not known to be complete")]
    IncompleteBasis,
    /// `TipReduce` and `CompletelyReduce` require the cached tip-reduced
    /// property to be known true.
    #[error("the Groebner basis has not been tip reduced")]
    BasisNotTipReduced,
    /// An internal relation/dictionary index violated QPA's object invariant.
    #[error("the Groebner basis has an invalid relation/dictionary state")]
    InvalidBasisState,
    /// QPA's dynamic suffix-tree method leaves a sparse relation list in this
    /// case and its final `Sort` raises an error.
    #[error("QPA TipReduceGroebnerBasis produced the pinned sparse-list error")]
    PinnedSparseTipReduction,
    /// The requested construction applies only to two-sided bases.
    #[error("the operation requires a two-sided Groebner basis")]
    NotTwoSidedBasis,
    /// QPA's direct construction of a non-tip-reduced right basis is not
    /// reachable through `RightGroebnerBasis` and retains source defects.
    #[error("QPA does not provide a sound non-tip-reduced right-basis constructor")]
    UnsupportedRightTipReductionConstruction,
    /// `RightGroebnerBasis` requires a finite-dimensional quotient.
    #[error("the quotient must admit finitely many nontips")]
    InfiniteNontips,
    /// Source-translated path-algebra arithmetic failed.
    #[error(transparent)]
    PathAlgebra(#[from] PathAlgebraError),
    /// QPA path construction failed.
    #[error(transparent)]
    Path(#[from] crate::PathError),
    /// QPA static-dictionary construction failed.
    #[error(transparent)]
    StaticDictionary(#[from] StaticDictionaryError),
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

    use crate::{ArrowIndex, ArrowSpec, PathAlgebra, PathAlgebraElement, Quiver, VertexIndex};

    use super::{
        GroebnerBasis, GroebnerError, leftmost_occurrence, make_uniform, remainder_of_division,
    };

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn arrow(index: usize) -> ArrowIndex {
        ArrowIndex::new(index).unwrap()
    }

    fn small_integer(value: &Integer) -> i64 {
        value
            .to_string()
            .parse()
            .expect("the Groebner fixture uses only small integers")
    }

    fn element_json<F: Field>(
        element: &PathAlgebraElement<F>,
        coefficient_json: impl Fn(&F::Element) -> Value,
    ) -> Value {
        Value::Array(
            element
                .terms()
                .iter()
                .map(|term| {
                    json!({
                        "path": term.path().to_string(),
                        "coefficient": coefficient_json(term.coefficient())
                    })
                })
                .collect(),
        )
    }

    fn rational_json(element: &PathAlgebraElement<RationalField>) -> Value {
        element_json(element, |coefficient| {
            json!([
                small_integer(coefficient.numerator_ref()),
                small_integer(coefficient.denominator_ref())
            ])
        })
    }

    fn zp_json(element: &PathAlgebraElement<Zp>, field: &Zp) -> Value {
        element_json(element, |coefficient| {
            small_integer(&field.to_integer(coefficient)).into()
        })
    }

    fn z2_json(element: &PathAlgebraElement<Z2>) -> Value {
        element_json(element, |coefficient| {
            small_integer(&Z2.to_integer(coefficient)).into()
        })
    }

    fn external_element_json<F: Field>(
        element: &PathAlgebraElement<F>,
        coefficient_json: &impl Fn(&F::Element) -> Value,
    ) -> Value {
        Value::Array(
            element
                .terms()
                .iter()
                .map(|term| {
                    json!({
                        "path": term.path().external_representation(),
                        "coefficient": coefficient_json(term.coefficient())
                    })
                })
                .collect(),
        )
    }

    fn external_element_list_json<F: Field>(
        elements: &[PathAlgebraElement<F>],
        coefficient_json: &impl Fn(&F::Element) -> Value,
    ) -> Value {
        Value::Array(
            elements
                .iter()
                .map(|element| external_element_json(element, coefficient_json))
                .collect(),
        )
    }

    fn rational_coefficient_json(
        coefficient: &<RationalField as numerica::domains::Set>::Element,
    ) -> Value {
        json!([
            small_integer(coefficient.numerator_ref()),
            small_integer(coefficient.denominator_ref())
        ])
    }

    fn loop_quiver() -> Quiver {
        Quiver::from_indexed(
            1,
            vec![
                ArrowSpec::named(vertex(1), vertex(1), "x"),
                ArrowSpec::named(vertex(1), vertex(1), "y"),
            ],
        )
        .unwrap()
    }

    fn single_loop_quiver() -> Quiver {
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap()
    }

    #[test]
    fn leftmost_occurrence_translates_qpa_inclusive_result() {
        assert_eq!(leftmost_occurrence(&[1, 2, 1, 2], &[1, 2]), Some((0, 1)));
        assert_eq!(leftmost_occurrence(&[1, 2], &[2, 1]), None);
        assert_eq!(leftmost_occurrence::<u8>(&[], &[]), None);
    }

    #[test]
    fn make_uniform_splits_source_target_components() {
        let quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(2), vertex(3), "b"),
            ],
        )
        .unwrap();
        let algebra = PathAlgebra::new(Q, quiver).unwrap();
        let a = algebra
            .element(algebra.quiver().arrow(arrow(1)).unwrap().path())
            .unwrap();
        let b = algebra
            .element(algebra.quiver().arrow(arrow(2)).unwrap().path())
            .unwrap();
        let components = make_uniform(&a.add(&b).unwrap(), &algebra).unwrap();
        assert_eq!(components, vec![a, b]);
    }

    #[test]
    fn division_and_overlap_completion_follow_qpa_high_level_algorithm() {
        let quiver = loop_quiver();
        let algebra = PathAlgebra::new(Q, quiver).unwrap();
        let x = algebra
            .element(algebra.quiver().arrow(arrow(1)).unwrap().path())
            .unwrap();
        let y = algebra
            .element(algebra.quiver().arrow(arrow(2)).unwrap().path())
            .unwrap();
        let xx = x.multiply(&x).unwrap();
        let xy_minus_yx = x
            .multiply(&y)
            .unwrap()
            .sub(&y.multiply(&x).unwrap())
            .unwrap();
        let yy = y.multiply(&y).unwrap();
        let basis =
            GroebnerBasis::high_level(&algebra, &[xx.clone(), xy_minus_yx, yy.clone()]).unwrap();

        let word = y.multiply(&x).unwrap().multiply(&y).unwrap();
        assert!(basis.remainder(word).unwrap().is_zero());
        assert!(
            remainder_of_division(xx, basis.relations(), &algebra)
                .unwrap()
                .is_zero()
        );
        assert_eq!(basis.relations().len(), 3);
        assert_eq!(Q.one(), algebra.field().one());
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn pinned_qpa_high_level_groebner_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/groebner_high_level.schema-v1.json"
        ))
        .unwrap();
        assert_eq!(
            expected["oracle"],
            json!({
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1"
            })
        );

        let quiver = loop_quiver();
        let algebra = PathAlgebra::new(Q, quiver.clone()).unwrap();
        let x = algebra
            .element(quiver.arrow(arrow(1)).unwrap().path())
            .unwrap();
        let y = algebra
            .element(quiver.arrow(arrow(2)).unwrap().path())
            .unwrap();
        let rational_basis = GroebnerBasis::high_level(
            &algebra,
            &[
                x.pow(2).unwrap(),
                x.multiply(&y)
                    .unwrap()
                    .sub(&y.multiply(&x).unwrap())
                    .unwrap(),
                y.pow(2).unwrap(),
            ],
        )
        .unwrap();
        let yx = y.multiply(&x).unwrap();
        let yxy_plus_x = yx.multiply(&y).unwrap().add(&x).unwrap();
        assert_eq!(
            json!({
                "basis": rational_basis.relations().iter().map(rational_json).collect::<Vec<_>>(),
                "remainder_yx": rational_json(&rational_basis.remainder(yx).unwrap()),
                "remainder_yxy_plus_x": rational_json(&rational_basis.remainder(yxy_plus_x).unwrap())
            }),
            expected["rational"]
        );

        let field5 = Zp::new(5);
        let algebra5 = PathAlgebra::new(field5.clone(), quiver.clone()).unwrap();
        let x5 = algebra5
            .element(quiver.arrow(arrow(1)).unwrap().path())
            .unwrap();
        let y5 = algebra5
            .element(quiver.arrow(arrow(2)).unwrap().path())
            .unwrap();
        let relation5 = x5
            .multiply(&y5)
            .unwrap()
            .sub(
                &y5.multiply(&x5)
                    .unwrap()
                    .scale_left(&field5.to_element(2))
                    .unwrap(),
            )
            .unwrap();
        let basis5 = GroebnerBasis::high_level(
            &algebra5,
            &[x5.pow(2).unwrap(), relation5, y5.pow(2).unwrap()],
        )
        .unwrap();
        assert_eq!(
            json!({
                "basis": basis5.relations().iter().map(|element| zp_json(element, &field5)).collect::<Vec<_>>(),
                "remainder_yx": zp_json(&basis5.remainder(y5.multiply(&x5).unwrap()).unwrap(), &field5)
            }),
            expected["f5"]
        );

        let algebra2 = PathAlgebra::new(Z2, quiver).unwrap();
        let x2 = algebra2
            .element(algebra2.quiver().arrow(arrow(1)).unwrap().path())
            .unwrap();
        let y2 = algebra2
            .element(algebra2.quiver().arrow(arrow(2)).unwrap().path())
            .unwrap();
        let relation2 = x2
            .multiply(&y2)
            .unwrap()
            .add(&y2.multiply(&x2).unwrap())
            .unwrap();
        let basis2 = GroebnerBasis::high_level(
            &algebra2,
            &[x2.pow(2).unwrap(), relation2, y2.pow(2).unwrap()],
        )
        .unwrap();
        assert_eq!(
            json!({
                "basis": basis2.relations().iter().map(z2_json).collect::<Vec<_>>(),
                "remainder_yx": z2_json(&basis2.remainder(y2.multiply(&x2).unwrap()).unwrap())
            }),
            expected["f2"]
        );

        let uniform_quiver = Quiver::from_indexed(
            2,
            vec![
                ArrowSpec::named(vertex(1), vertex(1), "a"),
                ArrowSpec::named(vertex(1), vertex(2), "b"),
                ArrowSpec::named(vertex(2), vertex(1), "c"),
            ],
        )
        .unwrap();
        let uniform_algebra = PathAlgebra::new(Q, uniform_quiver).unwrap();
        let mut sum = uniform_algebra.zero();
        for index in 1..=3 {
            sum = sum
                .add(
                    &uniform_algebra
                        .element(uniform_algebra.quiver().arrow(arrow(index)).unwrap().path())
                        .unwrap(),
                )
                .unwrap();
        }
        assert_eq!(
            json!(
                make_uniform(&sum, &uniform_algebra)
                    .unwrap()
                    .iter()
                    .map(rational_json)
                    .collect::<Vec<_>>()
            ),
            expected["uniform_components"]
        );
    }

    #[test]
    #[allow(clippy::similar_names, clippy::too_many_lines)]
    fn pinned_qpa_groebner_object_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/groebner_objects.schema-v1.json"
        ))
        .unwrap();
        assert_eq!(
            expected["oracle"],
            json!({
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1"
            })
        );

        let quiver = single_loop_quiver();
        let algebra = PathAlgebra::new(Q, quiver.clone()).unwrap();
        let generators = algebra.generators();
        let vertex = generators[0].clone();
        let x = generators[1].clone();
        let x_squared = x.multiply(&x).unwrap();
        let redundant = [x.clone(), x_squared.clone()];
        let basis = GroebnerBasis::from_complete_relations(&algebra, &redundant).unwrap();
        let initial = basis.relations().to_vec();
        let initial_flags = [
            basis.tip_reduced_property().is_some(),
            basis.completely_reduced_property().is_some(),
            basis.is_complete(),
        ];
        let tip_reduced = basis.tip_reduced_basis().unwrap();
        let x_cubed = x_squared.multiply(&x).unwrap();
        let tip_reduction = tip_reduced
            .tip_reduce(x_cubed.add(&vertex).unwrap())
            .unwrap();
        let complete_reduction = tip_reduced
            .completely_reduce(x_cubed.add(&x).unwrap().add(&vertex).unwrap())
            .unwrap();
        let completely_reduced = tip_reduced.completely_reduced_basis().unwrap();
        let final_flags = [
            completely_reduced.tip_reduced_property() == Some(true),
            completely_reduced.completely_reduced_property() == Some(true),
            completely_reduced.is_complete(),
        ];

        let cubic_basis =
            GroebnerBasis::from_complete_relations(&algebra, std::slice::from_ref(&x_cubed))
                .unwrap();
        let nontips = cubic_basis.nontips().unwrap().unwrap();
        let nontip_size = cubic_basis.nontip_size().unwrap().unwrap();
        let prefixes = [
            cubic_basis
                .is_prefix_of_tip_in_tip_ideal(&x_squared)
                .unwrap(),
            cubic_basis.is_prefix_of_tip_in_tip_ideal(&x_cubed).unwrap(),
            cubic_basis
                .is_prefix_of_tip_in_tip_ideal(&x_cubed.multiply(&x).unwrap())
                .unwrap(),
        ];
        let right_basis = cubic_basis.right_basis().unwrap();
        let right_reduction = right_basis
            .tip_reduce(x_cubed.multiply(&x).unwrap().add(&x).unwrap())
            .unwrap();

        let field5 = Zp::new(5);
        let algebra5 = PathAlgebra::new(field5.clone(), quiver.clone()).unwrap();
        let generators5 = algebra5.generators();
        let vertex5 = generators5[0].clone();
        let x5 = generators5[1].clone();
        let x5_squared = x5.multiply(&x5).unwrap();
        let basis5 =
            GroebnerBasis::from_complete_relations(&algebra5, &[x5.clone(), x5_squared.clone()])
                .unwrap()
                .tip_reduced_basis()
                .unwrap();
        let tip_reduction5 = basis5
            .tip_reduce(
                x5_squared
                    .scale_left(&field5.to_element(3))
                    .unwrap()
                    .add(&vertex5)
                    .unwrap(),
            )
            .unwrap();
        let f5_coefficient_json = |coefficient: &<Zp as numerica::domains::Set>::Element| {
            small_integer(&field5.to_integer(coefficient)).into()
        };

        let algebra2 = PathAlgebra::new(Z2, quiver).unwrap();
        let generators2 = algebra2.generators();
        let vertex2 = generators2[0].clone();
        let x2 = generators2[1].clone();
        let x2_squared = x2.multiply(&x2).unwrap();
        let basis2 =
            GroebnerBasis::from_complete_relations(&algebra2, &[x2.clone(), x2_squared.clone()])
                .unwrap()
                .tip_reduced_basis()
                .unwrap();
        let tip_reduction2 = basis2
            .tip_reduce(x2_squared.add(&vertex2).unwrap())
            .unwrap();
        let f2_coefficient_json = |coefficient: &<Z2 as numerica::domains::Set>::Element| {
            small_integer(&Z2.to_integer(coefficient)).into()
        };

        assert_eq!(
            json!({
                "schema_version": 1,
                "kind": "groebner_objects",
                "oracle": expected["oracle"].clone(),
                "rational": {
                    "initial": external_element_list_json(&initial, &rational_coefficient_json),
                    "initial_flags": initial_flags,
                    "tip_reduced": external_element_list_json(
                        tip_reduced.relations(),
                        &rational_coefficient_json
                    ),
                    "tip_reduction": external_element_json(
                        &tip_reduction,
                        &rational_coefficient_json
                    ),
                    "complete_reduction": external_element_json(
                        &complete_reduction,
                        &rational_coefficient_json
                    ),
                    "completely_reduced": external_element_list_json(
                        completely_reduced.relations(),
                        &rational_coefficient_json
                    ),
                    "final_flags": final_flags,
                    "nontips": nontips
                        .iter()
                        .map(crate::Path::external_representation)
                        .collect::<Vec<_>>(),
                    "nontip_size": small_integer(&nontip_size),
                    "prefix": prefixes,
                    "right_relations": external_element_list_json(
                        right_basis.relations(),
                        &rational_coefficient_json
                    ),
                    "right_reduction": external_element_json(
                        &right_reduction,
                        &rational_coefficient_json
                    )
                },
                "f5": {
                    "tip_reduced": external_element_list_json(
                        basis5.relations(),
                        &f5_coefficient_json
                    ),
                    "tip_reduction": external_element_json(
                        &tip_reduction5,
                        &f5_coefficient_json
                    )
                },
                "f2": {
                    "tip_reduced": external_element_list_json(
                        basis2.relations(),
                        &f2_coefficient_json
                    ),
                    "tip_reduction": external_element_json(
                        &tip_reduction2,
                        &f2_coefficient_json
                    )
                }
            }),
            expected
        );

        let broken = x_squared.sub(&x).unwrap();
        let pinned_bug = GroebnerBasis::from_complete_relations(&algebra, &[broken, x_squared])
            .unwrap()
            .tip_reduced_basis();
        assert_eq!(
            pinned_bug.unwrap_err(),
            GroebnerError::PinnedSparseTipReduction
        );
    }
}
