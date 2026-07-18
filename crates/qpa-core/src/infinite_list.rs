//! QPA's lazy, repeating, and doubly-infinite list machinery.
//!
//! The implementation follows `inflist.gd` and `inflist.gi`: a finite middle
//! part is joined to independently described positive and negative halves.
//! Half-lists may repeat, compute each value from the previous one, or compute
//! values directly from their positions.  In particular, `next/repeat` uses
//! QPA's same collision test and rewrites the owning infinite list when a
//! period is discovered.

use std::{
    fmt,
    sync::{Arc, Mutex, Weak},
};

use numerica::domains::integer::Integer;
use thiserror::Error;

/// QPA's two legal directions for a half-infinite list.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum InfiniteListDirection {
    /// Positions increase from the start position.
    Positive,
    /// Positions decrease from the start position.
    Negative,
}

impl InfiniteListDirection {
    const fn factor(self) -> i64 {
        match self {
            Self::Positive => 1,
            Self::Negative => -1,
        }
    }
}

impl fmt::Display for InfiniteListDirection {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Positive => "1",
            Self::Negative => "-1",
        })
    }
}

/// The normalized `InfListType` of a QPA half-infinite list.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum HalfInfiniteListType {
    /// A fixed finite list repeats forever.
    Repeat,
    /// Each value is computed from its predecessor.
    Next,
    /// Each value is computed directly from its position.
    Position,
}

/// A GAP integer enlarged by the two infinity values used by QPA.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub enum ExtendedInteger {
    /// GAP's `-infinity`.
    NegativeInfinity,
    /// An ordinary GAP integer.
    Finite(Integer),
    /// GAP's `infinity`.
    PositiveInfinity,
}

/// A function used by QPA's `"next"` and `"next/repeat"` half-lists.
pub type InfiniteListNextFunction<T, E> = Arc<dyn Fn(&T) -> Result<T, E> + Send + Sync + 'static>;

/// A function used by QPA's `"pos"` half-lists.
pub type InfiniteListPositionFunction<T, E> =
    Arc<dyn Fn(&Integer) -> Result<T, E> + Send + Sync + 'static>;

/// QPA's callback invoked after a lazy value is stored.
pub type InfiniteListNewValueCallback<E> = Arc<
    dyn Fn(&Integer, InfiniteListDirection, HalfInfiniteListType) -> Result<(), E>
        + Send
        + Sync
        + 'static,
>;

/// QPA's callback used to replace a detected `"next/repeat"` tail and period.
pub type InfiniteListRepeatifyCallback<T, E> =
    Arc<dyn Fn(Vec<T>, HalfInfiniteList<T, E>) -> Result<(), E> + Send + Sync + 'static>;

type InfiniteListMapFunction<T, U, E> = Arc<dyn Fn(&T) -> Result<U, E> + Send + Sync + 'static>;

/// The typed counterpart of QPA's `typeWithArgs` list.
pub enum HalfInfiniteListSpecification<T, E> {
    /// `[` `"repeat"`, repeating-list `]`.
    Repeat(Vec<T>),
    /// `[` `"next"`, function, initial-value `]`, optionally using QPA's
    /// `"next/repeat"` collision conversion.
    Next {
        /// Function applied successively, beginning with `initial_value`.
        function: InfiniteListNextFunction<T, E>,
        /// Value passed to the function to obtain the first stored entry.
        initial_value: T,
        /// Whether QPA's automatic repeating-list detection is enabled.
        repeatify: bool,
    },
    /// `[` `"pos"`, function, store-values `]`.
    Position {
        /// Function evaluated on an arbitrary GAP integer position.
        function: InfiniteListPositionFunction<T, E>,
        /// Whether values are computed consecutively and retained.
        store_values: bool,
    },
}

impl<T: Clone, E> Clone for HalfInfiniteListSpecification<T, E> {
    fn clone(&self) -> Self {
        match self {
            Self::Repeat(values) => Self::Repeat(values.clone()),
            Self::Next {
                function,
                initial_value,
                repeatify,
            } => Self::Next {
                function: Arc::clone(function),
                initial_value: initial_value.clone(),
                repeatify: *repeatify,
            },
            Self::Position {
                function,
                store_values,
            } => Self::Position {
                function: Arc::clone(function),
                store_values: *store_values,
            },
        }
    }
}

/// The callable value returned by QPA's `ElementFunction` accessor.
pub enum HalfInfiniteListElementFunction<T, E> {
    /// A predecessor-based element function.
    Next(InfiniteListNextFunction<T, E>),
    /// A position-based element function.
    Position(InfiniteListPositionFunction<T, E>),
}

impl<T, E> Clone for HalfInfiniteListElementFunction<T, E> {
    fn clone(&self) -> Self {
        match self {
            Self::Next(function) => Self::Next(Arc::clone(function)),
            Self::Position(function) => Self::Position(Arc::clone(function)),
        }
    }
}

enum HalfInfiniteListKind<T, E> {
    Repeat(Vec<T>),
    Next {
        function: InfiniteListNextFunction<T, E>,
        initial_value: T,
        repeatify: bool,
    },
    Position {
        function: InfiniteListPositionFunction<T, E>,
        store_values: bool,
    },
}

struct HalfInfiniteListData<T, E> {
    values: Vec<T>,
    kind: HalfInfiniteListKind<T, E>,
    callback: Option<InfiniteListNewValueCallback<E>>,
    repeatify_callback: Option<InfiniteListRepeatifyCallback<T, E>>,
}

/// QPA's `IsHalfInfList` object.
pub struct HalfInfiniteList<T, E = InfiniteListError> {
    start: Integer,
    direction: InfiniteListDirection,
    data: Arc<Mutex<HalfInfiniteListData<T, E>>>,
}

impl<T, E> Clone for HalfInfiniteList<T, E> {
    fn clone(&self) -> Self {
        Self {
            start: self.start.clone(),
            direction: self.direction,
            data: Arc::clone(&self.data),
        }
    }
}

impl<T, E> fmt::Debug for HalfInfiniteList<T, E> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("HalfInfiniteList")
            .field("start", &self.start)
            .field("direction", &self.direction)
            .finish_non_exhaustive()
    }
}

impl<T, E> HalfInfiniteList<T, E>
where
    T: Clone + PartialEq + Send + Sync + 'static,
    E: From<InfiniteListError> + Send + Sync + 'static,
{
    /// Translate `MakeHalfInfList` from a typed specification.
    #[must_use]
    pub fn new(
        start: Integer,
        direction: InfiniteListDirection,
        specification: HalfInfiniteListSpecification<T, E>,
        callback: Option<InfiniteListNewValueCallback<E>>,
        repeatify_callback: Option<InfiniteListRepeatifyCallback<T, E>>,
    ) -> Self {
        Self::new_with_repeatify(
            start,
            direction,
            specification,
            callback,
            repeatify_callback,
        )
    }

    fn new_with_repeatify(
        start: Integer,
        direction: InfiniteListDirection,
        specification: HalfInfiniteListSpecification<T, E>,
        callback: Option<InfiniteListNewValueCallback<E>>,
        repeatify_callback: Option<InfiniteListRepeatifyCallback<T, E>>,
    ) -> Self {
        let kind = match specification {
            HalfInfiniteListSpecification::Repeat(values) => HalfInfiniteListKind::Repeat(values),
            HalfInfiniteListSpecification::Next {
                function,
                initial_value,
                repeatify,
            } => HalfInfiniteListKind::Next {
                function,
                initial_value,
                repeatify,
            },
            HalfInfiniteListSpecification::Position {
                function,
                store_values,
            } => HalfInfiniteListKind::Position {
                function,
                store_values,
            },
        };
        Self {
            start,
            direction,
            data: Arc::new(Mutex::new(HalfInfiniteListData {
                values: Vec::new(),
                kind,
                callback,
                repeatify_callback,
            })),
        }
    }

    /// Construct QPA's repeating half-list.
    #[must_use]
    pub fn repeating(start: Integer, direction: InfiniteListDirection, values: Vec<T>) -> Self {
        Self::new(
            start,
            direction,
            HalfInfiniteListSpecification::Repeat(values),
            None,
            None,
        )
    }

    /// Construct QPA's `"next"` half-list.
    #[must_use]
    pub fn next(
        start: Integer,
        direction: InfiniteListDirection,
        initial_value: T,
        function: impl Fn(&T) -> Result<T, E> + Send + Sync + 'static,
        callback: Option<InfiniteListNewValueCallback<E>>,
    ) -> Self {
        Self::new(
            start,
            direction,
            HalfInfiniteListSpecification::Next {
                function: Arc::new(function),
                initial_value,
                repeatify: false,
            },
            callback,
            None,
        )
    }

    /// Construct QPA's `"next/repeat"` half-list.
    #[must_use]
    pub fn next_repeating(
        start: Integer,
        direction: InfiniteListDirection,
        initial_value: T,
        function: impl Fn(&T) -> Result<T, E> + Send + Sync + 'static,
        callback: Option<InfiniteListNewValueCallback<E>>,
        repeatify_callback: Option<InfiniteListRepeatifyCallback<T, E>>,
    ) -> Self {
        Self::new(
            start,
            direction,
            HalfInfiniteListSpecification::Next {
                function: Arc::new(function),
                initial_value,
                repeatify: true,
            },
            callback,
            repeatify_callback,
        )
    }

    /// Construct QPA's `"pos"` half-list.
    #[must_use]
    pub fn by_position(
        start: Integer,
        direction: InfiniteListDirection,
        function: impl Fn(&Integer) -> Result<T, E> + Send + Sync + 'static,
        store_values: bool,
        callback: Option<InfiniteListNewValueCallback<E>>,
    ) -> Self {
        Self::new(
            start,
            direction,
            HalfInfiniteListSpecification::Position {
                function: Arc::new(function),
                store_values,
            },
            callback,
            None,
        )
    }

    /// Translate `StartPosition`.
    #[must_use]
    pub fn start_position(&self) -> &Integer {
        &self.start
    }

    /// Translate `Direction`.
    #[must_use]
    pub const fn direction(&self) -> InfiniteListDirection {
        self.direction
    }

    /// Translate `InfListType`.
    ///
    /// # Errors
    ///
    /// Reports a poisoned lazy-list cache.
    pub fn list_type(&self) -> Result<HalfInfiniteListType, E> {
        let data = self.lock_data()?;
        Ok(match data.kind {
            HalfInfiniteListKind::Repeat(_) => HalfInfiniteListType::Repeat,
            HalfInfiniteListKind::Next { .. } => HalfInfiniteListType::Next,
            HalfInfiniteListKind::Position { .. } => HalfInfiniteListType::Position,
        })
    }

    /// Translate `RepeatingList`, returning `None` for GAP's `fail`.
    ///
    /// # Errors
    ///
    /// Reports a poisoned lazy-list cache.
    pub fn repeating_list(&self) -> Result<Option<Vec<T>>, E> {
        let data = self.lock_data()?;
        Ok(match &data.kind {
            HalfInfiniteListKind::Repeat(values) => Some(values.clone()),
            _ => None,
        })
    }

    /// Translate `ElementFunction`, returning `None` for GAP's `fail`.
    ///
    /// # Errors
    ///
    /// Reports a poisoned lazy-list cache.
    pub fn element_function(&self) -> Result<Option<HalfInfiniteListElementFunction<T, E>>, E> {
        let data = self.lock_data()?;
        Ok(match &data.kind {
            HalfInfiniteListKind::Repeat(_) => None,
            HalfInfiniteListKind::Next { function, .. } => {
                Some(HalfInfiniteListElementFunction::Next(Arc::clone(function)))
            }
            HalfInfiniteListKind::Position { function, .. } => Some(
                HalfInfiniteListElementFunction::Position(Arc::clone(function)),
            ),
        })
    }

    /// Translate `IsStoringValues`.
    ///
    /// # Errors
    ///
    /// Reports a poisoned lazy-list cache.
    pub fn is_storing_values(&self) -> Result<bool, E> {
        let data = self.lock_data()?;
        Ok(match data.kind {
            HalfInfiniteListKind::Repeat(_) => false,
            HalfInfiniteListKind::Next { .. } => true,
            HalfInfiniteListKind::Position { store_values, .. } => store_values,
        })
    }

    /// Translate `NewValueCallback`.
    ///
    /// # Errors
    ///
    /// Reports a poisoned lazy-list cache.
    pub fn new_value_callback(&self) -> Result<Option<InfiniteListNewValueCallback<E>>, E> {
        Ok(self.lock_data()?.callback.clone())
    }

    /// Translate `IsRepeating`.
    ///
    /// # Errors
    ///
    /// Reports a poisoned lazy-list cache.
    pub fn is_repeating(&self) -> Result<bool, E> {
        Ok(matches!(
            self.lock_data()?.kind,
            HalfInfiniteListKind::Repeat(_)
        ))
    }

    /// Translate `InitialValue`, returning `None` for GAP's `fail`.
    ///
    /// # Errors
    ///
    /// Reports a poisoned lazy-list cache.
    pub fn initial_value(&self) -> Result<Option<T>, E> {
        let data = self.lock_data()?;
        Ok(match &data.kind {
            HalfInfiniteListKind::Next { initial_value, .. } => Some(initial_value.clone()),
            _ => None,
        })
    }

    /// Return the values computed and stored by a lazy half-list.
    ///
    /// This exposes QPA's `list!.values` state for diagnostics and exact
    /// laziness tests.
    ///
    /// # Errors
    ///
    /// Reports a poisoned lazy-list cache.
    pub fn stored_values(&self) -> Result<Vec<T>, E> {
        Ok(self.lock_data()?.values.clone())
    }

    /// Return whether two handles denote the identical QPA half-list object.
    #[must_use]
    pub fn is_identical(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.data, &other.data)
    }

    /// Translate half-list exponentiation, `list ^ position`.
    ///
    /// # Errors
    ///
    /// Reports QPA's illegal-position and empty-period failures, a position
    /// outside Rust's address space, a poisoned cache, or an element/callback
    /// failure supplied by the caller.
    #[allow(clippy::too_many_lines)]
    pub fn value(&self, position: &Integer) -> Result<T, E> {
        let distance = (position - &self.start) * self.direction.factor();
        if distance < 0 {
            return Err(InfiniteListError::IllegalHalfPosition {
                position: position.clone(),
                start: self.start.clone(),
                direction: self.direction,
            }
            .into());
        }
        let index = usize::try_from(distance.clone()).map_err(|_| {
            E::from(InfiniteListError::PositionNotRepresentable {
                position: position.clone(),
            })
        })?;

        loop {
            enum Action<T, E> {
                Return(T),
                Next {
                    function: InfiniteListNextFunction<T, E>,
                    input: T,
                    computed_index: usize,
                },
                Position {
                    function: InfiniteListPositionFunction<T, E>,
                    position: Integer,
                    computed_index: Option<usize>,
                },
            }

            let action = {
                let data = self.lock_data()?;
                match &data.kind {
                    HalfInfiniteListKind::Repeat(repeating) => {
                        if repeating.is_empty() {
                            return Err(InfiniteListError::EmptyRepeatingList.into());
                        }
                        Action::Return(repeating[index % repeating.len()].clone())
                    }
                    HalfInfiniteListKind::Next {
                        function,
                        initial_value,
                        ..
                    } => {
                        if let Some(value) = data.values.get(index) {
                            Action::Return(value.clone())
                        } else {
                            Action::Next {
                                function: Arc::clone(function),
                                input: data
                                    .values
                                    .last()
                                    .cloned()
                                    .unwrap_or_else(|| initial_value.clone()),
                                computed_index: data.values.len(),
                            }
                        }
                    }
                    HalfInfiniteListKind::Position {
                        function,
                        store_values,
                    } => {
                        if !store_values {
                            Action::Position {
                                function: Arc::clone(function),
                                position: position.clone(),
                                computed_index: None,
                            }
                        } else if let Some(value) = data.values.get(index) {
                            Action::Return(value.clone())
                        } else {
                            let computed_index = data.values.len();
                            Action::Position {
                                function: Arc::clone(function),
                                position: &self.start
                                    + Integer::from(computed_index) * self.direction.factor(),
                                computed_index: Some(computed_index),
                            }
                        }
                    }
                }
            };

            match action {
                Action::Return(value) => return Ok(value),
                Action::Position {
                    function,
                    position: computed_position,
                    computed_index: None,
                } => return function(&computed_position),
                Action::Position {
                    function,
                    position: computed_position,
                    computed_index: Some(computed_index),
                } => {
                    let value = function(&computed_position)?;
                    let callback = {
                        let mut data = self.lock_data()?;
                        if data.values.len() != computed_index {
                            continue;
                        }
                        data.values.push(value);
                        data.callback.clone()
                    };
                    if let Some(callback) = callback {
                        callback(
                            &computed_position,
                            self.direction,
                            HalfInfiniteListType::Position,
                        )?;
                    }
                }
                Action::Next {
                    function,
                    input,
                    computed_index,
                } => {
                    let value = function(&input)?;
                    let (callback, should_repeatify) = {
                        let mut data = self.lock_data()?;
                        if data.values.len() != computed_index {
                            continue;
                        }
                        data.values.push(value);
                        let repeatify = matches!(
                            data.kind,
                            HalfInfiniteListKind::Next {
                                repeatify: true,
                                ..
                            }
                        );
                        (data.callback.clone(), repeatify)
                    };
                    let computed_position =
                        &self.start + Integer::from(computed_index) * self.direction.factor();
                    if let Some(callback) = callback {
                        callback(
                            &computed_position,
                            self.direction,
                            HalfInfiniteListType::Next,
                        )?;
                    }
                    if should_repeatify && computed_index >= 2 && computed_index % 2 == 0 {
                        let collision = {
                            let data = self.lock_data()?;
                            data.values[computed_index / 2] == data.values[computed_index]
                        };
                        if collision {
                            let (tail, repeating) = self.make_repeating_list_zero_based(
                                computed_index / 2,
                                computed_index,
                            )?;
                            let repeatify_callback = self
                                .lock_data()?
                                .repeatify_callback
                                .clone()
                                .ok_or_else(|| {
                                E::from(InfiniteListError::MissingRepeatifyCallback)
                            })?;
                            repeatify_callback(tail, repeating.clone())?;
                            if let HalfInfiniteListKind::Next { repeatify, .. } =
                                &mut self.lock_data()?.kind
                            {
                                *repeatify = false;
                            }
                            return repeating.value(position);
                        }
                    }
                }
            }
        }
    }

    /// Translate `MakeRepeatingList` using QPA's one-based collision indices.
    ///
    /// # Errors
    ///
    /// The indices must address stored values and QPA's repeated first value
    /// must occur again after the computed tail.
    pub fn make_repeating_list(
        &self,
        collision_index_one: usize,
        collision_index_two: usize,
    ) -> Result<(Vec<T>, Self), E> {
        if collision_index_one == 0 || collision_index_two == 0 {
            return Err(InfiniteListError::CollisionIndexNotPositive.into());
        }
        self.make_repeating_list_zero_based(collision_index_one - 1, collision_index_two - 1)
    }

    fn make_repeating_list_zero_based(
        &self,
        collision_index_one: usize,
        collision_index_two: usize,
    ) -> Result<(Vec<T>, Self), E> {
        let values = self.lock_data()?.values.clone();
        if collision_index_one >= values.len() {
            return Err(InfiniteListError::CollisionIndexOutsideStoredValues {
                index: collision_index_one + 1,
                stored: values.len(),
            }
            .into());
        }
        if collision_index_two >= values.len() {
            return Err(InfiniteListError::CollisionIndexOutsideStoredValues {
                index: collision_index_two + 1,
                stored: values.len(),
            }
            .into());
        }

        let mut first = collision_index_one + 1;
        let mut second = collision_index_two + 1;
        while first > 0 && second > 0 && values[first - 1] == values[second - 1] {
            first -= 1;
            second -= 1;
        }
        let repeat_start = first;
        let Some(repeat_value) = values.get(repeat_start) else {
            return Err(InfiniteListError::RepeatingCollisionNotFound.into());
        };
        let repeat_end = values
            .iter()
            .enumerate()
            .skip(repeat_start + 1)
            .find_map(|(index, value)| (value == repeat_value).then_some(index))
            .ok_or_else(|| E::from(InfiniteListError::RepeatingCollisionNotFound))?;

        let mut tail = values[..repeat_start].to_vec();
        if self.direction == InfiniteListDirection::Negative {
            tail.reverse();
        }
        let repeating_values = values[repeat_start..repeat_end].to_vec();
        let repeating_start = &self.start + Integer::from(repeat_start) * self.direction.factor();
        Ok((
            tail,
            Self::repeating(repeating_start, self.direction, repeating_values),
        ))
    }

    /// Translate half-list `LowestKnownPosition`.
    ///
    /// `None` is QPA's `"none"` sentinel.
    ///
    /// # Errors
    ///
    /// Reports a poisoned lazy-list cache.
    pub fn lowest_known_position(&self) -> Result<Option<ExtendedInteger>, E> {
        let data = self.lock_data()?;
        if self.direction == InfiniteListDirection::Positive {
            if data.values.is_empty() && !matches!(data.kind, HalfInfiniteListKind::Repeat(_)) {
                Ok(None)
            } else {
                Ok(Some(ExtendedInteger::Finite(self.start.clone())))
            }
        } else if matches!(data.kind, HalfInfiniteListKind::Repeat(_)) {
            Ok(Some(ExtendedInteger::NegativeInfinity))
        } else if data.values.is_empty() {
            Ok(None)
        } else {
            Ok(Some(ExtendedInteger::Finite(
                &self.start - Integer::from(data.values.len() - 1),
            )))
        }
    }

    /// Translate half-list `HighestKnownPosition`.
    ///
    /// `None` is QPA's `"none"` sentinel.
    ///
    /// # Errors
    ///
    /// Reports a poisoned lazy-list cache.
    pub fn highest_known_position(&self) -> Result<Option<ExtendedInteger>, E> {
        let data = self.lock_data()?;
        if self.direction == InfiniteListDirection::Negative {
            if data.values.is_empty() && !matches!(data.kind, HalfInfiniteListKind::Repeat(_)) {
                Ok(None)
            } else {
                Ok(Some(ExtendedInteger::Finite(self.start.clone())))
            }
        } else if matches!(data.kind, HalfInfiniteListKind::Repeat(_)) {
            Ok(Some(ExtendedInteger::PositiveInfinity))
        } else if data.values.is_empty() {
            Ok(None)
        } else {
            Ok(Some(ExtendedInteger::Finite(
                &self.start + Integer::from(data.values.len() - 1),
            )))
        }
    }

    /// Translate QPA's half-list `Shift` method.
    ///
    /// # Errors
    ///
    /// Reports a poisoned source cache.
    pub fn shift(&self, shift: &Integer) -> Result<Self, E> {
        let new_start = &self.start - shift;
        if let Some(repeating) = self.repeating_list()? {
            return Ok(Self::repeating(new_start, self.direction, repeating));
        }
        let original = self.clone();
        let shift = shift.clone();
        Ok(Self::by_position(
            new_start,
            self.direction,
            move |position| original.value(&(position + &shift)),
            false,
            None,
        ))
    }

    /// Translate QPA's `Cut` method.
    ///
    /// The return type is an `InfiniteList`, matching the source
    /// implementation (despite its half-list wording in the comment).
    ///
    /// # Errors
    ///
    /// Reports a position outside the half-list or propagates lazy access and
    /// list-construction failures.
    pub fn cut(&self, position: &Integer) -> Result<InfiniteList<T, E>, E> {
        let before_start = match self.direction {
            InfiniteListDirection::Positive => position < &self.start,
            InfiniteListDirection::Negative => position > &self.start,
        };
        if before_start {
            return Err(InfiniteListError::CutPositionOutsideRange {
                position: position.clone(),
                start: self.start.clone(),
                direction: self.direction,
            }
            .into());
        }

        let (middle, half) = if position == &self.start {
            (Vec::new(), self.clone())
        } else if let Some(repeating) = self.repeating_list()? {
            if repeating.is_empty() {
                return Err(InfiniteListError::EmptyRepeatingList.into());
            }
            let displacement = (&self.start - position) * self.direction.factor();
            let period = Integer::from(repeating.len());
            let mut remainder = displacement % &period;
            if remainder < 0 {
                remainder += &period;
            }
            let middle_length = usize::try_from(remainder).map_err(|_| {
                E::from(InfiniteListError::PositionNotRepresentable {
                    position: position.clone(),
                })
            })?;
            let middle_start = if self.direction == InfiniteListDirection::Positive {
                position.clone()
            } else {
                position - Integer::from(middle_length.saturating_sub(1))
            };
            let middle_end = if self.direction == InfiniteListDirection::Positive {
                position + Integer::from(middle_length.saturating_sub(1))
            } else {
                position.clone()
            };
            let middle = if middle_length == 0 {
                Vec::new()
            } else {
                integer_range_values(&middle_start, &middle_end, |index| self.value(index))?
            };
            let new_start = position + Integer::from(middle_length) * self.direction.factor();
            let half = self.shift(&(&self.start - &new_start))?;
            (middle, half)
        } else {
            let original = self.clone();
            (
                Vec::new(),
                Self::by_position(
                    position.clone(),
                    self.direction,
                    move |index| original.value(index),
                    false,
                    None,
                ),
            )
        };

        if self.direction == InfiniteListDirection::Positive {
            InfiniteList::from_parts(position.clone(), middle, Some(half), None)
        } else {
            InfiniteList::from_parts(
                position + Integer::from(1) - Integer::from(middle.len()),
                middle,
                None,
                Some(half),
            )
        }
    }

    /// Translate QPA's `HalfInfList(list, func)` mapping method.
    ///
    /// # Errors
    ///
    /// Reports a poisoned source cache or propagates eager mapping failures
    /// on a repeating period.
    pub fn map<U>(
        &self,
        function: InfiniteListMapFunction<T, U, E>,
    ) -> Result<HalfInfiniteList<U, E>, E>
    where
        U: Clone + PartialEq + Send + Sync + 'static,
    {
        if let Some(repeating) = self.repeating_list()? {
            let mapped = repeating
                .iter()
                .map(|value| function(value))
                .collect::<Result<Vec<_>, _>>()?;
            return Ok(HalfInfiniteList::repeating(
                self.start.clone(),
                self.direction,
                mapped,
            ));
        }
        let original = self.clone();
        Ok(HalfInfiniteList::by_position(
            self.start.clone(),
            self.direction,
            move |position| {
                let value = original.value(position)?;
                function(&value)
            },
            false,
            None,
        ))
    }

    fn lock_data(&self) -> Result<std::sync::MutexGuard<'_, HalfInfiniteListData<T, E>>, E> {
        self.data
            .lock()
            .map_err(|_| E::from(InfiniteListError::StatePoisoned))
    }
}

struct InfiniteListData<T, E> {
    middle: Vec<T>,
    positive: Option<HalfInfiniteList<T, E>>,
    negative: Option<HalfInfiniteList<T, E>>,
}

/// QPA's `IsInfList` object.
pub struct InfiniteList<T, E = InfiniteListError> {
    base_position: Integer,
    data: Arc<Mutex<InfiniteListData<T, E>>>,
}

impl<T, E> Clone for InfiniteList<T, E> {
    fn clone(&self) -> Self {
        Self {
            base_position: self.base_position.clone(),
            data: Arc::clone(&self.data),
        }
    }
}

impl<T, E> fmt::Debug for InfiniteList<T, E> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("InfiniteList")
            .field("base_position", &self.base_position)
            .finish_non_exhaustive()
    }
}

impl<T, E> InfiniteList<T, E>
where
    T: Clone + PartialEq + Send + Sync + 'static,
    E: From<InfiniteListError> + Send + Sync + 'static,
{
    /// Translate `MakeInfList`.
    #[must_use]
    pub fn new(
        base_position: Integer,
        middle: Vec<T>,
        positive: HalfInfiniteListSpecification<T, E>,
        negative: HalfInfiniteListSpecification<T, E>,
        callback: Option<InfiniteListNewValueCallback<E>>,
    ) -> Self {
        let base_for_data = base_position.clone();
        let data = Arc::new_cyclic(|weak: &Weak<Mutex<InfiniteListData<T, E>>>| {
            let positive_parent = Weak::clone(weak);
            let positive_repeatify: InfiniteListRepeatifyCallback<T, E> =
                Arc::new(move |tail, repeating| {
                    let parent = positive_parent
                        .upgrade()
                        .ok_or_else(|| E::from(InfiniteListError::ParentListDropped))?;
                    let mut parent = parent
                        .lock()
                        .map_err(|_| E::from(InfiniteListError::StatePoisoned))?;
                    parent.middle.extend(tail);
                    parent.positive = Some(repeating);
                    Ok(())
                });
            let negative_parent = Weak::clone(weak);
            let negative_repeatify: InfiniteListRepeatifyCallback<T, E> =
                Arc::new(move |tail, repeating| {
                    let parent = negative_parent
                        .upgrade()
                        .ok_or_else(|| E::from(InfiniteListError::ParentListDropped))?;
                    let mut parent = parent
                        .lock()
                        .map_err(|_| E::from(InfiniteListError::StatePoisoned))?;
                    let mut new_middle = tail;
                    new_middle.append(&mut parent.middle);
                    parent.middle = new_middle;
                    parent.negative = Some(repeating);
                    Ok(())
                });

            Mutex::new(InfiniteListData {
                positive: Some(HalfInfiniteList::new_with_repeatify(
                    &base_for_data + Integer::from(middle.len()),
                    InfiniteListDirection::Positive,
                    positive,
                    callback.clone(),
                    Some(positive_repeatify),
                )),
                negative: Some(HalfInfiniteList::new_with_repeatify(
                    &base_for_data - Integer::from(1),
                    InfiniteListDirection::Negative,
                    negative,
                    callback,
                    Some(negative_repeatify),
                )),
                middle,
            })
        });
        Self {
            base_position,
            data,
        }
    }

    /// Translate `MakeInfListFromHalfInfLists`.
    ///
    /// # Errors
    ///
    /// QPA requires the negative half to start immediately below the base
    /// position and the positive half immediately above the middle part.
    pub fn from_parts(
        base_position: Integer,
        middle: Vec<T>,
        positive: Option<HalfInfiniteList<T, E>>,
        negative: Option<HalfInfiniteList<T, E>>,
    ) -> Result<Self, E> {
        if let Some(negative) = &negative {
            let expected = negative.start_position() + Integer::from(1);
            if base_position != expected {
                return Err(InfiniteListError::NegativeStartMismatch {
                    base_position,
                    negative_start: negative.start_position().clone(),
                }
                .into());
            }
        }
        if let Some(positive) = &positive {
            let expected = &base_position + Integer::from(middle.len());
            if positive.start_position() != &expected {
                return Err(InfiniteListError::PositiveStartMismatch {
                    expected,
                    positive_start: positive.start_position().clone(),
                }
                .into());
            }
        }
        Ok(Self {
            base_position,
            data: Arc::new(Mutex::new(InfiniteListData {
                middle,
                positive,
                negative,
            })),
        })
    }

    /// Translate `FunctionInfList`.
    #[must_use]
    pub fn function(function: impl Fn(&Integer) -> Result<T, E> + Send + Sync + 'static) -> Self {
        let function: InfiniteListPositionFunction<T, E> = Arc::new(function);
        Self::new(
            Integer::from(0),
            Vec::new(),
            HalfInfiniteListSpecification::Position {
                function: Arc::clone(&function),
                store_values: false,
            },
            HalfInfiniteListSpecification::Position {
                function,
                store_values: false,
            },
            None,
        )
    }

    /// Translate `ConstantInfList`.
    #[must_use]
    pub fn constant(value: T) -> Self {
        Self::new(
            Integer::from(0),
            Vec::new(),
            HalfInfiniteListSpecification::Repeat(vec![value.clone()]),
            HalfInfiniteListSpecification::Repeat(vec![value]),
            None,
        )
    }

    /// Translate `FiniteInfList`.
    #[must_use]
    pub fn finite(base_position: Integer, values: Vec<T>) -> Self {
        Self {
            base_position,
            data: Arc::new(Mutex::new(InfiniteListData {
                middle: values,
                positive: None,
                negative: None,
            })),
        }
    }

    /// Translate `MiddleStart`.
    #[must_use]
    pub fn middle_start(&self) -> &Integer {
        &self.base_position
    }

    /// Translate `MiddleEnd`.
    ///
    /// # Errors
    ///
    /// Reports a poisoned list state.
    pub fn middle_end(&self) -> Result<Integer, E> {
        Ok(&self.base_position + Integer::from(self.lock_data()?.middle.len()) - Integer::from(1))
    }

    /// Translate `MiddlePart`.
    ///
    /// # Errors
    ///
    /// Reports a poisoned list state.
    pub fn middle_part(&self) -> Result<Vec<T>, E> {
        Ok(self.lock_data()?.middle.clone())
    }

    /// Translate `PositivePart`.
    ///
    /// # Errors
    ///
    /// Reports a poisoned list state.
    pub fn positive_part(&self) -> Result<Option<HalfInfiniteList<T, E>>, E> {
        Ok(self.lock_data()?.positive.clone())
    }

    /// Translate `NegativePart`.
    ///
    /// # Errors
    ///
    /// Reports a poisoned list state.
    pub fn negative_part(&self) -> Result<Option<HalfInfiniteList<T, E>>, E> {
        Ok(self.lock_data()?.negative.clone())
    }

    /// Return whether two handles denote the identical QPA infinite-list object.
    #[must_use]
    pub fn is_identical(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.data, &other.data)
    }

    /// Translate full-list `LowestKnownPosition`.
    ///
    /// `None` is QPA's `"none"` sentinel.
    ///
    /// # Errors
    ///
    /// A finite empty list reproduces QPA's missing-positive-part dispatch
    /// failure; other lazy and state failures propagate.
    pub fn lowest_known_position(&self) -> Result<Option<ExtendedInteger>, E> {
        let (middle_nonempty, positive, negative) = {
            let data = self.lock_data()?;
            (
                !data.middle.is_empty(),
                data.positive.clone(),
                data.negative.clone(),
            )
        };
        if let Some(negative) = negative
            && let Some(position) = negative.lowest_known_position()?
        {
            return Ok(Some(position));
        }
        if middle_nonempty {
            return Ok(Some(ExtendedInteger::Finite(self.base_position.clone())));
        }
        positive
            .ok_or_else(|| {
                E::from(InfiniteListError::MissingPositivePart {
                    operation: "LowestKnownPosition",
                })
            })?
            .lowest_known_position()
    }

    /// Translate full-list `HighestKnownPosition`.
    ///
    /// `None` is QPA's `"none"` sentinel.
    ///
    /// # Errors
    ///
    /// A finite empty list reproduces QPA's missing-negative-part dispatch
    /// failure; other lazy and state failures propagate.
    pub fn highest_known_position(&self) -> Result<Option<ExtendedInteger>, E> {
        let (middle_nonempty, positive, negative) = {
            let data = self.lock_data()?;
            (
                !data.middle.is_empty(),
                data.positive.clone(),
                data.negative.clone(),
            )
        };
        if let Some(positive) = positive
            && let Some(position) = positive.highest_known_position()?
        {
            return Ok(Some(position));
        }
        if middle_nonempty {
            return Ok(Some(ExtendedInteger::Finite(self.middle_end()?)));
        }
        negative
            .ok_or_else(|| {
                E::from(InfiniteListError::MissingNegativePart {
                    operation: "HighestKnownPosition",
                })
            })?
            .highest_known_position()
    }

    /// Translate full-list `UpperBound`.
    ///
    /// # Errors
    ///
    /// Reports a poisoned list state.
    pub fn upper_bound(&self) -> Result<ExtendedInteger, E> {
        let data = self.lock_data()?;
        if data.positive.is_some() {
            Ok(ExtendedInteger::PositiveInfinity)
        } else if !data.middle.is_empty() {
            Ok(ExtendedInteger::Finite(
                &self.base_position + Integer::from(data.middle.len()) - Integer::from(1),
            ))
        } else if let Some(negative) = &data.negative {
            Ok(ExtendedInteger::Finite(negative.start_position().clone()))
        } else {
            Ok(ExtendedInteger::NegativeInfinity)
        }
    }

    /// Translate full-list `LowerBound`.
    ///
    /// # Errors
    ///
    /// Reports a poisoned list state.
    pub fn lower_bound(&self) -> Result<ExtendedInteger, E> {
        let data = self.lock_data()?;
        if data.negative.is_some() {
            Ok(ExtendedInteger::NegativeInfinity)
        } else if !data.middle.is_empty() {
            Ok(ExtendedInteger::Finite(self.base_position.clone()))
        } else if let Some(positive) = &data.positive {
            Ok(ExtendedInteger::Finite(positive.start_position().clone()))
        } else {
            Ok(ExtendedInteger::PositiveInfinity)
        }
    }

    /// Translate full-list exponentiation, `list ^ position`.
    ///
    /// # Errors
    ///
    /// Reports a position outside a finite list or propagates half-list and
    /// state failures.
    pub fn value(&self, position: &Integer) -> Result<T, E> {
        let (middle, positive, negative) = {
            let data = self.lock_data()?;
            (
                data.middle.clone(),
                data.positive.clone(),
                data.negative.clone(),
            )
        };
        let middle_end = &self.base_position + Integer::from(middle.len()) - Integer::from(1);
        if position >= &self.base_position && position <= &middle_end {
            let index = usize::try_from(position - &self.base_position).map_err(|_| {
                E::from(InfiniteListError::PositionNotRepresentable {
                    position: position.clone(),
                })
            })?;
            return Ok(middle[index].clone());
        }
        if let Some(positive) = positive
            && position >= positive.start_position()
        {
            return positive.value(position);
        }
        if let Some(negative) = negative
            && position <= negative.start_position()
        {
            return negative.value(position);
        }
        Err(InfiniteListError::PositionOutsideInfiniteList {
            position: position.clone(),
        }
        .into())
    }

    /// Translate full-list `Shift`.
    ///
    /// QPA has no `Shift(fail, i)` method, so finite lists with either absent
    /// half produce the corresponding typed dispatch failure.
    ///
    /// # Errors
    ///
    /// Reports a missing half or propagates half-list and state failures.
    pub fn shift(&self, shift: &Integer) -> Result<Self, E> {
        let middle = self.middle_part()?;
        let positive = self.positive_part()?.ok_or_else(|| {
            E::from(InfiniteListError::MissingPositivePart { operation: "Shift" })
        })?;
        let negative = self.negative_part()?.ok_or_else(|| {
            E::from(InfiniteListError::MissingNegativePart { operation: "Shift" })
        })?;
        Self::from_parts(
            &self.base_position - shift,
            middle,
            Some(positive.shift(shift)?),
            Some(negative.shift(shift)?),
        )
    }

    /// Translate `FinitePartAsList`.
    ///
    /// # Errors
    ///
    /// Propagates value access failures.
    pub fn finite_part_as_list(
        &self,
        start_position: &Integer,
        end_position: &Integer,
    ) -> Result<Vec<T>, E> {
        integer_range_values(start_position, end_position, |position| {
            self.value(position)
        })
    }

    /// Translate `PositivePartFrom`.
    ///
    /// # Errors
    ///
    /// QPA dispatches on the positive half before considering the requested
    /// position, so an absent half is always an error. Other list failures
    /// propagate.
    pub fn positive_part_from(&self, position: &Integer) -> Result<Self, E> {
        let positive = self.positive_part()?.ok_or_else(|| {
            E::from(InfiniteListError::MissingPositivePart {
                operation: "PositivePartFrom",
            })
        })?;
        let positive_start = positive.start_position().clone();
        if position >= &positive_start {
            positive.cut(position)
        } else {
            let middle =
                self.finite_part_as_list(position, &(&positive_start - Integer::from(1)))?;
            Self::from_parts(position.clone(), middle, Some(positive), None)
        }
    }

    /// Translate `NegativePartFrom`.
    ///
    /// # Errors
    ///
    /// QPA dispatches on the negative half before considering the requested
    /// position, so an absent half is always an error. Other list failures
    /// propagate.
    pub fn negative_part_from(&self, position: &Integer) -> Result<Self, E> {
        let negative = self.negative_part()?.ok_or_else(|| {
            E::from(InfiniteListError::MissingNegativePart {
                operation: "NegativePartFrom",
            })
        })?;
        let negative_start = negative.start_position().clone();
        if position <= &negative_start {
            negative.cut(position)
        } else {
            let middle =
                self.finite_part_as_list(&(&negative_start + Integer::from(1)), position)?;
            Self::from_parts(
                &negative_start + Integer::from(1),
                middle,
                None,
                Some(negative),
            )
        }
    }

    /// Translate `Splice(positiveList, negativeList, joinPosition)`.
    ///
    /// # Errors
    ///
    /// Propagates the two source cuts and the checked join.
    pub fn splice(
        positive_list: &Self,
        negative_list: &Self,
        join_position: &Integer,
    ) -> Result<Self, E> {
        let positive_cut = positive_list.positive_part_from(&(join_position + Integer::from(1)))?;
        let negative_cut = negative_list.negative_part_from(join_position)?;
        let mut middle = negative_cut.middle_part()?;
        middle.extend(positive_cut.middle_part()?);
        Self::from_parts(
            negative_cut.middle_start().clone(),
            middle,
            positive_cut.positive_part()?,
            negative_cut.negative_part()?,
        )
    }

    /// Translate `InfConcatenation`.
    ///
    /// # Errors
    ///
    /// Propagates source metadata, half-list shifting, and checked join
    /// failures.
    pub fn concatenate(lists: &[Self]) -> Result<Self, E> {
        let Some(first) = lists.first() else {
            return Ok(Self::finite(Integer::from(0), Vec::new()));
        };
        let Some(last) = lists.last() else {
            return Ok(Self::finite(Integer::from(0), Vec::new()));
        };
        if lists.len() == 1 {
            return Ok(first.clone());
        }

        let mut middle = Vec::new();
        for list in lists.iter().rev() {
            middle.extend(list.middle_part()?);
        }
        let base_position = first.middle_end()? - Integer::from(middle.len()) + Integer::from(1);
        let positive = first.positive_part()?;
        let mut negative = last.negative_part()?;
        if let Some(negative_half) = negative {
            let shift = negative_half.start_position() - &base_position + Integer::from(1);
            negative = Some(negative_half.shift(&shift)?);
        }
        Self::from_parts(base_position, middle, positive, negative)
    }

    /// Translate `InfList(list, func)`.
    ///
    /// QPA eagerly dispatches `HalfInfList` on both halves, so a finite list
    /// with an absent half fails even when its middle could be mapped.
    ///
    /// # Errors
    ///
    /// Reports a missing half or propagates mapping and state failures.
    pub fn map<U>(
        &self,
        function: impl Fn(&T) -> Result<U, E> + Send + Sync + 'static,
    ) -> Result<InfiniteList<U, E>, E>
    where
        U: Clone + PartialEq + Send + Sync + 'static,
    {
        let function: InfiniteListMapFunction<T, U, E> = Arc::new(function);
        let middle = self
            .middle_part()?
            .iter()
            .map(|value| function(value))
            .collect::<Result<Vec<_>, _>>()?;
        let positive = self.positive_part()?.ok_or_else(|| {
            E::from(InfiniteListError::MissingPositivePart {
                operation: "InfList",
            })
        })?;
        let negative = self.negative_part()?.ok_or_else(|| {
            E::from(InfiniteListError::MissingNegativePart {
                operation: "InfList",
            })
        })?;
        InfiniteList::from_parts(
            self.base_position.clone(),
            middle,
            Some(positive.map(Arc::clone(&function))?),
            Some(negative.map(function)?),
        )
    }

    fn lock_data(&self) -> Result<std::sync::MutexGuard<'_, InfiniteListData<T, E>>, E> {
        self.data
            .lock()
            .map_err(|_| E::from(InfiniteListError::StatePoisoned))
    }
}

/// Translate QPA's declared but uninstalled `ShiftedSplice` operation.
///
/// The pinned QPA source declares this operation in `inflist.gd` but installs
/// no method in `inflist.gi`; every source call therefore reaches GAP's
/// `NoMethodFound`. Rust exposes that observable state as a typed failure
/// rather than inventing an implementation absent from QPA.
///
/// # Errors
///
/// Always returns [`InfiniteListError::ShiftedSpliceNotInstalled`].
pub fn shifted_splice<T, E>(
    _left: &InfiniteList<T, E>,
    _left_shift: &Integer,
    _right: &InfiniteList<T, E>,
    _right_shift: &Integer,
    _join_position: &Integer,
) -> Result<InfiniteList<T, E>, E>
where
    T: Clone + PartialEq + Send + Sync + 'static,
    E: From<InfiniteListError> + Send + Sync + 'static,
{
    Err(InfiniteListError::ShiftedSpliceNotInstalled.into())
}

/// Return QPA's global `IntegersList`.
#[must_use]
pub fn integers_list() -> InfiniteList<Integer> {
    InfiniteList::function(|position| Ok(position.clone()))
}

fn integer_range_values<T, E>(
    start: &Integer,
    end: &Integer,
    mut value: impl FnMut(&Integer) -> Result<T, E>,
) -> Result<Vec<T>, E> {
    let mut result = Vec::new();
    let mut position = start.clone();
    while position <= *end {
        result.push(value(&position)?);
        position += Integer::from(1);
    }
    Ok(result)
}

/// Failure in QPA's infinite-list machinery.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum InfiniteListError {
    /// A requested half-list position lies in the opposite direction.
    #[error(
        "illegal half-infinite-list position {position} from start {start} in direction {direction}"
    )]
    IllegalHalfPosition {
        /// Rejected position.
        position: Integer,
        /// Half-list start position.
        start: Integer,
        /// Half-list direction.
        direction: InfiniteListDirection,
    },
    /// A GAP integer distance cannot address a Rust vector.
    #[error("infinite-list position {position} is outside Rust's addressable range")]
    PositionNotRepresentable {
        /// Rejected arbitrary-precision position.
        position: Integer,
    },
    /// QPA permits construction of an empty repeating list but fails on use.
    #[error("a repeating half-infinite list cannot have an empty period")]
    EmptyRepeatingList,
    /// A mutex was poisoned by a panic while list state was borrowed.
    #[error("the infinite-list lazy cache is poisoned")]
    StatePoisoned,
    /// A cyclic construction callback outlived its parent list.
    #[error("the parent infinite list was dropped during repeat conversion")]
    ParentListDropped,
    /// A standalone `next/repeat` half-list reached QPA's attempted call of
    /// its non-function `false` repeatify callback.
    #[error("a next/repeat half-list has no repeatify callback")]
    MissingRepeatifyCallback,
    /// QPA declares collision positions to be positive integers.
    #[error("repeating-list collision indices must be positive")]
    CollisionIndexNotPositive,
    /// A collision index does not address a stored lazy value.
    #[error("collision index {index} is outside the {stored} stored values")]
    CollisionIndexOutsideStoredValues {
        /// Rejected one-based collision index.
        index: usize,
        /// Current number of stored values.
        stored: usize,
    },
    /// The entered collision data has no second occurrence delimiting a period.
    #[error("the stored values do not contain QPA's expected repeating collision")]
    RepeatingCollisionNotFound,
    /// `MakeInfListFromHalfInfLists` received a misaligned negative half.
    #[error(
        "negative half starts at {negative_start}, but the middle base position is {base_position}"
    )]
    NegativeStartMismatch {
        /// Entered middle base position.
        base_position: Integer,
        /// Entered negative-half start position.
        negative_start: Integer,
    },
    /// `MakeInfListFromHalfInfLists` received a misaligned positive half.
    #[error("positive half starts at {positive_start}, expected {expected}")]
    PositiveStartMismatch {
        /// Expected positive-half start.
        expected: Integer,
        /// Entered positive-half start.
        positive_start: Integer,
    },
    /// `Cut` received a position before the half-list start.
    #[error(
        "cut position {position} is outside the half-list starting at {start} in direction {direction}"
    )]
    CutPositionOutsideRange {
        /// Rejected cut position.
        position: Integer,
        /// Half-list start position.
        start: Integer,
        /// Half-list direction.
        direction: InfiniteListDirection,
    },
    /// A full-list position lies outside all three constituent parts.
    #[error("position {position} is outside the infinite-list index range")]
    PositionOutsideInfiniteList {
        /// Rejected position.
        position: Integer,
    },
    /// QPA attempted to dispatch an operation on an absent positive half.
    #[error("{operation} has no method for an absent positive half-list")]
    MissingPositivePart {
        /// QPA operation whose `fail` dispatch was reached.
        operation: &'static str,
    },
    /// QPA attempted to dispatch an operation on an absent negative half.
    #[error("{operation} has no method for an absent negative half-list")]
    MissingNegativePart {
        /// QPA operation whose `fail` dispatch was reached.
        operation: &'static str,
    },
    /// `ShiftedSplice` is declared but has no method in the pinned QPA source.
    #[error("the pinned QPA source installs no ShiftedSplice method")]
    ShiftedSpliceNotInstalled,
}

#[cfg(test)]
mod tests {
    use std::sync::{Arc, Mutex};

    use numerica::domains::integer::Integer;
    use serde_json::{Value, json};

    use super::{
        ExtendedInteger, HalfInfiniteList, HalfInfiniteListSpecification, HalfInfiniteListType,
        InfiniteList, InfiniteListDirection, InfiniteListError, InfiniteListNewValueCallback,
        InfiniteListNextFunction, integers_list, shifted_splice,
    };

    fn integer(value: i64) -> Integer {
        Integer::from(value)
    }

    fn values(list: &InfiniteList<i64>, positions: &[i64]) -> Result<Vec<i64>, InfiniteListError> {
        positions
            .iter()
            .map(|position| list.value(&integer(*position)))
            .collect()
    }

    fn extended_json(value: Option<ExtendedInteger>) -> Value {
        match value {
            None => json!("none"),
            Some(ExtendedInteger::NegativeInfinity) => json!("-infinity"),
            Some(ExtendedInteger::PositiveInfinity) => json!("infinity"),
            Some(ExtendedInteger::Finite(value)) => json!(i64::try_from(value).unwrap()),
        }
    }

    fn bound_json(value: ExtendedInteger) -> Value {
        extended_json(Some(value))
    }

    fn direction_value(direction: InfiniteListDirection) -> i64 {
        match direction {
            InfiniteListDirection::Positive => 1,
            InfiniteListDirection::Negative => -1,
        }
    }

    fn list_type_value(list_type: HalfInfiniteListType) -> &'static str {
        match list_type {
            HalfInfiniteListType::Repeat => "repeat",
            HalfInfiniteListType::Next => "next",
            HalfInfiniteListType::Position => "pos",
        }
    }

    fn callback_json(
        callbacks: &[(Integer, InfiniteListDirection, HalfInfiniteListType)],
    ) -> Value {
        json!(
            callbacks
                .iter()
                .map(|(position, direction, list_type)| json!([
                    i64::try_from(position.clone()).unwrap(),
                    direction_value(*direction),
                    list_type_value(*list_type),
                ]))
                .collect::<Vec<_>>()
        )
    }

    #[allow(clippy::too_many_lines)]
    #[test]
    fn pinned_qpa_half_infinite_list_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/inflist.schema-v1.json"
        ))
        .unwrap();

        let positive_repeat: HalfInfiniteList<i64> = HalfInfiniteList::repeating(
            integer(2),
            InfiniteListDirection::Positive,
            vec![10, 20, 30],
        );
        let negative_repeat: HalfInfiniteList<i64> =
            HalfInfiniteList::repeating(integer(-1), InfiniteListDirection::Negative, vec![7, 8]);
        assert_eq!(
            json!({
                "start": 2,
                "direction": direction_value(positive_repeat.direction()),
                "type": list_type_value(positive_repeat.list_type().unwrap()),
                "repeating": positive_repeat.repeating_list().unwrap().unwrap(),
                "storing": positive_repeat.is_storing_values().unwrap(),
                "initial_is_fail": positive_repeat.initial_value().unwrap().is_none(),
                "values": (2..=8).map(|position|
                    positive_repeat.value(&integer(position)).unwrap()).collect::<Vec<_>>(),
                "lowest": extended_json(positive_repeat.lowest_known_position().unwrap()),
                "highest": extended_json(positive_repeat.highest_known_position().unwrap()),
            }),
            expected["half_lists"]["positive_repeat"]
        );
        assert_eq!(
            json!({
                "values": (-6..=-1).map(|position|
                    negative_repeat.value(&integer(position)).unwrap()).collect::<Vec<_>>(),
                "lowest": extended_json(negative_repeat.lowest_known_position().unwrap()),
                "highest": extended_json(negative_repeat.highest_known_position().unwrap()),
            }),
            expected["half_lists"]["negative_repeat"]
        );

        let next_callbacks = Arc::new(Mutex::new(Vec::new()));
        let next_callback_log = Arc::clone(&next_callbacks);
        let next_callback: InfiniteListNewValueCallback<InfiniteListError> =
            Arc::new(move |position, direction, list_type| {
                next_callback_log
                    .lock()
                    .unwrap()
                    .push((position.clone(), direction, list_type));
                Ok(())
            });
        let repeat_data = Arc::new(Mutex::new(None));
        let repeat_data_log = Arc::clone(&repeat_data);
        let next_list: HalfInfiniteList<i64> = HalfInfiniteList::next_repeating(
            integer(0),
            InfiniteListDirection::Positive,
            -1,
            |value| Ok((value + 1).rem_euclid(4)),
            Some(next_callback),
            Some(Arc::new(move |tail, repeating| {
                *repeat_data_log.lock().unwrap() = Some((
                    tail,
                    repeating.start_position().clone(),
                    repeating.direction(),
                    repeating.repeating_list()?.unwrap(),
                ));
                Ok(())
            })),
        );
        let next_value = next_list.value(&integer(10)).unwrap();
        let (tail, repeat_start, repeat_direction, period) =
            repeat_data.lock().unwrap().clone().unwrap();
        assert_eq!(
            json!({
                "value_10": next_value,
                "stored": next_list.stored_values().unwrap(),
                "callbacks": callback_json(&next_callbacks.lock().unwrap()),
                "repeat_data": [
                    json!(tail),
                    json!(i64::try_from(repeat_start).unwrap()),
                    json!(direction_value(repeat_direction)),
                    json!(period),
                ],
            }),
            expected["half_lists"]["next_repeat"]
        );

        let position_callbacks = Arc::new(Mutex::new(Vec::new()));
        let position_callback_log = Arc::clone(&position_callbacks);
        let position_callback: InfiniteListNewValueCallback<InfiniteListError> =
            Arc::new(move |position, direction, list_type| {
                position_callback_log.lock().unwrap().push((
                    position.clone(),
                    direction,
                    list_type,
                ));
                Ok(())
            });
        let position_list: HalfInfiniteList<i64> = HalfInfiniteList::by_position(
            integer(-2),
            InfiniteListDirection::Positive,
            |position| Ok(i64::try_from(position.clone()).unwrap().pow(2)),
            true,
            Some(position_callback),
        );
        let lowest_before = position_list.lowest_known_position().unwrap();
        let value_two = position_list.value(&integer(2)).unwrap();
        let value_zero = position_list.value(&integer(0)).unwrap();
        assert_eq!(
            json!({
                "lowest_before": extended_json(lowest_before),
                "value_2": value_two,
                "value_0_again": value_zero,
                "stored": position_list.stored_values().unwrap(),
                "callbacks": callback_json(&position_callbacks.lock().unwrap()),
                "lowest_after": extended_json(position_list.lowest_known_position().unwrap()),
                "highest_after": extended_json(position_list.highest_known_position().unwrap()),
            }),
            expected["half_lists"]["position_storing"]
        );

        let nonstoring_callbacks = Arc::new(Mutex::new(Vec::new()));
        let nonstoring_callback_log = Arc::clone(&nonstoring_callbacks);
        let nonstoring_callback: InfiniteListNewValueCallback<InfiniteListError> =
            Arc::new(move |position, direction, list_type| {
                nonstoring_callback_log.lock().unwrap().push((
                    position.clone(),
                    direction,
                    list_type,
                ));
                Ok(())
            });
        let nonstoring: HalfInfiniteList<i64> = HalfInfiniteList::by_position(
            integer(5),
            InfiniteListDirection::Negative,
            |position| Ok(i64::try_from(position.clone()).unwrap() + 100),
            false,
            Some(nonstoring_callback),
        );
        assert_eq!(
            json!({
                "value_1": nonstoring.value(&integer(1)).unwrap(),
                "stored": nonstoring.stored_values().unwrap(),
                "callbacks": callback_json(&nonstoring_callbacks.lock().unwrap()),
                "lowest": extended_json(nonstoring.lowest_known_position().unwrap()),
            }),
            expected["half_lists"]["position_nonstoring"]
        );
    }

    #[allow(clippy::too_many_lines)]
    #[test]
    fn pinned_qpa_full_infinite_list_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/inflist.schema-v1.json"
        ))
        .unwrap();
        let full = InfiniteList::new(
            integer(0),
            vec![100_i64, 200],
            HalfInfiniteListSpecification::Repeat(vec![1, 2, 3]),
            HalfInfiniteListSpecification::Repeat(vec![7, 8]),
            None,
        );
        assert_eq!(
            json!({
                "middle_start": 0,
                "middle_end": i64::try_from(full.middle_end().unwrap()).unwrap(),
                "middle": full.middle_part().unwrap(),
                "positive_start": i64::try_from(full.positive_part().unwrap().unwrap()
                    .start_position().clone()).unwrap(),
                "negative_start": i64::try_from(full.negative_part().unwrap().unwrap()
                    .start_position().clone()).unwrap(),
                "values": values(&full, &[-5,-2,-1,0,1,2,5]).unwrap(),
                "lowest_known": extended_json(full.lowest_known_position().unwrap()),
                "highest_known": extended_json(full.highest_known_position().unwrap()),
                "lower_bound": bound_json(full.lower_bound().unwrap()),
                "upper_bound": bound_json(full.upper_bound().unwrap()),
            }),
            expected["full"]
        );

        let shifted = full.shift(&integer(3)).unwrap();
        let mapped = full.map(|value| Ok(2 * value)).unwrap();
        let positive_cut = full
            .positive_part()
            .unwrap()
            .unwrap()
            .cut(&integer(3))
            .unwrap();
        let negative_cut = full
            .negative_part()
            .unwrap()
            .unwrap()
            .cut(&integer(-2))
            .unwrap();
        let positive_from = full.positive_part_from(&integer(3)).unwrap();
        let negative_from = full.negative_part_from(&integer(-2)).unwrap();
        let right = InfiniteList::new(
            integer(0),
            vec![20_i64, 21],
            HalfInfiniteListSpecification::Repeat(vec![22, 23]),
            HalfInfiniteListSpecification::Repeat(vec![18, 19]),
            None,
        );
        let spliced = InfiniteList::splice(&full, &right, &integer(0)).unwrap();
        let concatenated = InfiniteList::concatenate(&[full.clone(), right]).unwrap();
        let positive_repeat: HalfInfiniteList<i64> = HalfInfiniteList::repeating(
            integer(2),
            InfiniteListDirection::Positive,
            vec![10, 20, 30],
        );
        let half_mapped = positive_repeat
            .map(Arc::new(|value| Ok(value + 1)))
            .unwrap();
        assert_eq!(
            json!({
                "shift": {
                    "middle_start": i64::try_from(shifted.middle_start().clone()).unwrap(),
                    "values": values(&shifted, &[-5,-3,-2,-1,2]).unwrap(),
                },
                "map": values(&mapped, &[-2,-1,0,1,2]).unwrap(),
                "positive_cut": {
                    "middle_start": i64::try_from(positive_cut.middle_start().clone()).unwrap(),
                    "middle": positive_cut.middle_part().unwrap(),
                    "positive_start": i64::try_from(positive_cut.positive_part().unwrap().unwrap()
                        .start_position().clone()).unwrap(),
                    "values": values(&positive_cut, &[3,4,5,6,7,8]).unwrap(),
                },
                "negative_cut": {
                    "middle_start": i64::try_from(negative_cut.middle_start().clone()).unwrap(),
                    "middle": negative_cut.middle_part().unwrap(),
                    "negative_start": i64::try_from(negative_cut.negative_part().unwrap().unwrap()
                        .start_position().clone()).unwrap(),
                    "values": values(&negative_cut, &[-6,-5,-4,-3,-2]).unwrap(),
                },
                "positive_from_values": values(&positive_from, &[3,4,5,6,7,8]).unwrap(),
                "negative_from_values": values(&negative_from, &[-6,-5,-4,-3,-2]).unwrap(),
                "splice": values(&spliced, &[-2,-1,0,1,2,3]).unwrap(),
                "concatenation": {
                    "middle_start": i64::try_from(concatenated.middle_start().clone()).unwrap(),
                    "middle": concatenated.middle_part().unwrap(),
                    "positive_start": i64::try_from(concatenated.positive_part().unwrap().unwrap()
                        .start_position().clone()).unwrap(),
                    "negative_start": i64::try_from(concatenated.negative_part().unwrap().unwrap()
                        .start_position().clone()).unwrap(),
                    "values": values(&concatenated, &[-4,-3,-2,-1,0,1,2,3]).unwrap(),
                },
                "half_map": (2..=7).map(|position|
                    half_mapped.value(&integer(position)).unwrap()).collect::<Vec<_>>(),
            }),
            expected["transformations"]
        );

        let finite = InfiniteList::finite(integer(-1), vec![4_i64, 5, 6]);
        let constant = InfiniteList::constant(9_i64);
        let integers = integers_list();
        let integer_values = [-3, 0, 4]
            .map(|position| i64::try_from(integers.value(&integer(position)).unwrap()).unwrap());
        assert_eq!(
            json!({
                "finite": {
                    "middle_start": -1,
                    "middle_end": i64::try_from(finite.middle_end().unwrap()).unwrap(),
                    "values": values(&finite, &[-1,0,1]).unwrap(),
                    "lower_bound": bound_json(finite.lower_bound().unwrap()),
                    "upper_bound": bound_json(finite.upper_bound().unwrap()),
                },
                "constant": values(&constant, &[-3,0,4]).unwrap(),
                "integers": integer_values,
            }),
            expected["constructors"]
        );
    }

    #[test]
    fn pinned_qpa_repeat_conversion_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/inflist.schema-v1.json"
        ))
        .unwrap();
        let tail_next: InfiniteListNextFunction<i64, InfiniteListError> = Arc::new(|value| {
            Ok(match value {
                -1 => 9,
                5 => 6,
                6 => 7,
                _ => 5,
            })
        });
        let positive_parent = InfiniteList::new(
            integer(0),
            vec![99_i64],
            HalfInfiniteListSpecification::Next {
                function: Arc::clone(&tail_next),
                initial_value: -1,
                repeatify: true,
            },
            HalfInfiniteListSpecification::Repeat(vec![8]),
            None,
        );
        let positive_parent_value = positive_parent.value(&integer(10)).unwrap();
        let negative_parent = InfiniteList::new(
            integer(0),
            vec![99_i64],
            HalfInfiniteListSpecification::Repeat(vec![8]),
            HalfInfiniteListSpecification::Next {
                function: tail_next,
                initial_value: -1,
                repeatify: true,
            },
            None,
        );
        let negative_parent_value = negative_parent.value(&integer(-10)).unwrap();
        assert_eq!(
            json!({
                "positive_parent": {
                    "value_10": positive_parent_value,
                    "middle_start": 0,
                    "middle": positive_parent.middle_part().unwrap(),
                    "positive_start": i64::try_from(positive_parent.positive_part().unwrap()
                        .unwrap().start_position().clone()).unwrap(),
                    "period": positive_parent.positive_part().unwrap().unwrap()
                        .repeating_list().unwrap().unwrap(),
                },
                "negative_parent": {
                    "value_minus_10": negative_parent_value,
                    "middle_start": 0,
                    "middle": negative_parent.middle_part().unwrap(),
                    "negative_start": i64::try_from(negative_parent.negative_part().unwrap()
                        .unwrap().start_position().clone()).unwrap(),
                    "period": negative_parent.negative_part().unwrap().unwrap()
                        .repeating_list().unwrap().unwrap(),
                },
            }),
            expected["repeat_conversion"]
        );
    }

    #[test]
    fn repeating_position_and_next_halves_follow_qpa_laziness() {
        let repeating: HalfInfiniteList<i64> = HalfInfiniteList::repeating(
            integer(2),
            InfiniteListDirection::Positive,
            vec![10, 20, 30],
        );
        assert_eq!(
            (2..=8)
                .map(|position| repeating.value(&integer(position)).unwrap())
                .collect::<Vec<_>>(),
            vec![10, 20, 30, 10, 20, 30, 10]
        );
        assert_eq!(repeating.list_type().unwrap(), HalfInfiniteListType::Repeat);
        assert_eq!(
            repeating.highest_known_position().unwrap(),
            Some(ExtendedInteger::PositiveInfinity)
        );

        let callbacks = Arc::new(Mutex::new(Vec::new()));
        let callback_log = Arc::clone(&callbacks);
        let callback: InfiniteListNewValueCallback<InfiniteListError> =
            Arc::new(move |position, direction, list_type| {
                callback_log
                    .lock()
                    .unwrap()
                    .push((position.clone(), direction, list_type));
                Ok(())
            });
        let positioned: HalfInfiniteList<i64> = HalfInfiniteList::by_position(
            integer(0),
            InfiniteListDirection::Positive,
            |position| Ok(i64::try_from(position.clone()).unwrap().pow(2)),
            true,
            Some(callback),
        );
        assert_eq!(positioned.value(&integer(3)).unwrap(), 9);
        assert_eq!(positioned.stored_values().unwrap(), vec![0, 1, 4, 9]);
        assert_eq!(callbacks.lock().unwrap().len(), 4);

        let next: HalfInfiniteList<i64> = HalfInfiniteList::next_repeating(
            integer(0),
            InfiniteListDirection::Positive,
            -1_i64,
            |value| Ok((value + 1).rem_euclid(4)),
            None,
            Some(Arc::new(|_, _| Ok(()))),
        );
        assert_eq!(next.value(&integer(10)).unwrap(), 2);
        assert_eq!(
            next.stored_values().unwrap(),
            vec![0, 1, 2, 3, 0, 1, 2, 3, 0]
        );
    }

    #[test]
    fn full_list_operations_follow_qpa_parts_and_orientation() {
        let list = InfiniteList::new(
            integer(0),
            vec![100, 200],
            HalfInfiniteListSpecification::Repeat(vec![1, 2, 3]),
            HalfInfiniteListSpecification::Repeat(vec![7, 8]),
            None,
        );
        assert_eq!(
            values(&list, &[-5, -2, -1, 0, 1, 2, 5]).unwrap(),
            vec![7, 8, 7, 100, 200, 1, 1]
        );
        assert_eq!(list.middle_start(), &integer(0));
        assert_eq!(list.middle_end().unwrap(), integer(1));
        assert_eq!(
            list.upper_bound().unwrap(),
            ExtendedInteger::PositiveInfinity
        );
        assert_eq!(
            list.lower_bound().unwrap(),
            ExtendedInteger::NegativeInfinity
        );

        let shifted = list.shift(&integer(3)).unwrap();
        assert_eq!(shifted.middle_start(), &integer(-3));
        assert_eq!(
            values(&shifted, &[-3, -2, -1, 2]).unwrap(),
            vec![100, 200, 1, 1]
        );

        let positive = list.positive_part_from(&integer(3)).unwrap();
        assert_eq!(values(&positive, &[3, 4, 5, 6]).unwrap(), vec![2, 3, 1, 2]);
        let negative = list.negative_part_from(&integer(-2)).unwrap();
        assert_eq!(
            values(&negative, &[-5, -4, -3, -2]).unwrap(),
            vec![7, 8, 7, 8]
        );
    }

    #[test]
    fn finite_dispatch_failures_and_declared_only_surface_are_typed() {
        let finite = InfiniteList::finite(integer(0), vec![1, 2]);
        assert!(matches!(
            finite.shift(&integer(1)),
            Err(InfiniteListError::MissingPositivePart { operation: "Shift" })
        ));
        assert!(matches!(
            finite.map(|value| Ok(value * 2)),
            Err(InfiniteListError::MissingPositivePart {
                operation: "InfList"
            })
        ));
        assert!(matches!(
            shifted_splice(&finite, &integer(0), &finite, &integer(0), &integer(0)),
            Err(InfiniteListError::ShiftedSpliceNotInstalled)
        ));
        assert_eq!(
            finite.value(&integer(2)),
            Err(InfiniteListError::PositionOutsideInfiniteList {
                position: integer(2)
            })
        );
    }

    #[test]
    fn pinned_qpa_infinite_list_error_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/inflist_errors.schema-v1.json"
        ))
        .unwrap();
        let finite: InfiniteList<i64> = InfiniteList::finite(integer(0), vec![1_i64, 2]);
        let full: InfiniteList<i64> = InfiniteList::constant(1_i64);
        let positive: HalfInfiniteList<i64> =
            HalfInfiniteList::repeating(integer(1), InfiniteListDirection::Positive, vec![1]);
        let bad_positive: HalfInfiniteList<i64> =
            HalfInfiniteList::repeating(integer(2), InfiniteListDirection::Positive, vec![1]);
        let negative: HalfInfiniteList<i64> =
            HalfInfiniteList::repeating(integer(-1), InfiniteListDirection::Negative, vec![1]);
        let empty_repeat: HalfInfiniteList<i64> =
            HalfInfiniteList::repeating(integer(0), InfiniteListDirection::Positive, Vec::new());
        let empty: InfiniteList<i64> = InfiniteList::finite(integer(0), Vec::new());
        assert_eq!(
            json!({
                "shift_finite": finite.shift(&integer(1)).is_ok(),
                "map_finite": finite.map(|value| Ok(2 * value)).is_ok(),
                "empty_known": empty.lowest_known_position().is_ok(),
                "shifted_splice": shifted_splice(
                    &full, &integer(0), &full, &integer(0), &integer(0)).is_ok(),
                "illegal_half_position": positive.value(&integer(0)).is_ok(),
                "empty_repeat_value": empty_repeat.value(&integer(0)).is_ok(),
                "cut_outside": positive.cut(&integer(0)).is_ok(),
                "finite_outside": finite.value(&integer(2)).is_ok(),
                "bad_positive_start": InfiniteList::from_parts(
                    integer(0), vec![1], Some(bad_positive), None).is_ok(),
                "bad_negative_start": InfiniteList::from_parts(
                    integer(1), vec![1], None, Some(negative)).is_ok(),
            }),
            expected["success"]
        );
    }

    #[test]
    fn constructors_mapping_splicing_and_concatenation_match_source_order() {
        let integers = integers_list();
        assert_eq!(
            [-3, 0, 4]
                .into_iter()
                .map(|position| integers.value(&integer(position)).unwrap())
                .collect::<Vec<_>>(),
            vec![integer(-3), integer(0), integer(4)]
        );
        let constant = InfiniteList::constant(9_i64);
        assert_eq!(values(&constant, &[-10, 0, 10]).unwrap(), vec![9, 9, 9]);

        let left = InfiniteList::new(
            integer(0),
            vec![10, 11],
            HalfInfiniteListSpecification::Repeat(vec![12, 13]),
            HalfInfiniteListSpecification::Repeat(vec![8, 9]),
            None,
        );
        let right = InfiniteList::new(
            integer(0),
            vec![20, 21],
            HalfInfiniteListSpecification::Repeat(vec![22, 23]),
            HalfInfiniteListSpecification::Repeat(vec![18, 19]),
            None,
        );
        let splice = InfiniteList::splice(&left, &right, &integer(0)).unwrap();
        assert_eq!(
            values(&splice, &[-2, -1, 0, 1, 2]).unwrap(),
            vec![19, 18, 20, 11, 12]
        );

        let concatenated = InfiniteList::concatenate(&[left.clone(), right.clone()]).unwrap();
        assert_eq!(concatenated.middle_part().unwrap(), vec![20, 21, 10, 11]);
        let mapped = left.map(|value| Ok(value * 2)).unwrap();
        assert_eq!(
            values(&mapped, &[-1, 0, 1, 2]).unwrap(),
            vec![16, 20, 22, 24]
        );
    }
}
