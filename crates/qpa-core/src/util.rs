//! General utilities translated from QPA's `lib/util.gi`.
//!
//! These operations are not specific to quivers.  QPA implements
//! `PositionsNonzero` through GAP's generic `IsZero` predicate and implements
//! `NullList` as the single row of a `1 x length` zero matrix.

use std::num::NonZeroUsize;

use numerica::domains::Field;

/// Translate QPA's `PositionsNonzero(list)` operation.
///
/// The returned positions are one-based, as in GAP.  The predicate is the
/// Rust analogue of GAP's generic `IsZero`, so this operation remains usable
/// for scalars, matrices, algebra elements, and future QPA object types.
pub fn positions_nonzero<T>(list: &[T], mut is_zero: impl FnMut(&T) -> bool) -> Vec<usize> {
    list.iter()
        .enumerate()
        .filter_map(|(index, value)| (!is_zero(value)).then_some(index + 1))
        .collect()
}

/// Apply [`positions_nonzero`] to elements of a Numerica field.
#[must_use]
pub fn positions_nonzero_in_field<F: Field>(list: &[F::Element], field: &F) -> Vec<usize> {
    positions_nonzero(list, |value| field.is_zero(value))
}

/// Translate QPA's `NullList(length, field)` operation.
///
/// QPA installs this method only for a positive integer length, represented
/// here by [`NonZeroUsize`].  Extracting the sole row of GAP's zero matrix is
/// observably the same ordered list of independent field-zero values.
#[must_use]
pub fn null_list<F: Field>(length: NonZeroUsize, field: &F) -> Vec<F::Element> {
    vec![field.zero(); length.get()]
}

#[cfg(test)]
mod tests {
    use std::num::NonZeroUsize;

    use numerica::domains::{Ring, finite_field::Z2, rational::Q};
    use serde_json::{Value, json};

    use super::{null_list, positions_nonzero, positions_nonzero_in_field};

    #[test]
    fn pinned_qpa_util_fixture_matches() {
        let expected: Value =
            serde_json::from_str(include_str!("../../../tests/fixtures/util.schema-v1.json"))
                .unwrap();

        let rational = vec![Q.zero(), 2.into(), Q.zero(), (-3).into(), Q.zero()];
        assert_eq!(
            json!(positions_nonzero_in_field(&rational, &Q)),
            expected["rational"]["positions"]
        );
        let rational_null = null_list(NonZeroUsize::new(4).unwrap(), &Q);
        assert_eq!(
            json!(
                rational_null
                    .iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>()
            ),
            expected["rational"]["null_list"]
        );

        let binary = Z2;
        let binary_values = vec![
            binary.zero(),
            binary.one(),
            binary.zero(),
            binary.one(),
            binary.zero(),
        ];
        assert_eq!(
            json!(positions_nonzero_in_field(&binary_values, &binary)),
            expected["binary"]["positions"]
        );
        assert_eq!(
            null_list(NonZeroUsize::new(4).unwrap(), &binary),
            vec![binary.zero(); 4]
        );

        let matrices = [[[0]], [[1]], [[0]], [[-2]]];
        assert_eq!(
            json!(positions_nonzero(&matrices, |matrix| matrix[0][0] == 0)),
            expected["generic_matrices"]["positions"]
        );
    }
}
