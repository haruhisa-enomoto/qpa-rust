//! QPA's combinatorial-map construction for gentle algebras.
//!
//! The algorithms here are direct translations of
//! `MaximalPathsOfGentleAlgebra` and `CombinatorialMapOfGentleAlgebra` in
//! `lib/combinatorialmap.gi`.  Ideal membership continues to use QPA's
//! configured GBNP-backed quotient machinery.

use numerica::domains::{Field, integer::Integer, rational::Rational};
use qpa_core::{
    Arrow, CombinatorialMap, CombinatorialMapError, Path, PathAlgebraError, PathError, Permutation,
    PermutationError, QpaCoefficientField, QuiverAlgebra, QuiverError, QuotientError,
};
use thiserror::Error;

use crate::pathalg::{PathalgGbnpError, path_algebra_is_gentle, quotient_is_gentle};

/// Translate QPA's `MaximalPathsOfGentleAlgebra` attribute.
///
/// Only positive-length paths are returned.  Results are cached on the
/// algebra handle, preserving GAP attribute identity and source arrow order.
///
/// # Errors
///
/// Returns an error when the algebra is not gentle and propagates QPA's
/// GBNP-backed gentle test, ideal membership, incidence, and path failures.
pub fn maximal_paths_of_gentle_algebra<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<&[Path], GentleCombinatorialMapError> {
    if let Some(paths) = algebra.cached_gentle_maximal_paths() {
        return Ok(paths);
    }
    require_gentle(algebra)?;

    let quiver = algebra.quiver();
    let mut maximal_paths = Vec::new();
    for arrow in quiver.arrows() {
        let mut walk = vec![arrow];
        loop {
            let first = walk
                .first()
                .cloned()
                .ok_or(GentleCombinatorialMapError::InvalidMaximalPathState)?;
            let last = walk
                .last()
                .cloned()
                .ok_or(GentleCombinatorialMapError::InvalidMaximalPathState)?;
            let mut finished = true;

            for previous in quiver.incoming_arrows(&first.source())? {
                if !arrow_product_is_zero(algebra, &previous, &first)? {
                    ensure_new_arrow(&walk, &previous)?;
                    walk.insert(0, previous);
                    finished = false;
                }
            }
            for next in quiver.outgoing_arrows(&last.target())? {
                if !arrow_product_is_zero(algebra, &last, &next)? {
                    ensure_new_arrow(&walk, &next)?;
                    walk.push(next);
                    finished = false;
                }
            }
            if finished {
                break;
            }
        }

        let path = Path::from_arrow_list(quiver, &walk)?;
        if !maximal_paths.contains(&path) {
            maximal_paths.push(path);
        }
    }
    Ok(algebra.cache_gentle_maximal_paths(maximal_paths))
}

/// Translate QPA's `CombinatorialMapOfGentleAlgebra` attribute.
///
/// QPA labels the two half-edges at vertex `i` by `2*i-1` and `2*i`, then
/// applies its transposition product for every maximal path.  Marked points,
/// including duplicates, retain source order.
///
/// # Errors
///
/// Propagates maximal-path and permutation construction failures.  In
/// particular, QPA's radical-square-zero loop attempts the invalid GAP cycle
/// `(1,1)`; Rust retains that failure as a typed permutation error.
pub fn combinatorial_map_of_gentle_algebra<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<&CombinatorialMap, GentleCombinatorialMapError> {
    if let Some(map) = algebra.cached_gentle_combinatorial_map() {
        return Ok(map);
    }
    let maximal_paths = maximal_paths_of_gentle_algebra(algebra)?;
    let quiver = algebra.quiver();
    let vertex_count = quiver.number_of_vertices();
    let size = vertex_count
        .checked_mul(2)
        .ok_or(GentleCombinatorialMapError::MapSizeOverflow)?;

    let pairing = Permutation::from_images(
        (1..=vertex_count)
            .flat_map(|index| [2 * index, 2 * index - 1])
            .collect(),
    )?;
    let mut ordering = Permutation::identity();
    let mut marked = Vec::new();

    for path in maximal_paths {
        let walk = path
            .walk()
            .ok_or(GentleCombinatorialMapError::InvalidMaximalPathState)?;
        let (&first_index, _) = walk
            .split_first()
            .ok_or(GentleCombinatorialMapError::InvalidMaximalPathState)?;
        let first_arrow = quiver
            .arrow(first_index)
            .ok_or(GentleCombinatorialMapError::InvalidMaximalPathState)?;
        let start_odd = 2 * first_arrow.source().index().get() - 1;
        let start_index = if ordering.image(start_odd) == start_odd {
            start_odd
        } else {
            start_odd + 1
        };

        for &arrow_index in walk {
            let arrow = quiver
                .arrow(arrow_index)
                .ok_or(GentleCombinatorialMapError::InvalidMaximalPathState)?;
            let target_odd = 2 * arrow.target().index().get() - 1;
            let target_index = if ordering.image(target_odd) == target_odd {
                target_odd
            } else {
                target_odd + 1
            };
            let cycle = [start_index, target_index];
            let transposition = Permutation::from_cycles(size, &[&cycle])?;
            ordering = ordering.gap_product(&transposition);
        }
        marked.push(ordering.inverse().image(start_index));
    }

    marked.extend((1..=size).filter(|&point| ordering.image(point) == point));
    let map = CombinatorialMap::new(size, ordering, pairing, marked)?;
    Ok(algebra.cache_gentle_combinatorial_map(map))
}

/// Translate QPA's derived-equivalence classifier for two gentle algebras.
///
/// The comparison follows the source branch order: coefficient field, genus,
/// boundary and vertex counts, boundary winding data, the genus-one gcd, and
/// finally the higher-genus parity/Arf invariant.  This includes QPA's
/// asymmetric higher-genus parity loops exactly as installed in version 1.37.
///
/// # Errors
///
/// Returns an error if either algebra is not gentle, if construction of its
/// combinatorial map or topological invariants fails, or if QPA's final
/// `RemInt` input would not be integral.
pub fn are_derived_equivalent<F: QpaCoefficientField, G: QpaCoefficientField>(
    left: &QuiverAlgebra<F>,
    right: &QuiverAlgebra<G>,
) -> Result<bool, GentleCombinatorialMapError> {
    require_gentle(left)?;
    require_gentle(right)?;
    if left.field().qpa_characteristic() != right.field().qpa_characteristic() {
        return Ok(false);
    }

    let left_map = combinatorial_map_of_gentle_algebra(left)?;
    let right_map = combinatorial_map_of_gentle_algebra(right)?;
    compare_gentle_maps(left_map, right_map)
}

fn compare_gentle_maps(
    left: &CombinatorialMap,
    right: &CombinatorialMap,
) -> Result<bool, GentleCombinatorialMapError> {
    if left.genus() != right.genus()
        || left.faces().len() != right.faces().len()
        || vertex_count(left) != vertex_count(right)
    {
        return Ok(false);
    }
    if left.genus() == &Rational::from(0) && left.faces().len() == 1 {
        return Ok(true);
    }

    let left_boundary = boundary_winding_data(left)?;
    let right_boundary = boundary_winding_data(right)?;
    if !right_boundary
        .iter()
        .all(|entry| left_boundary.contains(entry))
    {
        return Ok(false);
    }
    if left.genus() == &Rational::from(0) {
        return Ok(true);
    }

    let left_boundary_winding = left_boundary
        .iter()
        .map(|&(winding, _)| winding)
        .collect::<Vec<_>>();
    let right_boundary_winding = right_boundary
        .iter()
        .map(|&(winding, _)| winding)
        .collect::<Vec<_>>();
    let left_basis_winding = basis_winding_data(left)?;
    let right_basis_winding = basis_winding_data(right)?;

    if left.genus() == &Rational::from(1) {
        return Ok(genus_one_gcd(&left_boundary_winding, &left_basis_winding)
            == genus_one_gcd(&right_boundary_winding, &right_basis_winding));
    }

    let mut left_odd = left_boundary_winding.iter().any(|winding| winding % 2 != 0);
    let right_odd = right_boundary_winding
        .iter()
        .any(|winding| winding % 2 != 0);
    let left_mod_four = left_boundary_winding.iter().any(|winding| winding % 4 == 0);
    if left_basis_winding
        .iter()
        .any(|&(a, b)| a % 2 != 0 || b % 2 != 0)
    {
        left_odd = true;
    }

    // QPA 1.37 repeats the left-basis loop here instead of inspecting the
    // right basis.  Preserve that installed computation rather than repairing
    // the apparent typo.
    if left_basis_winding
        .iter()
        .any(|&(a, b)| a % 2 != 0 || b % 2 != 0)
    {
        left_odd = true;
    }
    if (left_odd && right_odd) || left_mod_four {
        return Ok(true);
    }

    let arf_difference = arf_sum(&left_basis_winding) - arf_sum(&right_basis_winding);
    if arf_difference.denominator_ref() != &Integer::one() {
        return Err(GentleCombinatorialMapError::NonIntegralArfDifference);
    }
    Ok(arf_difference.numerator_ref() % 2_i64 == Integer::zero())
}

fn vertex_count(map: &CombinatorialMap) -> usize {
    map.ordering()
        .orbits_on(&(1..=map.size()).collect::<Vec<_>>())
        .len()
}

fn boundary_winding_data(
    map: &CombinatorialMap,
) -> Result<Vec<(i64, usize)>, CombinatorialMapError> {
    map.boundary_curves()?
        .iter()
        .map(|(curve, marked)| Ok((-map.winding_number(curve)?, *marked)))
        .collect()
}

fn basis_winding_data(map: &CombinatorialMap) -> Result<Vec<(i64, i64)>, CombinatorialMapError> {
    map.homology_basis()?
        .iter()
        .map(|(a, b)| Ok((map.winding_number(a)?, map.winding_number(b)?)))
        .collect()
}

fn genus_one_gcd(boundary: &[i64], basis: &[(i64, i64)]) -> Integer {
    let mut gcd = Integer::zero();
    for winding in boundary {
        gcd = gcd.gcd(&(Integer::from(*winding) + Integer::from(2)));
    }
    for &(a, b) in basis {
        gcd = gcd.gcd(&Integer::from(a));
        gcd = gcd.gcd(&Integer::from(b));
    }
    gcd.abs()
}

fn arf_sum(basis: &[(i64, i64)]) -> Rational {
    basis.iter().fold(Rational::from(0), |sum, &(a, b)| {
        let left = Rational::from(a) / Rational::from(2) + Rational::from(1);
        let right = Rational::from(b) / Rational::from(2) + Rational::from(1);
        sum + left * right
    })
}

fn require_gentle<F: Field>(algebra: &QuiverAlgebra<F>) -> Result<(), GentleCombinatorialMapError> {
    let gentle = match algebra {
        QuiverAlgebra::Path(path_algebra) => path_algebra_is_gentle(path_algebra),
        QuiverAlgebra::Quotient(quotient) => quotient_is_gentle(quotient)?,
    };
    if gentle {
        Ok(())
    } else {
        Err(GentleCombinatorialMapError::NonGentleAlgebra)
    }
}

fn arrow_product_is_zero<F: Field>(
    algebra: &QuiverAlgebra<F>,
    left: &Arrow,
    right: &Arrow,
) -> Result<bool, GentleCombinatorialMapError> {
    let QuiverAlgebra::Quotient(quotient) = algebra else {
        return Ok(false);
    };
    let parent = quotient.original_path_algebra();
    let product = left.path().multiply(&right.path())?;
    Ok(quotient.ideal().contains(parent.element(product)?)?)
}

fn ensure_new_arrow(walk: &[Arrow], extension: &Arrow) -> Result<(), GentleCombinatorialMapError> {
    if walk.contains(extension) {
        Err(GentleCombinatorialMapError::InfiniteMaximalPath)
    } else {
        Ok(())
    }
}

/// Failure in QPA's combinatorial-map construction for a gentle algebra.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum GentleCombinatorialMapError {
    /// The GAP method has no applicable gentle-algebra dispatch.
    #[error("the entered algebra is not gentle")]
    NonGentleAlgebra,
    /// A repeated arrow would make QPA's extension loop nonterminating.
    #[error("QPA's maximal gentle path extension does not terminate")]
    InfiniteMaximalPath,
    /// A cached or constructed maximal path violated QPA's positive walk invariant.
    #[error("QPA's maximal gentle path has no valid positive arrow walk")]
    InvalidMaximalPathState,
    /// Twice the number of vertices does not fit Rust's index type.
    #[error("the gentle combinatorial-map size overflows usize")]
    MapSizeOverflow,
    /// QPA's final `RemInt` would receive a nonintegral Arf difference.
    #[error("the derived-equivalence Arf difference is not integral")]
    NonIntegralArfDifference,
    /// QPA's GBNP-backed algebra-property computation failed.
    #[error(transparent)]
    Pathalg(#[from] PathalgGbnpError),
    /// QPA path construction or multiplication failed.
    #[error(transparent)]
    Path(#[from] PathError),
    /// Parent path-algebra embedding failed.
    #[error(transparent)]
    PathAlgebra(#[from] PathAlgebraError),
    /// Quotient ideal membership failed.
    #[error(transparent)]
    Quotient(#[from] QuotientError),
    /// Quiver incidence lookup failed.
    #[error(transparent)]
    Quiver(#[from] QuiverError),
    /// A QPA-generated permutation was invalid.
    #[error(transparent)]
    Permutation(#[from] PermutationError),
    /// The final QPA combinatorial-map constructor failed.
    #[error(transparent)]
    CombinatorialMap(#[from] CombinatorialMapError),
}

#[cfg(test)]
mod tests {
    use numerica::domains::{Field, finite_field::Z2, rational::Q};
    use qpa_core::{
        ArrowSpec, CombinatorialMap, PathAlgebra, PermutationError, Quiver, QuiverAlgebra,
        VertexIndex,
    };
    use serde_json::{Value, json};

    use super::{
        GentleCombinatorialMapError, are_derived_equivalent, basis_winding_data,
        boundary_winding_data, combinatorial_map_of_gentle_algebra, genus_one_gcd,
        maximal_paths_of_gentle_algebra, vertex_count,
    };
    use crate::pathalg::quotient_by_relations;

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn map_json(map: &CombinatorialMap) -> Value {
        let genus = map.genus();
        json!({
            "size": map.size(),
            "ordering": map.ordering().images_on(map.size()).unwrap(),
            "pairing": map.pairing().images_on(map.size()).unwrap(),
            "marked": map.marked_edges(),
            "faces": map.faces(),
            "genus": [
                genus.numerator().to_string().parse::<i64>().unwrap(),
                genus.denominator().to_string().parse::<i64>().unwrap()
            ],
        })
    }

    fn paths_json<F: Field>(algebra: &QuiverAlgebra<F>) -> Value {
        json!(
            maximal_paths_of_gentle_algebra(algebra)
                .unwrap()
                .iter()
                .map(|path| {
                    path.walk()
                        .unwrap()
                        .iter()
                        .map(|&index| algebra.quiver().arrow(index).unwrap().name().to_owned())
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        )
    }

    fn case_json<F: Field>(algebra: &QuiverAlgebra<F>) -> Value {
        json!({
            "maximal_paths": paths_json(algebra),
            "map": map_json(combinatorial_map_of_gentle_algebra(algebra).unwrap()),
        })
    }

    fn field_cases<F: Field>(field: F) -> Value {
        let chain_quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(2), vertex(3), "b"),
            ],
        )
        .unwrap();
        let chain = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), chain_quiver).unwrap());

        let discrete = QuiverAlgebra::Path(
            PathAlgebra::new(field.clone(), Quiver::from_indexed(2, vec![]).unwrap()).unwrap(),
        );

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let loop_parent = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
        let x = loop_parent.generators()[1].clone();
        let square = QuiverAlgebra::Quotient(
            quotient_by_relations(loop_parent, &[x.pow(2).unwrap()]).unwrap(),
        );

        let branch_quiver = Quiver::from_indexed(
            4,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(2), vertex(3), "b"),
                ArrowSpec::named(vertex(2), vertex(4), "c"),
            ],
        )
        .unwrap();
        let branch_parent = PathAlgebra::new(field, branch_quiver).unwrap();
        let generators = branch_parent.generators();
        let relation = generators[4].multiply(&generators[5]).unwrap();
        let branch =
            QuiverAlgebra::Quotient(quotient_by_relations(branch_parent, &[relation]).unwrap());

        json!({
            "chain": case_json(&chain),
            "discrete": case_json(&discrete),
            "square_loop": {"maximal_paths": paths_json(&square)},
            "branch": case_json(&branch),
        })
    }

    fn quotient_case<F: Field>(
        field: F,
        vertices: usize,
        arrows: Vec<ArrowSpec>,
        relation_pairs: &[(usize, usize)],
    ) -> QuiverAlgebra<F> {
        let parent =
            PathAlgebra::new(field, Quiver::from_indexed(vertices, arrows).unwrap()).unwrap();
        let generators = parent.generators();
        let relations = relation_pairs
            .iter()
            .map(|&(left, right)| {
                generators[vertices + left]
                    .multiply(&generators[vertices + right])
                    .unwrap()
            })
            .collect::<Vec<_>>();
        QuiverAlgebra::Quotient(quotient_by_relations(parent, &relations).unwrap())
    }

    fn genus_one_a<F: Field>(field: F) -> QuiverAlgebra<F> {
        quotient_case(
            field,
            6,
            vec![
                ArrowSpec::named(vertex(2), vertex(1), "aa"),
                ArrowSpec::named(vertex(3), vertex(2), "bb"),
                ArrowSpec::named(vertex(1), vertex(4), "cc"),
                ArrowSpec::named(vertex(4), vertex(2), "dd"),
                ArrowSpec::named(vertex(3), vertex(5), "ee"),
                ArrowSpec::named(vertex(5), vertex(4), "ff"),
                ArrowSpec::named(vertex(4), vertex(6), "gg"),
                ArrowSpec::named(vertex(6), vertex(5), "hh"),
                ArrowSpec::named(vertex(5), vertex(6), "ii"),
            ],
            &[(0, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 8)],
        )
    }

    fn genus_one_b<F: Field>(field: F) -> QuiverAlgebra<F> {
        quotient_case(
            field,
            6,
            vec![
                ArrowSpec::named(vertex(2), vertex(1), "aaa"),
                ArrowSpec::named(vertex(1), vertex(3), "bbb"),
                ArrowSpec::named(vertex(3), vertex(2), "ccc"),
                ArrowSpec::named(vertex(1), vertex(5), "ddd"),
                ArrowSpec::named(vertex(2), vertex(6), "eee"),
                ArrowSpec::named(vertex(3), vertex(4), "fff"),
                ArrowSpec::named(vertex(4), vertex(5), "ggg"),
                ArrowSpec::named(vertex(5), vertex(6), "hhh"),
                ArrowSpec::named(vertex(6), vertex(4), "iii"),
            ],
            &[(0, 1), (1, 2), (2, 0), (5, 6), (6, 7), (7, 8)],
        )
    }

    fn chain<F: Field>(field: F) -> QuiverAlgebra<F> {
        QuiverAlgebra::Path(
            PathAlgebra::new(
                field,
                Quiver::from_indexed(
                    3,
                    vec![
                        ArrowSpec::named(vertex(1), vertex(2), "chain_a"),
                        ArrowSpec::named(vertex(2), vertex(3), "chain_b"),
                    ],
                )
                .unwrap(),
            )
            .unwrap(),
        )
    }

    fn branch<F: Field>(field: F) -> QuiverAlgebra<F> {
        quotient_case(
            field,
            4,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "branch_a"),
                ArrowSpec::named(vertex(2), vertex(3), "branch_b"),
                ArrowSpec::named(vertex(2), vertex(4), "branch_c"),
            ],
            &[(0, 1)],
        )
    }

    fn derived_data<F: Field>(algebra: &QuiverAlgebra<F>) -> Value {
        let map = combinatorial_map_of_gentle_algebra(algebra).unwrap();
        let boundary = boundary_winding_data(map).unwrap();
        let basis_winding = basis_winding_data(map).unwrap();
        let gcd = if map.genus() == &numerica::domains::rational::Rational::from(1) {
            Some(
                genus_one_gcd(
                    &boundary
                        .iter()
                        .map(|&(winding, _)| winding)
                        .collect::<Vec<_>>(),
                    &basis_winding,
                )
                .to_i64()
                .unwrap(),
            )
        } else {
            None
        };
        json!({
            "genus": [
                map.genus().numerator_ref().to_i64().unwrap(),
                map.genus().denominator_ref().to_i64().unwrap()
            ],
            "face_count": map.faces().len(),
            "vertex_count": vertex_count(map),
            "boundary": boundary,
            "basis_winding": basis_winding,
            "gcd": gcd,
        })
    }

    #[test]
    fn pinned_qpa_gentle_combinatorial_map_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/combinatorialmap_gentle.schema-v1.json"
        ))
        .unwrap();
        assert_eq!(field_cases(Q), expected["rational"]);
        assert_eq!(field_cases(Z2), expected["f2"]);
    }

    #[test]
    fn pinned_qpa_derived_equivalence_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/combinatorialmap_derived.schema-v1.json"
        ))
        .unwrap();
        let genus_one_a = genus_one_a(Q);
        let genus_one_b = genus_one_b(Q);
        let chain_q = chain(Q);
        let chain_q_copy = chain(Q);
        let chain_f2 = chain(Z2);
        let branch_q = branch(Q);

        combinatorial_map_of_gentle_algebra(&genus_one_a).expect("genus-one A map");
        combinatorial_map_of_gentle_algebra(&genus_one_b).expect("genus-one B map");
        assert_eq!(derived_data(&genus_one_a), expected["genus_one_a"]);
        assert_eq!(derived_data(&genus_one_b), expected["genus_one_b"]);
        assert_eq!(
            json!({
                "source_genus_one": are_derived_equivalent(&genus_one_a, &genus_one_b).unwrap(),
                "reflexive_genus_one": are_derived_equivalent(&genus_one_a, &genus_one_a).unwrap(),
                "genus_zero_disc": are_derived_equivalent(&chain_q, &chain_q_copy).unwrap(),
                "different_field": are_derived_equivalent(&chain_q, &chain_f2).unwrap(),
                "different_vertex_count": are_derived_equivalent(&chain_q, &branch_q).unwrap(),
                "different_genus": are_derived_equivalent(&chain_q, &genus_one_a).unwrap(),
            }),
            expected["results"]
        );
    }

    #[test]
    fn gentle_attributes_cache_identity_and_preserve_qpa_failures() {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
        let first_paths = maximal_paths_of_gentle_algebra(&algebra).unwrap();
        let second_paths = maximal_paths_of_gentle_algebra(&algebra).unwrap();
        assert!(std::ptr::eq(first_paths, second_paths));
        let first_map = combinatorial_map_of_gentle_algebra(&algebra).unwrap();
        assert!(first_map.is_identical(combinatorial_map_of_gentle_algebra(&algebra).unwrap()));

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let parent = PathAlgebra::new(Q, loop_quiver).unwrap();
        let x = parent.generators()[1].clone();
        let square = QuiverAlgebra::Quotient(
            quotient_by_relations(parent.clone(), &[x.pow(2).unwrap()]).unwrap(),
        );
        assert!(matches!(
            combinatorial_map_of_gentle_algebra(&square),
            Err(GentleCombinatorialMapError::Permutation(
                PermutationError::RepeatedCyclePoint(1)
            ))
        ));
        let cubic =
            QuiverAlgebra::Quotient(quotient_by_relations(parent, &[x.pow(3).unwrap()]).unwrap());
        assert!(matches!(
            are_derived_equivalent(&cubic, &algebra),
            Err(GentleCombinatorialMapError::NonGentleAlgebra)
        ));
        assert!(matches!(
            maximal_paths_of_gentle_algebra(&cubic),
            Err(GentleCombinatorialMapError::NonGentleAlgebra)
        ));
    }
}
