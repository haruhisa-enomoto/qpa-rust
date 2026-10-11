use numerica::domains::{
    Field,
    finite_field::{Z2, Zp},
    rational::Q,
};
use qpa_core::{ArrowSpec, PathAlgebra, Quiver, VertexIndex, make_uniform_list};

use super::{component_ids, precondition};
use crate::qpa_conversion::{cohen_to_path, gbnp_groebner_basis, path_to_cohen};
use crate::sgrobner;

fn compare<F: Field>(field: &F) {
    for swap in [false, true] {
        for shift in 0..5 {
            let left = VertexIndex::new(if swap { 2 } else { 1 }).unwrap();
            let right = VertexIndex::new(if swap { 1 } else { 2 }).unwrap();
            let mut arrows = vec![
                ArrowSpec::named(left, left, "p"),
                ArrowSpec::named(left, left, "q"),
                ArrowSpec::named(left, right, "bridge"),
                ArrowSpec::named(right, right, "r"),
                ArrowSpec::named(right, right, "s"),
            ];
            arrows.rotate_left(shift);
            let parent =
                PathAlgebra::new(field.clone(), Quiver::from_indexed(3, arrows).unwrap()).unwrap();
            let element = |name| {
                parent
                    .element(
                        parent
                            .quiver()
                            .arrows()
                            .find(|a| a.name() == name)
                            .unwrap()
                            .path(),
                    )
                    .unwrap()
            };
            let p = element("p");
            let q = element("q");
            let bridge = element("bridge");
            let r = element("r");
            let s = element("s");
            let two = field.add(&field.one(), &field.one());
            let mut relations = vec![
                // Nonuniform AND nonhomogeneous: split corners before completion.
                p.pow(2)
                    .unwrap()
                    .sub(&q.pow(3).unwrap())
                    .unwrap()
                    .add(&r.pow(2).unwrap())
                    .unwrap(),
                q.pow(4).unwrap(),
                s.pow(4).unwrap(),
                p.multiply(&q)
                    .unwrap()
                    .sub(&q.multiply(&p).unwrap().scale_left(&two).unwrap())
                    .unwrap(),
                r.multiply(&s)
                    .unwrap()
                    .sub(&s.multiply(&r).unwrap())
                    .unwrap(),
                p.multiply(&bridge)
                    .unwrap()
                    .sub(&bridge.multiply(&r).unwrap())
                    .unwrap(),
                q.multiply(&bridge)
                    .unwrap()
                    .sub(&bridge.multiply(&s).unwrap())
                    .unwrap(),
                parent.zero(),
            ];
            relations.rotate_left(shift);
            let uniform = make_uniform_list(&relations, &parent).unwrap();
            assert!(precondition(&uniform, &parent).unwrap().is_some());
            // Frozen original checked-adapter algorithm, not the new helper.
            let raw = path_to_cohen(&uniform).unwrap();
            let input = raw.iter().map(|p| p.clean().unwrap()).collect::<Vec<_>>();
            let expected = cohen_to_path(&sgrobner(&input).unwrap(), &parent).unwrap();
            assert_eq!(gbnp_groebner_basis(&relations, &parent).unwrap(), expected);
        }
    }
}

#[test]
fn structural_completion_matches_ordered_unpreconditioned_basis() {
    compare(&Q);
    compare(&Z2);
    // The p^2=q^3, pq=2qp critical pair changes in characteristic seven.
    compare(&Zp::new(7));
    compare(&Zp::new(101));
}

#[test]
fn components_use_reachability_not_labels_or_weak_connectivity() {
    let v = |i| VertexIndex::new(i).unwrap();
    let quiver = Quiver::from_indexed(
        6,
        vec![
            ArrowSpec::named(v(1), v(2), "a"),
            ArrowSpec::named(v(2), v(1), "b"),
            ArrowSpec::named(v(2), v(3), "c"),
            ArrowSpec::named(v(3), v(4), "d"),
            ArrowSpec::named(v(4), v(3), "e"),
            ArrowSpec::named(v(5), v(5), "f"),
        ],
    )
    .unwrap();
    let ids = component_ids(&quiver);
    assert_eq!(ids[0], ids[1]);
    assert_eq!(ids[2], ids[3]);
    assert_ne!(ids[1], ids[2]);
    assert_ne!(ids[4], ids[5]);
    assert_ne!(ids[3], ids[4]);
    let isolated = component_ids(&Quiver::from_indexed(2, vec![]).unwrap());
    assert_ne!(isolated[0], isolated[1]);
}

#[test]
fn no_split_preserves_existing_completion_path() {
    let v = |i| VertexIndex::new(i).unwrap();
    let quiver = Quiver::from_indexed(
        2,
        vec![
            ArrowSpec::named(v(1), v(2), "a"),
            ArrowSpec::named(v(2), v(1), "b"),
        ],
    )
    .unwrap();
    let parent = PathAlgebra::new(Q, quiver).unwrap();
    let arrows = parent
        .quiver()
        .arrows()
        .map(|a| parent.element(a.path()).unwrap())
        .collect::<Vec<_>>();
    assert!(precondition(&arrows, &parent).unwrap().is_none());
    assert!(precondition(&[], &parent).unwrap().is_none());
    let quiver = Quiver::from_indexed(2, vec![ArrowSpec::named(v(1), v(2), "renamed")]).unwrap();
    let parent = PathAlgebra::new(Q, quiver).unwrap();
    let relation = parent
        .element(parent.quiver().arrows().next().unwrap().path())
        .unwrap();
    assert!(precondition(&[relation], &parent).unwrap().is_none());
}
