use super::*;
use crate::{ArrowSpec, GaloisField, PathAlgebra, Quiver, QuotientPathAlgebra, VertexIndex};
use numerica::domains::{
    finite_field::{Z2, Zp},
    rational::Q,
};

fn compare<F: Field>(field: &F) {
    let vertex = |i| VertexIndex::new(i).unwrap();
    let parent = PathAlgebra::new(
        field.clone(),
        Quiver::from_indexed(
            2,
            vec![
                ArrowSpec::named(vertex(1), vertex(1), "x"),
                ArrowSpec::named(vertex(1), vertex(2), "y"),
                ArrowSpec::named(vertex(2), vertex(2), "z"),
            ],
        )
        .unwrap(),
    )
    .unwrap();
    let arrows = parent
        .quiver()
        .arrows()
        .map(|a| parent.element(a.path()).unwrap())
        .collect::<Vec<_>>();
    let x2 = arrows[0].multiply(&arrows[0]).unwrap();
    let x3 = x2.multiply(&arrows[0]).unwrap();
    let yz = arrows[1].multiply(&arrows[2]).unwrap();
    let relation = x2
        .multiply(&arrows[1])
        .unwrap()
        .sub(&arrows[0].multiply(&yz).unwrap())
        .unwrap();
    for nonuniform in [false, true] {
        let mut relations = vec![x3.clone(), relation.clone()];
        if nonuniform {
            relations.push(arrows[0].add(&arrows[1]).unwrap());
        }
        let algebra = QuiverAlgebra::Quotient(
            QuotientPathAlgebra::by_ideal(
                &crate::PathAlgebraIdeal::new(parent.clone(), relations).unwrap(),
            )
            .unwrap(),
        );
        for left in 0..=3 {
            for right in 0..=3 {
                for seed in 0..6 {
                    let dimensions = [left, right];
                    let matrices = [(left, left), (left, right), (right, right)]
                        .into_iter()
                        .enumerate()
                        .map(|(arrow, (rows, columns))| {
                            let mut matrix =
                                DenseMatrix::zero(rows.max(1), columns.max(1), field.clone());
                            for i in 0..rows {
                                for j in 0..columns {
                                    if (i + 3 * j + seed + arrow) % 3 == 0
                                        && (arrow == 1 || j > i || seed == 5)
                                    {
                                        matrix[(i, j)] = field.one();
                                    }
                                }
                            }
                            matrix
                        })
                        .collect::<Vec<_>>();
                    assert_eq!(
                        relations_are_satisfied(&algebra, &dimensions, &matrices, field),
                        relations_are_satisfied_uncached(&algebra, &dimensions, &matrices, field)
                    );
                }
            }
        }
    }
}

#[test]
fn bounded_relation_cache_matches_original_including_zero_spaces_and_nonuniform_relations() {
    compare(&Q);
    compare(&Z2);
    compare(&Zp::new(101));
    compare(&GaloisField::conway(2, 2).unwrap());
}

#[test]
fn zero_products_retain_nc_field_errors() {
    let field = Zp::new(7);
    let vertex = VertexIndex::new(1).unwrap();
    let parent = PathAlgebra::new(
        field.clone(),
        Quiver::from_indexed(
            1,
            vec![
                ArrowSpec::named(vertex, vertex, "x"),
                ArrowSpec::named(vertex, vertex, "y"),
            ],
        )
        .unwrap(),
    )
    .unwrap();
    let generators = parent.generators();
    let relation = generators[1].multiply(&generators[2]).unwrap();
    let algebra = QuiverAlgebra::Quotient(
        QuotientPathAlgebra::by_relations_high_level(parent, &[relation]).unwrap(),
    );
    let matrices = [
        DenseMatrix::zero(2, 2, field.clone()),
        DenseMatrix::zero(2, 2, Zp::new(11)),
    ];
    assert_eq!(
        relations_are_satisfied(&algebra, &[2], &matrices, &field),
        Err(MatrixError::DifferentFields)
    );
    assert_eq!(
        relations_are_satisfied_uncached(&algebra, &[2], &matrices, &field),
        Err(MatrixError::DifferentFields)
    );
}

#[test]
fn nonuniform_zero_space_direct_sum_keeps_legacy_shape_failure() {
    let vertex = |i| VertexIndex::new(i).unwrap();
    let parent = PathAlgebra::new(
        Q,
        Quiver::from_indexed(
            2,
            vec![
                ArrowSpec::named(vertex(1), vertex(1), "x"),
                ArrowSpec::named(vertex(2), vertex(2), "y"),
            ],
        )
        .unwrap(),
    )
    .unwrap();
    let generators = parent.generators();
    let relation = generators[2].sub(&generators[3]).unwrap();
    let algebra = QuiverAlgebra::Quotient(
        QuotientPathAlgebra::by_ideal(
            &crate::PathAlgebraIdeal::new(parent, vec![relation]).unwrap(),
        )
        .unwrap(),
    );
    // QPA's physical stand-ins make 0_(1x1) - 0_(1x1) valid here. The
    // direct sum has shapes 2x2 and 1x1, so it must NOT inherit certification.
    let module = right_module_over_path_algebra(algebra, vec![1, 0], &[])
        .unwrap()
        .unwrap();
    assert!(direct_sum_of_qpa_modules(&[module.clone(), module]).is_err());
}
