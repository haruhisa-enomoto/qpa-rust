use numerica::domains::{
    Field,
    finite_field::{Z2, Zp},
    rational::Q,
};
use qpa_core::{
    ArrowSpec, GaloisField, PathAlgebra, Quiver, QuiverAlgebra, VertexIndex,
    direct_sum_of_qpa_modules, hom_over_algebra, zero_module,
};

use super::{
    EnvelopingModuleError, ExtDimensionSource, ext_dimension_from_hom_dimensions, ext_over_algebra,
    indec_projective_modules, projective_hom_dimension_from_multiplicities, simple_modules,
};

#[test]
fn ext_dimension_arithmetic_is_checked_without_intermediate_sum_overflow() {
    assert_eq!(ext_dimension_from_hom_dimensions(5, 9, 6).unwrap(), 2);
    assert_eq!(
        ext_dimension_from_hom_dimensions(usize::MAX, usize::MAX, usize::MAX).unwrap(),
        usize::MAX
    );
    for (syzygy, projective, source) in [(1, 2, 3), (4, 9, 1)] {
        assert!(matches!(
            ext_dimension_from_hom_dimensions(syzygy, projective, source),
            Err(EnvelopingModuleError::InternalExtHomRejected)
        ));
    }
    assert_eq!(
        projective_hom_dimension_from_multiplicities(&[2, 0, 3], &[7, 99, 5]).unwrap(),
        29
    );
    assert_eq!(
        projective_hom_dimension_from_multiplicities(&[usize::MAX, 0], &[0, usize::MAX]).unwrap(),
        0
    );
    for (counts, dims) in [
        (vec![usize::MAX], vec![2]),
        (vec![usize::MAX, 1], vec![1, 1]),
    ] {
        assert!(matches!(
            projective_hom_dimension_from_multiplicities(&counts, &dims),
            Err(EnvelopingModuleError::Module(
                qpa_core::PathAlgebraModuleError::DimensionOverflow
            ))
        ));
    }
}

fn compare<F: Field>(field: F) {
    let vertex = |i| VertexIndex::new(i).unwrap();
    let quiver = Quiver::from_indexed(
        3,
        vec![
            ArrowSpec::named(vertex(1), vertex(2), "a"),
            ArrowSpec::named(vertex(2), vertex(3), "b"),
        ],
    )
    .unwrap();
    let parent = PathAlgebra::new(field, quiver).unwrap();
    let arrows = parent
        .quiver()
        .arrows()
        .map(|a| parent.element(a.path()).unwrap())
        .collect::<Vec<_>>();
    let algebra = QuiverAlgebra::Quotient(
        crate::quotient_by_relations(parent, &[arrows[0].multiply(&arrows[1]).unwrap()]).unwrap(),
    );
    let simples = simple_modules(&algebra).unwrap();
    let projectives = indec_projective_modules(&algebra).unwrap();
    let source = direct_sum_of_qpa_modules(&[
        simples[0].clone(),
        simples[1].clone(),
        simples[1].clone(),
        projectives[0].clone(),
        simples[0].clone(),
        simples[1].clone(),
    ])
    .unwrap()
    .unwrap();
    let prepared = ExtDimensionSource::new(&source).unwrap();
    assert_eq!(source.dimension_vector(), &[3, 4, 0]);
    // Counts are multiplicities of e_i A, NOT dimensions of M or P.
    assert_eq!(prepared.projective_multiplicities, vec![3, 3, 0]);
    assert_eq!(prepared.cover.source().dimension_vector(), &[3, 6, 3]);
    let mut targets = simples;
    targets.extend(projectives);
    targets.push(zero_module(&algebra).unwrap());
    targets.push(source.clone());
    for target in targets.iter().chain(targets.iter().rev()) {
        assert_eq!(
            projective_hom_dimension_from_multiplicities(
                &prepared.projective_multiplicities,
                target.dimension_vector()
            )
            .unwrap(),
            hom_over_algebra(prepared.cover.source(), target)
                .unwrap()
                .unwrap()
                .len()
        );
        assert_eq!(
            prepared.dimension(target).unwrap(),
            ext_over_algebra(&source, target).unwrap().dimension()
        );
    }
}

#[test]
fn ext_dimension_matches_basis_with_repeated_projectives_at_distinct_vertices() {
    compare(Q);
    compare(Z2);
    compare(Zp::new(101));
    compare(GaloisField::conway(2, 2).unwrap());
}
