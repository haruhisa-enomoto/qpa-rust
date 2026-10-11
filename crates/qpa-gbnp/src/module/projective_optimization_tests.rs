use super::*;
use qpa_core::{ArrowSpec, GaloisField, PathAlgebra, Quiver, VertexIndex};

fn reference_projective<F: Field>(
    algebra: &QuiverAlgebra<F>,
    data: &ProjectiveBasisData<F>,
    vertex: usize,
) -> PathAlgebraModule<F> {
    let intervals = &data.intervals[vertex];
    let dimensions = intervals.iter().map(Vec::len).collect::<Vec<_>>();
    let mut specification = Vec::new();
    for arrow in algebra.quiver().arrows() {
        let source = arrow.source().index().get() - 1;
        let target = arrow.target().index().get() - 1;
        if dimensions[source] == 0 || dimensions[target] == 0 {
            continue;
        }
        let element = include_in_path_algebra(&arrow.path(), algebra).unwrap();
        let rows = intervals[source]
            .iter()
            .map(|&index| {
                let image = data.basis[index].multiply(&element).unwrap();
                let coefficients =
                    canonical_basis_coefficients(algebra, &data.basis, &image).unwrap();
                intervals[target]
                    .iter()
                    .map(|&i| coefficients[i].clone())
                    .collect()
            })
            .collect();
        specification.push(NamedArrowMatrix::new(
            arrow.name(),
            DenseMatrix::from_rows(rows, algebra.field().clone()).unwrap(),
        ));
    }
    right_module_over_path_algebra(algebra.clone(), dimensions, &specification)
        .unwrap()
        .unwrap()
}

fn compare<F: Field>(field: F, scalar: &F::Element) {
    let vertex = |i| VertexIndex::new(i).unwrap();
    let parent = PathAlgebra::new(
        field,
        Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(1), "x"),
                ArrowSpec::named(vertex(1), vertex(1), "y"),
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(2), vertex(3), "b"),
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
    let relations = [
        arrows[0].multiply(&arrows[0]).unwrap(),
        arrows[1].multiply(&arrows[1]).unwrap(),
        arrows[1]
            .multiply(&arrows[0])
            .unwrap()
            .sub(
                &arrows[0]
                    .multiply(&arrows[1])
                    .unwrap()
                    .scale_left(scalar)
                    .unwrap(),
            )
            .unwrap(),
        arrows[2].multiply(&arrows[3]).unwrap(),
    ];
    let algebra =
        QuiverAlgebra::Quotient(crate::quotient_by_relations(parent, &relations).unwrap());
    let data = projective_basis_data(&algebra).unwrap().unwrap();
    let simple = simple_modules(&algebra).unwrap().remove(0);
    let generator = simple.basis_elements().remove(0);
    let selected = hom_from_projective(&generator, &simple).unwrap().unwrap();
    let held = selected.source().clone();
    assert!(algebra.cached_indec_projective_modules().is_none());
    let family = indec_projective_modules(&algebra).unwrap();
    assert!(held.is_identical(&family[0]));
    for (vertex, projective) in family.iter().enumerate() {
        let reference = reference_projective(&algebra, &data, vertex);
        assert_eq!(projective.dimension_vector(), reference.dimension_vector());
        assert_eq!(projective.arrow_matrices(), reference.arrow_matrices());
    }
    let source = direct_sum_of_qpa_modules(&[simple.clone(), family[1].clone(), simple])
        .unwrap()
        .unwrap();
    let cover = projective_cover(&source).unwrap();
    let maps = minimal_generating_set_of_module(source.clone())
        .unwrap()
        .iter()
        .map(|g| hom_from_projective(g, &source).unwrap().unwrap())
        .collect::<Vec<_>>();
    let sum =
        direct_sum_of_qpa_modules(&maps.iter().map(|m| m.source().clone()).collect::<Vec<_>>())
            .unwrap()
            .unwrap();
    let composites = direct_sum_projections(&sum)
        .unwrap()
        .iter()
        .zip(&maps)
        .map(|(p, m)| p.compose(m).unwrap())
        .collect::<Vec<_>>();
    let reference = composites
        .iter()
        .skip(1)
        .fold(composites[0].clone(), |a, b| a.add(b).unwrap());
    assert_eq!(cover.source().arrow_matrices(), sum.arrow_matrices());
    assert_eq!(cover.vertex_matrices(), reference.vertex_matrices());
    for (summand, map) in cover
        .source()
        .direct_sum_summands()
        .unwrap()
        .iter()
        .zip(&maps)
    {
        assert!(summand.is_identical(map.source()));
    }
    drop(family);
    let rebuilt = indec_projective_modules(&algebra).unwrap();
    assert!(held.is_identical(&rebuilt[0]));
}

#[test]
fn selective_projectives_and_cover_blocks_match_original_constructions() {
    compare(Q, &Q.one());
    compare(Z2, &Z2.one());
    let field = Zp::new(101);
    compare(field.clone(), &field.nth(Integer::from(3)));
    let extension = GaloisField::conway(2, 2).unwrap();
    let scalar = extension.primitive_element();
    compare(extension, &scalar);
}
