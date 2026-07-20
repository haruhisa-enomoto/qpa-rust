from __future__ import annotations

import pytest

import qpa


def test_exact_rational_and_prime_field_conversion() -> None:
    rationals = qpa.Field.rational()
    huge = 10**80 + 1
    value = rationals.element(huge, 3)
    assert value.numerator == str(huge)
    assert value.denominator == "3"
    assert repr(value) == f"{huge}/3"

    field5 = qpa.Field.prime(5)
    residue = field5.element(7, 3)
    assert residue.representative == 4
    assert repr(residue) == "4 mod 5"
    assert residue.field == field5

    with pytest.raises(qpa.InvalidInputError, match="not prime"):
        qpa.Field.prime(6)
    with pytest.raises(qpa.InvalidInputError, match=r"zero in GF\(2\)"):
        qpa.Field.prime(2).element(1, 2)


def test_quiver_identity_and_left_to_right_paths() -> None:
    quiver = qpa.Quiver(
        ["v1", "v2", "v3"],
        [("v1", "v2", "a"), ("v2", "v3", "b")],
    )
    path = quiver.path(["a", "b"])
    assert path.source == "v1"
    assert path.target == "v3"
    assert path.length == 2
    assert path.external_representation == [4, 5]
    assert quiver.arrow_path("a") * quiver.arrow_path("b") == path

    noncomposable = quiver.arrow_path("b") * quiver.arrow_path("a")
    assert noncomposable.length is None
    assert noncomposable.source is None

    reconstructed = qpa.Quiver(quiver.vertices, quiver.arrows)
    assert not quiver.is_identical(reconstructed)


def test_extended_quiver_domain_delegates_to_qpa_core() -> None:
    quiver = qpa.Quiver.indexed(
        3,
        [(1, 2, "a"), (2, 3, "b")],
    )
    assert quiver.number_of_vertices == 3
    assert quiver.number_of_arrows == 2
    assert quiver.adjacency_matrix == [[0, 1, 0], [0, 0, 1], [0, 0, 0]]
    assert quiver.is_acyclic and quiver.is_finite and quiver.is_connected
    assert quiver.is_unoriented_acyclic and quiver.is_tree
    assert quiver.cardinality == "7"
    assert quiver.topological_order == ["v1", "v2", "v3"]
    assert quiver.incoming_arrows("v2") == ["a"]
    assert quiver.outgoing_arrows("v2") == ["b"]
    assert quiver.neighbors("v2") == ["v3"]

    path = quiver.path(["a", "b"])
    assert path.quiver.is_identical(quiver)
    assert path.is_zero is False and path.is_vertex is False and path.is_arrow is False
    assert quiver.compare_paths(quiver.vertex_path("v1"), path) < 0
    assert quiver.path_at(6) == path
    assert quiver.path_from_external_representation([4, 5]) == path

    opposite = quiver.opposite()
    assert path.opposite().quiver.is_identical(opposite)
    assert opposite.opposite().is_identical(quiver)
    assert len(quiver.double().arrows) == 4
    assert quiver.separated().number_of_vertices == 6
    assert quiver.full_subquiver(["v1", "v2"]).number_of_arrows == 1
    assert len(quiver.connected_components()) == 1
    assert quiver.is_special_biserial

    dynkin = qpa.Quiver.dynkin("A", 3, ["r", "r"])
    assert dynkin is not None and dynkin.dynkin_type == "A_3"
    assert qpa.Quiver.from_adjacency([[0, 1], [0, 0]]).is_acyclic


def test_path_algebra_uses_exact_rust_operations() -> None:
    quiver = qpa.Quiver(
        ["v1", "v2", "v3"],
        [("v1", "v2", "a"), ("v2", "v3", "b")],
    )
    algebra = qpa.PathAlgebra(qpa.Field.rational(), quiver)
    assert algebra.quiver.is_identical(quiver)

    a = algebra.generator("a")
    b = algebra.generator("b")
    ab = a * b
    assert not ab.is_zero
    assert (b * a).is_zero

    half_a = algebra.term(quiver.arrow_path("a"), algebra.field.element(1, 2))
    assert "1/2" in repr(half_a)
    assert not (half_a + half_a).is_zero

    independent = qpa.PathAlgebra(qpa.Field.rational(), quiver)
    assert not algebra.is_identical(independent)
    with pytest.raises(qpa.InvalidInputError, match="different path algebras"):
        _ = a + independent.generator("a")


def test_quotient_runs_the_existing_high_level_groebner_route() -> None:
    quiver = qpa.Quiver(["v"], [("v", "v", "x")])
    algebra = qpa.PathAlgebra(qpa.Field.prime(5), quiver)
    x = algebra.generator("x")
    quotient = algebra.quotient([x * x])

    assert quotient.dimension() == "2"
    residue = quotient.generator("x")
    assert residue.is_normal_form
    assert (residue * residue).is_zero
    assert not quotient.one().is_zero


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)])
def test_ideal_groebner_and_quotient_domain_uses_qpa_algorithms(field: qpa.Field) -> None:
    quiver = qpa.Quiver(["v"], [("v", "v", "x")])
    algebra = qpa.PathAlgebra(field, quiver)
    vertex, x = algebra.generators
    x2 = x * x

    assert algebra.global_dimension == 1
    assert not algebra.has_finite_canonical_basis
    assert algebra.canonical_basis() is None
    assert algebra.relations == []
    assert algebra.is_monomial
    assert not algebra.is_special_biserial
    assert not algebra.is_string
    assert not algebra.is_gentle
    assert vertex.vertex_position == 1
    assert x.is_uniform and x.leading_monomial.length == 1
    assert x2.external_representation[0][0] == [2, 2]
    assert x2.leading_coefficient == field.element(1)
    assert x2 - x2 == algebra.zero()
    assert -(-x) == x
    unit = field.element(-1)
    assert x.scale_left(unit).divide_scalar(unit) == x
    assert algebra.nth_power_of_arrow_ideal(2) == [x2]
    assert algebra.make_uniform_on_right([x]) == [x]
    assert algebra.generators_times_arrows_on_right([x]) == [x2]

    ideal = algebra.ideal([x2])
    assert ideal.algebra.is_identical(algebra)
    assert ideal.attached_groebner_basis is None
    assert ideal.generators_are_monomial and ideal.is_monomial()
    assert ideal.contains(x2)
    assert ideal.attached_groebner_basis is not None
    assert ideal.is_admissible()
    assert ideal.minimal_generating_set() == [x2]

    basis = algebra.groebner_basis([x2]).completely_reduced_basis()
    gbnp_basis = algebra.gbnp_groebner_basis([x2]).completely_reduced_basis()
    assert basis.algebra.is_identical(algebra)
    assert basis.kind == "two-sided" and basis.is_complete
    assert basis.tip_reduced_property is True
    assert basis.completely_reduced_property is True
    assert len(basis) == 1
    assert gbnp_basis.relations == basis.relations
    assert basis.admits_finitely_many_nontips
    assert basis.nontip_size == "2"
    assert basis.nontips() == [quiver.vertex_path("v"), quiver.arrow_path("x")]
    assert basis.tip_reduce(x2 * x).is_zero
    assert basis.completely_reduce(x2 * x).is_zero
    assert basis.is_prefix_of_tip_in_tip_ideal(x2 * x)
    assert basis.right_basis().kind == "right"
    assert ideal.right_groebner_basis().kind == "right"
    assert ideal.product(ideal).contains(x2.pow(2))

    quotient = ideal.quotient()
    residue = quotient.generator("x")
    assert quotient.original_path_algebra.is_identical(algebra)
    assert quotient.ideal.is_identical(ideal)
    assert quotient.quiver.is_identical(quiver)
    assert quotient.has_normal_form_function
    assert quotient.is_finite_dimensional and quotient.dimension() == "2"
    assert quotient.is_admissible
    assert quotient.is_monomial
    assert quotient.is_special_biserial and quotient.is_string and quotient.is_gentle
    assert quotient.is_triangular_reduced
    assert not quotient.is_semisimple
    assert quotient.is_elementary and quotient.is_basic
    assert len(quotient.generators) == 2
    assert len(quotient.canonical_basis() or []) == 2
    assert residue.representative == x
    assert residue.leading_monomial == quiver.arrow_path("x")
    assert residue.is_uniform and (residue * residue).is_zero
    assert residue.inverse() is None and residue.pow_integer(-1) is None
    assert quotient.one().inverse() == quotient.one()
    assert residue.scale_right(field.element(2)) - residue.scale_left(field.element(2)) == quotient.zero()
    assert residue.opposite().algebra.opposite().is_identical(quotient)

    truncated = qpa.QuotientAlgebra.truncated(field, quiver, 2)
    assert truncated.dimension() == "2"
    assert (truncated.generator("x") * truncated.generator("x")).is_zero

    high_level = algebra.quotient_high_level([x2])
    assert high_level.dimension() == "2"


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)])
def test_checked_matrix_module_calls_qpa_core(field: qpa.Field) -> None:
    quiver = qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    algebra = qpa.PathAlgebra(field, quiver)
    module = algebra.module_from_integer_matrices([1, 1], [("a", [[1]])])

    assert module.dimension_vector == [1, 1]
    assert module.dimension == 2
    assert module.arrow_matrix_shapes == [(1, 1)]
    assert "dimension vector [ 1, 1 ]" in repr(module)
    assert module.is_identical(module)


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)])
def test_module_elements_homomorphisms_and_structure_follow_qpa(field: qpa.Field) -> None:
    quiver = qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    algebra = qpa.PathAlgebra(field, quiver)
    module = algebra.module_from_integer_matrices([1, 1], [("a", [[1]])])
    zero = field.element(0)
    one = field.element(1)

    assert module.field.characteristic == field.characteristic
    assert not module.is_zero
    assert len(module.arrow_matrices) == 1
    assert module.zero_element().is_zero

    basis = module.basis_elements()
    assert len(basis) == 2
    assert basis[0].support_vertices == [1]
    assert basis[1].support_vertices == [2]
    assert module.element([[one], [zero]]) == basis[0]
    assert basis[0].act_by_algebra_element(algebra.generator("a")) == basis[1]
    assert basis[0] + basis[1] - basis[1] == basis[0]
    assert basis[0].scale_left(one) == basis[0]

    identity = module.identity()
    entered_identity = module.homomorphism_to(module, [[[one]], [[one]]])
    assert identity == entered_identity
    assert identity.image(basis[0]) == basis[0]
    assert identity.preimage_representative(basis[1]) == basis[1]
    assert identity.is_injective and identity.is_surjective and identity.is_isomorphism
    assert identity.is_split_monomorphism and identity.is_split_epimorphism
    assert identity.inverse().compose(identity) == identity
    assert identity.right_inverse() == identity
    assert identity.left_inverse() == identity
    assert identity.kernel().is_zero
    assert identity.image_module() == module
    assert identity.cokernel().is_zero

    direct_sum = qpa.MatrixModule.direct_sum([module, module])
    assert direct_sum is not None
    assert direct_sum.dimension_vector == [2, 2]
    assert direct_sum.is_direct_sum
    assert len(direct_sum.direct_sum_summands() or []) == 2
    assert len(direct_sum.direct_sum_projections()) == 2
    assert len(direct_sum.direct_sum_inclusions()) == 2

    assert module.radical().dimension_vector == [0, 1]
    assert module.top().dimension_vector == [1, 0]
    assert module.socle().dimension_vector == [0, 1]
    assert not module.is_simple and not module.is_semisimple
    assert module.loewy_length == 2
    assert module.minimal_generating_set() == [basis[0]]

    loop = qpa.Quiver(["w"], [("w", "w", "x")])
    quotient = qpa.QuotientAlgebra.truncated(field, loop, 2)
    quotient_module = quotient.module_from_integer_matrices([2], [("x", [[0, 1], [0, 0]])])
    quotient_basis = quotient_module.basis_elements()
    assert quotient_basis[0].act_by_quotient_element(quotient.generator("x")) == quotient_basis[1]
    assert quotient_module.radical().dimension_vector == [1]
    assert quotient_module.top().dimension_vector == [1]
    assert quotient_module.socle().dimension_vector == [1]


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)])
def test_submodules_projectives_and_injectives_follow_qpa(field: qpa.Field) -> None:
    quiver = qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    algebra = qpa.PathAlgebra(field, quiver)
    module = algebra.module_from_integer_matrices([1, 1], [("a", [[1]])])
    basis = module.basis_elements()

    assert algebra.zero_module().dimension_vector == [0, 0]
    simples = algebra.simple_modules()
    projectives = algebra.indecomposable_projectives()
    injectives = algebra.indecomposable_injectives()
    assert [simple.dimension_vector for simple in simples] == [[1, 0], [0, 1]]
    assert [projective.dimension_vector for projective in projectives] == [[1, 1], [0, 1]]
    assert [injective.dimension_vector for injective in injectives] == [[1, 0], [1, 1]]
    assert module == projectives[0]
    assert module.is_projective and module.is_injective

    cover = simples[0].projective_cover()
    envelope = simples[1].injective_envelope()
    assert cover.source == projectives[0] and cover.range == simples[0]
    assert cover.is_surjective
    assert envelope.source == simples[1] and envelope.range == injectives[1]
    assert envelope.is_injective

    concrete = module.submodule([basis[0]])
    assert concrete.parent.is_identical(module)
    assert concrete.generators == [basis[0]]
    assert concrete.canonical_basis_elements() == basis
    assert concrete.as_module().dimension_vector == [1, 1]
    with pytest.raises(qpa.InvalidInputError, match="zero-space arrow"):
        module.submodule([basis[1]]).as_module()

    submodule = module.subrepresentation([basis[1]])
    inclusion = module.subrepresentation_inclusion([basis[1]])
    assert submodule.dimension_vector == [0, 1]
    assert inclusion.source == submodule and inclusion.range.is_identical(module)
    assert inclusion.is_injective

    radical = module.radical_inclusion()
    top = module.top_projection()
    socle = module.socle_inclusion()
    assert radical.source.dimension_vector == [0, 1]
    assert top.range.dimension_vector == [1, 0]
    assert socle.source.dimension_vector == [0, 1]

    identity = module.identity()
    assert identity.kernel_inclusion().source.is_zero
    image_projection = identity.image_projection()
    image_inclusion = identity.image_inclusion()
    factorization = identity.image_projection_inclusion()
    assert image_projection.compose(image_inclusion) == identity
    assert factorization[0].compose(factorization[1]) == identity
    assert identity.cokernel_projection().range.is_zero
    assert identity.top().source.dimension_vector == [1, 0]

    intersection = radical.intersection(socle)
    total = radical.sum(socle)
    assert intersection[2].source.dimension_vector == [0, 1]
    assert total[0].source.dimension_vector == [0, 1]

    loop = qpa.Quiver(["w"], [("w", "w", "x")])
    quotient = qpa.QuotientAlgebra.truncated(field, loop, 2)
    assert quotient.zero_module().dimension_vector == [0]
    assert [module.dimension_vector for module in quotient.simple_modules()] == [[1]]
    assert [module.dimension_vector for module in quotient.indecomposable_projectives()] == [[2]]
    assert [module.dimension_vector for module in quotient.indecomposable_injectives()] == [[2]]


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)])
def test_lazy_resolutions_and_syzygies_follow_qpa(field: qpa.Field) -> None:
    quiver = qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    algebra = qpa.PathAlgebra(field, quiver)
    first_simple, second_simple = algebra.simple_modules()

    assert first_simple.nth_syzygy(0).is_identical(first_simple)
    assert first_simple.first_syzygy().dimension_vector == [0, 1]
    assert first_simple.nth_syzygy(1).dimension_vector == [0, 1]
    assert first_simple.nth_syzygy(2).is_zero
    with pytest.raises(qpa.InvalidInputError, match="negative"):
        first_simple.nth_syzygy(-1)

    projective = first_simple.projective_resolution()
    assert projective.object(-1).is_identical(first_simple)
    assert projective.object(0).dimension_vector == [1, 1]
    assert projective.object(1).dimension_vector == [0, 1]
    assert projective.object(2).is_zero
    assert projective.object(5).is_zero
    differential_zero = projective.differential(0)
    differential_one = projective.differential(1)
    assert differential_zero.range.is_identical(first_simple)
    assert differential_one.compose(differential_zero).is_zero

    injective = second_simple.injective_resolution()
    assert injective.object(0).is_identical(second_simple)
    assert injective.differential(0).source.is_identical(second_simple)
    assert injective.differential(0).range.dimension_vector == [1, 1]
    assert injective.differential(0).compose(injective.differential(-1)).is_zero

    loop = qpa.Quiver(["w"], [("w", "w", "x")])
    quotient = qpa.QuotientAlgebra.truncated(field, loop, 2)
    simple = quotient.simple_modules()[0]
    periodic = simple.projective_resolution()
    assert periodic.object(0).dimension_vector == [2]
    assert periodic.object(1).dimension_vector == [2]


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)])
def test_module_complex_surface_follows_qpa(field: qpa.Field) -> None:
    quiver = qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    algebra = qpa.PathAlgebra(field, quiver)
    module = algebra.indecomposable_projectives()[0]
    first_simple, _second_simple = algebra.simple_modules()

    zero = module.zero_complex()
    assert zero.is_zero and zero.object(100).is_zero
    assert zero.upper_bound == "-infinity" and zero.lower_bound == "infinity"
    assert zero.length == "0" and zero.is_finite
    assert zero.is_exact_sequence

    stalk = first_simple.stalk_complex(0)
    assert not stalk.is_zero
    assert stalk.object(0).is_identical(first_simple)
    assert stalk.homology(0).dimension_vector == [1, 0]
    assert stalk.homology(1).is_zero
    assert stalk.is_finite and stalk.length == "1"
    assert stalk.shift(2).object(-2).dimension_vector == [1, 0]
    assert stalk.shift_unsigned(-1).object(1).dimension_vector == [1, 0]

    radical = module.radical_inclusion()
    top = module.top_projection()
    short_exact = module.finite_complex(0, [top, radical])
    checked_short_exact = radical.short_exact_sequence_with(top)
    assert checked_short_exact.is_short_exact_sequence
    assert short_exact.is_short_exact_sequence
    assert short_exact.is_exact_in_degree(0)
    assert short_exact.differential(0).source.is_identical(module)
    assert short_exact.cycles(0).dimension_vector == [0, 1]
    assert short_exact.boundaries(0).dimension_vector == [0, 1]
    assert short_exact.homology(0).is_zero
    assert short_exact.highest_known_degree is not None
    assert short_exact.lowest_known_degree is not None
    assert short_exact.positive_repeat_degrees() == ["3"]
    assert short_exact.negative_repeat_degrees() == ["-2"]
    assert short_exact.brutal_truncation_above(0).object(1).is_zero
    assert short_exact.brutal_truncation_below(0).object(-1).is_zero
    assert short_exact.good_truncation_above(0).object(1).is_zero
    assert short_exact.good_truncation_below(0).object(-1).is_zero
    assert short_exact.syzygy_truncation(0).object(1).dimension_vector == [0, 1]
    assert short_exact.cosyzygy_truncation(0).object(-2).is_zero
    good = short_exact.good_truncation(0, 0)
    assert good.object(0).is_zero and good.object(1).is_zero
    assert short_exact.brutal_truncation(0, 0).object(0).dimension_vector == [1, 1]
    both = short_exact.syzygy_cosyzygy_truncation(0, 0)
    assert both.object(-1).dimension_vector == [1, 0]
    assert both.object(1).dimension_vector == [0, 1]
    assert short_exact.is_identical(short_exact)
    assert " -> " in repr(short_exact)

    projective = first_simple.projective_resolution_complex()
    assert projective.object(-1).is_identical(first_simple)
    assert projective.object(0).dimension_vector == [1, 1]
    assert projective.object(1).dimension_vector == [0, 1]
    assert projective.object(2).is_zero
    assert projective.object(5).is_zero
    assert projective.is_positive_repeating
    assert projective.is_exact_in_degree(0)
    assert projective.cut_above().upper_bound == "1"

    assert first_simple.injective_resolution_complex().cut_below().lower_bound == "-1"
    injective = first_simple.injective_resolution_complex()
    assert injective.object(0).is_identical(first_simple)
    assert injective.object(-1).dimension_vector == [1, 0]
    assert injective.object(-2).is_zero
    assert injective.object(-3).is_zero
    assert injective.is_negative_repeating


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)])
def test_module_chain_maps_and_mapping_cones_follow_qpa(field: qpa.Field) -> None:
    quiver = qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    algebra = qpa.PathAlgebra(field, quiver)
    module = algebra.indecomposable_projectives()[0]
    complex = module.finite_complex(0, [module.top_projection(), module.radical_inclusion()])

    components = [complex.object(degree).identity() for degree in (-1, 0, 1)]
    identity = complex.finite_chain_map_to(complex, -1, components)
    assert identity.source.is_identical(complex)
    assert identity.range.is_identical(complex)
    assert identity.morphism(0).is_isomorphism
    assert identity.highest_known_degree == "1"
    assert identity.lowest_known_degree == "-1"
    assert identity.is_identical(identity)
    assert repr(identity) == "<chain map>"

    zero = complex.zero_chain_map_to(complex)
    assert zero.morphism(50).is_zero
    assert zero.morphism(-50).is_zero

    cone, inclusion, projection = identity.mapping_cone()
    assert cone.object(0).dimension_vector == [2, 1]
    assert inclusion.source.is_identical(complex)
    assert inclusion.range.is_identical(cone)
    assert projection.source.is_identical(cone)
    assert projection.range.object(0).dimension_vector == complex.object(-1).dimension_vector
    assert inclusion.morphism(0).range.is_identical(cone.object(0))
    assert projection.morphism(0).source.is_identical(cone.object(0))

    simple = algebra.simple_modules()[0]
    module_identity = simple.identity()
    projective = simple.projective_resolution_complex()
    lifting = module_identity.comparison_lifting(projective, projective)
    assert lifting.morphism(-1) == module_identity
    canonical_lifting = module_identity.comparison_lifting_to_projective_resolution()
    assert canonical_lifting.morphism(-1) == module_identity


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)])
def test_yoneda_product_follows_qpa_splicing_convention(field: qpa.Field) -> None:
    algebra = qpa.PathAlgebra(field, qpa.Quiver(["v"], []))
    one = algebra.module_from_integer_matrices([1], [])
    two = algebra.module_from_integer_matrices([2], [])
    zero = field.element(0)
    unit = field.element(1)
    inclusion = one.homomorphism_to(two, [[[unit, zero]]])
    projection = two.homomorphism_to(one, [[[zero], [unit]]])
    extension = inclusion.short_exact_sequence_with(projection)

    product = extension.yoneda_product(extension)
    assert product.lower_bound == "-1"
    assert product.upper_bound == "2"
    assert [product.object(i).dimension_vector for i in range(-1, 3)] == [
        [1],
        [2],
        [2],
        [1],
    ]


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)])
def test_derived_category_foundations_follow_qpa(field: qpa.Field) -> None:
    algebra = qpa.PathAlgebra(field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")]))
    projectives = algebra.indecomposable_projectives()
    injectives = algebra.indecomposable_injectives()
    simples = algebra.simple_modules()
    projective_sum = qpa.MatrixModule.direct_sum([projectives[0], projectives[1], projectives[0]])
    injective_sum = qpa.MatrixModule.direct_sum([injectives[1], injectives[0], injectives[1]])
    assert projective_sum is not None and injective_sum is not None

    assert projectives[0].compare_with_indecomposable_projective() == 1
    assert projective_sum.compare_with_indecomposable_projective() is None
    assert simples[0].compare_with_indecomposable_projective() is None
    assert injectives[1].compare_with_indecomposable_injective() == 2
    assert injective_sum.compare_with_indecomposable_injective() is None
    assert projective_sum.projective_description() == [1, 2, 1]
    assert algebra.zero_module().projective_description() == []

    projective_complex = projective_sum.stalk_complex(1)
    injective_complex = injective_sum.stalk_complex(-1)
    assert projective_complex.is_projective_complex
    assert injective_complex.is_injective_complex
    assert not simples[0].stalk_complex(0).is_projective_complex
    assert not simples[1].stalk_complex(0).is_injective_complex
    assert projective_complex.projective_description_in_degree(1) == [2, 1]
    assert projective_complex.description_in_degree(1, True) == [2, 1]
    assert injective_complex.injective_description_in_degree(-1) == [1, 2]
    assert injective_complex.description_in_degree(-1, False) == [1, 2]
    assert projective_complex.finite_projective_description() == [[2, 1]]
    assert projective_complex.finite_description(True) == [[2, 1]]
    assert injective_complex.finite_injective_description() == [[1, 2]]
    assert injective_complex.finite_description(False) == [[1, 2]]

    identity = projective_sum.identity()
    components = identity.all_components()
    assert len(components) == 3 and all(len(row) == 3 for row in components)
    assert components[0][0].is_isomorphism and components[0][1].is_zero
    assert [len(matrix) for matrix in identity.dual_vertex_matrices] == [2, 3]
    inclusion, projection = projective_sum.direct_sum_minus_summands([2])
    assert inclusion.source.dimension_vector == [2, 2]
    assert projection.range.dimension_vector == [2, 2]


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)])
def test_derived_category_transforms_follow_qpa(field: qpa.Field) -> None:
    algebra = qpa.PathAlgebra(field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")]))
    projectives = algebra.indecomposable_projectives()
    simple = algebra.simple_modules()[0]

    indecomposable_star = projectives[0].identity().star_between_indecomposable_projectives(1, [1])
    assert indecomposable_star.is_isomorphism
    projective_sum = qpa.MatrixModule.direct_sum(projectives)
    assert projective_sum is not None
    description = projective_sum.projective_description()
    identity = projective_sum.identity()
    decomposable_star = identity.star_between_decomposable_projectives(description, description)
    dispatched_star = identity.star_between_projectives(description, description)
    assert decomposable_star.is_isomorphism and dispatched_star == decomposable_star

    projective_stalk = projectives[0].stalk_complex(2)
    assert projective_stalk.projective_resolution().is_identical(projective_stalk)
    resolved_stalk = simple.stalk_complex(2).projective_resolution()
    assert resolved_stalk.object(2).dimension_vector == [1, 1]
    assert resolved_stalk.object(3).dimension_vector == [0, 1]
    assert resolved_stalk.object(4).is_zero
    assert resolved_stalk.is_projective_complex

    finite_injective = projective_stalk.projective_to_injective_finite()
    general_injective = projective_stalk.projective_to_injective()
    assert finite_injective.object(2).dimension_vector == [1, 0]
    assert finite_injective.is_injective_complex and general_injective.is_injective_complex

    semisimple = qpa.PathAlgebra(field, qpa.Quiver(["v"], []))
    tau = semisimple.indecomposable_projectives()[0].stalk_complex(0).tau()
    assert tau.lower_bound == "-1" and tau.upper_bound == "-1"
    assert tau.object(-1).dimension_vector == [1]
    assert tau.is_projective_complex


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)])
def test_derived_category_reductions_follow_qpa(field: qpa.Field) -> None:
    algebra = qpa.PathAlgebra(field, qpa.Quiver(["v"], []))
    projective = algebra.indecomposable_projectives()[0]
    source = qpa.MatrixModule.direct_sum([projective, projective])
    range_module = qpa.MatrixModule.direct_sum([projective, projective])
    assert source is not None and range_module is not None
    zero = field.element(0)
    unit = field.element(1)
    differential = source.homomorphism_to(range_module, [[[unit, unit], [unit, zero]]])
    zero_map = source.homomorphism_to(range_module, [[[zero, zero], [zero, zero]]])
    assert zero_map.check_for_projective_isomorphism() is None

    cancellation = differential.check_for_projective_isomorphism()
    assert cancellation is not None
    assert cancellation.vertex == 1
    assert cancellation.source_position == 1 and cancellation.range_position == 1
    assert cancellation.map.is_isomorphism
    assert "source summand 1" in repr(cancellation)

    source_inclusion, _source_projection = source.direct_sum_minus_summands([1])
    range_inclusion, _range_projection = range_module.direct_sum_minus_summands([1])
    reduced = differential.reduce_differential(
        cancellation, source_inclusion.source, range_inclusion.source
    )
    assert reduced.source.dimension_vector == [1]
    assert reduced.range.dimension_vector == [1]
    assert reduced.is_isomorphism

    t_map = source.direct_sum_projections()[0]
    reduced_t = t_map.reduce_t_map(differential, cancellation, source_inclusion)
    assert reduced_t.source.dimension_vector == [1] and reduced_t.is_isomorphism
    range_projections = range_module.direct_sum_projections()
    previous = range_projections[0].add(range_projections[1])
    assert previous.reduce_previous_differential(range_inclusion).is_isomorphism
    source_projections = source.direct_sum_projections()
    previous_t = source_projections[0].add(source_projections[1])
    assert previous_t.reduce_previous_t_map(source_inclusion).is_isomorphism


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)])
def test_homological_dimensions_and_periods_follow_qpa(field: qpa.Field) -> None:
    algebra = qpa.PathAlgebra(field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")]))
    first_simple, second_simple = algebra.simple_modules()
    assert algebra.global_dimension_bounded(0) == "1"
    assert algebra.cached_gorenstein_dimension is None
    assert algebra.gorenstein_dimension_bounded(0) == 1
    assert algebra.cached_gorenstein_dimension == 1
    assert algebra.dominant_dimension_bounded(3) == "1"

    assert first_simple.cached_projective_dimension is None
    assert first_simple.cached_injective_dimension is None
    assert first_simple.projective_dimension(3) == 1
    assert first_simple.injective_dimension(3) == 0
    assert first_simple.cached_projective_dimension == 1
    assert first_simple.cached_injective_dimension == 0
    assert first_simple.dominant_dimension(3) == "0"
    assert first_simple.is_nth_syzygy(0)
    assert not first_simple.is_nth_syzygy(1)
    assert first_simple.omega_period(3) is None
    assert first_simple.tau_period(3) is None
    assert second_simple.projective_dimension(3) == 0
    assert second_simple.injective_dimension(3) == 1

    dual_numbers = qpa.QuotientAlgebra.truncated(
        field, qpa.Quiver(["v"], [("v", "v", "x")]), 2
    )
    periodic_simple = dual_numbers.simple_modules()[0]
    assert dual_numbers.global_dimension_bounded(3) == "infinity"
    assert dual_numbers.gorenstein_dimension_bounded(3) == 0
    assert dual_numbers.cached_gorenstein_dimension == 0
    assert dual_numbers.dominant_dimension_bounded(3) == "infinity"
    assert periodic_simple.projective_dimension(3) is None
    assert periodic_simple.injective_dimension(3) is None
    assert periodic_simple.dominant_dimension(3) == "exceeds_bound"
    assert periodic_simple.is_nth_syzygy(1)
    assert periodic_simple.omega_period(3) == "1"
    assert periodic_simple.tau_period(3) == "1"


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)])
def test_ext_results_and_generators_follow_qpa(field: qpa.Field) -> None:
    algebra = qpa.PathAlgebra(field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")]))
    first_simple, second_simple = algebra.simple_modules()
    ext = first_simple.ext1(second_simple)
    assert ext.dimension == 1
    assert ext.syzygy_inclusion.source.dimension_vector == [0, 1]
    assert len(ext.basis) == 1
    coordinates = ext.coefficients(ext.basis[0])
    assert coordinates is not None and len(coordinates) == 1
    assert coordinates[0].field.characteristic == field.characteristic
    zero_ext = second_simple.ext1(first_simple)
    assert zero_ext.dimension == 0 and zero_ext.basis == []

    dual_numbers = qpa.QuotientAlgebra.truncated(
        field, qpa.Quiver(["v"], [("v", "v", "x")]), 2
    )
    generators = dual_numbers.simple_modules()[0].ext_algebra_generators(3)
    assert generators.dimensions == [1, 1, 1, 1]
    assert generators.generator_counts == [1, 1, 0, 0]
    assert [len(degree) for degree in generators.generators] == [1, 1, 0, 0]


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)])
def test_module_dualities_and_stable_functors_follow_qpa(field: qpa.Field) -> None:
    algebra = qpa.PathAlgebra(field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")]))
    first_simple, second_simple = algebra.simple_modules()
    projective = algebra.indecomposable_projectives()[0]

    assert first_simple.dual().dimension_vector == [1, 0]
    assert first_simple.transpose().dimension_vector == [0, 1]
    assert first_simple.d_tr().dimension_vector == [0, 1]
    assert first_simple.tr_d().is_zero
    assert first_simple.d_tr(0).is_identical(first_simple)
    assert first_simple.d_tr(-1).is_zero
    assert first_simple.star().is_zero
    assert first_simple.nakayama().is_zero
    assert first_simple.opposite_nakayama().dimension_vector == [1, 1]
    assert second_simple.identity_to_double_star().is_isomorphism

    identity = projective.identity()
    assert identity.dual().is_isomorphism
    assert identity.transpose().is_isomorphism
    assert identity.d_tr().is_isomorphism
    assert identity.tr_d().is_isomorphism
    assert identity.d_tr(0).is_identical(identity)
    assert identity.d_tr(-1).is_isomorphism
    assert identity.star().source.dimension_vector == [1, 0]
    assert identity.nakayama().source.dimension_vector == [1, 0]
    assert identity.opposite_nakayama().source.dimension_vector == [0, 1]


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)])
def test_rigidity_predicates_follow_qpa(field: qpa.Field) -> None:
    algebra = qpa.PathAlgebra(field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")]))
    simple = algebra.simple_modules()[0]
    assert simple.is_rigid
    assert simple.is_tau_rigid
    assert simple.is_exceptional is True
    assert simple.is_n_rigid(2)

    dual_numbers = qpa.QuotientAlgebra.truncated(
        field, qpa.Quiver(["v"], [("v", "v", "x")]), 2
    )
    periodic_simple = dual_numbers.simple_modules()[0]
    assert not periodic_simple.is_rigid
    assert not periodic_simple.is_tau_rigid
    assert periodic_simple.is_exceptional is False
    assert periodic_simple.is_n_rigid(0)
    assert not periodic_simple.is_n_rigid(1)


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)])
def test_tilting_and_cotilting_results_follow_qpa(field: qpa.Field) -> None:
    algebra = qpa.PathAlgebra(field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")]))
    projective_generator = qpa.MatrixModule.direct_sum(algebra.indecomposable_projectives())
    injective_generator = qpa.MatrixModule.direct_sum(algebra.indecomposable_injectives())
    assert projective_generator is not None and injective_generator is not None
    assert projective_generator.cached_is_tilting is None
    assert injective_generator.cached_is_cotilting is None

    coresolution = projective_generator.finite_coresolution_in_add(projective_generator, 0)
    assert coresolution is not None
    assert coresolution.base_degree == "0"
    assert len(coresolution.differentials) == 1
    assert coresolution.differential(0).is_isomorphism
    assert coresolution.object(0).dimension_vector == [1, 2]

    resolution = injective_generator.finite_resolution_in_add(injective_generator, 0)
    assert resolution is not None and resolution.is_in_add
    assert resolution.differentials == []
    assert resolution.differential(0) is None and resolution.object(0) is None

    tilting = projective_generator.tilting(2)
    assert tilting is not None and tilting.projective_dimension == 0
    assert len(tilting.coresolutions) == 2
    assert projective_generator.cached_is_tilting is True
    cotilting = injective_generator.cotilting(2)
    assert cotilting is not None and cotilting.injective_dimension == 0
    assert len(cotilting.resolutions) == 2
    assert all(result.is_in_add for result in cotilting.resolutions)
    assert injective_generator.cached_is_cotilting is True


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)])
def test_approximation_operations_follow_qpa(field: qpa.Field) -> None:
    algebra = qpa.PathAlgebra(field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")]))
    projective = algebra.indecomposable_projectives()[0]
    injectives = algebra.indecomposable_injectives()
    simples = algebra.simple_modules()
    simple = simples[0]

    right = projective.right_approximation_by_add(simple)
    assert right.source.dimension_vector == [1, 1]
    assert right.range.dimension_vector == [1, 0]
    right_list = qpa.MatrixModule.right_approximation_by_add_list(
        [projective, simple], simple
    )
    assert right_list.source.dimension_vector == [2, 1]
    assert right_list.range.dimension_vector == [1, 0]
    minimal_right = projective.minimal_right_approximation(simple)
    assert minimal_right.source.dimension_vector == [1, 1]
    assert projective.minimal_right_add_approximation(simple) == minimal_right

    left = projective.left_approximation_by_add(simple)
    assert left.source.dimension_vector == [1, 1]
    assert left.range.dimension_vector == [1, 0]
    minimal_left = projective.minimal_left_approximation(simple)
    assert minimal_left.range.dimension_vector == [1, 0]
    assert projective.minimal_left_add_approximation(simple) == minimal_left

    right_fac = projective.right_fac_approximation(simple)
    assert right_fac.source.dimension_vector == [1, 0]
    assert right_fac.range.dimension_vector == [1, 0]
    assert projective.right_fac_m_approximation(simple) == right_fac
    assert projective.minimal_right_fac_m_approximation(simple) == right_fac
    left_fac = projective.left_fac_approximation(simple)
    assert left_fac.source.dimension_vector == [1, 1]
    assert left_fac.range.dimension_vector == [1, 0]
    assert projective.minimal_left_fac_approximation(simple) == left_fac

    right_sub = projective.right_sub_approximation(simple)
    assert right_sub.source.dimension_vector == [1, 1]
    assert right_sub.range.dimension_vector == [1, 0]
    assert projective.minimal_right_sub_approximation(simple) == right_sub
    left_sub = projective.left_sub_approximation(simple)
    assert left_sub.source.dimension_vector == [1, 1]
    assert left_sub.range.dimension_vector == [1, 0]
    assert projective.minimal_left_sub_approximation(simple) == left_sub

    radical = qpa.MatrixModule.radical_right_approximation_by_add(
        [projective, simple], 2
    )
    assert radical.source.dimension_vector == [1, 1]
    assert radical.range.dimension_vector == [1, 0]

    cotilting = qpa.MatrixModule.direct_sum([injectives[1], simples[1]])
    assert cotilting is not None and cotilting.cotilting(1) is not None
    perpendicular = cotilting.right_approximation_by_perp(simple)
    assert perpendicular.source.dimension_vector == [1, 1]
    assert perpendicular.range.dimension_vector == [1, 0]
    t_hat = cotilting.left_approximation_by_add_t_hat(simple)
    assert t_hat.source.dimension_vector == [1, 0]
    assert t_hat.range.dimension_vector == [1, 0]


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)])
def test_complements_and_mutations_follow_qpa(field: qpa.Field) -> None:
    algebra = qpa.PathAlgebra(field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")]))
    first_projective, second_projective = algebra.indecomposable_projectives()
    first_simple = algebra.simple_modules()[0]

    coresolution_side = first_projective.all_complements(second_projective)
    assert coresolution_side.resolution is None
    assert coresolution_side.coresolution is not None
    assert coresolution_side.coresolution.differential(-1).source.dimension_vector == [1, 1]
    resolution_side = first_projective.all_complements(first_simple)
    assert resolution_side.resolution is not None
    assert resolution_side.coresolution is None
    assert resolution_side.resolution.differential(1).source.dimension_vector == [1, 1]
    cotilting_alias = first_projective.all_cotilting_complements(second_projective)
    assert cotilting_alias.resolution is None and cotilting_alias.coresolution is not None

    assert first_projective.cached_faithful_dimension is None
    assert first_projective.faithful_dimension() == "1"
    assert first_projective.cached_faithful_dimension == "1"
    assert first_projective.number_of_tilting_complements() == "2"
    assert first_projective.number_of_cotilting_complements() == "2"
    regular = qpa.MatrixModule.direct_sum([first_projective, second_projective])
    assert regular is not None
    assert regular.number_of_tilting_complements() == "infinity"
    assert regular.number_of_cotilting_complements() == "infinity"

    left = first_projective.left_tilting_complement_mutation(second_projective)
    assert left is not None and left.dimension_vector == [1, 0]
    assert first_projective.right_tilting_complement_mutation(second_projective) is None
    left_cotilting = first_projective.left_cotilting_complement_mutation(second_projective)
    assert left_cotilting is not None and left_cotilting.dimension_vector == [1, 0]
    assert first_projective.right_cotilting_complement_mutation(second_projective) is None


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)])
def test_almost_split_sequences_follow_qpa(field: qpa.Field) -> None:
    algebra = qpa.PathAlgebra(field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")]))
    projectives = algebra.indecomposable_projectives()
    injectives = algebra.indecomposable_injectives()
    simples = algebra.simple_modules()

    ending = simples[0].almost_split_sequence()
    assert ending is not None
    monomorphism, epimorphism = ending
    assert monomorphism.source.dimension_vector == [0, 1]
    assert monomorphism.range.dimension_vector == [1, 1]
    assert epimorphism.range.dimension_vector == [1, 0]
    assert simples[0].almost_split_sequence("r") is not None
    assert simples[0].almost_split_sequence("starting") is None
    assert simples[1].almost_split_sequence("l") is not None
    assert projectives[0].almost_split_sequence() is None

    cotilting = qpa.MatrixModule.direct_sum([injectives[1], simples[1]])
    assert cotilting is not None and cotilting.cotilting(1) is not None
    perpendicular = cotilting.almost_split_sequence_in_perp(simples[0])
    assert perpendicular is not None
    assert perpendicular[0].source.dimension_vector == [0, 1]
    assert perpendicular[1].range.dimension_vector == [1, 0]


@pytest.mark.parametrize("field", [qpa.Field.prime(2), qpa.Field.prime(5)])
def test_irreducible_morphisms_and_predecessors_follow_qpa(field: qpa.Field) -> None:
    algebra = qpa.PathAlgebra(field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")]))
    simples = algebra.simple_modules()
    ending = simples[0].irreducible_morphisms_ending_in()
    starting = simples[1].irreducible_morphisms_starting_in()
    assert len(ending) == 1 and ending[0].source.dimension_vector == [1, 1]
    assert len(starting) == 1 and starting[0].range.dimension_vector == [1, 1]

    first = simples[0].predecessors(1)
    assert first.is_module and first.module == simples[0]
    assert first.layers is None and first.valuations is None
    second = simples[0].predecessors(2)
    assert not second.is_module and second.module is None
    assert second.layers is not None
    assert [len(layer) for layer in second.layers] == [1, 1, 1]
    assert second.valuations is not None
    valuation = second.valuations[0][0]
    assert valuation.predecessor_position == 1
    assert valuation.successor_position == 1
    assert valuation.pair_multiplicity == (1, None)
    assert valuation.scalar_multiplicity is None


def test_irreducible_morphisms_preserve_qpa_finite_field_requirement() -> None:
    algebra = qpa.PathAlgebra(
        qpa.Field.rational(), qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    with pytest.raises(qpa.InvalidInputError, match="finite field"):
        algebra.simple_modules()[0].irreducible_morphisms_ending_in()


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)])
def test_module_comparison_operations_follow_qpa(field: qpa.Field) -> None:
    algebra = qpa.PathAlgebra(field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")]))
    projective = algebra.indecomposable_projectives()[0]
    first_simple, second_simple = algebra.simple_modules()
    mixed = qpa.MatrixModule.direct_sum([first_simple, projective])
    repeated = qpa.MatrixModule.direct_sum([first_simple, first_simple])
    assert mixed is not None and repeated is not None

    assert first_simple.isomorphic_to(first_simple)
    assert not first_simple.isomorphic_to(second_simple)
    isomorphism = first_simple.isomorphism_to(first_simple)
    assert isomorphism is not None and isomorphism.is_isomorphism
    assert first_simple.isomorphism_to(second_simple) is None
    assert first_simple.is_direct_summand_of(mixed)
    assert first_simple.is_in_additive_closure_of(repeated) is True
    assert first_simple.is_indecomposable is True
    if field.characteristic is None:
        assert mixed.is_indecomposable is None
        with pytest.raises(qpa.InvalidInputError, match="finite field"):
            mixed.decompose()
    else:
        assert mixed.is_indecomposable is False


@pytest.mark.parametrize("field", [qpa.Field.prime(2), qpa.Field.prime(5)])
def test_module_decomposition_operations_follow_qpa(field: qpa.Field) -> None:
    algebra = qpa.PathAlgebra(field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")]))
    projective = algebra.indecomposable_projectives()[0]
    simple = algebra.simple_modules()[0]
    mixed = qpa.MatrixModule.direct_sum([simple, projective])
    repeated = qpa.MatrixModule.direct_sum([simple, simple])
    assert mixed is not None and repeated is not None

    assert [module.dimension_vector for module in mixed.decompose()] == [[1, 0], [1, 1]]
    inclusions = mixed.decompose_with_inclusions()
    assert len(inclusions) == 2
    assert all(inclusion.range == mixed for inclusion in inclusions)
    multiplicities = repeated.decompose_with_multiplicities()
    assert multiplicities is not None
    assert [module.dimension_vector for module in multiplicities.basic_summands] == [[1, 0]]
    assert multiplicities.multiplicities == [2]
    assert repeated.basic_version().dimension_vector == [1, 0]

    idempotents = mixed.block_splitting_idempotents()
    blocks = mixed.block_decomposition()
    assert len(idempotents) == 2
    assert [module.dimension_vector for module in blocks] == [[1, 1], [1, 0]]
    assert [module.dimension_vector for module in qpa.MatrixModule.make_basic([simple, simple])] == [
        [1, 0]
    ]
    summand_data = mixed.nonisomorphic_direct_summand_data()
    assert summand_data.count == 2
    assert summand_data.block_dimensions == [1, 1]
    assert algebra.zero_module().decompose_with_multiplicities() is None


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)])
def test_hom_spaces_categorical_maps_and_minimal_versions_follow_qpa(
    field: qpa.Field,
) -> None:
    algebra = qpa.PathAlgebra(field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")]))
    projective = algebra.indecomposable_projectives()[0]
    simple = algebra.simple_modules()[0]

    hom_basis = projective.hom_basis_to(simple)
    assert hom_basis is not None and len(hom_basis) == 1
    factoring = projective.hom_factoring_through_projectives_to(simple)
    assert len(factoring) == 1 and factoring[0].range == simple
    assert projective.zero_map_to(simple).is_zero
    identity = simple.homomorphism_from_images(simple, simple.basis_elements())
    assert identity.is_isomorphism
    assert identity.images(simple.basis_elements()) == simple.basis_elements()

    cover = simple.projective_cover()
    pushout = cover.push_out_with(projective.identity())
    assert pushout is not None
    assert pushout[0].source.dimension_vector == [1, 1]
    assert pushout[1].source.dimension_vector == [1, 0]
    pullback = cover.pull_back_with(simple.identity())
    assert pullback is not None
    assert pullback[0].range.dimension_vector == [1, 0]
    assert pullback[1].range.dimension_vector == [1, 1]

    repeated = qpa.MatrixModule.direct_sum([simple, simple])
    assert repeated is not None
    inclusions = repeated.direct_sum_inclusions()
    projections = repeated.direct_sum_projections()
    inclusion_lift = inclusions[0].lift_inclusion(inclusions[0])
    assert inclusion_lift is not None and inclusion_lift.is_isomorphism
    projective_lift = cover.lift_from_projective(cover)
    assert projective_lift is not None and projective_lift.is_isomorphism

    right_nonminimal = projections[0].add(projections[1])
    assert not right_nonminimal.is_right_minimal
    right_minimal, removed_sources = right_nonminimal.right_minimal_version()
    assert right_minimal.is_right_minimal
    assert right_minimal.source.dimension_vector == [1, 0]
    assert len(removed_sources) == 1
    left_nonminimal = inclusions[0].add(inclusions[1])
    assert not left_nonminimal.is_left_minimal
    left_minimal, removed_ranges = left_nonminimal.left_minimal_version()
    assert left_minimal.is_left_minimal
    assert left_minimal.range.dimension_vector == [1, 0]
    assert len(removed_ranges) == 1


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)])
def test_structured_homomorphism_operations_follow_qpa(field: qpa.Field) -> None:
    algebra = qpa.PathAlgebra(field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")]))
    first_simple, second_simple = algebra.simple_modules()
    repeated = qpa.MatrixModule.direct_sum([first_simple, first_simple])
    assert repeated is not None
    inclusions = repeated.direct_sum_inclusions()
    projections = repeated.direct_sum_projections()

    right_reduction = projections[0].add(projections[1]).more_right_minimal_version()
    assert not right_reduction.is_unchanged
    assert right_reduction.minimal_part.source.dimension_vector == [1, 0]
    assert right_reduction.removed_part is not None
    assert right_reduction.removed_part.source.dimension_vector == [1, 0]
    left_reduction = inclusions[0].add(inclusions[1]).more_left_minimal_version()
    assert not left_reduction.is_unchanged
    assert left_reduction.minimal_part.range.dimension_vector == [1, 0]
    assert left_reduction.removed_part is not None
    assert left_reduction.removed_part.range.dimension_vector == [1, 0]
    assert first_simple.identity().more_right_minimal_version().is_unchanged
    assert first_simple.identity().more_left_minimal_version().is_unchanged

    common = repeated.common_direct_summand(first_simple)
    assert common.status == "found"
    assert common.source_common is not None
    assert common.source_common.dimension_vector == [1, 0]
    assert common.source_complement is not None
    assert common.source_complement.dimension_vector == [1, 0]
    assert common.range_common is not None
    assert common.range_common.dimension_vector == [1, 0]
    assert common.range_complement is not None
    assert common.range_complement.dimension_vector == [0, 0]
    assert first_simple.common_direct_summand(second_simple).status == "none"

    maximal = repeated.maximal_common_direct_summand(first_simple)
    assert maximal is not None
    assert [summand.dimension_vector for summand in maximal.common_summands] == [[1, 0]]
    assert maximal.source_remainder.dimension_vector == [1, 0]
    assert maximal.range_remainder.dimension_vector == [0, 0]
    assert first_simple.maximal_common_direct_summand(second_simple) is None

    cover = first_simple.projective_cover()
    induced_maps = [
        cover.morphism_on_kernel(cover, cover.source.identity(), cover.range.identity()),
        cover.morphism_on_image(cover, cover.source.identity(), cover.range.identity()),
        cover.morphism_on_cokernel(cover, cover.source.identity(), cover.range.identity()),
    ]
    assert all(induced is not None and induced.is_isomorphism for induced in induced_maps)
    assert cover.morphism_on_kernel(
        cover, cover.source.zero_map_to(cover.source), cover.range.identity()
    ) is None

    identity = first_simple.identity()
    zero = first_simple.zero_map_to(first_simple)
    block_identity = qpa.ModuleHomomorphism.multiply_lists(
        projections, [[identity, zero], [zero, identity]], inclusions
    )
    assert block_identity == repeated.identity()


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)])
def test_endomorphism_algebras_and_stable_quotients_follow_qpa(field: qpa.Field) -> None:
    algebra = qpa.PathAlgebra(field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")]))
    projective = algebra.indecomposable_projectives()[0]
    first_simple, second_simple = algebra.simple_modules()
    zero = field.element(0)
    one = field.element(1)

    hom = projective.hom_basis_with_coefficients_to(first_simple)
    assert hom is not None and hom.dimension == 1
    assert hom.coefficients(hom.basis[0]) == [one]
    zero_hom = first_simple.hom_basis_with_coefficients_to(second_simple)
    assert zero_hom is not None and zero_hom.dimension == 0
    assert zero_hom.coefficients(first_simple.zero_map_to(second_simple)) == [zero]

    module = algebra.module_from_integer_matrices(
        [2, 2], [("a", [[1, 0], [0, 0]])]
    )
    endomorphisms = module.endomorphism_algebra()
    assert endomorphisms.field.characteristic == field.characteristic
    assert endomorphisms.dimension == 5
    assert len(endomorphisms.basis_matrices) == 5
    assert endomorphisms.is_identical(module.endomorphism_algebra())
    structure_constants = endomorphisms.structure_constants_algebra()
    assert structure_constants.dimension == 5
    assert len(structure_constants.multiplication_tensor) == 5
    one_coordinates = endomorphisms.coefficients(endomorphisms.one)
    assert one_coordinates is not None and len(one_coordinates) == 5
    assert endomorphisms.matrix_from_coefficients(one_coordinates) == endomorphisms.one
    assert (
        endomorphisms.multiply(endomorphisms.one, endomorphisms.basis_matrices[0])
        == endomorphisms.basis_matrices[0]
    )
    assert all(
        module.endomorphism_from_matrix(matrix).endomorphism_matrix() == matrix
        for matrix in endomorphisms.basis_matrices
    )
    assert first_simple.identity().endomorphism_matrix() == [[one]]
    assert first_simple.identity().underlying_linear_matrix() == [
        [one, zero],
        [zero, zero],
    ]

    natural = module.stable_endomorphism_quotient()
    stable = natural.range
    assert natural.source.is_identical(endomorphisms)
    assert stable.field.characteristic == field.characteristic
    assert stable.dimension == 1
    assert len(natural.kernel_basis_matrices) == 4
    assert [natural.image(matrix).coordinates for matrix in endomorphisms.basis_matrices] == [
        [zero],
        [one],
        [zero],
        [zero],
        [zero],
    ]
    stable_one = stable.one()
    assert stable_one is not None
    assert stable_one.algebra.is_identical(stable)
    assert stable.multiply(stable_one, stable_one).coordinates == [one]
    assert stable.zero().is_zero

    projective_stable = second_simple.stable_endomorphism_quotient().range
    assert projective_stable.dimension == 0
    assert projective_stable.basis_elements == []
    assert projective_stable.one() is None

    with pytest.raises(qpa.InvalidInputError, match="direct sum"):
        algebra.zero_module().endomorphism_algebra()
    with pytest.raises(qpa.InvalidInputError, match="zero module"):
        algebra.zero_module().stable_endomorphism_quotient()


@pytest.mark.parametrize("characteristic", [2, 5])
def test_endomorphism_as_quiver_algebra_follows_qpa(characteristic: int) -> None:
    field = qpa.Field.prime(characteristic)
    dual_numbers = qpa.QuotientAlgebra.truncated(
        field, qpa.Quiver(["v"], [("v", "v", "x")]), 2
    )
    regular = dual_numbers.indecomposable_projectives()[0]

    result = regular.endomorphism_as_quiver_algebra()
    assert result.endomorphism_algebra.dimension == 2
    assert result.algebra.is_quotient
    assert result.algebra.quiver.adjacency_matrix == [[1]]
    assert result.adjacency_matrix == [[1]]
    assert result.quiver.adjacency_matrix == [[1]]
    assert result.is_quotient
    assert result.path_algebra is None
    assert result.quotient_algebra is not None
    assert result.quotient_algebra.dimension() == "2"

    repeated_simple = qpa.PathAlgebra(
        field, qpa.Quiver(["v"], [])
    ).module_from_integer_matrices([2], [])
    with pytest.raises(qpa.InvalidInputError, match="not basic"):
        repeated_simple.endomorphism_as_quiver_algebra()


@pytest.mark.parametrize(
    "field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)]
)
def test_basic_endomorphism_as_quiver_algebra_follows_qpa(field: qpa.Field) -> None:
    algebra = qpa.PathAlgebra(
        field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    first_simple, second_simple = algebra.simple_modules()
    first_projective = algebra.indecomposable_projectives()[0]
    modules = [first_simple, first_projective, second_simple]

    result = qpa.MatrixModule.endomorphism_as_quiver_algebra_of_basic(modules)
    assert result.search_bound == 2
    assert result.algebra.is_quotient
    assert result.final_search_length == 2
    assert not result.relations_may_be_incomplete
    assert result.quiver.adjacency_matrix == [
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
    ]
    assert result.is_quotient
    assert result.path_algebra is None
    assert result.quotient_algebra is not None
    assert result.quotient_algebra.dimension() == "5"

    bounded = qpa.MatrixModule.endomorphism_as_quiver_algebra_of_basic(modules, 1)
    assert bounded.search_bound == 1
    assert bounded.final_search_length == 2
    assert bounded.relations_may_be_incomplete
    assert bounded.quiver.adjacency_matrix == result.quiver.adjacency_matrix
    assert bounded.quotient_algebra is not None
    assert bounded.quotient_algebra.dimension() == "5"

    with pytest.raises(qpa.InvalidInputError, match="positive"):
        qpa.MatrixModule.endomorphism_as_quiver_algebra_of_basic(modules, 0)


def test_endomorphism_quiver_algebra_validation_follows_qpa() -> None:
    with pytest.raises(qpa.InvalidInputError, match="nonempty"):
        qpa.MatrixModule.endomorphism_as_quiver_algebra_of_basic([])

    field = qpa.Field.rational()
    dual_numbers = qpa.QuotientAlgebra.truncated(
        field, qpa.Quiver(["v"], [("v", "v", "x")]), 2
    )
    with pytest.raises(qpa.InvalidInputError, match="finite field"):
        dual_numbers.indecomposable_projectives()[
            0
        ].endomorphism_as_quiver_algebra()


def _structure_constants_tensor(
    field: qpa.Field,
    dimension: int,
    nonzero_products: list[tuple[int, int, int]],
) -> list[list[list[qpa.Scalar]]]:
    zero = field.element(0)
    one = field.element(1)
    products = [
        [[zero for _ in range(dimension)] for _ in range(dimension)]
        for _ in range(dimension)
    ]
    for left, right, output in nonzero_products:
        products[left][right][output] = one
    return products


def _representative_rows(
    rows: list[list[qpa.Scalar]],
) -> list[list[int]]:
    return [[coefficient.representative for coefficient in row] for row in rows]


def test_structure_constants_algebra_operations_follow_qpa() -> None:
    field2 = qpa.Field.prime(2)
    upper = qpa.StructureConstantsAlgebra(
        field2,
        _structure_constants_tensor(
            field2,
            3,
            [(0, 0, 0), (0, 1, 1), (1, 2, 1), (2, 2, 2)],
        ),
    )

    basis = upper.basis_elements
    assert _representative_rows([element.coordinates for element in basis]) == [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
    ]
    assert [coefficient.representative for coefficient in upper.zero().coordinates] == [
        0,
        0,
        0,
    ]
    assert basis[0] * basis[1] == basis[1]
    assert basis[1] * basis[2] == basis[1]

    identity = upper.find_multiplicative_identity()
    assert identity is not None
    assert [coefficient.representative for coefficient in identity.coordinates] == [
        1,
        0,
        1,
    ]
    center = upper.center()
    radical = upper.radical()
    assert _representative_rows(center.basis) == [[1, 0, 1]]
    assert _representative_rows(radical.basis) == [[0, 1, 0]]
    assert center.ambient.dimension == 3
    assert center.algebra.dimension == 1
    assert [
        coefficient.representative
        for coefficient in center.lift_coordinates([field2.element(1)])
    ] == [1, 0, 1]
    assert [
        _representative_rows(term) for term in upper.radical_series
    ] == [
        [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        [[0, 1, 0]],
        [],
    ]
    assert _representative_rows(upper.power_subalgebra(1).basis) == [[1, 0, 1]]
    assert [_representative_rows(block.basis) for block in upper.blocks()] == [
        [[1, 0, 1], [0, 0, 1], [0, 1, 0]]
    ]
    assert not upper.is_semisimple
    assert upper.is_radical_square_zero
    assert upper.is_basic is True
    assert upper.is_elementary is True

    zero_algebra = qpa.StructureConstantsAlgebra(
        field2, _structure_constants_tensor(field2, 1, [])
    )
    assert zero_algebra.find_multiplicative_identity() is None

    rational = qpa.Field.rational()
    rational_dual = qpa.StructureConstantsAlgebra(
        rational,
        _structure_constants_tensor(
            rational, 2, [(0, 0, 0), (0, 1, 1), (1, 0, 1)]
        ),
    )
    assert rational_dual.is_basic is None
    assert rational_dual.is_elementary is None

    field5 = qpa.Field.prime(5)
    shifted_products = _structure_constants_tensor(field5, 3, [])
    shifted_products[0][0][0] = field5.element(1)
    shifted_products[0][1][2] = field5.element(2)
    shifted_products[0][2][2] = field5.element(1)
    shifted_products[1][1][1] = field5.element(1)
    shifted_products[2][1][2] = field5.element(1)
    shifted = qpa.StructureConstantsAlgebra(field5, shifted_products)
    assert _representative_rows(
        [element.coordinates for element in shifted.idempotents_for_decomposition()]
    ) == [[0, 1, 0], [1, 0, 3]]

    matrix_products = _structure_constants_tensor(field2, 4, [])
    for left_row in range(2):
        for left_column in range(2):
            for right_row in range(2):
                for right_column in range(2):
                    if left_column == right_row:
                        left = 2 * left_row + left_column
                        right = 2 * right_row + right_column
                        output = 2 * left_row + right_column
                        matrix_products[left][right][output] = field2.element(1)
    matrix_algebra = qpa.StructureConstantsAlgebra(field2, matrix_products)
    primitive = matrix_algebra.primitive_idempotents()
    primitive_rows = _representative_rows(
        [element.coordinates for element in primitive]
    )
    assert len(primitive) == 2
    assert [sum(row[column] for row in primitive_rows) % 2 for column in range(4)] == [
        1,
        0,
        0,
        1,
    ]
    assert all(element * element == element for element in primitive)
    assert primitive[0] * primitive[1] == matrix_algebra.zero()
    assert primitive[1] * primitive[0] == matrix_algebra.zero()


def test_opposite_structure_constants_algebra_follows_qpa() -> None:
    field5 = qpa.Field.prime(5)
    upper = qpa.StructureConstantsAlgebra(
        field5,
        _structure_constants_tensor(
            field5,
            3,
            [(0, 0, 0), (0, 1, 1), (1, 2, 1), (2, 2, 2)],
        ),
    )
    underlying_basis = upper.basis_elements
    upper.bind_named_element("y", underlying_basis[1])

    opposite = upper.opposite()
    another_opposite = upper.opposite()
    assert isinstance(opposite, qpa.OppositeStructureConstantsAlgebra)
    assert opposite.field == field5
    assert opposite.dimension == 3
    assert opposite.is_opposite_algebra
    assert not opposite.is_identical_to_underlying
    assert opposite.is_finite_dimensional
    assert opposite.underlying_algebra.dimension == upper.dimension
    assert opposite.multiplication_algebra.dimension == upper.dimension

    basis = opposite.basis_elements
    canonical = opposite.canonical_basis_elements
    assert _representative_rows([element.coordinates for element in basis]) == [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
    ]
    assert _representative_rows([element.coordinates for element in canonical]) == [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
    ]
    x, y, z = basis
    assert all(element.is_opposite_algebra_element for element in basis)
    assert _representative_rows(
        [(x * y).coordinates, (y * x).coordinates, (y * z).coordinates, (z * y).coordinates]
    ) == [[0, 0, 0], [0, 1, 0], [0, 0, 0], [0, 1, 0]]
    assert [coefficient.representative for coefficient in (y + z).coordinates] == [
        0,
        1,
        1,
    ]
    assert [coefficient.representative for coefficient in y.scale(field5.element(2)).coordinates] == [
        0,
        2,
        0,
    ]
    assert [coefficient.representative for coefficient in (-y).coordinates] == [0, 4, 0]
    assert [
        coefficient.representative for coefficient in y.zero_same_parent().coordinates
    ] == [0, 0, 0]
    assert [coefficient.representative for coefficient in opposite.coefficients(y)] == [
        0,
        1,
        0,
    ]
    assert y.external_representation == y.coordinates
    assert y.underlying_element == underlying_basis[1]
    assert opposite.element_from_underlying(underlying_basis[1]) == y
    assert opposite.element([field5.element(0), field5.element(1), field5.element(0)]) == y
    assert opposite.element_by_name("y") == y
    assert opposite.element_by_name("missing") is None
    assert y in opposite
    assert not another_opposite.contains(y)
    assert not (x < y)
    assert y < x
    assert y.compare(x) == -1
    identity = opposite.one()
    assert identity is not None
    assert [coefficient.representative for coefficient in identity.coordinates] == [1, 0, 1]
    assert repr(y) == "op( v.2 )"
    assert repr(y + z) == "op( v.2+v.3 )"
    assert repr(y.scale(field5.element(2))) == "op( (Z(5))*v.2 )"
    assert repr(-y) == "op( (Z(5)^2)*v.2 )"
    assert repr(opposite.zero()) == "op( 0*v.1 )"

    rational = qpa.Field.rational()
    rational_upper = qpa.StructureConstantsAlgebra(
        rational,
        _structure_constants_tensor(
            rational,
            3,
            [(0, 0, 0), (0, 1, 1), (1, 2, 1), (2, 2, 2)],
        ),
    ).opposite()
    rational_x, rational_y, _ = rational_upper.basis_elements
    assert repr(rational_x + rational_y) == "op( v.1+v.2 )"
    assert repr(rational_y.scale(rational.element(2))) == "op( (2)*v.2 )"
    assert repr(rational_y.scale(rational.element(1, 2))) == "op( (1/2)*v.2 )"
    assert repr(-rational_y) == "op( (-1)*v.2 )"
    assert repr(rational_x + rational_y.scale(rational.element(-2))) == (
        "op( v.1+(-2)*v.2 )"
    )

    dual_numbers = qpa.StructureConstantsAlgebra(
        field5,
        _structure_constants_tensor(
            field5, 2, [(0, 0, 0), (0, 1, 1), (1, 0, 1)]
        ),
    ).opposite()
    assert dual_numbers.is_identical_to_underlying
    assert not dual_numbers.is_opposite_algebra
    assert not dual_numbers.basis_elements[0].is_opposite_algebra_element


def test_algebra_as_quiver_algebra_follows_qpa() -> None:
    field2 = qpa.Field.prime(2)
    upper_tensor = _structure_constants_tensor(
        field2,
        3,
        [(0, 0, 0), (0, 1, 1), (1, 2, 1), (2, 2, 2)],
    )
    upper = qpa.StructureConstantsAlgebra(field2, upper_tensor)
    assert upper.field == field2
    assert upper.dimension == 3
    assert upper.multiplication_tensor == upper_tensor

    converted = upper.as_quiver_algebra()
    assert isinstance(converted, qpa.AlgebraAsQuiverAlgebraResult)
    assert converted.loewy_length == 2
    assert converted.algebra.quiver.adjacency_matrix == [[0, 0], [1, 0]]
    assert not converted.algebra.is_quotient
    assert converted.algebra.quotient_algebra is None
    assert converted.algebra.path_algebra is not None
    assert len(converted.algebra.path_algebra.canonical_basis()) == 3
    assert [
        [coefficient.representative for coefficient in row]
        for row in converted.generator_images
    ] == [[0, 0, 1], [1, 0, 0], [0, 1, 0]]

    field5 = qpa.Field.prime(5)
    dual_numbers = qpa.StructureConstantsAlgebra(
        field5,
        _structure_constants_tensor(
            field5, 2, [(0, 0, 0), (0, 1, 1), (1, 0, 1)]
        ),
    )
    dual_converted = dual_numbers.as_quiver_algebra()
    assert dual_converted.loewy_length == 2
    assert dual_converted.algebra.quiver.adjacency_matrix == [[1]]
    assert dual_converted.algebra.is_quotient
    assert dual_converted.algebra.path_algebra is None
    assert dual_converted.algebra.quotient_algebra is not None
    assert dual_converted.algebra.quotient_algebra.dimension() == "2"
    assert [
        [coefficient.representative for coefficient in row]
        for row in dual_converted.generator_images
    ] == [[1, 0], [0, 1]]

    rational = qpa.Field.rational()
    rational_dual = qpa.StructureConstantsAlgebra(
        rational,
        _structure_constants_tensor(
            rational, 2, [(0, 0, 0), (0, 1, 1), (1, 0, 1)]
        ),
    )
    with pytest.raises(qpa.InvalidInputError, match="finite field"):
        rational_dual.as_quiver_algebra()

    product = qpa.StructureConstantsAlgebra(
        field2,
        _structure_constants_tensor(field2, 2, [(0, 0, 0), (1, 1, 1)]),
    )
    with pytest.raises(qpa.InvalidInputError, match="indecomposable"):
        product.as_quiver_algebra()

    with pytest.raises(qpa.InvalidInputError, match="tensor"):
        qpa.StructureConstantsAlgebra(field2, upper_tensor[:-1])


@pytest.mark.parametrize(
    "field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)]
)
def test_adr_algebra_of_algebra_follows_qpa(field: qpa.Field) -> None:
    algebra = qpa.PathAlgebra(
        field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    result = algebra.adr_algebra()
    assert isinstance(result, qpa.AdrAlgebraOfAlgebraResult)
    assert result.is_finite_global_dimension_algebra
    assert result.algebra.is_quotient
    assert result.algebra.quiver.adjacency_matrix == [
        [0, 0, 1],
        [1, 0, 0],
        [0, 0, 0],
    ]
    assert result.algebra.quotient_algebra is not None
    assert result.algebra.quotient_algebra.dimension() == "5"

    dual_numbers = qpa.QuotientAlgebra.truncated(
        field, qpa.Quiver(["v"], [("v", "v", "x")]), 2
    )
    dual_result = dual_numbers.adr_algebra()
    assert dual_result.is_finite_global_dimension_algebra
    assert dual_result.algebra.quiver.adjacency_matrix == [[0, 1], [1, 0]]
    assert dual_result.algebra.quotient_algebra is not None
    assert dual_result.algebra.quotient_algebra.dimension() == "5"


@pytest.mark.parametrize("characteristic", [2, 5])
def test_adr_algebra_of_basic_module_follows_qpa(characteristic: int) -> None:
    field = qpa.Field.prime(characteristic)
    algebra = qpa.PathAlgebra(
        field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    first_projective = algebra.indecomposable_projectives()[0]
    result = qpa.MatrixModule.adr_algebra_of_basic([first_projective])
    assert isinstance(result, qpa.QuiverAlgebra)
    assert result.is_quotient
    assert result.quiver.adjacency_matrix == [[0, 1], [0, 0]]
    assert result.quotient_algebra is not None
    assert result.quotient_algebra.dimension() == "3"


def test_adr_algebra_validation_follows_qpa() -> None:
    with pytest.raises(qpa.InvalidInputError, match="empty"):
        qpa.MatrixModule.adr_algebra_of_basic([])

    rational_algebra = qpa.PathAlgebra(
        qpa.Field.rational(),
        qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")]),
    )
    with pytest.raises(qpa.InvalidInputError, match="finite field"):
        qpa.MatrixModule.adr_algebra_of_basic(
            [rational_algebra.indecomposable_projectives()[0]]
        )

    loop_algebra = qpa.PathAlgebra(
        qpa.Field.prime(2), qpa.Quiver(["v"], [("v", "v", "x")])
    )
    with pytest.raises(qpa.InvalidInputError, match="finite-dimensional"):
        loop_algebra.adr_algebra()
    with pytest.raises(qpa.InvalidInputError, match="finite-dimensional"):
        qpa.MatrixModule.adr_algebra_of_basic([loop_algebra.simple_modules()[0]])


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)])
def test_projective_matrix_conversions_follow_qpa(field: qpa.Field) -> None:
    algebra = qpa.PathAlgebra(field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")]))
    projectives = qpa.MatrixModule.direct_sum(algebra.indecomposable_projectives())
    assert projectives is not None
    matrix = projectives.identity().projective_matrix()
    assert len(matrix) == 2 and all(len(row) == 2 for row in matrix)
    assert matrix[0][0].external_representation[0][0] == [1]
    assert matrix[0][1].is_zero
    assert matrix[1][0].is_zero
    assert matrix[1][1].external_representation[0][0] == [2]
    roundtrip = algebra.homomorphism_of_projectives_from_matrix(
        matrix, [1, 2], [1, 2]
    )
    assert roundtrip == projectives.identity()

    quotient = algebra.quotient([])
    quotient_projectives = qpa.MatrixModule.direct_sum(
        quotient.indecomposable_projectives()
    )
    assert quotient_projectives is not None
    quotient_matrix = quotient_projectives.identity().projective_matrix()
    quotient_entries = [
        [quotient.element(element) for element in row] for row in quotient_matrix
    ]
    quotient_roundtrip = quotient.homomorphism_of_projectives_from_matrix(
        quotient_entries, [1, 2], [1, 2]
    )
    assert quotient_roundtrip == quotient_projectives.identity()


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)])
def test_projective_generators_rejects_and_map_provenance_follow_qpa(
    field: qpa.Field,
) -> None:
    algebra = qpa.PathAlgebra(field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")]))
    projective = algebra.indecomposable_projectives()[0]
    first_simple, second_simple = algebra.simple_modules()
    cover = first_simple.projective_cover()

    kernel = cover.kernel_inclusion()
    image = cover.image_inclusion()
    cokernel = cover.cokernel_projection()
    assert kernel.kernel_of == cover
    assert image.image_of == cover
    assert cokernel.cokernel_of == cover
    assert cover.kernel_of is None
    assert cover.image_of is None
    assert cover.cokernel_of is None

    projective_basis = projective.basis_elements()
    generator_map = projective_basis[0].hom_from_projective()
    assert generator_map is not None and generator_map.is_isomorphism
    assert generator_map.source.dimension_vector == [1, 1]
    assert generator_map.range == projective
    assert projective_basis[0].add(projective_basis[1]).hom_from_projective() is None

    disjoint_reject = first_simple.reject_of(second_simple)
    assert disjoint_reject.is_isomorphism
    projective_reject = projective.reject_of(first_simple)
    assert projective_reject.source.dimension_vector == [0, 1]
    assert projective_reject.range == projective

    repeated = qpa.MatrixModule.direct_sum([first_simple, first_simple])
    assert repeated is not None
    if field.characteristic is None:
        with pytest.raises(qpa.InvalidInputError, match="finite field"):
            repeated.simple_submodule_inclusions()
    else:
        expected_count = field.characteristic + 1
        simple_inclusions = repeated.simple_submodule_inclusions()
        assert len(simple_inclusions) == expected_count
        assert all(inclusion.source.dimension_vector == [1, 0] for inclusion in simple_inclusions)
        assert all(inclusion.range == repeated for inclusion in simple_inclusions)


@pytest.mark.parametrize("characteristic", [2, 5])
def test_all_submodules_follow_qpa(characteristic: int) -> None:
    field = qpa.Field.prime(characteristic)
    algebra = qpa.PathAlgebra(field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")]))
    zero_result = algebra.zero_module().all_submodules()
    assert zero_result.is_flat
    assert zero_result.flat is not None and len(zero_result.flat) == 1
    assert zero_result.by_dimension is None

    simple_result = algebra.simple_modules()[0].all_submodules()
    assert simple_result.is_flat
    assert simple_result.flat is not None and len(simple_result.flat) == 2
    assert [inclusion.source.dimension for inclusion in simple_result.flattened] == [0, 1]

    semisimple = algebra.module_from_integer_matrices(
        [2, 1], [("a", [[0], [0]])]
    )
    result = semisimple.all_submodules()
    assert not result.is_flat and result.flat is None
    assert result.by_dimension is not None
    assert [len(level) for level in result.by_dimension] == [
        1,
        characteristic + 2,
        characteristic + 2,
        1,
    ]
    assert len(result.flattened) == 2 * characteristic + 6
    assert all(inclusion.range == semisimple for inclusion in result.flattened)

    rational = qpa.PathAlgebra(
        qpa.Field.rational(), qpa.Quiver(["v"], [])
    ).simple_modules()[0]
    with pytest.raises(qpa.InvalidInputError, match="finite field"):
        rational.all_submodules()


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)])
def test_matrix_submodule_natural_mapping_follows_qpa(field: qpa.Field) -> None:
    algebra = qpa.PathAlgebra(field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")]))
    module = algebra.module_from_integer_matrices(
        [3, 3],
        [("a", [[1, 0, 0], [0, 1, 0], [0, 0, 1]])],
    )
    basis = module.basis_elements()
    submodule = module.submodule([basis[0].add(basis[4]), module.zero_element()])
    mapping = submodule.natural_quotient_mapping()

    assert mapping.source == module
    assert mapping.range.dimension_vector == [2, 1]
    assert len(mapping.basis_images) == module.dimension
    assert mapping.is_surjective
    assert mapping.is_identical(mapping)
    assert [mapping.image(element) for element in basis] == mapping.basis_images
    # QPA's mutable `base` object leaks the last nonzero vertex-1 coordinates
    # into later basis images; retaining this verifies the literal source method.
    assert (
        mapping.basis_images[2].external_representation
        == mapping.basis_images[3].external_representation
    )

    zero_submodule = module.submodule([])
    with pytest.raises(qpa.InvalidInputError, match="zero submodule"):
        zero_submodule.natural_quotient_mapping()


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)])
def test_restriction_via_algebra_homomorphism_follows_qpa(field: qpa.Field) -> None:
    algebra = qpa.PathAlgebra(field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")]))
    quotient = algebra.quotient([])
    quotient_projection = qpa.AlgebraHomomorphism.path_to_quotient(
        algebra,
        quotient,
        [quotient.element(generator) for generator in algebra.generators],
    )
    assert quotient_projection.source_kind == "path"
    assert quotient_projection.range_kind == "quotient"
    assert quotient_projection.generator_count == 3
    assert quotient_projection.uses_generators_with_one
    assert quotient_projection.source_path_algebra is not None
    assert quotient_projection.source_path_algebra.is_identical(algebra)
    assert quotient_projection.range_quotient_algebra is not None
    assert quotient_projection.range_quotient_algebra.is_identical(quotient)

    quotient_module = quotient.module_from_integer_matrices([1, 1], [("a", [[1]])])
    expected = algebra.module_from_integer_matrices([1, 1], [("a", [[1]])])
    restricted = quotient_module.restrict_scalars(quotient_projection)
    restricted_identity = quotient_module.identity().restrict_scalars(quotient_projection)
    assert restricted == expected
    assert restricted_identity == expected.identity()
    assert restricted_identity.source == restricted

    path_identity = qpa.AlgebraHomomorphism.path_to_path(
        algebra, algebra, algebra.generators
    )
    quotient_identity = qpa.AlgebraHomomorphism.quotient_to_quotient(
        quotient, quotient, [quotient.one(), *quotient.generators]
    )
    quotient_identity_with_one = qpa.AlgebraHomomorphism.quotient_to_quotient(
        quotient, quotient, quotient.generators, generators_with_one=True
    )
    quotient_lift = qpa.AlgebraHomomorphism.quotient_to_path(
        quotient, algebra, [algebra.one(), *algebra.generators]
    )
    assert path_identity.source_kind == path_identity.range_kind == "path"
    assert quotient_identity.generator_count == 4
    assert not quotient_identity.uses_generators_with_one
    assert quotient_identity_with_one.generator_count == 3
    assert quotient_identity_with_one.uses_generators_with_one
    assert quotient_lift.source_kind == "quotient"
    assert quotient_lift.range_kind == "path"
    assert expected.restrict_scalars(quotient_lift).dimension_vector == [1, 1]

    with pytest.raises(qpa.InvalidInputError, match="3 images for 4"):
        qpa.AlgebraHomomorphism.quotient_to_quotient(
            quotient, quotient, quotient.generators
        )
