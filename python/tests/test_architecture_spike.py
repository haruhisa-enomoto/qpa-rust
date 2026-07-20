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
