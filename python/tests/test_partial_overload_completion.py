from __future__ import annotations

from fractions import Fraction

import pytest

import qpa


def _value(scalar: qpa.Scalar) -> Fraction | int:
    if scalar.representative is not None:
        return scalar.representative
    assert scalar.numerator is not None and scalar.denominator is not None
    return Fraction(int(scalar.numerator), int(scalar.denominator))


def _matrices(
    matrices: list[list[list[qpa.Scalar]]],
) -> list[list[list[Fraction | int]]]:
    return [
        [[_value(entry) for entry in row] for row in matrix]
        for matrix in matrices
    ]


def _matrix_support(matrix: list[list[qpa.Scalar]]) -> list[int]:
    return [
        position
        for position, entry in enumerate(
            (entry for row in matrix for entry in row), start=1
        )
        if _value(entry) != 0
    ]


@pytest.mark.parametrize(
    "field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)]
)
def test_checked_module_constructor_overloads_follow_qpa(field: qpa.Field) -> None:
    quiver = qpa.Quiver(
        ["v1", "v2", "v3"],
        [("v1", "v2", "a"), ("v2", "v3", "b")],
    )
    path = qpa.PathAlgebra(field, quiver)
    algebra = qpa.QuiverAlgebra.from_path_algebra(path)
    a_matrix = [[field.element(1), field.element(2)]]
    b_matrix = [[field.element(3)], [field.element(4)]]

    modules = [
        algebra.right_module_over_path_algebra_with_dimensions(
            [1, 2, 1], [("a", a_matrix), ("b", b_matrix)]
        ),
        algebra.right_module_over_path_algebra([a_matrix, b_matrix]),
        algebra.right_module_over_path_algebra_named(
            [("b", b_matrix), ("a", a_matrix)]
        ),
        algebra.right_module_over_path_algebra_indexed(
            [(quiver.arrow("b"), b_matrix), (quiver.arrow("a"), a_matrix)]
        ),
    ]
    expected = [
        [[_value(field.element(1)), _value(field.element(2))]],
        [[_value(field.element(3))], [_value(field.element(4))]],
    ]
    for module in modules:
        assert module is not None
        assert module.dimension_vector == [1, 2, 1]
        assert _matrices(module.arrow_matrices) == expected

    zero_source = algebra.right_module_over_path_algebra_named(
        [("a", [0, 2]), ("b", [2, 0])]
    )
    assert zero_source is not None
    assert zero_source.dimension_vector == [0, 2, 0]

    isolated_quiver = qpa.Quiver(
        ["i1", "i2", "i3"], [("i1", "i2", "c")]
    )
    isolated = qpa.QuiverAlgebra.from_path_algebra(
        qpa.PathAlgebra(field, isolated_quiver)
    ).right_module_over_path_algebra([[[field.element(1)]]])
    assert isolated is not None
    assert isolated.dimension_vector == [1, 1, 0]


def test_checked_constructor_rejects_foreign_field_and_returns_qpa_fail() -> None:
    field3 = qpa.Field.prime(3)
    field5 = qpa.Field.prime(5)
    path = qpa.PathAlgebra(
        field3, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    algebra = qpa.QuiverAlgebra.from_path_algebra(path)
    with pytest.raises(qpa.InvalidInputError, match=r"GF\(3\)"):
        algebra.right_module_over_path_algebra([[[field5.element(1)]]])

    bound_path = qpa.PathAlgebra(
        field3,
        qpa.Quiver(
            ["q1", "q2", "q3"],
            [("q1", "q2", "x"), ("q2", "q3", "y")],
        ),
    )
    quotient = bound_path.quotient_high_level(
        [bound_path.generator("x") * bound_path.generator("y")]
    )
    bound = qpa.QuiverAlgebra.from_quotient_algebra(quotient)
    one = field3.element(1)
    assert bound.right_module_over_path_algebra([[[one]], [[one]]]) is None


@pytest.mark.parametrize(
    "field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)]
)
def test_list_intersection_and_sum_follow_qpa_folds(field: qpa.Field) -> None:
    path = qpa.PathAlgebra(
        field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    simple1, simple2 = path.simple_modules()
    target = qpa.MatrixModule.direct_sum([simple1, simple1, simple2])
    assert target is not None
    inclusions = target.direct_sum_inclusions()

    intersection = qpa.ModuleHomomorphism.intersection_list(inclusions)
    binary_intersection = inclusions[0].intersection(inclusions[1])[2].intersection(
        inclusions[2]
    )[2]
    assert intersection == binary_intersection
    assert intersection.source.is_zero
    assert intersection.range == target
    assert qpa.ModuleHomomorphism.intersection_list(
        [inclusions[0]]
    ).is_identical(inclusions[0])

    total = qpa.ModuleHomomorphism.sum_list(inclusions)
    binary_total = inclusions[0].sum(inclusions[1])[0].sum(inclusions[2])[0]
    assert total == binary_total
    assert total.source.dimension_vector == target.dimension_vector
    assert total.range == target
    assert qpa.ModuleHomomorphism.sum_list([inclusions[0]]).is_identical(
        inclusions[0]
    )

    with pytest.raises(qpa.InvalidInputError, match="empty"):
        qpa.ModuleHomomorphism.intersection_list([])
    with pytest.raises(qpa.InvalidInputError, match="empty"):
        qpa.ModuleHomomorphism.sum_list([])


def test_product_decomposition_supports_qpa_one_based_indexing() -> None:
    left = qpa.Quiver(["u1", "u2"], [("u1", "u2", "a")])
    right = qpa.Quiver(["v1", "v2"], [("v1", "v2", "b")])
    decomposition = left.product(right).decomposition
    assert decomposition[1].is_identical(left)
    assert decomposition[2].is_identical(right)
    with pytest.raises(qpa.InvalidInputError, match="outside 1..=2"):
        _ = decomposition[0]


def test_display_overloads_preserve_qpa_rendering() -> None:
    field = qpa.Field.rational()
    path = qpa.PathAlgebra(
        field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    module = path.module_from_integer_matrices([1, 1], [("a", [[1]])])
    module_print = (
        "<Module over <Rationals[<quiver with 2 vertices and 1 arrows>]>"
        " with dimension vector [ 1, 1 ]>"
    )
    assert module.view() == "<[ 1, 1 ]>"
    assert module.print_obj() == module_print
    assert module.display() == (
        f"{module_print} and linear maps given by\n"
        "for arrow a:\n"
        "[ [  1 ] ]\n"
    )

    identity = module.identity()
    map_print = f"<{module_print} ---> {module_print}>"
    assert identity.view() == "<<[ 1, 1 ]> ---> <[ 1, 1 ]>>"
    assert identity.print_obj() == map_print
    assert identity.display() == (
        f"{map_print}\n"
        "with linear map for vertex number 1:\n"
        "[ [  1 ] ]\n"
        "linear map for vertex number 2:\n"
        "[ [  1 ] ]\n"
    )

    decomposition = path.quiver.product(path.quiver).decomposition
    assert decomposition.view() == "<quiver product decomposition>"
    assert decomposition.print_obj() == "<quiver product decomposition>"
    assert decomposition.display() == "<quiver product decomposition>"


def test_print_obj_and_view_obj_overloads_preserve_qpa_rendering() -> None:
    field = qpa.Field.rational()
    quiver = qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    path = qpa.PathAlgebra(field, quiver)
    arrow = path.generator("a")

    assert quiver.print_obj() == 'Quiver( ["v1","v2"], [["v1","v2","a"]] )'
    assert quiver.view() == "<quiver with 2 vertices and 1 arrows>"
    assert quiver.vertex("v1").print_obj() == "v1"
    assert quiver.arrow("a").print_obj() == "a"
    assert quiver.arrow_path("a").print_obj() == "a"
    assert quiver.zero_path().print_obj() == "0"
    assert quiver.ordering.print_obj() == "<length left lexicographic ordering>"
    assert quiver.ordering.view() == "<length left lexicographic ordering>"

    path_view = "<Rationals[<quiver with 2 vertices and 1 arrows>]>"
    assert path.view() == path_view
    assert path.print_obj() == (
        "<Path algebra of the quiver <quiver with 2 vertices and 1 arrows> "
        "over the field Rationals>"
    )
    assert arrow.print_obj() == "(1)*a"
    assert path.ideal([arrow]).view() == (
        f"<two-sided ideal in {path_view}, (1 generator)>"
    )

    quotient = path.quotient([])
    ideal_view = f"<two-sided ideal in {path_view}, (0 generators)>"
    assert quotient.view() == f"<Rationals[<quiver with 2 vertices and 1 arrows>]/{ideal_view}>"
    assert quotient.print_obj() == (
        f"<A quotient of the path algebra {path_view} modulo the ideal {ideal_view}>"
    )
    assert quotient.generator("a").print_obj() == "[(1)*a]"

    loop = qpa.PathAlgebra(
        field, qpa.Quiver(["w"], [("w", "w", "x")])
    )
    x = loop.generator("x")
    basis = loop.groebner_basis([x * x]).completely_reduced_basis()
    assert basis.view() == "<complete two-sided Groebner basis containing 1 elements>"
    assert basis.print_obj() == (
        "GroebnerBasis( <Path algebra of the quiver <quiver with 1 vertices "
        "and 1 arrows> over the field Rationals>, [ (1)*x^2 ] )"
    )

    assert qpa.StaticDictionary([[1], [2]]).view() == (
        "<static dictionary containing 2 patterns>"
    )
    assert qpa.QuiverStaticDictionary(quiver, [quiver.arrow_path("a")]).view() == (
        "<static dictionary for quivers containing 1 path>"
    )
    assert qpa.UnitForm([[2, -1], [-1, 2]]).view() == (
        "[ [ 2, -1 ], [ -1, 2 ] ]"
    )
    assert qpa.Poset(["a", "b"], [["a", "b"]]).view() == (
        "<A poset on 2 points>"
    )
    ar_quiver = qpa.ArQuiverNumerical.predefined("A3 zero")
    ar_view = "<ARQuiverNumerical with 5 indecomposables and 3 projectives>"
    assert ar_quiver.print_obj() == ar_view
    assert ar_quiver.view() == ar_view

    module = path.module_from_integer_matrices([1, 1], [("a", [[1]])])
    assert module.basis_elements()[0].print_obj() == "[ [ 1 ], [ 0 ] ]"
    category = qpa.ModuleCategory.from_path_algebra(path)
    assert category.print_obj() == "<cat: right modules over algebra>"
    assert category.view() == "<cat: right modules over algebra>"

    callback_category = qpa.cat(
        "toy",
        {
            "name": "toy category",
            "zeroObj": 0,
            "isZeroObj": lambda value: value == 0,
            "zeroMap": lambda source, range: (source, range, 0),
            "isZeroMapping": lambda value: value[2] == 0,
            "composeMaps": lambda left, right: (right[0], left[1], 0),
            "ker": lambda value: value[0],
            "im": lambda value: value[1],
            "isExact": lambda left, right: True,
            "objStr": str,
        },
    )
    assert callback_category.print_obj() == "<cat: toy category>"
    assert callback_category.view() == "<cat: toy category>"


@pytest.mark.parametrize(
    "field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)]
)
def test_path_vector_print_obj_overloads_cover_every_python_family(
    field: qpa.Field,
) -> None:
    path = qpa.PathAlgebra(
        field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    v2 = path.generator("v2")
    arrow = path.generator("a")

    uniform = path.right_projective_module_from_uniform_elements([arrow]).generators[0]
    assert uniform.print_obj() == "[ (1)*a ]"
    assert repr(uniform) == uniform.print_obj()

    projective = path.right_projective_module([1, 2])
    vector = projective.vector([arrow, v2])
    assert vector.print_obj() == "[ (1)*a, (1)*v2 ]"
    assert repr(vector) == vector.print_obj()

    module_basis = projective.right_groebner_basis()
    assert module_basis.view() == repr(module_basis)
    assert module_basis.view_without_basis() == (
        f"<Groebner basis of <right-module over <{field}[<quiver with 2 vertices "
        "and 1 arrows>]>>, ...>"
    )

    fp_module = projective.natural_quotient(
        projective.submodule([projective.zero()])
    ).range
    fp_vector = fp_module.generators[0]
    assert fp_vector.print_obj() == "[ [ (1)*v1, 0 ] ]"
    assert repr(fp_vector) == fp_vector.print_obj()

    quotient = path.quotient([])
    quotient_vector = quotient.right_projective_module([1, 2]).vector(
        [quotient.generator("a"), quotient.generator("v2")]
    )
    assert quotient_vector.print_obj() == "[ [(1)*a], [(1)*v2] ]"
    assert repr(quotient_vector) == quotient_vector.print_obj()


@pytest.mark.parametrize(
    "field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)]
)
def test_complex_print_obj_overloads_cover_qpa_dispatch_priority(
    field: qpa.Field,
) -> None:
    path = qpa.PathAlgebra(
        field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    projective = path.indecomposable_projectives()[0]
    first_simple = path.simple_modules()[0]

    zero = projective.zero_complex()
    projective_stalk = projective.stalk_complex(0)
    injective_stalk = first_simple.stalk_complex(0)
    generic = projective.finite_complex(
        0, [projective.top_projection(), projective.radical_inclusion()]
    )
    chain_map = generic.zero_chain_map_to(generic)

    assert zero.print_obj() == "0 -> 0"
    assert projective_stalk.print_obj() == "0 -> 0: P1 -> 0"
    assert injective_stalk.print_obj() == "0 -> 0: I1 -> 0"
    assert generic.print_obj() == (
        "0 -> 1:(0,1) -> 0:(1,1) -> -1:(1,0) -> 0"
    )
    assert chain_map.print_obj() == "<chain map>"
    assert repr(chain_map) == chain_map.print_obj()


@pytest.mark.parametrize(
    "field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)]
)
def test_loewy_length_algebra_overload_follows_qpa(field: qpa.Field) -> None:
    path = qpa.PathAlgebra(
        field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    assert path.loewy_length == 2
    assert qpa.QuiverAlgebra.from_path_algebra(path).loewy_length == 2

    quotient = qpa.QuotientAlgebra.truncated(
        field, qpa.Quiver(["w"], [("w", "w", "x")]), 2
    )
    assert quotient.loewy_length == 2
    assert qpa.QuiverAlgebra.from_quotient_algebra(quotient).loewy_length == 2


@pytest.mark.parametrize(
    "field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)]
)
def test_finite_type_decisions_preserve_qpa_values_and_messages(
    field: qpa.Field,
) -> None:
    path = qpa.PathAlgebra(
        field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    path_decision = path.finite_type_decision()
    assert path.is_finite_type_algebra
    assert path_decision.value is True
    assert path_decision.reason == "DynkinComponents"
    assert path_decision.messages == [
        "A_2\n",
        "Finite type!\nQuiver is a (union of) Dynkin quiver(s).\n",
    ]

    quotient = qpa.QuotientAlgebra.truncated(
        field, qpa.Quiver(["w"], [("w", "w", "x")]), 2
    )
    quotient_decision = quotient.finite_type_decision()
    assert quotient.is_finite_type_algebra is True
    assert quotient_decision.value is True
    assert quotient_decision.reason == "RadicalSquareZero"
    assert quotient_decision.messages == [
        "Quiver contains an (un)oriented cycle.\n",
        "A_2\n",
    ]


@pytest.mark.parametrize(
    "field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)]
)
def test_uniform_element_projective_family_preserves_qpa_behavior(
    field: qpa.Field,
) -> None:
    path = qpa.PathAlgebra(
        field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    v1, v2, arrow = path.generators
    module = path.right_projective_module_from_uniform_elements([arrow])
    assert isinstance(module, qpa.RightUniformPathAlgebraModule)
    assert module.algebra.is_identical(path)
    assert module.elements == [arrow]
    assert module.is_vertex_projective

    generator = module.generators[0]
    zero = module.zero()
    assert isinstance(generator, qpa.RightUniformPathAlgebraVector)
    assert generator.module.is_identical(module)
    assert generator.components == [arrow]
    assert generator.leading_position == 1
    assert generator.is_normal_form and not generator.is_zero
    assert generator.leading_term == generator
    assert generator.leading_coefficient == field.element(1)
    assert generator.leading_component == arrow
    assert generator.source_vertex == v1.leading_monomial
    assert generator.target_vertex == v2.leading_monomial
    assert generator.is_left_uniform
    assert generator.is_right_uniform
    assert generator.is_uniform
    assert generator.compare(generator) == 0
    assert generator.is_left_divisible(generator)

    assert zero.leading_position == 2 and zero.is_zero
    assert generator.act(v1).is_zero
    negative = -generator
    assert generator + negative == zero
    assert generator.scale_left(field.element(0)) == zero
    assert generator.scale_right(field.element(1)) == generator
    expected_two = 0 if field.characteristic == 2 else 2
    assert generator.scale_left(field.element(2)).components == [
        arrow.scale_left(field.element(expected_two))
    ]

    with pytest.raises(qpa.InvalidInputError, match="no vertices component"):
        generator.act(v2)
    with pytest.raises(qpa.InvalidInputError, match="no vertices component"):
        module.vector([arrow])
    with pytest.raises(qpa.InvalidInputError, match="no vertices component"):
        module.canonical_basis()
    with pytest.raises(qpa.InvalidInputError, match="no vertices component"):
        module.uniform_generators()

    duplicate = qpa.RightUniformPathAlgebraModule.from_elements(path, [arrow])
    assert duplicate.elements == module.elements
    assert not duplicate.is_identical(module)
    with pytest.raises(qpa.InvalidInputError, match="not uniform"):
        path.right_projective_module_from_uniform_elements([v1 + v2])
    with pytest.raises(qpa.InvalidInputError, match="cannot extract coefficients"):
        path.is_left_uniform_list_at_path([arrow], v1.leading_monomial)
    with pytest.raises(qpa.InvalidInputError, match="cannot extract coefficients"):
        path.is_right_uniform_list_at_path([arrow], v1.leading_monomial)


@pytest.mark.parametrize(
    "field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)]
)
def test_quotient_canonical_basis_object_and_coefficients_follow_qpa(
    field: qpa.Field,
) -> None:
    quotient = qpa.QuotientAlgebra.truncated(
        field, qpa.Quiver(["w"], [("w", "w", "x")]), 2
    )
    basis = quotient.canonical_basis_object()
    assert isinstance(basis, qpa.QuotientCanonicalBasis)
    assert basis.quotient.is_identical(quotient)
    assert basis.has_basis_vectors
    assert basis.basis_vectors == quotient.canonical_basis()
    assert basis.is_identical(quotient.canonical_basis_object())

    element = quotient.one().scale_left(field.element(2)) + quotient.generator(
        "x"
    ).scale_left(field.element(3))
    assert [_value(value) for value in basis.coefficients(element)] == [
        _value(field.element(2)),
        _value(field.element(3)),
    ]


@pytest.mark.parametrize(
    "field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)]
)
def test_generic_path_ring_hom_and_end_follow_qpa_hopf_oracle(
    field: qpa.Field,
) -> None:
    path = qpa.PathAlgebra(
        field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    module = path.module_from_integer_matrices(
        [2, 2], [("a", [[1, 0], [0, 1]])]
    )
    hom = path.hom(module, module)
    end = path.end(module)

    expected = [[1, 11], [2, 12], [5, 15], [6, 16]]
    assert isinstance(hom, qpa.GenericModuleHomSpace)
    assert isinstance(end, qpa.GenericModuleEndomorphismAlgebra)
    assert hom.field == field and end.field == field
    assert hom.source.is_identical(module) and hom.range.is_identical(module)
    assert end.module.is_identical(module)
    assert hom.dimension == end.dimension == 4
    assert [_matrix_support(mapping.matrix) for mapping in hom.generators] == expected
    assert [_matrix_support(mapping.matrix) for mapping in end.generators] == expected
    assert all(
        isinstance(mapping, qpa.GenericModuleHomomorphism)
        and mapping.is_algebra_module_homomorphism
        for mapping in hom.generators
    )

    coefficients = hom.coefficients(hom.generators[0].matrix)
    assert coefficients is not None
    assert [_value(value) for value in coefficients] == [1, 0, 0, 0]
    assert _matrix_support(end.one.matrix) == [1, 6, 11, 16]
    assert _matrix_support(end.zero.matrix) == []
    assert end.multiply(end.generators[1], end.generators[2]) == end.generators[0]
    assert end.multiply(end.generators[2], end.generators[1]) == end.generators[3]

    basis = module.basis_elements()
    submodule = module.submodule([basis[0] + basis[3]])
    canonical = submodule.canonical_basis_elements()
    sub_hom = path.hom(submodule, submodule)
    sub_end = path.end(submodule)
    to_submodule = path.hom(module, submodule)
    from_submodule = path.hom(submodule, module)
    assert sub_hom.source.is_identical(submodule)
    assert sub_hom.range.is_identical(submodule)
    assert sub_end.module.is_identical(submodule)
    assert sub_hom.dimension == sub_end.dimension == 3
    assert [_matrix_support(mapping.matrix) for mapping in sub_hom.generators] == [
        [1, 5],
        [8],
        [9],
    ]
    assert [_matrix_support(mapping.matrix) for mapping in to_submodule.generators] == [
        [1, 8],
        [4, 11],
    ]
    assert [
        _matrix_support(mapping.matrix) for mapping in from_submodule.generators
    ] == [[1, 7], [2, 8], [11], [12]]
    assert sub_hom.generators[0](canonical[0]) == canonical[0]
    assert sub_hom.generators[0](canonical[1]) == canonical[1]
    assert sub_hom.generators[0](canonical[2]).is_zero


def test_generic_path_ring_hom_and_end_preserve_qpa_failures() -> None:
    field = qpa.Field.rational()
    path = qpa.PathAlgebra(
        field, qpa.Quiver(["x", "y"], [("x", "y", "a")])
    )
    module = path.module_from_integer_matrices([1, 1], [("a", [[1]])])
    zero = path.zero_module()
    zero_submodule = module.submodule([])
    for operation in (
        lambda: path.hom(zero, module),
        lambda: path.hom(module, zero),
        lambda: path.end(zero),
        lambda: path.hom(zero_submodule, zero_submodule),
        lambda: path.end(zero_submodule),
    ):
        with pytest.raises(qpa.InvalidInputError):
            operation()

    arrowless = qpa.PathAlgebra(field, qpa.Quiver(["s"], []))
    arrowless_module = arrowless.module_from_integer_matrices([2], [])
    with pytest.raises(qpa.InvalidInputError):
        arrowless.hom(arrowless_module, arrowless_module)
    with pytest.raises(qpa.InvalidInputError):
        arrowless.end(arrowless_module)

    other_path = qpa.PathAlgebra(
        field, qpa.Quiver(["x", "y"], [("x", "y", "a")])
    )
    other_module = other_path.module_from_integer_matrices(
        [1, 1], [("a", [[1]])]
    )
    with pytest.raises(qpa.InvalidInputError, match="right acting domain"):
        path.hom(module, other_module)
    with pytest.raises(TypeError, match="MatrixModule or ModuleSubmodule"):
        path.hom(object(), module)


@pytest.mark.parametrize(
    "field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)]
)
def test_complex_and_chain_map_callback_half_specs_follow_qpa(
    field: qpa.Field,
) -> None:
    path = qpa.PathAlgebra(field, qpa.Quiver(["v"], []))
    module = path.module_from_integer_matrices([1], [])
    category = qpa.ModuleCategory.from_path_algebra(path)
    zero = module.zero_map_to(module)
    identity = module.identity()

    complex_position_calls: list[tuple[qpa.ModuleComplex, int]] = []

    def complex_position(
        complex: qpa.ModuleComplex, degree: int
    ) -> qpa.ModuleHomomorphism:
        complex_position_calls.append((complex, degree))
        return zero

    positioned = qpa.ModuleComplex.from_parts(
        category,
        0,
        [zero],
        positive=["pos", complex_position, False],
        negative=["pos", complex_position, False],
    )
    assert positioned.differential(5) == zero
    assert positioned.differential(-5) == zero
    assert [degree for _complex, degree in complex_position_calls] == [5, -5]
    assert all(
        complex.is_identical(positioned)
        for complex, _degree in complex_position_calls
    )

    next_calls: list[qpa.ModuleHomomorphism] = []

    def next_zero(
        previous: qpa.ModuleHomomorphism,
    ) -> qpa.ModuleHomomorphism:
        next_calls.append(previous)
        return previous

    generated = qpa.ModuleComplex.from_parts(
        category,
        0,
        [zero],
        positive=["next/repeat", next_zero, zero],
        negative=[zero],
    )
    assert generated.differential(10).is_zero
    assert next_calls

    periodic = qpa.ModuleComplex.from_parts(
        category, 0, [zero], positive=[zero], negative=[zero]
    )
    chain_position_calls: list[tuple[qpa.ModuleChainMap, int]] = []

    def chain_position(
        chain_map: qpa.ModuleChainMap, degree: int
    ) -> qpa.ModuleHomomorphism:
        chain_position_calls.append((chain_map, degree))
        return identity

    positioned_map = qpa.ModuleChainMap.from_parts(
        periodic,
        periodic,
        0,
        [identity],
        positive=["pos", chain_position, False],
        negative=["pos", chain_position, False],
    )
    assert positioned_map.morphism(7).is_isomorphism
    assert positioned_map.morphism(-7).is_isomorphism
    assert [degree for _map, degree in chain_position_calls] == [7, -7]
    assert all(
        chain_map.is_identical(positioned_map)
        for chain_map, _degree in chain_position_calls
    )

    next_map = qpa.ModuleChainMap.from_parts(
        periodic,
        periodic,
        0,
        [identity],
        positive=["next/repeat", lambda previous: previous, identity],
        negative=[identity],
    )
    assert next_map.morphism(10).is_isomorphism


@pytest.mark.parametrize(
    "field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)]
)
def test_combined_complex_callback_half_specs_follow_qpa(
    field: qpa.Field,
) -> None:
    path = qpa.PathAlgebra(field, qpa.Quiver(["v"], []))
    module = path.module_from_integer_matrices([1], [])
    category = qpa.ModuleCategory.from_path_algebra(path)
    zero = module.zero_map_to(module)
    identity = module.identity()
    periodic = qpa.ModuleComplex.from_parts(
        category, 0, [zero], positive=[zero], negative=[zero]
    )
    position_calls: list[int] = []

    def position(
        complex: qpa.ModuleComplex,
        incoming: qpa.ModuleChainMap,
        outgoing: qpa.ModuleChainMap,
        degree: int,
    ) -> list[qpa.ModuleHomomorphism]:
        assert isinstance(complex, qpa.ModuleComplex)
        assert isinstance(incoming, qpa.ModuleChainMap)
        assert isinstance(outgoing, qpa.ModuleChainMap)
        position_calls.append(degree)
        return [zero, identity, identity]

    combined, incoming, outgoing = qpa.complex_and_chain_maps(
        [periodic],
        [periodic],
        0,
        [[zero, identity, identity]],
        positive=["pos", position, False],
        negative=["pos", position, False],
    )
    assert combined.differential(6).is_zero
    assert incoming[0].morphism(6).is_isomorphism
    assert outgoing[0].morphism(6).is_isomorphism
    assert combined.differential(-6).is_zero
    assert incoming[0].morphism(-6).is_isomorphism
    assert outgoing[0].morphism(-6).is_isomorphism
    assert position_calls == list(range(1, 7)) * 3 + list(range(-1, -7, -1)) * 3

    initial = [zero, identity, identity]
    generated, generated_in, generated_out = qpa.complex_and_chain_maps(
        [periodic],
        [periodic],
        0,
        [initial],
        positive=["next/repeat", lambda previous: previous, initial],
        negative=[initial],
    )
    assert generated.differential(9).is_zero
    assert generated_in[0].morphism(9).is_isomorphism
    assert generated_out[0].morphism(9).is_isomorphism


def test_complex_callback_specs_reject_invalid_python_results() -> None:
    path = qpa.PathAlgebra(qpa.Field.rational(), qpa.Quiver(["v"], []))
    module = path.module_from_integer_matrices([1], [])
    zero = module.zero_map_to(module)
    complex = qpa.ModuleComplex.from_parts(
        qpa.ModuleCategory.from_path_algebra(path),
        0,
        [zero],
        positive=["pos", lambda _complex, _degree: object(), False],
        negative=[zero],
    )
    with pytest.raises(qpa.InvalidInputError, match="ModuleHomomorphism"):
        complex.differential(5)
