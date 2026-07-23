from __future__ import annotations

from pathlib import Path

import pytest

import qpa


@pytest.mark.parametrize(
    ("prime", "degree", "modulus"),
    [
        (2, 2, [1, 1, 1]),
        (2, 3, [1, 1, 0, 1]),
        (3, 2, [2, 2, 1]),
    ],
)
def test_gap_conway_field_and_scalar_round_trip(
    prime: int, degree: int, modulus: list[int]
) -> None:
    field = qpa.Field.extension(prime, degree)
    generator = field.primitive_element()

    assert field.characteristic == prime
    assert field.degree == degree
    assert field.conway_polynomial == modulus
    assert generator.field == field
    assert hash(generator.field) == hash(field)
    assert generator.coefficients == [0, 1] + [0] * (degree - 2)
    assert field.element_from_coefficients(generator.coefficients) == generator
    assert hash(field.element_from_coefficients(generator.coefficients)) == hash(generator)
    assert field.element(1).coefficients == [1] + [0] * (degree - 1)
    assert all(value.is_zero for value in field.null_list(3))


def test_gf4_arithmetic_flows_through_path_algebra_elements() -> None:
    field = qpa.Field.extension(2, 2)
    generator = field.primitive_element()
    quiver = qpa.Quiver(["v"], [("v", "v", "x")])
    algebra = qpa.PathAlgebra(field, quiver)
    scalar_element = algebra.term(quiver.vertex_path("v"), generator)

    squared = scalar_element * scalar_element
    assert squared.terms[0][1].coefficients == [1, 1]
    assert squared + scalar_element == algebra.one()


def test_frobenius_and_canonical_subfield_embedding_follow_gap() -> None:
    field4 = qpa.Field.extension(2, 2)
    field16 = qpa.Field.extension(2, 4)
    generator4 = field4.primitive_element()
    generator16 = field16.primitive_element()

    assert field4.frobenius(generator4, 1).coefficients == [1, 1]
    assert field4.frobenius(generator4, -1).coefficients == [1, 1]
    assert field16.frobenius(generator16, 1).coefficients == [0, 0, 1, 0]
    assert field16.embed(generator4).coefficients == [0, 1, 1, 0]
    assert field16.embed(qpa.Field.prime(2).element(1)).coefficients == [1, 0, 0, 0]

    with pytest.raises(qpa.InvalidInputError, match="not a canonical subfield"):
        field16.embed(qpa.Field.extension(3, 2).primitive_element())


@pytest.mark.parametrize(("prime", "degree", "order"), [(2, 2, 4), (2, 3, 8), (3, 2, 9)])
def test_module_homological_enumeration_and_endomorphism_apis(
    prime: int, degree: int, order: int
) -> None:
    field = qpa.Field.extension(prime, degree)
    algebra = qpa.PathAlgebra(
        field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    simple_one, simple_two = algebra.simple_modules()
    projective_one = algebra.indecomposable_projectives()[0]

    assert simple_one.first_syzygy().dimension_vector == [0, 1]
    assert simple_one.projective_dimension(3) == 1

    mixed = qpa.MatrixModule.direct_sum([simple_one, projective_one])
    repeated = qpa.MatrixModule.direct_sum([simple_one, simple_one])
    assert mixed is not None and repeated is not None
    assert [module.dimension_vector for module in mixed.decompose()] == [
        [1, 0],
        [1, 1],
    ]
    multiplicities = repeated.decompose_with_multiplicities()
    assert multiplicities is not None and multiplicities.multiplicities == [2]
    assert len(repeated.simple_submodule_inclusions()) == order + 1

    semisimple = algebra.module_from_integer_matrices(
        [2, 1], [("a", [[0], [0]])]
    )
    submodules = semisimple.all_submodules()
    assert submodules.by_dimension is not None
    assert [len(level) for level in submodules.by_dimension] == [
        1,
        order + 2,
        order + 2,
        1,
    ]
    assert len(submodules.flattened) == 2 * order + 6

    basic_endomorphisms = qpa.MatrixModule.endomorphism_as_quiver_algebra_of_basic(
        [simple_one, projective_one, simple_two]
    )
    assert basic_endomorphisms.quiver.adjacency_matrix == [
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
    ]
    assert basic_endomorphisms.quotient_algebra is not None
    assert basic_endomorphisms.quotient_algebra.dimension() == "5"

    dual_numbers = qpa.QuotientAlgebra.truncated(
        field, qpa.Quiver(["w"], [("w", "w", "x")]), 2
    )
    regular_endomorphisms = (
        dual_numbers.indecomposable_projectives()[0].endomorphism_as_quiver_algebra()
    )
    assert regular_endomorphisms.endomorphism_algebra.dimension == 2
    assert regular_endomorphisms.quiver.adjacency_matrix == [[1]]


def test_extension_field_constructor_rejects_invalid_descriptors() -> None:
    with pytest.raises(qpa.InvalidInputError, match="not a prime"):
        qpa.Field.extension(4, 2)
    with pytest.raises(qpa.InvalidInputError, match="at least two"):
        qpa.Field.extension(2, 1)
    with pytest.raises(qpa.InvalidInputError, match="coordinates"):
        qpa.Field.extension(2, 3).element_from_coefficients([1, 0])


def test_guarded_python_dispatches_accept_extension_objects() -> None:
    field = qpa.Field.extension(2, 2)
    path = qpa.PathAlgebra(
        field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    module = path.module_from_integer_matrices([1, 1], [("a", [[1]])])
    basis = module.basis_elements()

    assert basis[0].act_by_algebra_element(path.generator("a")) == basis[1]
    space = module.element_space(basis)
    assert space.contains(basis[0])
    assert space.basis().coefficients(basis[0]) is not None

    hom = path.hom(module, module)
    end = path.end(module)
    assert hom.dimension == end.dimension == 1
    assert hom.generators[0](basis[0]) == basis[0]
    assert end.multiply(end.one, end.one) == end.one

    category = qpa.ModuleCategory.from_path_algebra(path)
    stalk = module.stalk_complex(0)
    rebuilt = qpa.ModuleComplex.from_differential_list(category, stalk.differentials)
    assert rebuilt.object(0).is_identical(module)
    identity = stalk.finite_chain_map_to(stalk, 0, [module.identity()])
    assert identity.morphism(0).is_isomorphism


def test_qpa_text_persistence_round_trips_extension_fields(tmp_path: Path) -> None:
    field = qpa.Field.extension(2, 2)
    parent = qpa.PathAlgebra(
        field,
        qpa.Quiver(
            ["v1", "v2", "v3"],
            [("v1", "v2", "a"), ("v2", "v3", "b"), ("v1", "v3", "c")],
        ),
    )
    a, b, c = (parent.generator(name) for name in ("a", "b", "c"))
    quotient = parent.quotient([a * b + c.scale_left(field.primitive_element())])
    algebra = qpa.QuiverAlgebra.from_quotient_algebra(quotient)
    target = tmp_path / "gf4-algebra.txt"

    text = algebra.qpa_text()
    assert "Field:\nGF(4)\n" in text and "z^" in text
    assert algebra.save_algebra(str(target), "delete")
    loaded = qpa.read_algebra(str(target))
    assert loaded.field == field
    assert loaded.quotient_algebra is not None
    assert loaded.quotient_algebra.dimension() == quotient.dimension()
