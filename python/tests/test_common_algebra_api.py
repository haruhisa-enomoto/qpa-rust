"""Common algebra operations retain concrete parents and quotient relations."""

import pytest

import qpa


@pytest.mark.parametrize(
    "field",
    [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5), qpa.Field.extension(2, 2)],
)
@pytest.mark.parametrize("quotient", [False, True])
def test_common_modules_keep_parent_and_vertex_order(field, quotient):
    path = qpa.PathAlgebra(field, qpa.Quiver.indexed(3, [(1, 2, "a"), (2, 3, "b")]))
    concrete = path.quotient([path.generator("a") * path.generator("b")]) if quotient else path
    algebra = concrete.as_quiver_algebra()
    if quotient:
        assert algebra.quotient_algebra.is_identical(concrete)
    else:
        assert algebra.path_algebra.is_identical(concrete)
    assert [m.dimension_vector for m in algebra.simple_modules()] == [
        [1, 0, 0], [0, 1, 0], [0, 0, 1]
    ]
    projectives = [[1, 1, 0], [0, 1, 1], [0, 0, 1]] if quotient else [
        [1, 1, 1], [0, 1, 1], [0, 0, 1]
    ]
    injectives = [[1, 0, 0], [1, 1, 0], [0, 1, 1]] if quotient else [
        [1, 0, 0], [1, 1, 0], [1, 1, 1]
    ]
    for method, expected in [
        ("indecomposable_projectives", projectives),
        ("indecomposable_injectives", injectives),
    ]:
        modules = getattr(algebra, method)()
        assert [m.dimension_vector for m in modules] == expected
        assert all(m.right_acting_algebra.is_identical(algebra) for m in modules)
    zero = algebra.zero_module()
    assert zero.dimension_vector == [0, 0, 0]
    assert zero.right_acting_algebra.is_identical(algebra)
    module = algebra.module_from_integer_matrices([1, 1, 1], [("a", [[1]]), ("b", [[0]])])
    assert module.right_acting_algebra.is_identical(algebra)
    assert module.dimension_vector == [1, 1, 1]
    if quotient:
        with pytest.raises(qpa.RelationNotSatisfiedError, match="do not satisfy relations"):
            algebra.module_from_integer_matrices([1, 1, 1], [("a", [[1]]), ("b", [[1]])])


@pytest.mark.parametrize("quotient", [False, True])
def test_common_generators_preserve_concrete_order_and_relations(quotient):
    path = qpa.PathAlgebra(qpa.Field.prime(2), qpa.Quiver.indexed(3, [(1, 2, "a"), (2, 3, "b")]))
    concrete = path.quotient([path.generator("a") * path.generator("b")]) if quotient else path
    algebra = concrete.as_quiver_algebra()
    names = ["v1", "v2", "v3", "a", "b"]
    assert len(algebra.generators) == len(concrete.generators) == len(names)
    for name, generator in zip(names, algebra.generators):
        assert isinstance(generator, qpa.QuiverAlgebraElement)
        assert generator.algebra.is_identical(algebra)
        assert (generator - algebra.generator(name)).is_zero
        concrete_generator = generator.quotient_element if quotient else generator.path_element
        assert concrete_generator == concrete.generator(name)
    assert algebra.zero().is_zero
    assert (algebra.one() - sum_vertices(algebra)).is_zero
    product = algebra.generator("a") * algebra.generator("b")
    assert product.is_zero == quotient
    with pytest.raises(qpa.InvalidInputError, match="unknown quiver generator"):
        algebra.generator("missing")


def sum_vertices(algebra):
    return algebra.generator("v1") + algebra.generator("v2") + algebra.generator("v3")


@pytest.mark.parametrize("left_kind", ["path", "quotient", "wrapper"])
@pytest.mark.parametrize("right_kind", ["path", "quotient", "wrapper"])
def test_tensor_products_accept_algebra_variants(left_kind, right_kind):
    field = qpa.Field.prime(2)
    path = qpa.PathAlgebra(field, qpa.Quiver.indexed(2, [(1, 2, "a")]))
    left = path.quotient([]) if left_kind == "quotient" else path
    left_wrapper = left.as_quiver_algebra()
    if left_kind == "wrapper":
        left = left_wrapper
    right = path.opposite()
    if right_kind == "quotient":
        right = right.quotient([])
    right_wrapper = right.as_quiver_algebra()
    if right_kind == "wrapper":
        right = right_wrapper
    tensor = left.tensor_product(right)
    assert tensor.quotient.dimension() == "9"
    assert tensor.factor(1).is_identical(left_wrapper)
    assert tensor.factor(2).is_identical(right_wrapper)
    assert len(tensor.algebra.indecomposable_projectives()) == 4
    enveloping = left.enveloping_algebra()
    assert enveloping.is_identical(left_wrapper.enveloping_algebra())
    assert len(enveloping.algebra.simple_modules()) == 4


@pytest.mark.parametrize("kind", ["path", "quotient", "wrapper"])
def test_tensor_errors_are_explicit(kind):
    path = qpa.PathAlgebra(qpa.Field.prime(2), qpa.Quiver.indexed(1, []))
    if kind == "quotient":
        path = path.quotient([])
    elif kind == "wrapper":
        path = path.as_quiver_algebra()
    other = qpa.PathAlgebra(qpa.Field.rational(), qpa.Quiver.indexed(1, []))
    with pytest.raises(TypeError, match="PathAlgebra, QuotientAlgebra, or QuiverAlgebra"):
        path.tensor_product(1)
    with pytest.raises(qpa.InvalidInputError, match="different coefficient fields"):
        path.tensor_product(other)


@pytest.mark.parametrize("kind", ["path", "quotient", "wrapper"])
def test_exact_module_matrices_accept_mixed_rational_entries(kind):
    from fractions import Fraction

    field = qpa.Field.rational()
    path = qpa.PathAlgebra(field, qpa.Quiver.indexed(2, [(1, 2, "a")]))
    algebra = path.quotient([]) if kind == "quotient" else path
    if kind == "wrapper":
        algebra = path.as_quiver_algebra()
    module = algebra.module_from_matrices(
        [1, 3], [("a", [[Fraction(1, 2), field.element(2, 3), 10**30]])]
    )
    assert module.arrow_matrices == [[[field.element(1, 2), field.element(2, 3), field.element(10**30)]]]


@pytest.mark.parametrize("kind", ["path", "quotient", "wrapper"])
def test_exact_module_matrices_preserve_extension_coefficients(kind):
    field = qpa.Field.extension(2, 2)
    path = qpa.PathAlgebra(field, qpa.Quiver.indexed(2, [(1, 2, "a")]))
    algebra = path.quotient([]) if kind == "quotient" else path
    if kind == "wrapper":
        algebra = path.as_quiver_algebra()
    alpha = field.primitive_element()
    module = algebra.module_from_matrices([1, 2], [("a", [[alpha, 1]])])
    assert module.arrow_matrices == [[[alpha, field.element(1)]]]


def test_exact_module_matrices_retain_relations_and_reject_inexact_entries():
    from fractions import Fraction

    field = qpa.Field.prime(2)
    path = qpa.PathAlgebra(field, qpa.Quiver.indexed(1, [(1, 1, "x")]))
    quotient = path.quotient([path.generator("x") * path.generator("x")])
    for algebra in [quotient, quotient.as_quiver_algebra()]:
        with pytest.raises(qpa.RelationNotSatisfiedError, match="do not satisfy relations"):
            algebra.module_from_matrices([1], [("x", [[1]])])
    for algebra in [path, quotient, quotient.as_quiver_algebra()]:
        with pytest.raises(TypeError):
            algebra.module_from_matrices([1], [("x", [[0.5]])])
        with pytest.raises(qpa.InvalidInputError, match="denominator is zero"):
            algebra.module_from_matrices([1], [("x", [[Fraction(1, 2)]])])
        with pytest.raises(qpa.InvalidInputError):
            algebra.module_from_matrices([1], [("x", [[qpa.Field.prime(5).element(1)]])])


@pytest.mark.parametrize("quotient", [False, True])
def test_module_element_act_by_accepts_concrete_and_common_generators(quotient):
    field = qpa.Field.prime(5)
    path = qpa.PathAlgebra(field, qpa.Quiver.indexed(2, [(1, 2, "a")]))
    concrete = path.quotient([]) if quotient else path
    algebra = concrete.as_quiver_algebra()
    module = algebra.module_from_integer_matrices([1, 1], [("a", [[2]])])
    vector = module.element([[field.element(1)], [field.element(0)]])
    expected = module.element([[field.element(0)], [field.element(2)]])
    assert vector.act_by(concrete.generator("a")) == expected
    assert vector.act_by(algebra.generator("a")) == expected
    assert vector.act_by(algebra.one()) == vector
    assert vector.act_by(algebra.zero()).is_zero
    same_looking_path = qpa.PathAlgebra(field, qpa.Quiver.indexed(2, [(1, 2, "a")]))
    unrelated = same_looking_path.quotient([]) if quotient else same_looking_path
    for generator in [unrelated.generator("a"), unrelated.as_quiver_algebra().generator("a")]:
        with pytest.raises(qpa.InvalidInputError):
            vector.act_by(generator)
