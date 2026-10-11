"""Ordinary quotient/Homological APIs need no supplied relation partition."""

import pytest
import qpa


@pytest.mark.parametrize('field', [qpa.Field.rational(), qpa.Field.prime(2),
                                  qpa.Field.prime(5), qpa.Field.extension(2, 2)])
def test_ordinary_quotient_on_renamed_interleaved_component_relations(field):
    arrows = [('south', 'south', 'delta'), ('north', 'south', 'bridge'),
              ('north', 'north', 'alpha'), ('south', 'south', 'gamma'),
              ('north', 'north', 'beta')]
    parent = qpa.PathAlgebra(field, qpa.Quiver(['north', 'south'], arrows))
    p, q, h, r, s = [parent.generator(name)
                     for name in ['alpha', 'beta', 'bridge', 'gamma', 'delta']]
    scalar = field.primitive_element() if (field.degree or 1) > 1 else field.element(1)
    # Deliberately interleave crossing and internal relations. Each diagonal
    # block has dimension 9 and the connecting bimodule also has dimension 9.
    relations = [p*h - (h*r).scale_left(scalar), s**3, p*q - q*p,
                 q*h - h*s, p**3, r*s - s*r, q**3, r**3]
    algebra = parent.quotient(relations)
    assert int(algebra.dimension()) == 27
    assert len(algebra.canonical_basis()) == 27
    assert algebra.is_admissible
    assert all(algebra.groebner_basis.remainder(relation).is_zero for relation in relations)
    first, second = algebra.simple_modules()
    assert first.ext1_dimension(second) == 1
    assert second.ext1_dimension(first) == 0
