"""A noncommutative Gröbner basis and normal-form reduction."""

import qpa


field = qpa.Field.prime(5)
quiver = qpa.Quiver(["v"], [("v", "v", "x"), ("v", "v", "y")])
path_algebra = qpa.PathAlgebra(field, quiver)
x = path_algebra.generator("x")
y = path_algebra.generator("y")

relations = [x.pow(2), x * y - y * x, y.pow(2)]
basis = path_algebra.groebner_basis(relations)
remainder = path_algebra.remainder_of_division(y * x, basis.relations)

assert len(basis.relations) == 3
assert remainder == x * y

print("basis:")
for relation in basis:
    print(f"  {relation.print_obj()}")
print(f"remainder of y*x: {remainder.print_obj()}")
