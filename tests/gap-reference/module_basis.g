if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

SizeScreen([100000, 100000]);

a2 := PathAlgebra(Rationals,
    Quiver(["v1", "v2"], [["v1", "v2", "a"]]));;
module := RightModuleOverPathAlgebra(a2, [2, 1],
    [["a", [[1], [2]]]]);;
outer_standard := GeneratorsOfLeftModule(module);;
standard := List(outer_standard, ExtRepOfObj);;
entered := [ExtRepOfObj(Zero(module)), standard[2], standard[1], standard[3]];;
space := LeftModuleByGenerators(Rationals, entered);;
space_generators := GeneratorsOfLeftModule(space);;
v := 2 * outer_standard[1] - 3 * outer_standard[2] + 5 * outer_standard[3];;
raw_v := ExtRepOfObj(v);;

new_basis := NewBasis(space, space_generators);;
unchecked := BasisNC(space, space_generators);;
checked := Basis(space, space_generators);;
domain_basis := BasisOfDomain(space);;
default_basis := Basis(space);;
canonical := CanonicalBasis(module);;

Print("NEW_VECTORS=", BasisVectors(new_basis), "\n");
Print("UNCHECKED_VECTORS=", BasisVectors(unchecked), "\n");
Print("CHECKED_VECTORS=", BasisVectors(checked), "\n");
Print("CHECKED_COEFFICIENTS=", Coefficients(checked, raw_v), "\n");
Print("DOMAIN_VECTORS=", BasisVectors(domain_basis), "\n");
Print("DEFAULT_VECTORS=", BasisVectors(default_basis), "\n");
Print("CANONICAL_VECTORS=", BasisVectors(canonical), "\n");
Print("CANONICAL_PROPERTY=", IsCanonicalBasis(canonical), "\n");

prime_a2 := PathAlgebra(GF(2),
    Quiver(["p1", "p2"], [["p1", "p2", "b"]]));;
prime_module := RightModuleOverPathAlgebra(prime_a2, [2, 1],
    [["b", [[One(GF(2))], [One(GF(2))]]]]);;
prime_outer_standard := GeneratorsOfLeftModule(prime_module);;
prime_standard := List(prime_outer_standard, ExtRepOfObj);;
prime_space := LeftModuleByGenerators(GF(2),
    [prime_standard[3], prime_standard[1], prime_standard[2]]);;
prime_checked := Basis(prime_space, GeneratorsOfLeftModule(prime_space));;
Print("GF2_COEFFICIENTS=",
    Coefficients(prime_checked,
        ExtRepOfObj(prime_outer_standard[1] + prime_outer_standard[3])), "\n");

QUIT_GAP(0);
