if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

BuildSubmodule := function(module, generators)
    return SubAlgebraModule(module, generators);
end;;
BuildCanonicalBasis := function(module)
    return CanonicalBasis(module);
end;;
ConvertSubmodule := function(module)
    return SubmoduleAsModule(module);
end;;
BuildNaturalMap := function(module, submodule)
    return NaturalHomomorphismBySubAlgebraModule(module, submodule);
end;;

algebra := PathAlgebra(Rationals,
    Quiver(["v1", "v2"], [["v1", "v2", "a"]]));;
module := RightModuleOverPathAlgebra(algebra, [3, 3],
    [["a", [[1, 0, 0], [0, 1, 0], [0, 0, 1]]]]);;
basis := BasisVectors(Basis(module));;
generator := basis[1] + basis[5];;

empty := CALL_WITH_CATCH(BuildSubmodule, [module, []]);;
zero_submodule := SubAlgebraModule(module, [Zero(module)]);;
zero_canonical := CALL_WITH_CATCH(BuildCanonicalBasis, [zero_submodule]);;
zero_conversion := CALL_WITH_CATCH(ConvertSubmodule, [zero_submodule]);;
zero_natural := CALL_WITH_CATCH(BuildNaturalMap, [module, zero_submodule]);;
dependent_submodule := SubAlgebraModule(module,
    [generator, basis[1]]);;
dependent_canonical := CALL_WITH_CATCH(BuildCanonicalBasis,
    [dependent_submodule]);;

small := RightModuleOverPathAlgebra(algebra, [1, 1],
    [["a", [[1]]]]);;
small_basis := BasisVectors(Basis(small));;
target_only := SubAlgebraModule(small, [small_basis[2]]);;
zero_space_conversion := CALL_WITH_CATCH(ConvertSubmodule, [target_only]);;

medium := RightModuleOverPathAlgebra(algebra, [2, 2],
    [["a", [[1, 0], [0, 1]]]]);;
medium_basis := BasisVectors(Basis(medium));;
zero_space_submodule := SubAlgebraModule(medium,
    [medium_basis[1] + medium_basis[4]]);;
zero_space_natural := CALL_WITH_CATCH(BuildNaturalMap,
    [medium, zero_space_submodule]);;
foreign_parent := RightModuleOverPathAlgebra(algebra, [2, 2],
    [["a", [[1, 0], [0, 1]]]]);;
wrong_parent := CALL_WITH_CATCH(BuildNaturalMap,
    [foreign_parent, zero_space_submodule]);;

loop_parent := PathAlgebra(Rationals,
    Quiver(["w"], [["w", "w", "x"]]));;
dual := loop_parent / [loop_parent.x^2];;
dual_module := RightModuleOverPathAlgebra(dual, [2],
    [["x", [[0, 1], [0, 0]]]]);;
quotient_submodule := CALL_WITH_CATCH(BuildSubmodule,
    [dual_module, [BasisVectors(Basis(dual_module))[2]]]);;

Print("{\"schema_version\":1,\"kind\":\"module_submodule_errors\",\"success\":{");
Print("\"empty\":", empty[1], ",");
Print("\"zero_canonical\":", zero_canonical[1], ",");
Print("\"zero_conversion\":", zero_conversion[1], ",");
Print("\"zero_natural\":", zero_natural[1], ",");
Print("\"dependent_canonical\":", dependent_canonical[1], ",");
Print("\"zero_space_conversion\":", zero_space_conversion[1], ",");
Print("\"zero_space_natural\":", zero_space_natural[1], ",");
Print("\"wrong_parent\":", wrong_parent[1], ",");
Print("\"quotient_submodule\":", quotient_submodule[1], "}}\n");

QUIT_GAP(0);
