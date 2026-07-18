if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

BuildHom := function(ring, source, range)
    return Hom(ring, source, range);
end;;
BuildEnd := function(ring, module)
    return End(ring, module);
end;;
GetBasisVectors := function(space)
    return BasisVectors(Basis(space));
end;;

algebra := PathAlgebra(Rationals,
    Quiver(["x", "y"], [["x", "y", "a"]]));;
module := RightModuleOverPathAlgebra(algebra, [1, 1],
    [["a", [[1]]]]);;
homspace := Hom(algebra, module, module);;
zero := ZeroModule(algebra);;
zero_submodule := SubAlgebraModule(module, []);;

arrowless := PathAlgebra(Rationals, Quiver(["s"], []));;
arrowless_module := RightModuleOverPathAlgebra(arrowless, [2], []);;

other_algebra := PathAlgebra(Rationals,
    Quiver(["x", "y"], [["x", "y", "a"]]));;
other_module := RightModuleOverPathAlgebra(other_algebra, [1, 1],
    [["a", [[1]]]]);;

loop_parent := PathAlgebra(Rationals,
    Quiver(["w"], [["w", "w", "t"]]));;
dual_numbers := loop_parent / [loop_parent.t^2];;
dual_module := RightModuleOverPathAlgebra(dual_numbers, [2],
    [["t", [[0, 1], [0, 0]]]]);;

valid_hom := CALL_WITH_CATCH(BuildHom,
    [algebra, module, module]);;
valid_end := CALL_WITH_CATCH(BuildEnd,
    [algebra, module]);;
basis_vectors := CALL_WITH_CATCH(GetBasisVectors, [homspace]);;
zero_source := CALL_WITH_CATCH(BuildHom,
    [algebra, zero, module]);;
zero_range := CALL_WITH_CATCH(BuildHom,
    [algebra, module, zero]);;
zero_end := CALL_WITH_CATCH(BuildEnd, [algebra, zero]);;
zero_sub_hom := CALL_WITH_CATCH(BuildHom,
    [algebra, zero_submodule, zero_submodule]);;
zero_sub_end := CALL_WITH_CATCH(BuildEnd,
    [algebra, zero_submodule]);;
arrowless_hom := CALL_WITH_CATCH(BuildHom,
    [arrowless, arrowless_module, arrowless_module]);;
arrowless_end := CALL_WITH_CATCH(BuildEnd,
    [arrowless, arrowless_module]);;
wrong_algebra := CALL_WITH_CATCH(BuildHom,
    [algebra, module, other_module]);;
quotient_hom := CALL_WITH_CATCH(BuildHom,
    [dual_numbers, dual_module, dual_module]);;
quotient_end := CALL_WITH_CATCH(BuildEnd,
    [dual_numbers, dual_module]);;

Print("{\"schema_version\":1,\"kind\":\"module_generic_hom_errors\",\"success\":{");
Print("\"valid_hom\":", valid_hom[1], ",");
Print("\"valid_end\":", valid_end[1], ",");
Print("\"basis_vectors\":", basis_vectors[1], ",");
Print("\"zero_source\":", zero_source[1], ",");
Print("\"zero_range\":", zero_range[1], ",");
Print("\"zero_end\":", zero_end[1], ",");
Print("\"zero_sub_hom\":", zero_sub_hom[1], ",");
Print("\"zero_sub_end\":", zero_sub_end[1], ",");
Print("\"arrowless_hom\":", arrowless_hom[1], ",");
Print("\"arrowless_end\":", arrowless_end[1], ",");
Print("\"wrong_algebra\":", wrong_algebra[1], ",");
Print("\"quotient_hom\":", quotient_hom[1], ",");
Print("\"quotient_end\":", quotient_end[1], "}}\n");

QUIT_GAP(0);
