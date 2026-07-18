if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

BuildChecked := function(algebra, actions)
    return RightModuleOverPathAlgebra(algebra, actions);
end;;
BuildUnchecked := function(algebra, actions)
    return RightModuleOverPathAlgebraNC(algebra, actions);
end;;

chain := PathAlgebra(Rationals,
    Quiver(["v1", "v2", "v3"],
        [["v1", "v2", "a"], ["v2", "v3", "b"]]));;
isolated := PathAlgebra(Rationals,
    Quiver(["i1", "i2", "i3"], [["i1", "i2", "c"]]));;
arrowless := PathAlgebra(Rationals, Quiver(["e"], []));;

too_few := CALL_WITH_CATCH(BuildChecked, [chain, [[[1]]]]);;
mismatch := CALL_WITH_CATCH(BuildChecked,
    [chain, [[[1, 2]], [[1], [2], [3]]]]);;
invalid_zero_pair := CALL_WITH_CATCH(BuildChecked,
    [chain, [["a", [2, 2]], ["b", [[1], [2]]]]]);;
wrong_field := CALL_WITH_CATCH(BuildChecked,
    [chain, [[[One(GF(2))]], [[One(GF(2))]]]]);;
unchecked_wrong_field := CALL_WITH_CATCH(BuildUnchecked,
    [chain, [[[One(GF(2))]], [[One(GF(2))]]]]);;
duplicate := CALL_WITH_CATCH(BuildChecked,
    [chain, [["a", [[1]]], ["a", [[1]]]]]);;
unchecked_isolated := CALL_WITH_CATCH(BuildUnchecked,
    [isolated, [[[1]]]]);;
empty_arrowless := CALL_WITH_CATCH(BuildChecked, [arrowless, []]);;

Print("{\"schema_version\":1,\"kind\":\"module_inferred_constructor_errors\",\"success\":{");
Print("\"too_few\":", too_few[1], ",");
Print("\"mismatch\":", mismatch[1], ",");
Print("\"invalid_zero_pair\":", invalid_zero_pair[1], ",");
Print("\"wrong_field\":", wrong_field[1], ",");
Print("\"unchecked_wrong_field\":", unchecked_wrong_field[1], ",");
Print("\"duplicate\":", duplicate[1], ",");
Print("\"unchecked_isolated\":", unchecked_isolated[1], ",");
Print("\"empty_arrowless\":", empty_arrowless[1], "}}\n");

QUIT_GAP(0);
