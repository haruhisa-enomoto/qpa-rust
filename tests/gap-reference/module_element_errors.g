if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

BuildPathModuleElem := function(family, coordinates)
    return PathModuleElem(family, coordinates);
end;;
AddRaw := function(left, right)
    return left + right;
end;;
ActRaw := function(element, path)
    return element^path;
end;;

a2 := PathAlgebra(Rationals,
    Quiver(["v1", "v2"], [["v1", "v2", "a"]]));;
module := RightModuleOverPathAlgebra(a2, [2, 1],
    [["a", [[1], [2]]]]);;
raw := ExtRepOfObj(BasisVectors(Basis(module))[1]);;
family := FamilyObj(raw);;

foreign_a2 := PathAlgebra(Rationals,
    Quiver(["w1", "w2"], [["w1", "w2", "b"]]));;
foreign_module := RightModuleOverPathAlgebra(foreign_a2, [2, 1],
    [["b", [[1], [2]]]]);;
foreign_raw := ExtRepOfObj(BasisVectors(Basis(foreign_module))[1]);;

wrong_vertices := CALL_WITH_CATCH(BuildPathModuleElem,
    [family, [[1, 0]]]);;
wrong_dimension := CALL_WITH_CATCH(BuildPathModuleElem,
    [family, [[1], [0]]]);;
different_families := CALL_WITH_CATCH(AddRaw, [raw, foreign_raw]);;
foreign_action := CALL_WITH_CATCH(ActRaw, [raw, foreign_a2.b]);;

Print("{\"schema_version\":1,\"kind\":\"module_element_errors\",\"success\":{");
Print("\"wrong_vertices\":", wrong_vertices[1], ",");
Print("\"wrong_dimension\":", wrong_dimension[1], ",");
Print("\"different_families\":", different_families[1], ",");
Print("\"foreign_action\":", foreign_action[1], "}}\n");

QUIT_GAP(0);
