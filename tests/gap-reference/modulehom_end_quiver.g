if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintIntegerList := function(values)
    Print("[", JoinStringsWithSeparator(List(values, String), ","), "]");
end;

PrintIntegerMatrix := function(matrix)
    local first, row;
    Print("["); first := true;
    for row in matrix do
        if not first then Print(","); fi;
        PrintIntegerList(row); first := false;
    od;
    Print("]");
end;

PrintResult := function(result)
    Print("{\"endomorphism_dimension\":", Dimension(result[1]));
    Print(",\"adjacency\":"); PrintIntegerMatrix(result[2]);
    Print(",\"quiver_algebra_dimension\":", Dimension(result[3]), "}");
end;

PrintCase := function(field)
    local loop_quiver, loop_parent, dual_numbers, regular, regular_result;
    loop_quiver := Quiver(["w"], [["w", "w", "x"]]);
    loop_parent := PathAlgebra(field, loop_quiver);
    dual_numbers := loop_parent / [loop_parent.x^2];
    regular := RightModuleOverPathAlgebra(dual_numbers, [2],
        [["x", [[0, 1], [0, 0]] * One(field)]]);
    regular_result := EndOfModuleAsQuiverAlgebra(regular);

    Print("{\"regular_dual_numbers\":"); PrintResult(regular_result);
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"modulehom_end_quiver\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"binary\":"); PrintCase(GF(2));
Print(",\"quinary\":"); PrintCase(GF(5));
Print("}\n");

QUIT_GAP(0);
