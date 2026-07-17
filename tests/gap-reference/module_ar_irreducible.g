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
    if Length(values) = 0 then
        Print("[]");
    else
        Print("[", JoinStringsWithSeparator(List(values, String), ","), "]");
    fi;
end;

PrintMatrix := function(matrix)
    local first, row;
    Print("["); first := true;
    for row in matrix do
        if not first then Print(","); fi;
        PrintIntegerList(List(row, IntFFE)); first := false;
    od;
    Print("]");
end;

PrintMap := function(map)
    local first, matrix;
    Print("{\"source_dimension\":");
    PrintIntegerList(DimensionVector(Source(map)));
    Print(",\"range_dimension\":");
    PrintIntegerList(DimensionVector(Range(map)));
    Print(",\"matrices\":["); first := true;
    for matrix in MatricesOfPathAlgebraMatModuleHomomorphism(map) do
        if not first then Print(","); fi;
        PrintMatrix(matrix); first := false;
    od;
    Print("]}");
end;

PrintMaps := function(maps)
    local first, map;
    Print("["); first := true;
    for map in maps do
        if not first then Print(","); fi;
        PrintMap(map); first := false;
    od;
    Print("]");
end;

PrintCase := function(field)
    local a2, simple_one, simple_two, projective,
          loop_parent, dual, dual_simple;
    a2 := PathAlgebra(field,
        Quiver(["v1", "v2"], [["v1", "v2", "a"]]));
    simple_one := RightModuleOverPathAlgebra(a2, [1, 0], []);
    simple_two := RightModuleOverPathAlgebra(a2, [0, 1], []);
    projective := RightModuleOverPathAlgebra(a2, [1, 1],
        [["a", [[1]] * One(field)]]);

    loop_parent := PathAlgebra(field,
        Quiver(["w"], [["w", "w", "x"]]));
    dual := loop_parent / [loop_parent.x^2];
    dual_simple := RightModuleOverPathAlgebra(dual, [1],
        [["x", [[0]] * One(field)]]);

    Print("{\"ending_simple_one\":");
    PrintMaps(IrreducibleMorphismsEndingIn(simple_one));
    Print(",\"ending_projective\":");
    PrintMaps(IrreducibleMorphismsEndingIn(projective));
    Print(",\"ending_dual_simple\":");
    PrintMaps(IrreducibleMorphismsEndingIn(dual_simple));
    Print(",\"starting_simple_two\":");
    PrintMaps(IrreducibleMorphismsStartingIn(simple_two));
    Print(",\"starting_dual_simple\":");
    PrintMaps(IrreducibleMorphismsStartingIn(dual_simple));
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"module_ar_irreducible\",\"oracle\":{");
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
