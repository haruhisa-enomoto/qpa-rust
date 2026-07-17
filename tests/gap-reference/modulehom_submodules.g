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
        PrintIntegerList(List(row, Int)); first := false;
    od;
    Print("]");
end;

PrintMap := function(map)
    local first, matrix;
    Print("{\"source_dimension\":");
    PrintIntegerList(DimensionVector(Source(map)));
    Print(",\"matrices\":["); first := true;
    for matrix in MatricesOfPathAlgebraMatModuleHomomorphism(map) do
        if not first then Print(","); fi;
        PrintIntegerMatrix(matrix); first := false;
    od;
    Print("]}");
end;

PrintMapList := function(maps)
    local first, map;
    Print("["); first := true;
    for map in maps do
        if not first then Print(","); fi;
        PrintMap(map); first := false;
    od;
    Print("]");
end;

PrintMapLevels := function(levels)
    local first, level;
    Print("["); first := true;
    for level in levels do
        if not first then Print(","); fi;
        PrintMapList(level); first := false;
    od;
    Print("]");
end;

PrintCase := function(field)
    local quiver, algebra, zero, simple, semisimple, projective;
    quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
    algebra := PathAlgebra(field, quiver);
    zero := RightModuleOverPathAlgebra(algebra, [0, 0], []);
    simple := RightModuleOverPathAlgebra(algebra, [1, 0], []);
    semisimple := RightModuleOverPathAlgebra(algebra, [2, 1],
        [["a", [[0], [0]] * One(field)]]);
    projective := RightModuleOverPathAlgebra(algebra, [1, 1],
        [["a", [[1]] * One(field)]]);

    Print("{\"all_simple\":{");
    Print("\"zero\":"); PrintMapList(AllSimpleSubmodulesOfModule(zero));
    Print(",\"simple\":"); PrintMapList(AllSimpleSubmodulesOfModule(simple));
    Print(",\"semisimple\":");
    PrintMapList(AllSimpleSubmodulesOfModule(semisimple));
    Print(",\"projective\":");
    PrintMapList(AllSimpleSubmodulesOfModule(projective));
    Print("},\"all_submodules\":{");
    Print("\"zero\":{\"shape\":\"flat\",\"maps\":");
    PrintMapList(AllSubmodulesOfModule(zero)); Print("}");
    Print(",\"simple\":{\"shape\":\"flat\",\"maps\":");
    PrintMapList(AllSubmodulesOfModule(simple)); Print("}");
    Print(",\"semisimple\":{\"shape\":\"levels\",\"levels\":");
    PrintMapLevels(AllSubmodulesOfModule(semisimple)); Print("}");
    Print(",\"projective\":{\"shape\":\"levels\",\"levels\":");
    PrintMapLevels(AllSubmodulesOfModule(projective)); Print("}");
    Print("}}");
end;

Print("{\"schema_version\":1,\"kind\":\"modulehom_submodules\",\"oracle\":{");
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
