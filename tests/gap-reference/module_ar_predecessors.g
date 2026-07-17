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

PrintModule := function(module)
    local first, matrix;
    Print("{\"dimension\":");
    PrintIntegerList(DimensionVector(module));
    Print(",\"matrices\":["); first := true;
    for matrix in MatricesOfPathAlgebraModule(module) do
        if not first then Print(","); fi;
        PrintMatrix(matrix); first := false;
    od;
    Print("]}");
end;

PrintModuleList := function(modules)
    local first, module;
    Print("["); first := true;
    for module in modules do
        if not first then Print(","); fi;
        PrintModule(module); first := false;
    od;
    Print("]");
end;

PrintLayers := function(layers)
    local first, layer;
    Print("["); first := true;
    for layer in layers do
        if not first then Print(","); fi;
        PrintModuleList(layer); first := false;
    od;
    Print("]");
end;

PrintScalarLists := function(value)
    local first, entry;
    if IsList(value) then
        Print("["); first := true;
        for entry in value do
            if not first then Print(","); fi;
            PrintScalarLists(entry); first := false;
        od;
        Print("]");
    elif value = false then
        Print("false");
    else
        Print(value);
    fi;
end;

PrintPredecessors := function(value, is_module)
    if is_module then
        Print("{\"kind\":\"module\",\"module\":");
        PrintModule(value);
        Print("}");
    else
        Print("{\"kind\":\"layers\",\"layers\":");
        PrintLayers(value[1]);
        Print(",\"valuations\":");
        PrintScalarLists(value[2]);
        Print("}");
    fi;
end;

PrintCase := function(field)
    local a2, simple_one, loop_parent, dual, dual_simple;
    a2 := PathAlgebra(field,
        Quiver(["v1", "v2"], [["v1", "v2", "a"]]));
    simple_one := RightModuleOverPathAlgebra(a2, [1, 0], []);

    loop_parent := PathAlgebra(field,
        Quiver(["w"], [["w", "w", "x"]]));
    dual := loop_parent / [loop_parent.x^2];
    dual_simple := RightModuleOverPathAlgebra(dual, [1],
        [["x", [[0]] * One(field)]]);

    Print("{\"a2_n1\":"); PrintPredecessors(PredecessorsOfModule(simple_one, 1), true);
    Print(",\"a2_n2\":"); PrintPredecessors(PredecessorsOfModule(simple_one, 2), false);
    Print(",\"a2_n3\":"); PrintPredecessors(PredecessorsOfModule(simple_one, 3), false);
    Print(",\"a2_n4\":"); PrintPredecessors(PredecessorsOfModule(simple_one, 4), false);
    Print(",\"dual_n2\":"); PrintPredecessors(PredecessorsOfModule(dual_simple, 2), false);
    Print(",\"dual_n3\":"); PrintPredecessors(PredecessorsOfModule(dual_simple, 3), false);
    Print(",\"dual_n4\":"); PrintPredecessors(PredecessorsOfModule(dual_simple, 4), false);
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"module_ar_predecessors\",\"oracle\":{");
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
