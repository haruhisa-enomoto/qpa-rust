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

PrintModule := function(module)
    local first, matrix;
    Print("{\"dimension\":"); PrintIntegerList(DimensionVector(module));
    Print(",\"matrices\":["); first := true;
    for matrix in MatricesOfPathAlgebraModule(module) do
        if not first then Print(","); fi;
        PrintIntegerMatrix(matrix); first := false;
    od;
    Print("]}");
end;

PrintModules := function(modules)
    local first, module;
    Print("["); first := true;
    for module in modules do
        if not first then Print(","); fi;
        PrintModule(module); first := false;
    od;
    Print("]");
end;

PrintCase := function(field)
    local quiver, algebra, simple_one, simple_two, projective,
          projective_plus_simple;
    quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
    algebra := PathAlgebra(field, quiver);
    simple_one := RightModuleOverPathAlgebra(algebra, [1, 0], []);
    simple_two := RightModuleOverPathAlgebra(algebra, [0, 1], []);
    projective := RightModuleOverPathAlgebra(algebra, [1, 1],
        [["a", [[1]] * One(field)]]);
    projective_plus_simple :=
        DirectSumOfQPAModules([projective, simple_one]);

    Print("{\"ordered_removal\":");
    PrintModules(MakeBasicListOfModules(
        [projective, simple_one, simple_two, projective,
         projective_plus_simple]));
    Print(",\"early_break\":");
    PrintModules(MakeBasicListOfModules(
        [projective, projective, simple_one, simple_one]));
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"modulehom_make_basic\",\"oracle\":{");
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
