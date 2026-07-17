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

PrintPrimeRow := function(row)
    PrintIntegerList(List(row, IntFFE));
end;

PrintMatrix := function(matrix)
    local first, row;
    if not IsList(matrix) then
        Print("[]");
        return;
    fi;
    Print("["); first := true;
    for row in matrix do
        if not first then Print(","); fi;
        PrintPrimeRow(row); first := false;
    od;
    Print("]");
end;

PrintModule := function(module)
    local first, matrix;
    Print("{\"dimension\":"); PrintIntegerList(DimensionVector(module));
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

PrintCase := function(field)
    local a2, simples, plus, at_most, indec, loop_parent, dual,
          dual_simples, dual_length_two, dual_length_three;

    a2 := PathAlgebra(field,
        Quiver(["v1", "v2"], [["v1", "v2", "a"]]));
    simples := SimpleModules(a2);
    plus := AllModulesOfLengthPlusOne(simples);
    at_most := AllModulesOfLengthAtMost(a2, 2);
    indec := AllIndecModulesOfLengthAtMost(a2, 2);

    loop_parent := PathAlgebra(field,
        Quiver(["w"], [["w", "w", "x"]]));
    dual := loop_parent / [loop_parent.x^2];
    dual_simples := SimpleModules(dual);
    dual_length_two := AllModulesOfLengthPlusOne(dual_simples);
    dual_length_three := AllModulesOfLengthPlusOne(dual_length_two);

    Print("{\"a2_plus\":"); PrintModuleList(plus);
    Print(",\"a2_at_most_two\":"); PrintModuleList(at_most);
    Print(",\"a2_indec_two\":"); PrintModuleList(indec);
    Print(",\"a2_negative\":");
    PrintModuleList(AllModulesOfLengthAtMost(a2, -1));
    Print(",\"a2_zero\":");
    PrintModuleList(AllModulesOfLengthAtMost(a2, 0));
    Print(",\"dual_length_two\":"); PrintModuleList(dual_length_two);
    Print(",\"dual_length_three\":"); PrintModuleList(dual_length_three);
    Print("}");
end;

PrintRationalIndecomposable := function()
    local a2, zero, projective;
    a2 := PathAlgebra(Rationals,
        Quiver(["v1", "v2"], [["v1", "v2", "a"]]));
    zero := ZeroModule(a2);
    projective := RightModuleOverPathAlgebra(a2, [1, 1],
        [["a", [[1]]]]);
    Print("{\"zero\":", LowercaseString(String(IsIndecomposableModule(zero))));
    Print(",\"projective\":",
        LowercaseString(String(IsIndecomposableModule(projective))), "}");
end;

Print("{\"schema_version\":1,\"kind\":\"modulehom_length\",\"oracle\":{");
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
Print(",\"rational_indecomposable\":"); PrintRationalIndecomposable();
Print("}\n");

QUIT_GAP(0);
