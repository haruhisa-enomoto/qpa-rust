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
    local first, value;
    Print("["); first := true;
    for value in values do
        if not first then Print(","); fi;
        Print(Int(value)); first := false;
    od;
    Print("]");
end;

PrintMatrixList := function(matrices)
    local first, matrix, firstrow, row;
    Print("["); first := true;
    for matrix in matrices do
        if not first then Print(","); fi;
        Print("["); firstrow := true;
        for row in matrix do
            if not firstrow then Print(","); fi;
            PrintIntegerList(row); firstrow := false;
        od;
        Print("]"); first := false;
    od;
    Print("]");
end;

PrintModule := function(module)
    Print("{\"dimension\":"); PrintIntegerList(DimensionVector(module));
    Print(",\"matrices\":");
    PrintMatrixList(MatricesOfPathAlgebraModule(module));
    Print("}");
end;

PrintCase := function(field)
    local a2, a2module, vectoralgebra, vector, parent, dual, x;
    a2 := PathAlgebra(field,
        Quiver(["v1", "v2"], [["v1", "v2", "a"]]));
    a2module := RightModuleOverPathAlgebra(a2, [2,1],
        [["a", [[1],[0]] * One(field)]]);
    vectoralgebra := PathAlgebra(field, Quiver(["w"], []));
    vector := RightModuleOverPathAlgebra(vectoralgebra, [2], []);
    parent := PathAlgebra(field, Quiver(["u"], [["u", "u", "x"]]));
    dual := parent / [parent.x^2];
    x := [[0,1],[0,0]] * One(field);
    dual := RightModuleOverPathAlgebra(dual, [2], [["x", x]]);

    Print("{\"a2_partial\":"); PrintModule(PartialIyamaGenerator(a2module));
    Print(",\"a2_iyama\":"); PrintModule(IyamaGenerator(a2module));
    Print(",\"vector_partial\":"); PrintModule(PartialIyamaGenerator(vector));
    Print(",\"vector_iyama\":"); PrintModule(IyamaGenerator(vector));
    Print(",\"dual_partial\":"); PrintModule(PartialIyamaGenerator(dual));
    Print(",\"dual_iyama\":"); PrintModule(IyamaGenerator(dual));
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"modulehomalg_iyama\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"rational\":"); PrintCase(Rationals);
Print(",\"binary\":"); PrintCase(GF(2));
Print("}\n");

QUIT_GAP(0);
