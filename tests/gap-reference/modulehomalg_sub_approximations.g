if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

SizeScreen([100000, 100000]);

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

PrintRationalRow := function(row)
    local first, entry;
    Print("["); first := true;
    for entry in row do
        if not first then Print(","); fi;
        Print("[", NumeratorRat(entry), ",", DenominatorRat(entry), "]");
        first := false;
    od;
    Print("]");
end;

PrintMatrix := function(matrix, row_printer)
    local first, row;
    Print("["); first := true;
    for row in matrix do
        if not first then Print(","); fi;
        row_printer(row); first := false;
    od;
    Print("]");
end;

PrintMap := function(map, row_printer)
    local first, matrix;
    Print("{\"source_dimension\":");
    PrintIntegerList(DimensionVector(Source(map)));
    Print(",\"range_dimension\":");
    PrintIntegerList(DimensionVector(Range(map)));
    Print(",\"matrices\":["); first := true;
    for matrix in MatricesOfPathAlgebraMatModuleHomomorphism(map) do
        if not first then Print(","); fi;
        PrintMatrix(matrix, row_printer); first := false;
    od;
    Print("]}");
end;

PrintApproximations := function(source, module, row_printer)
    Print("{\"right\":");
    PrintMap(RightSubMApproximation(module, source), row_printer);
    Print(",\"minimal_right\":");
    PrintMap(MinimalRightSubMApproximation(module, source), row_printer);
    Print(",\"left\":");
    PrintMap(LeftSubMApproximation(source, module), row_printer);
    Print(",\"minimal_left\":");
    PrintMap(MinimalLeftSubMApproximation(source, module), row_printer);
    Print("}");
end;

PrintCase := function(field, row_printer)
    local a2, zero, simple_one, projective,
          loop_parent, dual, dual_simple, dual_regular;
    a2 := PathAlgebra(field,
        Quiver(["v1", "v2"], [["v1", "v2", "a"]]));
    zero := ZeroModule(a2);
    simple_one := RightModuleOverPathAlgebra(a2, [1, 0], []);
    projective := RightModuleOverPathAlgebra(a2, [1, 1],
        [["a", [[1]] * One(field)]]);

    loop_parent := PathAlgebra(field,
        Quiver(["w"], [["w", "w", "x"]]));
    dual := loop_parent / [loop_parent.x^2];
    dual_simple := RightModuleOverPathAlgebra(dual, [1],
        [["x", [[0]] * One(field)]]);
    dual_regular := RightModuleOverPathAlgebra(dual, [2],
        [["x", [[0, 1], [0, 0]] * One(field)]]);

    Print("{\"zero_and_simple\":");
    PrintApproximations(zero, simple_one, row_printer);
    Print(",\"simple_and_zero\":");
    PrintApproximations(simple_one, zero, row_printer);
    Print(",\"projective_and_simple\":");
    PrintApproximations(projective, simple_one, row_printer);
    Print(",\"dual_regular_and_simple\":");
    PrintApproximations(dual_regular, dual_simple, row_printer);
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"modulehomalg_sub_approximations\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\n\"binary\":"); PrintCase(GF(2), PrintPrimeRow);
Print(",\n\"quinary\":"); PrintCase(GF(5), PrintPrimeRow);
Print(",\n\"rational\":"); PrintCase(Rationals, PrintRationalRow);
Print("}\n");

QUIT_GAP(0);
