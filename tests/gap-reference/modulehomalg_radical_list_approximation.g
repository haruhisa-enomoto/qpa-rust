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

PrintCase := function(field, row_printer)
    local a2, simples, projectives,
          loop_parent, dual, dual_simple, dual_regular;
    a2 := PathAlgebra(field,
        Quiver(["v1", "v2"], [["v1", "v2", "a"]]));
    simples := SimpleModules(a2);
    projectives := IndecProjectiveModules(a2);

    loop_parent := PathAlgebra(field,
        Quiver(["w"], [["w", "w", "x"]]));
    dual := loop_parent / [loop_parent.x^2];
    dual_simple := SimpleModules(dual)[1];
    dual_regular := IndecProjectiveModules(dual)[1];

    Print("{\"projectives_at_one\":");
    PrintMap(RadicalRightApproximationByAddM(projectives, 1), row_printer);
    Print(",\n\"s1_p1_at_one\":");
    PrintMap(RadicalRightApproximationByAddM([simples[1], projectives[1]], 1), row_printer);
    Print(",\n\"p1_s1_at_two\":");
    PrintMap(RadicalRightApproximationByAddM([projectives[1], simples[1]], 2), row_printer);
    Print(",\n\"single_s1\":");
    PrintMap(RadicalRightApproximationByAddM([simples[1]], 1), row_printer);
    Print(",\n\"repeated_s1\":");
    PrintMap(RadicalRightApproximationByAddM([simples[1], simples[1]], 1), row_printer);
    Print(",\n\"single_dual_simple\":");
    PrintMap(RadicalRightApproximationByAddM([dual_simple], 1), row_printer);
    Print(",\n\"single_dual_regular\":");
    PrintMap(RadicalRightApproximationByAddM([dual_regular], 1), row_printer);
    Print(",\n\"dual_regular_simple_at_one\":");
    PrintMap(RadicalRightApproximationByAddM([dual_regular, dual_simple], 1), row_printer);
    Print(",\n\"dual_regular_simple_at_two\":");
    PrintMap(RadicalRightApproximationByAddM([dual_regular, dual_simple], 2), row_printer);
    Print(",\n\"repeated_dual_regular\":");
    PrintMap(RadicalRightApproximationByAddM([dual_regular, dual_regular], 1), row_printer);
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"modulehomalg_radical_list_approximation\",\"oracle\":{");
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
