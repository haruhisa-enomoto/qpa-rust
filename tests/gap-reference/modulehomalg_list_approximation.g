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
    local a2, simples, projectives, zero,
          loop_parent, dual, dual_simple, dual_regular;
    a2 := PathAlgebra(field,
        Quiver(["v1", "v2"], [["v1", "v2", "a"]]));
    simples := SimpleModules(a2);
    projectives := IndecProjectiveModules(a2);
    zero := ZeroModule(a2);

    loop_parent := PathAlgebra(field,
        Quiver(["w"], [["w", "w", "x"]]));
    dual := loop_parent / [loop_parent.x^2];
    dual_simple := SimpleModules(dual)[1];
    dual_regular := IndecProjectiveModules(dual)[1];

    Print("{\"empty_to_simple_one\":");
    PrintMap(RightApproximationByAddM([], simples[1]), row_printer);
    Print(",\n\"p1_p2_to_simple_one\":");
    PrintMap(RightApproximationByAddM(projectives, simples[1]), row_printer);
    Print(",\n\"s1_p1_to_simple_one\":");
    PrintMap(RightApproximationByAddM([simples[1], projectives[1]], simples[1]), row_printer);
    Print(",\n\"p1_s1_to_simple_one\":");
    PrintMap(RightApproximationByAddM([projectives[1], simples[1]], simples[1]), row_printer);
    Print(",\n\"repeated_s1_to_simple_one\":");
    PrintMap(RightApproximationByAddM([simples[1], simples[1]], simples[1]), row_printer);
    Print(",\n\"modules_to_zero\":");
    PrintMap(RightApproximationByAddM(projectives, zero), row_printer);
    Print(",\n\"dual_simple_regular_to_simple\":");
    PrintMap(RightApproximationByAddM([dual_simple, dual_regular], dual_simple), row_printer);
    Print(",\n\"dual_regular_simple_to_simple\":");
    PrintMap(RightApproximationByAddM([dual_regular, dual_simple], dual_simple), row_printer);
    Print(",\n\"dual_regular_to_simple\":");
    PrintMap(RightApproximationByAddM([dual_regular], dual_simple), row_printer);
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"modulehomalg_list_approximation\",\"oracle\":{");
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
