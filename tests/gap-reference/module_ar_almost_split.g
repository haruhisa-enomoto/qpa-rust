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

PrintSequence := function(sequence, row_printer)
    if sequence = fail then
        Print("null");
        return;
    fi;
    Print("{\"monomorphism\":");
    PrintMap(sequence[1], row_printer);
    Print(",\"epimorphism\":");
    PrintMap(sequence[2], row_printer);
    Print("}");
end;

PrintCase := function(field, row_printer)
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
    PrintSequence(AlmostSplitSequence(simple_one), row_printer);
    Print(",\"ending_simple_two\":");
    PrintSequence(AlmostSplitSequence(simple_two), row_printer);
    Print(",\"ending_projective\":");
    PrintSequence(AlmostSplitSequence(projective), row_printer);
    Print(",\"ending_dual_simple\":");
    PrintSequence(AlmostSplitSequence(dual_simple), row_printer);
    Print(",\"starting_simple_one\":");
    PrintSequence(AlmostSplitSequence(simple_one, "l"), row_printer);
    Print(",\"starting_simple_two\":");
    PrintSequence(AlmostSplitSequence(simple_two, "l"), row_printer);
    Print(",\"starting_dual_simple\":");
    PrintSequence(AlmostSplitSequence(dual_simple, "l"), row_printer);
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"module_ar_almost_split\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"binary\":"); PrintCase(GF(2), PrintPrimeRow);
Print(",\"quinary\":"); PrintCase(GF(5), PrintPrimeRow);
Print(",\"rational\":"); PrintCase(Rationals, PrintRationalRow);
Print("}\n");

QUIT_GAP(0);
