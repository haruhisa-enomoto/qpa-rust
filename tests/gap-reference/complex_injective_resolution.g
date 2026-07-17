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

PrintMatrixList := function(matrices, row_printer)
    local first, matrix;
    Print("["); first := true;
    for matrix in matrices do
        if not first then Print(","); fi;
        PrintMatrix(matrix, row_printer); first := false;
    od;
    Print("]");
end;

PrintMap := function(map, row_printer)
    Print("{\"source_dimension\":");
    PrintIntegerList(DimensionVector(Source(map)));
    Print(",\"range_dimension\":");
    PrintIntegerList(DimensionVector(Range(map)));
    Print(",\"matrices\":");
    PrintMatrixList(MatricesOfPathAlgebraMatModuleHomomorphism(map), row_printer);
    Print("}");
end;

PrintResolution := function(module, degrees, row_printer)
    local resolution, first, degree;
    resolution := InjectiveResolution(module);
    Print("{\"objects\":["); first := true;
    for degree in degrees do
        if not first then Print(","); fi;
        PrintIntegerList(DimensionVector(ObjectOfComplex(resolution, degree)));
        first := false;
    od;
    Print("],\"differentials\":["); first := true;
    for degree in degrees do
        if not first then Print(","); fi;
        PrintMap(DifferentialOfComplex(resolution, degree), row_printer);
        first := false;
    od;
    Print("]}");
end;

PrintCase := function(field, row_printer)
    local a2, a2_zero, simple_two, zero, loop_parent, dual, dual_simple;
    a2 := PathAlgebra(field,
        Quiver(["t1", "t2"], [["t1", "t2", "c"]]));
    simple_two := RightModuleOverPathAlgebra(a2, [0, 1], []);
    a2_zero := PathAlgebra(field,
        Quiver(["z1", "z2"], [["z1", "z2", "z"]]));
    zero := ZeroModule(a2_zero);

    loop_parent := PathAlgebra(field,
        Quiver(["w"], [["w", "w", "x"]]));
    dual := loop_parent / [loop_parent.x^2];
    dual_simple := RightModuleOverPathAlgebra(dual, [1],
        [["x", [[0]] * One(field)]]);

    Print("{\"dual_simple_resolution\":");
    PrintResolution(dual_simple, [-3,-2,-1,0,1,2], row_printer);
    Print(",\n\"a2_simple_two_resolution\":");
    PrintResolution(simple_two, [-3,-2,-1,0,1,2], row_printer);
    Print(",\n\"zero_resolution\":");
    PrintResolution(zero, [-3,-2,-1,0,1,2], row_printer);
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"complex_injective_resolution\",\"oracle\":{");
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
