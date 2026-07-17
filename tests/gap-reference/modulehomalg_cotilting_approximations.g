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

PrintCheck := function(cotilting, module, row_printer)
    Print("{\"right\":");
    PrintMap(RightApproximationByPerpT(cotilting, module), row_printer);
    Print(",\"left\":");
    PrintMap(LeftApproximationByAddTHat(cotilting, module), row_printer);
    Print("}");
end;

PrintCase := function(field, row_printer)
    local a2, injectives, simples, projectives, cotilting, zero;
    a2 := PathAlgebra(field,
        Quiver(["v1", "v2"], [["v1", "v2", "a"]]));
    injectives := IndecInjectiveModules(a2);
    simples := SimpleModules(a2);
    projectives := IndecProjectiveModules(a2);
    cotilting := DirectSumOfQPAModules([injectives[2], simples[2]]);
    zero := ZeroModule(a2);
    if CotiltingModule(cotilting, 1) = false then
        Error("failed to install the A2 cotilting attribute");
    fi;

    Print("{\"zero\":"); PrintCheck(cotilting, zero, row_printer);
    Print(",\n\"simple_one\":"); PrintCheck(cotilting, simples[1], row_printer);
    Print(",\n\"simple_two\":"); PrintCheck(cotilting, simples[2], row_printer);
    Print(",\n\"projective_one\":"); PrintCheck(cotilting, projectives[1], row_printer);
    Print(",\n\"projective_two\":"); PrintCheck(cotilting, projectives[2], row_printer);
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"modulehomalg_cotilting_approximations\",\"oracle\":{");
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
