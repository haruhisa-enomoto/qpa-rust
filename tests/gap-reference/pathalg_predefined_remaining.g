if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintRationalElement := function(element)
    local data, i;
    data := CoefficientsAndMagmaElements(element);
    Print("[");
    for i in [1, 3 .. Length(data) - 1] do
        if i > 1 then Print(","); fi;
        Print("{\"path\":", ExtRepOfObj(data[i]));
        Print(",\"coefficient\":[", NumeratorRat(data[i + 1]), ",",
              DenominatorRat(data[i + 1]), "]}");
    od;
    Print("]");
end;

PrintPrimeElement := function(element)
    local data, i;
    data := CoefficientsAndMagmaElements(element);
    Print("[");
    for i in [1, 3 .. Length(data) - 1] do
        if i > 1 then Print(","); fi;
        Print("{\"path\":", ExtRepOfObj(data[i]));
        Print(",\"coefficient\":", IntFFE(data[i + 1]), "}");
    od;
    Print("]");
end;

PrintElementList := function(elements, printer)
    local first, element;
    Print("[");
    first := true;
    for element in elements do
        if not first then Print(","); fi;
        printer(element);
        first := false;
    od;
    Print("]");
end;

PrintAlgebra := function(algebra, printer)
    local quiver, first, vertex, arrow;
    quiver := QuiverOfPathAlgebra(algebra);
    Print("{\"path_algebra\":", IsPathAlgebra(algebra));
    Print(",\"dimension\":", Dimension(algebra));
    Print(",\"vertices\":[");
    first := true;
    for vertex in VerticesOfQuiver(quiver) do
        if not first then Print(","); fi;
        PrintJsonString(String(vertex));
        first := false;
    od;
    Print("],\"arrows\":[");
    first := true;
    for arrow in ArrowsOfQuiver(quiver) do
        if not first then Print(","); fi;
        Print("{\"name\":");
        PrintJsonString(String(arrow));
        Print(",\"source\":", VertexIndex(SourceOfPath(arrow)));
        Print(",\"target\":", VertexIndex(TargetOfPath(arrow)), "}");
        first := false;
    od;
    Print("],\"relations\":");
    PrintElementList(RelationsOfAlgebra(algebra), printer);
    Print("}");
end;

a2 := PathAlgebra(Rationals, Quiver(2, [[1, 2, "x"]]));;
preprojective_a2 := PreprojectiveAlgebra(a2);;
green_2 := FiniteGlobalDimGreen(Rationals, 2);;
kirk_kuz_1 := FiniteGlobalDimKirkKuz(Rationals, 1);;

F2 := GF(2);;
a3_f2 := PathAlgebra(F2,
    Quiver(3, [[1, 2, "u"], [2, 3, "v"]]));;
preprojective_a3_f2 := PreprojectiveAlgebra(a3_f2);;
green_2_f2 := FiniteGlobalDimGreen(F2, 2);;
kirk_kuz_2_f2 := FiniteGlobalDimKirkKuz(F2, 2);;

Print("{\"schema_version\":1,\"kind\":\"pathalg_predefined_remaining\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":");
PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"sequences\":{");
Print("\"one_by_three\":", AdmissibleSequenceGenerator(1, 3));
Print(",\"three_by_three\":", AdmissibleSequenceGenerator(3, 3));
Print(",\"four_by_two\":", AdmissibleSequenceGenerator(4, 2));
Print(",\"one_by_one\":", AdmissibleSequenceGenerator(1, 1));
Print(",\"three_by_one\":", AdmissibleSequenceGenerator(3, 1), "}");
Print(",\"rational\":{");
Print("\"preprojective_a2\":");
PrintAlgebra(preprojective_a2, PrintRationalElement);
Print(",\"green_2\":");
PrintAlgebra(green_2, PrintRationalElement);
Print(",\"kirk_kuz_1\":");
PrintAlgebra(kirk_kuz_1, PrintRationalElement);
Print("},\"f2\":{");
Print("\"preprojective_a3\":");
PrintAlgebra(preprojective_a3_f2, PrintPrimeElement);
Print(",\"green_2\":");
PrintAlgebra(green_2_f2, PrintPrimeElement);
Print(",\"kirk_kuz_2\":");
PrintAlgebra(kirk_kuz_2_f2, PrintPrimeElement);
Print("}}\n");

QUIT_GAP(0);
