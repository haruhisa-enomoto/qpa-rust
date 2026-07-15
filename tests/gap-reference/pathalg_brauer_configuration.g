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

PrintOrientations := function(orientations)
    local first_row, first_name, row, name;
    Print("[");
    first_row := true;
    for row in orientations do
        if not first_row then Print(","); fi;
        Print("[");
        first_name := true;
        for name in row do
            if not first_name then Print(","); fi;
            PrintJsonString(name);
            first_name := false;
        od;
        Print("]");
        first_row := false;
    od;
    Print("]");
end;

trivial_input := [
    [["v", 1]],
    [["e", "v"]],
    [["e"]]
];;
trivial := BrauerConfigurationAlgebra(Rationals, trivial_input);;

documented_input := [
    [["v1", 1], ["v2", 1], ["v3", 2]],
    [["e1", "v1", "v2"], ["e2", "v2", "v3"]],
    [["e1"], ["e2", "e1"], ["e2"]]
];;
documented := BrauerConfigurationAlgebra(Rationals, documented_input);;

F2 := GF(2);;
repeated_input := [
    [["a", 1], ["b", 1], ["c", 1]],
    [["1", "a", "a"], ["2", "a", "b"],
     ["3", "b", "a"], ["4", "b", "c"]],
    [["1", "2", "3"], ["2", "4", "3"], ["4"]]
];;
repeated := BrauerConfigurationAlgebra(F2, repeated_input);;

Print("{\"schema_version\":1,\"kind\":\"pathalg_brauer_configuration\",\"oracle\":{");
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
Print("},\"rational\":{");
Print("\"trivial\":");
PrintAlgebra(trivial, PrintRationalElement);
Print(",\"documented\":");
PrintAlgebra(documented, PrintRationalElement);
Print(",\"orientations_after\":");
PrintOrientations(documented_input[3]);
Print(",\"properties\":[", IsSymmetricAlgebra(trivial), ",",
      IsSymmetricAlgebra(documented), "]}");
Print(",\"f2\":{\"repeated\":");
PrintAlgebra(repeated, PrintPrimeElement);
Print(",\"orientations_after\":");
PrintOrientations(repeated_input[3]);
Print(",\"properties\":[", IsSymmetricAlgebra(repeated), "]}}");
Print("\n");

QUIT_GAP(0);
