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

singleton_poset := Poset(["only"], []);;
singleton := PosetAlgebra(Rationals, singleton_poset);;
chain_poset := Poset(["c", "a", "b"], [["a", "b"], ["b", "c"]]);;
chain := PosetAlgebra(Rationals, chain_poset);;
diamond_poset := Poset([4, 2, 1, 3], [[1, 2, 3], [2, 4], [3, 4]]);;
diamond := PosetAlgebra(Rationals, diamond_poset);;

F2 := GF(2);;
parallel_poset := Poset([1, 2, 3, 4, 5],
    [[1, 2, 3, 4], [2, 5], [3, 5], [4, 5]]);;
parallel := PosetAlgebra(F2, parallel_poset);;
non_graded_poset := Poset([1, 2, 3, 4, 5],
    [[1, 2, 3], [2, 5], [3, 4], [4, 5]]);;
non_graded := PosetAlgebra(F2, non_graded_poset);;

Print("{\"schema_version\":1,\"kind\":\"pathalg_poset_algebra\",\"oracle\":{");
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
Print("\"singleton\":");
PrintAlgebra(singleton, PrintRationalElement);
Print(",\"chain\":");
PrintAlgebra(chain, PrintRationalElement);
Print(",\"diamond\":");
PrintAlgebra(diamond, PrintRationalElement);
Print(",\"attribute_identity\":[",
      IsIdenticalObj(PosetOfPosetAlgebra(singleton), singleton_poset), ",",
      IsIdenticalObj(PosetOfPosetAlgebra(chain), chain_poset), ",",
      IsIdenticalObj(PosetOfPosetAlgebra(diamond), diamond_poset), "]");
Print(",\"properties\":[", IsPosetAlgebra(diamond), ",",
      IsFiniteGlobalDimensionAlgebra(diamond), "]}");
Print(",\"f2\":{");
Print("\"parallel\":");
PrintAlgebra(parallel, PrintPrimeElement);
Print(",\"non_graded\":");
PrintAlgebra(non_graded, PrintPrimeElement);
Print(",\"attribute_identity\":[",
      IsIdenticalObj(PosetOfPosetAlgebra(parallel), parallel_poset), ",",
      IsIdenticalObj(PosetOfPosetAlgebra(non_graded), non_graded_poset), "]");
Print(",\"properties\":[", IsPosetAlgebra(non_graded), ",",
      IsFiniteGlobalDimensionAlgebra(non_graded), "]}}");
Print("\n");

QUIT_GAP(0);
