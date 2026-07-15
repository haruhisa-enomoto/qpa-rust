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
    Print(",\"admissible_quotient\":");
    if IsPathAlgebra(algebra) then
        Print("false");
    else
        Print(IsAdmissibleQuotientOfPathAlgebra(algebra));
    fi;
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

n_singleton := NakayamaAlgebra(Rationals, [1]);;
n_singleton_quotient := NakayamaAlgebra(Rationals, [2]);;
n_linear_path := NakayamaAlgebra(Rationals, [3, 2, 1]);;
n_linear_quotient := NakayamaAlgebra(Rationals, [3, 3, 2, 1]);;
n_cyclic := NakayamaAlgebra(Rationals, [3, 3, 3]);;
n_long_cyclic := NakayamaAlgebra(Rationals, [7, 6, 7, 7]);;
canonical_two := CanonicalAlgebra(Rationals, [2, 3]);;
canonical_four := CanonicalAlgebra(Rationals, [2, 2, 3, 4], [2, -3]);;

F2 := GF(2);;
n_cyclic_f2 := NakayamaAlgebra(F2, [3, 3, 3]);;
canonical_f2 := CanonicalAlgebra(F2, [2, 3, 4], [One(F2)]);;
kronecker_f2 := KroneckerAlgebra(F2, 3);;

Print("{\"schema_version\":1,\"kind\":\"pathalg_predefined\",\"oracle\":{");
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
Print("\"nakayama_singleton\":");
PrintAlgebra(n_singleton, PrintRationalElement);
Print(",\"nakayama_singleton_quotient\":");
PrintAlgebra(n_singleton_quotient, PrintRationalElement);
Print(",\"nakayama_linear_path\":");
PrintAlgebra(n_linear_path, PrintRationalElement);
Print(",\"nakayama_linear_quotient\":");
PrintAlgebra(n_linear_quotient, PrintRationalElement);
Print(",\"nakayama_cyclic\":");
PrintAlgebra(n_cyclic, PrintRationalElement);
Print(",\"nakayama_long_cyclic\":");
PrintAlgebra(n_long_cyclic, PrintRationalElement);
Print(",\"canonical_two\":");
PrintAlgebra(canonical_two, PrintRationalElement);
Print(",\"canonical_four\":");
PrintAlgebra(canonical_four, PrintRationalElement);
Print(",\"properties\":{");
Print("\"nakayama\":[", IsNakayamaAlgebra(n_cyclic), ",",
      IsFiniteTypeAlgebra(n_cyclic), "]");
Print(",\"canonical\":[", IsCanonicalAlgebra(canonical_four), ",",
      IsFiniteGlobalDimensionAlgebra(canonical_four), "]}}");
Print(",\"f2\":{");
Print("\"nakayama_cyclic\":");
PrintAlgebra(n_cyclic_f2, PrintPrimeElement);
Print(",\"canonical\":");
PrintAlgebra(canonical_f2, PrintPrimeElement);
Print(",\"kronecker\":");
PrintAlgebra(kronecker_f2, PrintPrimeElement);
Print(",\"properties\":{");
Print("\"nakayama\":[", IsNakayamaAlgebra(n_cyclic_f2), ",",
      IsFiniteTypeAlgebra(n_cyclic_f2), "]");
Print(",\"canonical\":[", IsCanonicalAlgebra(canonical_f2), ",",
      IsFiniteGlobalDimensionAlgebra(canonical_f2), "]");
Print(",\"kronecker\":[", IsKroneckerAlgebra(kronecker_f2), ",",
      IsFiniteGlobalDimensionAlgebra(kronecker_f2), ",",
      IsHereditaryAlgebra(kronecker_f2), "]}}");
Print("}\n");

QUIT_GAP(0);
