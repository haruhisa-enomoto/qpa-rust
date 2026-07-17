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

PrintIntegerMatrix := function(matrix)
    local first, row;
    Print("["); first := true;
    for row in matrix do
        if not first then Print(","); fi;
        PrintIntegerList(row); first := false;
    od;
    Print("]");
end;

PrintPrimeElement := function(element)
    local data, first, i;
    data := CoefficientsAndMagmaElements(element);
    Print("["); first := true;
    for i in [1,3..Length(data)-1] do
        if not first then Print(","); fi;
        Print("{\"walk\":");
        PrintIntegerList(List(WalkOfPath(data[i]), ArrowIndex));
        Print(",\"coefficient\":", IntFFE(data[i+1]), "}");
        first := false;
    od;
    Print("]");
end;

PrintRationalElement := function(element)
    local data, first, i;
    data := CoefficientsAndMagmaElements(element);
    Print("["); first := true;
    for i in [1,3..Length(data)-1] do
        if not first then Print(","); fi;
        Print("{\"walk\":");
        PrintIntegerList(List(WalkOfPath(data[i]), ArrowIndex));
        Print(",\"coefficient\":[", NumeratorRat(data[i+1]), ",",
              DenominatorRat(data[i+1]), "]}");
        first := false;
    od;
    Print("]");
end;

PrintAlgebra := function(algebra, element_printer)
    local quiver, first, arrow, relation;
    quiver := QuiverOfPathAlgebra(algebra);
    Print("{\"dimension\":", Dimension(algebra));
    Print(",\"adjacency\":");
    PrintIntegerMatrix(AdjacencyMatrixOfQuiver(quiver));
    Print(",\"arrows\":["); first := true;
    for arrow in ArrowsOfQuiver(quiver) do
        if not first then Print(","); fi;
        Print("[", VertexIndex(SourceOfPath(arrow)), ",",
              VertexIndex(TargetOfPath(arrow)), "]");
        first := false;
    od;
    Print("],\"relations\":["); first := true;
    for relation in RelationsOfAlgebra(algebra) do
        if not first then Print(","); fi;
        element_printer(relation); first := false;
    od;
    Print("]}");
end;

PrintAlgebraCase := function(field, element_printer)
    local quiver, algebra, adr, loop_quiver, loop_parent, dual_numbers,
          dual_adr;
    quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
    algebra := PathAlgebra(field, quiver);
    adr := ADRAlgebraOfAlgebra(algebra);
    loop_quiver := Quiver(["w"], [["w", "w", "x"]]);
    loop_parent := PathAlgebra(field, loop_quiver);
    dual_numbers := loop_parent / [loop_parent.x^2];
    dual_adr := ADRAlgebraOfAlgebra(dual_numbers);

    Print("{\"a2\":"); PrintAlgebra(adr, element_printer);
    Print(",\"a2_finite_global_dimension_cached\":",
        HasIsFiniteGlobalDimensionAlgebra(adr) and
        IsFiniteGlobalDimensionAlgebra(adr));
    Print(",\"dual_numbers\":"); PrintAlgebra(dual_adr, element_printer);
    Print(",\"dual_numbers_finite_global_dimension_cached\":",
        HasIsFiniteGlobalDimensionAlgebra(dual_adr) and
        IsFiniteGlobalDimensionAlgebra(dual_adr));
    Print("}");
end;

PrintBasicCase := function(field, element_printer)
    local quiver, algebra, projective, adr;
    quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
    algebra := PathAlgebra(field, quiver);
    projective := RightModuleOverPathAlgebra(algebra, [1, 1],
        [["a", [[1]] * One(field)]]);
    adr := ADRAlgebraOfBasicModule([projective]);
    PrintAlgebra(adr, element_printer);
end;

Print("{\"schema_version\":1,\"kind\":\"modulehom_adr\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"binary\":{\"algebra\":");
PrintAlgebraCase(GF(2), PrintPrimeElement);
Print(",\"basic_module\":"); PrintBasicCase(GF(2), PrintPrimeElement);
Print("},\"quinary\":{\"algebra\":");
PrintAlgebraCase(GF(5), PrintPrimeElement);
Print(",\"basic_module\":"); PrintBasicCase(GF(5), PrintPrimeElement);
Print("},\"rational\":{\"algebra\":");
PrintAlgebraCase(Rationals, PrintRationalElement);
Print("}}\n");

QUIT_GAP(0);
