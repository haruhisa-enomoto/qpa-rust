if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

UnderlyingRepresentative := function(element)
    if IsElementOfQuotientOfPathAlgebra(element) then
        return element![1];
    fi;
    return element;
end;

PrintRationalElement := function(element)
    local data, index;
    data := CoefficientsAndMagmaElements(UnderlyingRepresentative(element));
    Print("[");
    for index in [1, 3 .. Length(data) - 1] do
        if index > 1 then Print(","); fi;
        Print("{\"path\":", ExtRepOfObj(data[index]));
        Print(",\"coefficient\":[", NumeratorRat(data[index + 1]), ",",
              DenominatorRat(data[index + 1]), "]}");
    od;
    Print("]");
end;

PrintPrimeElement := function(element)
    local data, index;
    data := CoefficientsAndMagmaElements(UnderlyingRepresentative(element));
    Print("[");
    for index in [1, 3 .. Length(data) - 1] do
        if index > 1 then Print(","); fi;
        Print("{\"path\":", ExtRepOfObj(data[index]));
        Print(",\"coefficient\":", IntFFE(data[index + 1]), "}");
    od;
    Print("]");
end;

PrintElementList := function(elements, printer)
    local index;
    Print("[");
    for index in [1..Length(elements)] do
        if index > 1 then Print(","); fi;
        printer(elements[index]);
    od;
    Print("]");
end;

PrintQuiver := function(Q)
    local vertices, arrows, index;
    vertices := VerticesOfQuiver(Q);
    arrows := ArrowsOfQuiver(Q);
    Print("{\"vertices\":[");
    for index in [1..Length(vertices)] do
        if index > 1 then Print(","); fi;
        PrintJsonString(String(vertices[index]));
    od;
    Print("],\"arrows\":[");
    for index in [1..Length(arrows)] do
        if index > 1 then Print(","); fi;
        Print("{\"source\":");
        PrintJsonString(String(SourceOfPath(arrows[index])));
        Print(",\"target\":");
        PrintJsonString(String(TargetOfPath(arrows[index])));
        Print(",\"name\":");
        PrintJsonString(String(arrows[index]));
        Print("}");
    od;
    Print("]}");
end;

Q := Quiver(4, [[1, 2, "a"], [2, 3, "b"], [1, 3, "c"],
                 [3, 4, "d"], [2, 4, "e"], [1, 4, "f"]]);;
A := PathAlgebra(Rationals, Q);;
AQ := A/[A.a*A.b-A.c, A.b*A.d-A.e, A.c*A.d-A.f];;
B2 := QuiverAlgebraOfAmodAeA(AQ, [2]);;
B3 := QuiverAlgebraOfAmodAeA(AQ, [3]);;
Bnone := QuiverAlgebraOfAmodAeA(AQ, []);;
Bignored := QuiverAlgebraOfAmodAeA(AQ, [0, 5, 99]);;
Bpath := QuiverAlgebraOfAmodAeA(A, [2]);;

A2 := PathAlgebra(GF(2), Q);;
AQ2 := A2/[A2.a*A2.b+A2.c, A2.b*A2.d+A2.e,
            A2.c*A2.d+A2.f];;
B22 := QuiverAlgebraOfAmodAeA(AQ2, [2]);;

Print("{\"schema_version\":1,\"kind\":\"pathalg_vertex_quotient\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},\"rational\":{\"remove_two\":{\"quiver\":");
PrintQuiver(QuiverOfPathAlgebra(B2));
Print(",\"relations\":");
PrintElementList(RelationsOfAlgebra(B2), PrintRationalElement);
Print(",\"dimension\":", Dimension(B2), "}");
Print(",\"remove_three\":{\"quiver\":");
PrintQuiver(QuiverOfPathAlgebra(B3));
Print(",\"relations\":");
PrintElementList(RelationsOfAlgebra(B3), PrintRationalElement);
Print(",\"dimension\":", Dimension(B3), "}");
Print(",\"none\":{\"quiver\":");
PrintQuiver(QuiverOfPathAlgebra(Bnone));
Print(",\"relations\":");
PrintElementList(RelationsOfAlgebra(Bnone), PrintRationalElement);
Print(",\"dimension\":", Dimension(Bnone), "}");
Print(",\"ignored_same_quiver\":",
      String(QuiverOfPathAlgebra(Bignored)) = String(QuiverOfPathAlgebra(Bnone)));
Print(",\"path_algebra_input\":{\"is_quotient\":",
      IsQuotientOfPathAlgebra(Bpath), ",\"relations\":");
PrintElementList(RelationsOfAlgebra(Bpath), PrintRationalElement);
Print("}}");
Print(",\"f2\":{\"remove_two\":{\"quiver\":");
PrintQuiver(QuiverOfPathAlgebra(B22));
Print(",\"relations\":");
PrintElementList(RelationsOfAlgebra(B22), PrintPrimeElement);
Print(",\"dimension\":", Dimension(B22), "}}}");
Print("\n");

QUIT_GAP(0);
