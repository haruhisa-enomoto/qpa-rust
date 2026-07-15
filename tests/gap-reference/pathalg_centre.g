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

PrintRationalList := function(elements)
    local index;
    Print("[");
    for index in [1..Length(elements)] do
        if index > 1 then Print(","); fi;
        PrintRationalElement(elements[index]);
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

Q := Quiver(2, [[1, 2, "a"], [2, 1, "b"]]);;
P := PathAlgebra(Rationals, Q);;
CP := Centre(P);;

A := TruncatedPathAlgebra(Rationals, Q, 4);;
C := Centre(A);;
CB := BasisVectors(Basis(C));;
central := A.a*A.b+A.b*A.a;;

Qacyclic := Quiver(3, [[1, 2, "u"], [2, 3, "v"]]);;
Pacyclic := PathAlgebra(Rationals, Qacyclic);;
Cacyclic := Centre(Pacyclic);;

Qdiscrete := Quiver(2, []);;
Pdiscrete := PathAlgebra(Rationals, Qdiscrete);;
Cdiscrete := Centre(Pdiscrete);;

A2 := TruncatedPathAlgebra(GF(2), Q, 4);;
C2 := Centre(A2);;
CB2 := BasisVectors(Basis(C2));;

Print("{\"schema_version\":1,\"kind\":\"pathalg_centre\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},\"parent_cycle\":{\"generators\":");
PrintRationalList(GeneratorsOfAlgebraWithOne(CP));
Print("},\"parent_acyclic\":{\"generators\":");
PrintRationalList(GeneratorsOfAlgebraWithOne(Cacyclic));
Print("},\"parent_discrete\":{\"generators\":");
PrintRationalList(GeneratorsOfAlgebraWithOne(Cdiscrete));
Print("},\"quotient\":{\"ambient_dimension\":", Dimension(A));
Print(",\"basis\":");
PrintRationalList(CB);
Print(",\"dimension\":", Dimension(C));
Print(",\"central_member\":", central in C);
Print(",\"arrow_member\":", A.a in C);
Print("},\"f2\":{\"basis\":[");
for index in [1..Length(CB2)] do
    if index > 1 then Print(","); fi;
    PrintPrimeElement(CB2[index]);
od;
Print("],\"dimension\":", Dimension(C2));
Print("}}");
Print("\n");

QUIT_GAP(0);
