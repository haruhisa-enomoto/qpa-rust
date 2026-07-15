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
    local data, index;
    data := CoefficientsAndMagmaElements(element);
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
    data := CoefficientsAndMagmaElements(element);
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
    for index in [1 .. Length(elements)] do
        if index > 1 then Print(","); fi;
        printer(elements[index]);
    od;
    Print("]");
end;

Qfinite := Quiver(2, [[1, 2, "a"]]);;
P := PathAlgebra(Rationals, Qfinite);;
P2 := PathAlgebra(GF(2), Qfinite);;

Q := Quiver(1, [[1, 1, "x"], [1, 1, "y"]]);;
A := PathAlgebra(Rationals, Q);;
relations := [A.x^2 + A.y^3];;
AddNthPowerToRelations(A, relations, 4);;
B := A/relations;;
M := AssociatedMonomialAlgebra(B);;
G := AssociatedGradedAlgebra(B);;
T := TruncatedPathAlgebra(Rationals, Q, 3);;
TM := AssociatedMonomialAlgebra(T);;

A2 := PathAlgebra(GF(2), Q);;
relations2 := [A2.x^2 + A2.y^3];;
AddNthPowerToRelations(A2, relations2, 4);;
B2 := A2/relations2;;
M2 := AssociatedMonomialAlgebra(B2);;
G2 := AssociatedGradedAlgebra(B2);;
T2 := TruncatedPathAlgebra(GF(2), Q, 3);;
TM2 := AssociatedMonomialAlgebra(T2);;

Print("{\"schema_version\":1,\"kind\":\"pathalg_associated\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},\"rational\":{");
Print("\"finite_path_is_monomial\":", IsMonomialAlgebra(P));
Print(",\"finite_path_identity\":", IsIdenticalObj(P, AssociatedMonomialAlgebra(P)));
Print(",\"finite_path_graded_identity\":", IsIdenticalObj(P, AssociatedGradedAlgebra(P)));
Print(",\"nonmonomial\":", IsMonomialAlgebra(B));
Print(",\"monomial_identity\":", IsIdenticalObj(B, M));
Print(",\"monomial_result_is_monomial\":", IsMonomialAlgebra(M));
Print(",\"monomial_dimension\":", Dimension(M));
Print(",\"monomial_relations\":");
PrintElementList(RelationsOfAlgebra(M), PrintRationalElement);
Print(",\"graded_identity\":", IsIdenticalObj(B, G));
Print(",\"graded_dimension\":", Dimension(G));
Print(",\"graded_relations\":");
PrintElementList(RelationsOfAlgebra(G), PrintRationalElement);
Print(",\"already_monomial_identity\":", IsIdenticalObj(T, TM));
Print("},\"f2\":{");
Print("\"finite_path_is_monomial\":", IsMonomialAlgebra(P2));
Print(",\"finite_path_identity\":", IsIdenticalObj(P2, AssociatedMonomialAlgebra(P2)));
Print(",\"finite_path_graded_identity\":", IsIdenticalObj(P2, AssociatedGradedAlgebra(P2)));
Print(",\"nonmonomial\":", IsMonomialAlgebra(B2));
Print(",\"monomial_identity\":", IsIdenticalObj(B2, M2));
Print(",\"monomial_result_is_monomial\":", IsMonomialAlgebra(M2));
Print(",\"monomial_dimension\":", Dimension(M2));
Print(",\"monomial_relations\":");
PrintElementList(RelationsOfAlgebra(M2), PrintPrimeElement);
Print(",\"graded_identity\":", IsIdenticalObj(B2, G2));
Print(",\"graded_dimension\":", Dimension(G2));
Print(",\"graded_relations\":");
PrintElementList(RelationsOfAlgebra(G2), PrintPrimeElement);
Print(",\"already_monomial_identity\":", IsIdenticalObj(T2, TM2));
Print("}}");
Print("\n");

QUIT_GAP(0);
