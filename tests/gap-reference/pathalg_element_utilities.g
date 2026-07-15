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

Q := Quiver(4, [[1, 2, "a"], [2, 3, "b"], [1, 3, "c"],
                 [3, 4, "d"], [2, 4, "e"]]);;
A := PathAlgebra(Rationals, Q);;
x := 2*A.a*A.b + 3*A.c + 5*A.b*A.d + 7*A.d;;
removed_none := PathRemoval(x, []);;
removed_two := PathRemoval(x, [2]);;
removed_four := PathRemoval(x, [4]);;
removed_ignored := PathRemoval(x, [0, 5, 99]);;

m := A.v1 + 2*A.a + 3*A.a*A.b + 5*A.b + 7*A.d;;
left := LeftSupportOfQuiverAlgebraElement(A, m);;
right := RightSupportOfQuiverAlgebraElement(A, m);;
support := SupportOfQuiverAlgebraElement(A, m);;
zero_support := SupportOfQuiverAlgebraElement(A, Zero(A));;

AQ := A/[A.a*A.b-A.c, A.b*A.d];;
qm := AQ.v1 + 2*AQ.a + 3*AQ.b + 5*AQ.d;;
qmix := AQ.a*AQ.b + AQ.c;;
quotient_support := SupportOfQuiverAlgebraElement(AQ, qm);;
quotient_mix_support := SupportOfQuiverAlgebraElement(AQ, qmix);;
quotient_cancelled_support := SupportOfQuiverAlgebraElement(
    AQ, AQ.a*AQ.b-AQ.c);;

A2 := PathAlgebra(GF(2), Q);;
x2 := A2.a + A2.b*A2.d + A2.d;;
removed_two2 := PathRemoval(x2, [2]);;
removed_four2 := PathRemoval(x2, [4]);;
support2 := SupportOfQuiverAlgebraElement(A2, x2);;

Print("{\"schema_version\":1,\"kind\":\"pathalg_element_utilities\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},\"rational\":{\"path_removal\":{\"none\":");
PrintRationalElement(removed_none);
Print(",\"two\":");
PrintRationalElement(removed_two);
Print(",\"four\":");
PrintRationalElement(removed_four);
Print(",\"ignored\":");
PrintRationalElement(removed_ignored);
Print("},\"support\":{\"left\":", left, ",\"right\":", right,
      ",\"both\":", support, ",\"zero\":", zero_support, "}");
Print(",\"quotient\":{\"support\":", quotient_support,
      ",\"mix_support\":", quotient_mix_support,
      ",\"cancelled_support\":", quotient_cancelled_support, "}");
Print("},\"f2\":{\"removed_two\":");
PrintPrimeElement(removed_two2);
Print(",\"removed_four\":");
PrintPrimeElement(removed_four2);
Print(",\"support\":", support2, "}}");
Print("\n");

QUIT_GAP(0);
