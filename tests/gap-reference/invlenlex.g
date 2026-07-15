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

PrintRationalElementList := function(elements)
    local i;
    Print("[");
    for i in [1 .. Length(elements)] do
        if i > 1 then Print(","); fi;
        PrintRationalElement(elements[i]);
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

PrintPrimeElementList := function(elements)
    local i;
    Print("[");
    for i in [1 .. Length(elements)] do
        if i > 1 then Print(","); fi;
        PrintPrimeElement(elements[i]);
    od;
    Print("]");
end;

Q := Quiver(1, [[1, 1, "x"]]);;
KQ := PathAlgebra(Rationals, Q);;
v := ElementOfPathAlgebra(KQ, Q.v1);;
x := ElementOfPathAlgebra(KQ, Q.x);;
x2 := x * x;;
x3 := x2 * x;;
mixed := v + 2 * x + 3 * x2 + 4 * x3;;
positive := 2 * x + 3 * x2;;
inverse_remainder := RemainderOfDivisionInvLenInvLexGroebnerBasis(x3 + x, [x2], KQ);;
inverse_reduced := ReducedListInvLenInvLexQPA([x2, x3], KQ);;
inverse_tip_reduced := TipReducedListInvLenInvLex([x2, x3], KQ);;
inverse_basis := GroebnerBasisInvLenInvLex([x3], KQ);;

F5 := GF(5);;
KQ5 := PathAlgebra(F5, Q);;
v5 := ElementOfPathAlgebra(KQ5, Q.v1);;
x5 := ElementOfPathAlgebra(KQ5, Q.x);;
x25 := x5 * x5;;
x35 := x25 * x5;;
mixed5 := v5 + 2 * x5 + 3 * x25 + 4 * x35;;
inverse_basis5 := GroebnerBasisInvLenInvLex([x35], KQ5);;

F2 := GF(2);;
KQ2 := PathAlgebra(F2, Q);;
v2 := ElementOfPathAlgebra(KQ2, Q.v1);;
x_2 := ElementOfPathAlgebra(KQ2, Q.x);;
x22 := x_2 * x_2;;
x32 := x22 * x_2;;
mixed2 := v2 + x_2 + x22 + x32;;
inverse_basis2 := GroebnerBasisInvLenInvLex([x32], KQ2);;

Print("{\"schema_version\":1,\"kind\":\"invlenlex\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},\"rational\":{\"tip\":");
PrintRationalElement(TipInvLenInvLex(mixed));
Print(",\"tip_coefficient\":[", NumeratorRat(TipCoefficientInvLenInvLex(mixed)), ",",
      DenominatorRat(TipCoefficientInvLenInvLex(mixed)), "]");
Print(",\"tip_monomial\":", ExtRepOfObj(TipMonomialInvLenInvLex(mixed)));
Print(",\"tip_walk\":", List(TipWalkInvLenInvLex(mixed), ExtRepOfObj));
Print(",\"truncate_1\":");
PrintRationalElement(TruncateElement(mixed, 1));
Print(",\"truncate_2\":");
PrintRationalElement(TruncateElement(mixed, 2));
Print(",\"biggest_homogeneous\":");
PrintRationalElement(BiggestHomogeneousPartInvLenInvLex(positive));
Print(",\"biggest_homogeneous_list\":");
PrintRationalElementList(BiggestHomogeneousPartInvLenInvLex([mixed, positive]));
Print(",\"remainder\":");
PrintRationalElement(inverse_remainder);
Print(",\"reduced_list\":");
PrintRationalElementList(inverse_reduced);
Print(",\"tip_reduced_list\":");
PrintRationalElementList(inverse_tip_reduced);
Print(",\"leftmost_occurrence\":", LeftmostOccurrenceInvLenInvLex([1, 2, 1, 2], [1, 2]));
Print(",\"basis\":");
PrintRationalElementList(inverse_basis);
Print("},\"f5\":{\"tip\":");
PrintPrimeElement(TipInvLenInvLex(mixed5));
Print(",\"truncate_2\":");
PrintPrimeElement(TruncateElement(mixed5, 2));
Print(",\"basis\":");
PrintPrimeElementList(inverse_basis5);
Print("},\"f2\":{\"tip\":");
PrintPrimeElement(TipInvLenInvLex(mixed2));
Print(",\"truncate_2\":");
PrintPrimeElement(TruncateElement(mixed2, 2));
Print(",\"basis\":");
PrintPrimeElementList(inverse_basis2);
Print("}}\n");

QUIT_GAP(0);
