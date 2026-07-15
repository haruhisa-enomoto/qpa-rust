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

PrintPathList := function(paths)
    local i;
    Print("[");
    for i in [1 .. Length(paths)] do
        if i > 1 then Print(","); fi;
        Print(ExtRepOfObj(paths[i]));
    od;
    Print("]");
end;

Q := Quiver(1, [[1, 1, "x"]]);;
KQ := PathAlgebra(Rationals, Q);;
v := ElementOfPathAlgebra(KQ, Q.v1);;
x := ElementOfPathAlgebra(KQ, Q.x);;
redundant := [x, x * x];;
I := Ideal(KQ, redundant);;
gb := GroebnerBasis(I, redundant);;
initial := ShallowCopy(AsList(gb));;
initial_flags := [HasIsTipReducedGroebnerBasis(gb),
                  HasIsCompletelyReducedGroebnerBasis(gb),
                  IsCompleteGroebnerBasis(gb)];;
gb := TipReduceGroebnerBasis(gb);;
tip_reduced := ShallowCopy(AsList(gb));;
tip_reduction := TipReduce(gb, x * x * x + v);;
complete_reduction := CompletelyReduce(gb, x * x * x + x + v);;
gb := CompletelyReduceGroebnerBasis(gb);;
completely_reduced := ShallowCopy(AsList(gb));;
final_flags := [IsTipReducedGroebnerBasis(gb),
                IsCompletelyReducedGroebnerBasis(gb),
                IsCompleteGroebnerBasis(gb)];;

cubic := x * x * x;;
J := Ideal(KQ, [cubic]);;
cubic_gb := GroebnerBasis(J, [cubic]);;
nontips := Nontips(cubic_gb);;
nontip_size := NontipSize(cubic_gb);;
prefix_short := IsPrefixOfTipInTipIdeal(cubic_gb, x * x);;
prefix_equal := IsPrefixOfTipInTipIdeal(cubic_gb, cubic);;
prefix_long := IsPrefixOfTipInTipIdeal(cubic_gb, cubic * x);;
right_gb := RightGroebnerBasis(J);;
right_gb := TipReduceGroebnerBasis(right_gb);;
right_relations := ShallowCopy(AsList(right_gb));;
right_reduction := TipReduce(right_gb, cubic * x + x);;

F5 := GF(5);;
KQ5 := PathAlgebra(F5, Q);;
v5 := ElementOfPathAlgebra(KQ5, Q.v1);;
x5 := ElementOfPathAlgebra(KQ5, Q.x);;
redundant5 := [x5, x5 * x5];;
I5 := Ideal(KQ5, redundant5);;
gb5 := GroebnerBasis(I5, redundant5);;
gb5 := TipReduceGroebnerBasis(gb5);;
tip_reduction5 := TipReduce(gb5, 3 * x5 * x5 + v5);;

F2 := GF(2);;
KQ2 := PathAlgebra(F2, Q);;
v2 := ElementOfPathAlgebra(KQ2, Q.v1);;
x2 := ElementOfPathAlgebra(KQ2, Q.x);;
redundant2 := [x2, x2 * x2];;
I2 := Ideal(KQ2, redundant2);;
gb2 := GroebnerBasis(I2, redundant2);;
gb2 := TipReduceGroebnerBasis(gb2);;
tip_reduction2 := TipReduce(gb2, x2 * x2 + v2);;

Print("{\"schema_version\":1,\"kind\":\"groebner_objects\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},\"rational\":{\"initial\":");
PrintRationalElementList(initial);
Print(",\"initial_flags\":", initial_flags);
Print(",\"tip_reduced\":");
PrintRationalElementList(tip_reduced);
Print(",\"tip_reduction\":");
PrintRationalElement(tip_reduction);
Print(",\"complete_reduction\":");
PrintRationalElement(complete_reduction);
Print(",\"completely_reduced\":");
PrintRationalElementList(completely_reduced);
Print(",\"final_flags\":", final_flags);
Print(",\"nontips\":");
PrintPathList(nontips);
Print(",\"nontip_size\":", nontip_size);
Print(",\"prefix\":[", prefix_short, ",", prefix_equal, ",", prefix_long, "]");
Print(",\"right_relations\":");
PrintRationalElementList(right_relations);
Print(",\"right_reduction\":");
PrintRationalElement(right_reduction);
Print("},\"f5\":{\"tip_reduced\":");
PrintPrimeElementList(AsList(gb5));
Print(",\"tip_reduction\":");
PrintPrimeElement(tip_reduction5);
Print("},\"f2\":{\"tip_reduced\":");
PrintPrimeElementList(AsList(gb2));
Print(",\"tip_reduction\":");
PrintPrimeElement(tip_reduction2);
Print("}}\n");

QUIT_GAP(0);
