if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintJsonStringList := function(values)
    local first, value;
    Print("[");
    first := true;
    for value in values do
        if not first then
            Print(",");
        fi;
        PrintJsonString(value);
        first := false;
    od;
    Print("]");
end;

PrintRationalElement := function(element)
    local data, i;
    data := CoefficientsAndMagmaElements(element);
    Print("[");
    for i in [1, 3 .. Length(data) - 1] do
        if i > 1 then
            Print(",");
        fi;
        Print("{\"path\":");
        PrintJsonString(StringPrint(data[i]));
        Print(",\"path_external_representation\":", ExtRepOfObj(data[i]));
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
        if i > 1 then
            Print(",");
        fi;
        Print("{\"path\":");
        PrintJsonString(StringPrint(data[i]));
        Print(",\"path_external_representation\":", ExtRepOfObj(data[i]));
        Print(",\"coefficient\":", IntFFE(data[i + 1]), "}");
    od;
    Print("]");
end;

Q := Quiver(2, [[1, 1, "a"], [1, 2, "b"]]);;
a := Q.a;;
b := Q.b;;
aa := a * a;;

KQ := PathAlgebra(Rationals, Q);;
v1q := ElementOfPathAlgebra(KQ, Q.v1);;
aq := ElementOfPathAlgebra(KQ, a);;
bq := ElementOfPathAlgebra(KQ, b);;
xq := v1q + aq;;
squareq := xq * xq;;
mixedq := 5 * aq + (3 / 2) * ElementOfPathAlgebra(KQ, aa);;
uniformq := aq + bq;;
zeroq := Zero(KQ);;

KQ5 := PathAlgebra(GF(5), Q);;
v1f5 := ElementOfPathAlgebra(KQ5, Q.v1);;
af5 := ElementOfPathAlgebra(KQ5, a);;
f5square := (v1f5 + af5) * (v1f5 + af5);;
f5reduced := 4 * ElementOfPathAlgebra(KQ5, aa)
             + 3 * ElementOfPathAlgebra(KQ5, aa);;
f5cancelled := 2 * af5 + 3 * af5;;

KQ2 := PathAlgebra(GF(2), Q);;
v1f2 := ElementOfPathAlgebra(KQ2, Q.v1);;
af2 := ElementOfPathAlgebra(KQ2, a);;
f2square := (v1f2 + af2) * (v1f2 + af2);;
f2cancelled := af2 + af2;;

Line := Quiver(3, [[1, 2, "c"], [2, 3, "d"]]);;
KLine := PathAlgebra(Rationals, Line);;

Print("{\"schema_version\":1,\"kind\":\"path_algebra\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},\"rational\":{");
Print("\"zero\":");
PrintRationalElement(zeroq);
Print(",\"one\":");
PrintRationalElement(One(KQ));
Print(",\"x\":");
PrintRationalElement(xq);
Print(",\"x_squared\":");
PrintRationalElement(squareq);
Print(",\"mixed\":");
PrintRationalElement(mixedq);
Print(",\"zero_path_embedding\":");
PrintRationalElement(ElementOfPathAlgebra(KQ, Zero(Q)));
Print(",\"noncomposable_product\":");
PrintRationalElement(bq * aq);
Print(",\"leading\":{");
Print("\"term\":");
PrintRationalElement(LeadingTerm(mixedq));
Print(",\"coefficient\":[", NumeratorRat(LeadingCoefficient(mixedq)), ",",
      DenominatorRat(LeadingCoefficient(mixedq)), "]");
Print(",\"monomial\":");
PrintJsonString(StringPrint(LeadingMonomial(mixedq)));
Print(",\"zero_coefficient\":[", NumeratorRat(LeadingCoefficient(zeroq)), ",",
      DenominatorRat(LeadingCoefficient(zeroq)), "]");
Print(",\"zero_monomial\":");
PrintJsonString(StringPrint(LeadingMonomial(zeroq)));
Print("},\"uniformity\":{");
Print("\"zero\":[", IsLeftUniform(zeroq), ",", IsRightUniform(zeroq), ",",
      IsUniform(zeroq), "]");
Print(",\"a_plus_b\":[", IsLeftUniform(uniformq), ",",
      IsRightUniform(uniformq), ",", IsUniform(uniformq), "]}");
Print(",\"generator_terms\":[");
for i in [1 .. Length(GeneratorsOfAlgebra(KQ))] do
    if i > 1 then
        Print(",");
    fi;
    PrintRationalElement(GeneratorsOfAlgebra(KQ)[i]);
od;
Print("]");
Print(",\"global_dimension\":", GlobalDimension(KQ));
Print(",\"is_hereditary\":", IsHereditaryAlgebra(KQ));
Print("},\"prime_fields\":{");
Print("\"f5_square\":");
PrintPrimeElement(f5square);
Print(",\"f5_reduced\":");
PrintPrimeElement(f5reduced);
Print(",\"f5_cancelled\":");
PrintPrimeElement(f5cancelled);
Print(",\"f2_square\":");
PrintPrimeElement(f2square);
Print(",\"f2_cancelled\":");
PrintPrimeElement(f2cancelled);
Print("},\"acyclic_basis\":{");
Print("\"paths\":");
PrintJsonStringList(List(BasisVectors(CanonicalBasis(KLine)),
                         element -> StringPrint(LeadingMonomial(element))));
Print(",\"dimension\":", Dimension(KLine));
Print(",\"global_dimension\":", GlobalDimension(KLine));
Print("}}\n");

QUIT_GAP(0);
