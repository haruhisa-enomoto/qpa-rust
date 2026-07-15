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
        if i > 1 then
            Print(",");
        fi;
        Print("{\"path\":");
        PrintJsonString(StringPrint(data[i]));
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
        Print(",\"coefficient\":", IntFFE(data[i + 1]), "}");
    od;
    Print("]");
end;

PrintElementList := function(elements, printer)
    local first, element;
    Print("[");
    first := true;
    for element in elements do
        if not first then
            Print(",");
        fi;
        printer(element);
        first := false;
    od;
    Print("]");
end;

Q := Quiver(1, [[1, 1, "x"], [1, 1, "y"]]);;

KQ := PathAlgebra(Rationals, Q, HighLevelGroebnerBasis);;
xq := ElementOfPathAlgebra(KQ, Q.x);;
yq := ElementOfPathAlgebra(KQ, Q.y);;
relationsq := [xq^2, xq*yq-yq*xq, yq^2];;
Iq := Ideal(KQ, relationsq);;
Aq := KQ / relationsq;;
xaq := Aq.x;;
yaq := Aq.y;;

KQ5 := PathAlgebra(GF(5), Q, HighLevelGroebnerBasis);;
x5 := ElementOfPathAlgebra(KQ5, Q.x);;
y5 := ElementOfPathAlgebra(KQ5, Q.y);;
A5 := KQ5 / [x5^2, x5*y5-2*y5*x5, y5^2];;
xa5 := A5.x;;
ya5 := A5.y;;

KQ2 := PathAlgebra(GF(2), Q, HighLevelGroebnerBasis);;
x2 := ElementOfPathAlgebra(KQ2, Q.x);;
y2 := ElementOfPathAlgebra(KQ2, Q.y);;
A2 := KQ2 / [x2^2, x2*y2+y2*x2, y2^2];;
xa2 := A2.x;;
ya2 := A2.y;;

RegressionQ := Quiver(3, [[1, 2, "a"], [1, 2, "b"], [2, 2, "c"],
                          [2, 3, "d"], [3, 1, "e"]]);;
RegressionKQ := PathAlgebra(Rationals, RegressionQ, HighLevelGroebnerBasis);;
ra := ElementOfPathAlgebra(RegressionKQ, RegressionQ.a);;
rb := ElementOfPathAlgebra(RegressionKQ, RegressionQ.b);;
rc := ElementOfPathAlgebra(RegressionKQ, RegressionQ.c);;
rd := ElementOfPathAlgebra(RegressionKQ, RegressionQ.d);;
re := ElementOfPathAlgebra(RegressionKQ, RegressionQ.e);;
RegressionA := RegressionKQ / [rd*re, rc^2, ra*rc*rd-rb*rd, re*ra];;

Print("{\"schema_version\":1,\"kind\":\"quotient_path_algebra\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},\"rational\":{");
Print("\"relations\":");
PrintElementList(RelationsOfAlgebra(Aq), PrintRationalElement);
Print(",\"zero\":");
PrintRationalElement(Zero(Aq)![1]);
Print(",\"one\":");
PrintRationalElement(One(Aq)![1]);
Print(",\"yx\":");
PrintRationalElement((yaq*xaq)![1]);
Print(",\"x_squared\":");
PrintRationalElement((xaq^2)![1]);
Print(",\"yxy_plus_x\":");
PrintRationalElement((yaq*xaq*yaq+xaq)![1]);
Print(",\"one_plus_x_times_one_minus_x\":");
PrintRationalElement(((One(Aq)+xaq)*(One(Aq)-xaq))![1]);
Print(",\"finite_dimensional\":", IsFiniteDimensional(Aq));
Print(",\"dimension\":", Dimension(Aq));
Print(",\"canonical_basis_paths\":[");
for i in [1 .. Length(BasisVectors(CanonicalBasis(Aq)))] do
    if i > 1 then
        Print(",");
    fi;
    PrintJsonString(StringPrint(LeadingMonomial(
        BasisVectors(CanonicalBasis(Aq))[i]![1])));
od;
Print("]");
Print(",\"ideal_membership\":{");
Print("\"zero\":", Zero(KQ) in Iq);
Print(",\"x\":", xq in Iq);
Print(",\"x_cubed\":", xq^3 in Iq, "}");
Print("},\"f5\":{");
Print("\"yx\":");
PrintPrimeElement((ya5*xa5)![1]);
Print(",\"x_squared\":");
PrintPrimeElement((xa5^2)![1]);
Print("},\"f2\":{");
Print("\"yx\":");
PrintPrimeElement((ya2*xa2)![1]);
Print(",\"x_squared\":");
PrintPrimeElement((xa2^2)![1]);
Print("},\"qpa_regression_dimension_17\":{");
Print("\"finite_dimensional\":", IsFiniteDimensional(RegressionA));
Print(",\"dimension\":", Dimension(RegressionA));
Print(",\"canonical_basis_paths\":[");
for i in [1 .. Length(BasisVectors(CanonicalBasis(RegressionA)))] do
    if i > 1 then
        Print(",");
    fi;
    PrintJsonString(StringPrint(LeadingMonomial(
        BasisVectors(CanonicalBasis(RegressionA))[i]![1])));
od;
Print("]}}\n");

QUIT_GAP(0);
