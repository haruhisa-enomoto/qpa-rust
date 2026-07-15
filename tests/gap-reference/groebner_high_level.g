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
basisq := HighLevelGroebnerBasis(relationsq, KQ);;

KQ5 := PathAlgebra(GF(5), Q, HighLevelGroebnerBasis);;
x5 := ElementOfPathAlgebra(KQ5, Q.x);;
y5 := ElementOfPathAlgebra(KQ5, Q.y);;
relations5 := [x5^2, x5*y5-2*y5*x5, y5^2];;
basis5 := HighLevelGroebnerBasis(relations5, KQ5);;

KQ2 := PathAlgebra(GF(2), Q, HighLevelGroebnerBasis);;
x2 := ElementOfPathAlgebra(KQ2, Q.x);;
y2 := ElementOfPathAlgebra(KQ2, Q.y);;
relations2 := [x2^2, x2*y2+y2*x2, y2^2];;
basis2 := HighLevelGroebnerBasis(relations2, KQ2);;

UQ := Quiver(2, [[1, 1, "a"], [1, 2, "b"], [2, 1, "c"]]);;
UKQ := PathAlgebra(Rationals, UQ);;
au := ElementOfPathAlgebra(UKQ, UQ.a);;
bu := ElementOfPathAlgebra(UKQ, UQ.b);;
cu := ElementOfPathAlgebra(UKQ, UQ.c);;

Print("{\"schema_version\":1,\"kind\":\"groebner_high_level\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},\"rational\":{");
Print("\"basis\":");
PrintElementList(basisq, PrintRationalElement);
Print(",\"remainder_yx\":");
PrintRationalElement(RemainderOfDivision(yq*xq, basisq, KQ));
Print(",\"remainder_yxy_plus_x\":");
PrintRationalElement(RemainderOfDivision(yq*xq*yq+xq, basisq, KQ));
Print("},\"f5\":{");
Print("\"basis\":");
PrintElementList(basis5, PrintPrimeElement);
Print(",\"remainder_yx\":");
PrintPrimeElement(RemainderOfDivision(y5*x5, basis5, KQ5));
Print("},\"f2\":{");
Print("\"basis\":");
PrintElementList(basis2, PrintPrimeElement);
Print(",\"remainder_yx\":");
PrintPrimeElement(RemainderOfDivision(y2*x2, basis2, KQ2));
Print("},\"uniform_components\":");
PrintElementList(MakeUniform(au+bu+cu), PrintRationalElement);
Print("}\n");

QUIT_GAP(0);
