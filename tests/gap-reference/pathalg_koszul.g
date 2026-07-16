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

PrintList := function(elements, printer)
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

empty_quadratic := IsQuadraticIdeal([]);;

Q := Quiver(1, [[1, 1, "x"], [1, 1, "y"]]);;
KQ := PathAlgebra(Rationals, Q);;
x := ElementOfPathAlgebra(KQ, Q.x);;
y := ElementOfPathAlgebra(KQ, Q.y);;
relations := [x^2+y^2, x*y-2*y*x];;
perpendicular := QuadraticPerpOfPathAlgebraIdeal(relations);;
zero_perpendicular := QuadraticPerpOfPathAlgebraIdeal([Zero(KQ)]);;
full_perpendicular := QuadraticPerpOfPathAlgebraIdeal([x^2, x*y, y*x, y^2]);;

Qm := Quiver(2, [[1, 2, "a"], [2, 1, "b"], [1, 1, "c"]]);;

KQ2 := PathAlgebra(GF(2), Q);;
x2 := ElementOfPathAlgebra(KQ2, Q.x);;
y2 := ElementOfPathAlgebra(KQ2, Q.y);;
perpendicular2 := QuadraticPerpOfPathAlgebraIdeal([x2^2+y2^2, x2*y2]);;

KQ5 := PathAlgebra(GF(5), Q);;
x5 := ElementOfPathAlgebra(KQ5, Q.x);;
y5 := ElementOfPathAlgebra(KQ5, Q.y);;
perpendicular5 := QuadraticPerpOfPathAlgebraIdeal([x5^2+y5^2, x5*y5-2*y5*x5]);;

Print("{\"schema_version\":1,\"kind\":\"pathalg_koszul\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},\"paths_loop\":");
PrintList(List(PathsOfLengthTwo(Q), StringPrint), PrintJsonString);
Print(",\"paths_multi\":");
PrintList(List(PathsOfLengthTwo(Qm), StringPrint), PrintJsonString);
Print(",\"quadratic\":{");
Print("\"empty\":", empty_quadratic);
Print(",\"relations\":", IsQuadraticIdeal(relations));
Print(",\"zero\":", IsQuadraticIdeal([Zero(KQ)]));
Print(",\"linear\":", IsQuadraticIdeal([x]));
Print("},\"opposite_quiver\":{");
Print("\"vertices\":");
PrintList(List(VerticesOfQuiver(QuiverOfPathAlgebra(perpendicular[1])), StringPrint), PrintJsonString);
Print(",\"arrows\":");
PrintList(List(ArrowsOfQuiver(QuiverOfPathAlgebra(perpendicular[1])), StringPrint), PrintJsonString);
Print("},\"rational_perpendicular\":");
PrintList(perpendicular[2], PrintRationalElement);
Print(",\"zero_perpendicular\":");
PrintList(zero_perpendicular[2], PrintRationalElement);
Print(",\"full_perpendicular\":");
PrintList(full_perpendicular[2], PrintRationalElement);
Print(",\"f2_perpendicular\":");
PrintList(perpendicular2[2], PrintPrimeElement);
Print(",\"f5_perpendicular\":");
PrintList(perpendicular5[2], PrintPrimeElement);
Print("}\n");

QUIT_GAP(0);
