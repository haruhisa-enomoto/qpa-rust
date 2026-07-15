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

PrintRationalRows := function(rows)
    local row, column;
    Print("[");
    for row in [1..Length(rows)] do
        if row > 1 then Print(","); fi;
        Print("[");
        for column in [1..Length(rows[row])] do
            if column > 1 then Print(","); fi;
            PrintRationalElement(rows[row][column]);
        od;
        Print("]");
    od;
    Print("]");
end;

PrintRationalVector := function(vector)
    local index;
    Print("[");
    for index in [1..Length(vector)] do
        if index > 1 then Print(","); fi;
        Print("[", NumeratorRat(vector[index]), ",", DenominatorRat(vector[index]), "]");
    od;
    Print("]");
end;

PrintPrimeRows := function(rows)
    local row, column, data, index;
    Print("[");
    for row in [1..Length(rows)] do
        if row > 1 then Print(","); fi;
        Print("[");
        for column in [1..Length(rows[row])] do
            if column > 1 then Print(","); fi;
            data := CoefficientsAndMagmaElements(rows[row][column]);
            Print("[");
            for index in [1, 3 .. Length(data) - 1] do
                if index > 1 then Print(","); fi;
                Print("{\"path\":", ExtRepOfObj(data[index]));
                Print(",\"coefficient\":", IntFFE(data[index + 1]), "}");
            od;
            Print("]");
        od;
        Print("]");
    od;
    Print("]");
end;

Q := Quiver(2, [[1, 2, "a"]]);;
A := PathAlgebra(Rationals, Q);;
tip := TipMonomialandCoefficientOfVector(A, [A.v1+A.a, 2*A.a]);;
zero_tip := TipMonomialandCoefficientOfVector(A, [Zero(A), Zero(A)]);;
rows := [[A.a, Zero(A)], [2*A.a, A.v2]];;
same_rows := IsIdenticalObj(rows, TipReduceVectors(A, rows));;
basis := [[A.a, Zero(A)], [Zero(A), A.v2]];;
target := [3*A.a, 5*A.v2];;
coefficients := CoefficientsOfVectors(A, basis, target);;
not_span := [A.v1, Zero(A)];;
failure := CoefficientsOfVectors(A, basis, not_span);;

AQ := A/[];;
quotient_basis := [[AQ.a, Zero(AQ)], [Zero(AQ), AQ.v2]];;
quotient_target := [3*AQ.a, 5*AQ.v2];;
quotient_coefficients := CoefficientsOfVectors(AQ, quotient_basis, quotient_target);;

A2 := PathAlgebra(GF(2), Q);;
rows2 := [[A2.a, Zero(A2)], [A2.a, A2.v2]];;
TipReduceVectors(A2, rows2);;
basis2 := [[A2.a, Zero(A2)], [Zero(A2), A2.v2]];;
coefficients2 := CoefficientsOfVectors(A2, basis2, [A2.a, A2.v2]);;

Print("{\"schema_version\":1,\"kind\":\"pathalg_vector_reduction\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},\"rational\":{\"tip\":{\"position\":", tip[1]);
Print(",\"path\":", ExtRepOfObj(tip[2]));
Print(",\"coefficient\":[", NumeratorRat(tip[3]), ",", DenominatorRat(tip[3]), "]}");
Print(",\"zero_tip\":{\"position\":", zero_tip[1]);
Print(",\"path\":", ExtRepOfObj(zero_tip[2]));
Print(",\"coefficient\":[", NumeratorRat(zero_tip[3]), ",", DenominatorRat(zero_tip[3]), "]}");
Print(",\"same_rows\":", same_rows);
Print(",\"reduced\":");
PrintRationalRows(rows);
Print(",\"coefficients\":");
PrintRationalVector(coefficients);
Print(",\"failure_witness\":", failure = [basis, not_span]);
Print(",\"quotient_coefficients\":");
PrintRationalVector(quotient_coefficients);
Print("},\"f2\":{\"reduced\":");
PrintPrimeRows(rows2);
Print(",\"coefficients\":[");
for index in [1..Length(coefficients2)] do
    if index > 1 then Print(","); fi;
    Print(IntFFE(coefficients2[index]));
od;
Print("]}}");
Print("\n");

QUIT_GAP(0);
