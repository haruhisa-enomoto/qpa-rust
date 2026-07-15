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

Q := Quiver(2, [[1, 1, "x"], [2, 2, "y"], [1, 2, "a"]]);;
A := TruncatedPathAlgebra(Rationals, Q, 2);;
unit := 2*A.v1 + 3*A.v2 + A.x + A.y + A.a;;
nonunit := A.v1 + A.a;;
inverse := Inverse(unit);;
negative_two := unit^(-2);;
nonunit_inverse := Inverse(nonunit);;
nonunit_negative := nonunit^(-1);;
zero_power := nonunit^0;;

A2 := TruncatedPathAlgebra(GF(2), Q, 2);;
unit2 := A2.v1 + A2.v2 + A2.x + A2.y + A2.a;;
nonunit2 := A2.x;;
inverse2 := Inverse(unit2);;
negative_two2 := unit2^(-2);;
zero_power2 := nonunit2^0;;

Print("{\"schema_version\":1,\"kind\":\"pathalg_inverse\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},\"rational\":{\"inverse\":");
PrintRationalElement(inverse);
Print(",\"product\":");
PrintRationalElement(unit*inverse);
Print(",\"negative_two\":");
PrintRationalElement(negative_two);
Print(",\"nonunit_inverse_fail\":", nonunit_inverse = fail);
Print(",\"nonunit_negative_fail\":", nonunit_negative = fail);
Print(",\"zero_power\":");
PrintRationalElement(zero_power);
Print("},\"f2\":{\"inverse\":");
PrintPrimeElement(inverse2);
Print(",\"product\":");
PrintPrimeElement(unit2*inverse2);
Print(",\"negative_two\":");
PrintPrimeElement(negative_two2);
Print(",\"nonunit_inverse_fail\":", Inverse(nonunit2) = fail);
Print(",\"zero_power\":");
PrintPrimeElement(zero_power2);
Print("}}");
Print("\n");

QUIT_GAP(0);
