if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

SizeScreen([100000, 100000]);

PrintFiniteElement := function(field, value)
    if Size(field) = 2 then
        if value = Zero(field) then
            Print("0");
        else
            Print("1");
        fi;
    else
        Print(IntFFE(value));
    fi;
end;

PrintUnitElements := function(algebra, elements)
    local basis, element, coefficients, i;
    basis := CanonicalBasis(algebra);
    Print("[");
    for element in elements do
        if element <> elements[1] then
            Print(",");
        fi;
        coefficients := Coefficients(basis, element);
        Print("[");
        for i in [1 .. Length(coefficients)] do
            if i > 1 then
                Print(",");
            fi;
            PrintFiniteElement(LeftActingDomain(algebra), coefficients[i]);
        od;
        Print("]");
    od;
    Print("]");
end;

PrintTable := function(elements)
    local i, j;
    Print("[");
    for i in [1 .. Length(elements)] do
        if i > 1 then
            Print(",");
        fi;
        Print("[");
        for j in [1 .. Length(elements)] do
            if j > 1 then
                Print(",");
            fi;
            Print(Position(elements, elements[i] * elements[j]));
        od;
        Print("]");
    od;
    Print("]");
end;

q2 := Quiver(2, [[1, 2, "a"]]);;
a2 := PathAlgebra(GF(2), q2);;
a2_had_elements := HasElementsOfUnits(a2);;
a2_group := Units(a2);;
a2_units := ElementsOfUnits(a2);;

qd := Quiver(1, [[1, 1, "x"]]);;
pd := PathAlgebra(GF(3), qd);;
dual := pd / [pd.x^2];;
dual_group := Units(dual);;
dual_units := ElementsOfUnits(dual);;

pc := PathAlgebra(GF(2), qd);;
cubic := pc / [pc.x^3];;
cubic_group := Units(cubic);;
cubic_units := ElementsOfUnits(cubic);;

qs := Quiver(["u", "v"], []);;
semisimple := PathAlgebra(GF(3), qs);;
semisimple_group := Units(semisimple);;
semisimple_units := ElementsOfUnits(semisimple);;

q5 := Quiver(1, []);;
a1_5 := PathAlgebra(GF(5), q5);;
a1_5_group := Units(a1_5);;
a1_5_units := ElementsOfUnits(a1_5);;

loop := PathAlgebra(GF(2), Quiver(1, [[1, 1, "z"]]));;
infinite_result := CALL_WITH_CATCH(Units, [loop]);;
rational := PathAlgebra(Rationals, Quiver(1, []));;
rational_result := CALL_WITH_CATCH(Units, [rational]);;

Print("{\"schema_version\":1,\"kind\":\"pathalg_units\",\"oracle\":{");
Print("\"gap_version\":\""); Print(GAPInfo.Version); Print("\"");
Print(",\"gap_revision\":\"");
Print(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION); Print("\"");
Print(",\"qpa_version\":\""); Print(PackageInfo("qpa")[1].Version); Print("\"");
Print(",\"qpa_revision\":\"");
Print(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION); Print("\"");
Print(",\"gbnp_version\":\""); Print(PackageInfo("gbnp")[1].Version); Print("\"");
Print(",\"gbnp_revision\":\"");
Print(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION); Print("\"}");

Print(",\"binary_a2\":{\"had_elements\":", a2_had_elements);
Print(",\"has_elements\":", HasElementsOfUnits(a2));
Print(",\"is_group\":", IsGroup(a2_group));
Print(",\"size\":", Size(a2_group), ",\"elements\":");
PrintUnitElements(a2, a2_units);
Print(",\"table\":"); PrintTable(a2_units); Print("}");

Print(",\"ternary_dual_numbers\":{\"is_group\":", IsGroup(dual_group));
Print(",\"size\":", Size(dual_group), ",\"elements\":");
PrintUnitElements(dual, dual_units);
Print(",\"table\":"); PrintTable(dual_units); Print("}");

Print(",\"binary_cubic_numbers\":{\"is_group\":", IsGroup(cubic_group));
Print(",\"size\":", Size(cubic_group), ",\"elements\":");
PrintUnitElements(cubic, cubic_units);
Print(",\"table\":"); PrintTable(cubic_units); Print("}");

Print(",\"ternary_semisimple\":{\"is_group\":", IsGroup(semisimple_group));
Print(",\"size\":", Size(semisimple_group), ",\"elements\":");
PrintUnitElements(semisimple, semisimple_units);
Print(",\"table\":"); PrintTable(semisimple_units); Print("}");

Print(",\"quinary_a1\":{\"is_group\":", IsGroup(a1_5_group));
Print(",\"size\":", Size(a1_5_group), ",\"elements\":");
PrintUnitElements(a1_5, a1_5_units);
Print(",\"table\":"); PrintTable(a1_5_units); Print("}");

Print(",\"errors\":{\"infinite_success\":", infinite_result[1]);
Print(",\"rational_success\":", rational_result[1], "}}\n");

QUIT_GAP(0);
