if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintGlobalDimension := function(value)
    if value = infinity then
        PrintJsonString("infinity");
    else
        Print(value);
    fi;
end;

PrintDimensionList := function(values)
    local position;
    Print("[");
    for position in [1..Length(values)] do
        if position > 1 then
            Print(",");
        fi;
        PrintGlobalDimension(values[position]);
    od;
    Print("]");
end;

GlobalDimensionCases := function(field)
    local discrete_quiver, discrete, a2_quiver, a2, a3_quiver, parent,
          quotient, loop_quiver, loop_parent, dual_numbers;
    discrete_quiver := Quiver(["v"], []);
    discrete := PathAlgebra(field, discrete_quiver);

    a2_quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
    a2 := PathAlgebra(field, a2_quiver);

    a3_quiver := Quiver(["u1", "u2", "u3"],
        [["u1", "u2", "x"], ["u2", "u3", "y"]]);
    parent := PathAlgebra(field, a3_quiver);
    quotient := parent / [parent.x * parent.y];

    loop_quiver := Quiver(["w"], [["w", "w", "t"]]);
    loop_parent := PathAlgebra(field, loop_quiver);
    dual_numbers := loop_parent / [loop_parent.t * loop_parent.t];

    Print("{\"discrete\":");
    PrintDimensionList([GlobalDimensionOfAlgebra(discrete, 0)]);
    Print(",\"path_a2\":");
    PrintDimensionList([GlobalDimensionOfAlgebra(a2, 0),
                        GlobalDimensionOfAlgebra(a2, 3)]);
    Print(",\"radical_square_zero_a3\":");
    PrintDimensionList([GlobalDimensionOfAlgebra(quotient, 0),
                        GlobalDimensionOfAlgebra(quotient, 1),
                        GlobalDimensionOfAlgebra(quotient, 2),
                        GlobalDimensionOfAlgebra(quotient, 0)]);
    Print(",\"dual_numbers\":");
    PrintDimensionList([GlobalDimensionOfAlgebra(dual_numbers, 0)]);
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"modulehom_global_dimension\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":");
PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"f2\":");
GlobalDimensionCases(GF(2));
Print("}\n");

QUIT_GAP(0);
