if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintResult := function(value)
    if value = false then
        Print("false");
    elif value = infinity then
        Print("\"infinity\"");
    else
        Print(value);
    fi;
end;

PrintCase := function(field)
    local a2, simple_one, simple_two, projective, zero,
          loop_parent, dual, dual_simple, discrete;
    a2 := PathAlgebra(field,
        Quiver(["v1", "v2"], [["v1", "v2", "a"]]));
    simple_one := RightModuleOverPathAlgebra(a2, [1, 0], []);
    simple_two := RightModuleOverPathAlgebra(a2, [0, 1], []);
    projective := RightModuleOverPathAlgebra(a2, [1, 1],
        [["a", [[1]] * One(field)]]);
    zero := ZeroModule(a2);

    loop_parent := PathAlgebra(field,
        Quiver(["w"], [["w", "w", "x"]]));
    dual := loop_parent / [loop_parent.x^2];
    dual_simple := RightModuleOverPathAlgebra(dual, [1],
        [["x", [[0]] * One(field)]]);
    discrete := PathAlgebra(field, Quiver(["d"], []));

    Print("{\"simple_one_bound_0\":");
    PrintResult(DominantDimensionOfModule(simple_one, 0));
    Print(",\"simple_one_bound_2\":");
    PrintResult(DominantDimensionOfModule(simple_one, 2));
    Print(",\"simple_one_bound_negative\":");
    PrintResult(DominantDimensionOfModule(simple_one, -1));
    Print(",\"simple_two_bound_0\":");
    PrintResult(DominantDimensionOfModule(simple_two, 0));
    Print(",\"simple_two_bound_1\":");
    PrintResult(DominantDimensionOfModule(simple_two, 1));
    Print(",\"simple_two_bound_2\":");
    PrintResult(DominantDimensionOfModule(simple_two, 2));
    Print(",\"simple_two_bound_negative\":");
    PrintResult(DominantDimensionOfModule(simple_two, -1));
    Print(",\"projective_injective\":");
    PrintResult(DominantDimensionOfModule(projective, 0));
    Print(",\"zero\":");
    PrintResult(DominantDimensionOfModule(zero, 0));
    Print(",\"dual_simple_bound_2\":");
    PrintResult(DominantDimensionOfModule(dual_simple, 2));
    Print(",\"a2_bound_0\":");
    PrintResult(DominantDimensionOfAlgebra(a2, 0));
    Print(",\"a2_bound_1\":");
    PrintResult(DominantDimensionOfAlgebra(a2, 1));
    Print(",\"a2_bound_2\":");
    PrintResult(DominantDimensionOfAlgebra(a2, 2));
    Print(",\"dual_algebra\":");
    PrintResult(DominantDimensionOfAlgebra(dual, 2));
    Print(",\"discrete_algebra\":");
    PrintResult(DominantDimensionOfAlgebra(discrete, 0));
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"modulehomalg_dominant_dimension\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"binary\":"); PrintCase(GF(2));
Print(",\"quinary\":"); PrintCase(GF(5));
Print(",\"rational\":"); PrintCase(Rationals);
Print("}\n");

QUIT_GAP(0);
