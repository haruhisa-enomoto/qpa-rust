if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

SizeScreen([100000, 100000]);

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintIntegerList := function(values)
    Print("[", JoinStringsWithSeparator(List(values, String), ","), "]");
end;

PrintResult := function(result)
    local first, module;
    Print("{\"projective_dimension\":");
    first := result[1];
    if IsInt(first) then
        Print(first);
    else
        PrintJsonString(first);
    fi;
    Print(",\"syzygy_dimension_vectors\":[");
    first := true;
    for module in result[2] do
        if not first then Print(","); fi;
        PrintIntegerList(DimensionVector(module));
        first := false;
    od;
    Print("]}");
end;

PrintCase := function(field)
    local a2, simples, projectives,
          loop_parent, dual, dual_simple, dual_regular;
    a2 := PathAlgebra(field,
        Quiver(["v1", "v2"], [["v1", "v2", "a"]]));
    simples := SimpleModules(a2);
    projectives := IndecProjectiveModules(a2);

    loop_parent := PathAlgebra(field,
        Quiver(["w"], [["w", "w", "x"]]));
    dual := loop_parent / [loop_parent.x^2];
    dual_simple := SimpleModules(dual)[1];
    dual_regular := IndecProjectiveModules(dual)[1];

    Print("{\"single_s1_length_3\":");
    PrintResult(ProjectiveResolutionOfSimpleModuleOverEndo([simples[1]], 1, 3));
    Print(",\n\"projectives_at_one_length_3\":");
    PrintResult(ProjectiveResolutionOfSimpleModuleOverEndo(projectives, 1, 3));
    Print(",\n\"dual_regular_simple_at_two_length_3\":");
    PrintResult(ProjectiveResolutionOfSimpleModuleOverEndo(
        [dual_regular, dual_simple], 2, 3));
    Print(",\n\"dual_regular_simple_at_two_length_0\":");
    PrintResult(ProjectiveResolutionOfSimpleModuleOverEndo(
        [dual_regular, dual_simple], 2, 0));
    Print(",\n\"single_dual_regular_length_0\":");
    PrintResult(ProjectiveResolutionOfSimpleModuleOverEndo([dual_regular], 1, 0));
    Print(",\n\"single_dual_regular_length_2\":");
    PrintResult(ProjectiveResolutionOfSimpleModuleOverEndo([dual_regular], 1, 2));
    Print(",\n\"single_dual_regular_length_3\":");
    PrintResult(ProjectiveResolutionOfSimpleModuleOverEndo([dual_regular], 1, 3));
    Print(",\n\"single_dual_regular_length_4\":");
    PrintResult(ProjectiveResolutionOfSimpleModuleOverEndo([dual_regular], 1, 4));
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"modulehomalg_projective_resolution_simple_endo\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\n\"binary\":"); PrintCase(GF(2));
Print(",\n\"quinary\":"); PrintCase(GF(5));
Print(",\n\"rational\":"); PrintCase(Rationals);
Print("}\n");

QUIT_GAP(0);
