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
    else
        Print(value);
    fi;
end;

FreshA2Module := function(field, dimensions)
    local algebra;
    algebra := PathAlgebra(field,
        Quiver(["v1", "v2"], [["v1", "v2", "a"]]));
    return RightModuleOverPathAlgebra(algebra, dimensions, []);
end;

FreshA2Zero := function(field)
    local algebra;
    algebra := PathAlgebra(field,
        Quiver(["z1", "z2"], [["z1", "z2", "z"]]));
    return ZeroModule(algebra);
end;

FreshDualSimple := function(field)
    local parent, algebra;
    parent := PathAlgebra(field,
        Quiver(["w"], [["w", "w", "x"]]));
    algebra := parent / [parent.x^2];
    return RightModuleOverPathAlgebra(algebra, [1],
        [["x", [[0]] * One(field)]]);
end;

PrintCase := function(field)
    local pd_one, pd_two, pd_dual, pd_negative, pd_zero,
          id_two, id_one, id_dual, id_negative, id_zero,
          shortcut, shortcut_dual, value;

    pd_one := FreshA2Module(field, [1, 0]);
    Print("{\"pd_one_initially_bound\":", HasProjDimension(pd_one));
    Print(",\"pd_one_bound_0\":");
    PrintResult(ProjDimensionOfModule(pd_one, 0));
    Print(",\"pd_one_bound_after_0\":", HasProjDimension(pd_one));
    Print(",\"pd_one_bound_1\":");
    PrintResult(ProjDimensionOfModule(pd_one, 1));
    Print(",\"pd_one_attribute\":", ProjDimension(pd_one));
    Print(",\"pd_one_cached_bound_0\":");
    PrintResult(ProjDimensionOfModule(pd_one, 0));
    Print(",\"pd_one_cached_bound_2\":");
    PrintResult(ProjDimensionOfModule(pd_one, 2));

    pd_two := FreshA2Module(field, [0, 1]);
    Print(",\"pd_two_bound_0\":");
    PrintResult(ProjDimensionOfModule(pd_two, 0));
    pd_dual := FreshDualSimple(field);
    Print(",\"pd_dual_bound_2\":");
    PrintResult(ProjDimensionOfModule(pd_dual, 2));
    pd_negative := FreshA2Module(field, [0, 1]);
    Print(",\"pd_negative_bound\":");
    PrintResult(ProjDimensionOfModule(pd_negative, -1));
    Print(",\"pd_negative_bound_after\":", HasProjDimension(pd_negative));
    pd_zero := FreshA2Zero(field);
    Print(",\"pd_zero_bound_0\":");
    PrintResult(ProjDimensionOfModule(pd_zero, 0));

    id_two := FreshA2Module(field, [0, 1]);
    Print(",\"id_two_initially_bound\":", HasInjDimension(id_two));
    Print(",\"id_two_bound_0\":");
    PrintResult(InjDimensionOfModule(id_two, 0));
    Print(",\"id_two_bound_after_0\":", HasInjDimension(id_two));
    Print(",\"id_two_bound_1\":");
    PrintResult(InjDimensionOfModule(id_two, 1));
    Print(",\"id_two_attribute\":", InjDimension(id_two));
    Print(",\"id_two_cached_bound_0\":");
    PrintResult(InjDimensionOfModule(id_two, 0));
    Print(",\"id_two_cached_bound_2\":");
    PrintResult(InjDimensionOfModule(id_two, 2));

    id_one := FreshA2Module(field, [1, 0]);
    Print(",\"id_one_bound_0\":");
    PrintResult(InjDimensionOfModule(id_one, 0));
    id_dual := FreshDualSimple(field);
    Print(",\"id_dual_bound_2\":");
    PrintResult(InjDimensionOfModule(id_dual, 2));
    id_negative := FreshA2Module(field, [1, 0]);
    Print(",\"id_negative_bound\":");
    PrintResult(InjDimensionOfModule(id_negative, -1));
    Print(",\"id_negative_bound_after\":", HasInjDimension(id_negative));
    id_zero := FreshA2Zero(field);
    Print(",\"id_zero_bound_0\":");
    PrintResult(InjDimensionOfModule(id_zero, 0));

    shortcut := FreshA2Module(field, [0, 1]);
    shortcut_dual := DualOfModule(shortcut);
    value := ProjDimensionOfModule(shortcut_dual, 1);
    Print(",\"id_dual_cache_setup\":"); PrintResult(value);
    Print(",\"id_dual_cache_shortcut_bound_0\":");
    PrintResult(InjDimensionOfModule(shortcut, 0));
    Print(",\"id_dual_cache_shortcut_bound_after\":",
        HasInjDimension(shortcut));
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"modulehomalg_dimensions\",\"oracle\":{");
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
