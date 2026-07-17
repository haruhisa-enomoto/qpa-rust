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

PrintCoresolutionShape := function(complex)
    local degree, first, differential;
    Print("{"); first := true;
    for degree in [-2..0] do
        if not first then Print(","); fi;
        differential := DifferentialOfComplex(complex, degree);
        Print("\"", degree, "\":{");
        Print("\"source\":"); PrintIntegerList(DimensionVector(Source(differential)));
        Print(",\"range\":"); PrintIntegerList(DimensionVector(Range(differential)));
        Print("}"); first := false;
    od;
    Print("}");
end;

PrintTiltingValue := function(value)
    local first, coresolution;
    if IsBool(value) and value = false then
        Print("false");
        return;
    fi;
    Print("{\"projective_dimension\":", value[1]);
    Print(",\"coresolution_count\":", Length(value[2]));
    Print(",\"coresolutions\":["); first := true;
    for coresolution in value[2] do
        if not first then Print(","); fi;
        PrintCoresolutionShape(coresolution); first := false;
    od;
    Print("]}");
end;

PrintCheck := function(module, bound)
    local before, value, after;
    before := HasIsTiltingModule(module);
    value := TiltingModule(module, bound);
    after := HasIsTiltingModule(module);
    Print("{\"has_before\":", before, ",\"value\":");
    PrintTiltingValue(value);
    Print(",\"has_after\":", after, "}");
end;

PrintCase := function(field)
    local a2, projectives, regular, pd_one, simple,
          loop_parent, dual, dual_regular, dual_simple;
    a2 := PathAlgebra(field,
        Quiver(["v1", "v2"], [["v1", "v2", "a"]]));
    projectives := IndecProjectiveModules(a2);
    regular := DirectSumOfQPAModules(projectives);
    simple := SimpleModules(a2)[1];
    pd_one := DirectSumOfQPAModules([projectives[1], simple]);

    loop_parent := PathAlgebra(field,
        Quiver(["w"], [["w", "w", "x"]]));
    dual := loop_parent / [loop_parent.x^2];
    dual_regular := DirectSumOfQPAModules(IndecProjectiveModules(dual));
    dual_simple := SimpleModules(dual)[1];

    Print("{\"regular_bound_0\":"); PrintCheck(regular, 0);
    Print(",\n\"regular_cached_negative\":"); PrintCheck(regular, -1);
    Print(",\n\"pd_one_bound_0\":"); PrintCheck(pd_one, 0);
    Print(",\n\"pd_one_bound_1\":"); PrintCheck(pd_one, 1);
    Print(",\n\"simple_bound_1\":"); PrintCheck(simple, 1);
    Print(",\n\"dual_regular_bound_0\":"); PrintCheck(dual_regular, 0);
    Print(",\n\"dual_simple_bound_1\":"); PrintCheck(dual_simple, 1);
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"modulehomalg_tilting\",\"oracle\":{");
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
