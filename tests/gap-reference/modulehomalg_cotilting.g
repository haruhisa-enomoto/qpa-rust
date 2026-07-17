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

PrintResolutionShape := function(resolution)
    local degree, first, differential;
    if IsBool(resolution) then
        Print(resolution);
        return;
    fi;
    Print("{"); first := true;
    for degree in [0..3] do
        if not first then Print(","); fi;
        differential := DifferentialOfComplex(resolution, degree);
        Print("\"", degree, "\":{");
        Print("\"source\":"); PrintIntegerList(DimensionVector(Source(differential)));
        Print(",\"range\":"); PrintIntegerList(DimensionVector(Range(differential)));
        Print("}"); first := false;
    od;
    Print("}");
end;

PrintCotiltingValue := function(value)
    local first, resolution;
    if IsBool(value) and value = false then
        Print("false");
        return;
    fi;
    Print("{\"injective_dimension\":", value[1]);
    Print(",\"resolution_count\":", Length(value[2]));
    Print(",\"resolutions\":["); first := true;
    for resolution in value[2] do
        if not first then Print(","); fi;
        PrintResolutionShape(resolution); first := false;
    od;
    Print("]}");
end;

PrintCheck := function(module, bound)
    local before, value, after, cached_value;
    before := HasIsCotiltingModule(module);
    value := CotiltingModule(module, bound);
    after := HasIsCotiltingModule(module);
    Print("{\"has_before\":", before, ",\"value\":");
    PrintCotiltingValue(value);
    Print(",\"has_after\":", after);
    if after then
        cached_value := IsCotiltingModule(module);
        Print(",\"cached_value\":", cached_value);
    fi;
    Print("}");
end;

PrintCase := function(field)
    local a2, injectives, cogenerator, simples, id_one, semisimple,
          loop_parent, dual, dual_regular, dual_simple;
    a2 := PathAlgebra(field,
        Quiver(["v1", "v2"], [["v1", "v2", "a"]]));
    injectives := IndecInjectiveModules(a2);
    cogenerator := DirectSumOfQPAModules(injectives);
    simples := SimpleModules(a2);
    id_one := DirectSumOfQPAModules([injectives[2], simples[2]]);
    semisimple := DirectSumOfQPAModules(simples);

    loop_parent := PathAlgebra(field,
        Quiver(["w"], [["w", "w", "x"]]));
    dual := loop_parent / [loop_parent.x^2];
    dual_regular := DirectSumOfQPAModules(IndecInjectiveModules(dual));
    dual_simple := SimpleModules(dual)[1];

    Print("{\"cogenerator_bound_0\":"); PrintCheck(cogenerator, 0);
    Print(",\n\"cogenerator_cached_negative\":"); PrintCheck(cogenerator, -1);
    Print(",\n\"id_one_bound_0\":"); PrintCheck(id_one, 0);
    Print(",\n\"id_one_bound_1\":"); PrintCheck(id_one, 1);
    Print(",\n\"semisimple_bound_1\":"); PrintCheck(semisimple, 1);
    Print(",\n\"dual_regular_bound_0\":"); PrintCheck(dual_regular, 0);
    Print(",\n\"dual_simple_bound_1\":"); PrintCheck(dual_simple, 1);
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"modulehomalg_cotilting\",\"oracle\":{");
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
