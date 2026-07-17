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

PrintValue := function(value)
    if value = infinity then
        Print("\"infinity\"");
    else
        Print(value);
    fi;
end;

PrintCheck := function(module)
    local before, first, after_first, second, after_second;
    before := HasFaithfulDimension(module);
    first := FaithfulDimension(module);
    after_first := HasFaithfulDimension(module);
    second := FaithfulDimension(module);
    after_second := HasFaithfulDimension(module);
    Print("{\"has_before\":", before, ",\"first\":"); PrintValue(first);
    Print(",\"has_after_first\":", after_first, ",\"second\":"); PrintValue(second);
    Print(",\"has_after_second\":", after_second, "}");
end;

PrintCase := function(field)
    local a2, projectives, regular, simples, pd_one, zero,
          loop_parent, dual, dual_regular, dual_simple;
    a2 := PathAlgebra(field,
        Quiver(["v1", "v2"], [["v1", "v2", "a"]]));
    projectives := IndecProjectiveModules(a2);
    regular := DirectSumOfQPAModules(projectives);
    simples := SimpleModules(a2);
    pd_one := DirectSumOfQPAModules([projectives[1], simples[1]]);
    zero := ZeroModule(a2);

    loop_parent := PathAlgebra(field,
        Quiver(["w"], [["w", "w", "x"]]));
    dual := loop_parent / [loop_parent.x^2];
    dual_regular := IndecProjectiveModules(dual)[1];
    dual_simple := SimpleModules(dual)[1];

    Print("{\"regular\":"); PrintCheck(regular);
    Print(",\n\"projective_1\":"); PrintCheck(projectives[1]);
    Print(",\n\"simple_1\":"); PrintCheck(simples[1]);
    Print(",\n\"simple_2\":"); PrintCheck(simples[2]);
    Print(",\n\"pd_one\":"); PrintCheck(pd_one);
    Print(",\n\"zero\":"); PrintCheck(zero);
    Print(",\n\"dual_regular\":"); PrintCheck(dual_regular);
    Print(",\n\"dual_simple\":"); PrintCheck(dual_simple);
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"modulehomalg_faithful_dimension\",\"oracle\":{");
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
