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

PrintCase := function(field)
    local a2, simples, projectives, zero,
          loop_parent, dual, dual_simple, dual_regular;
    a2 := PathAlgebra(field,
        Quiver(["v1", "v2"], [["v1", "v2", "a"]]));
    simples := SimpleModules(a2);
    projectives := IndecProjectiveModules(a2);
    zero := ZeroModule(a2);

    loop_parent := PathAlgebra(field,
        Quiver(["w"], [["w", "w", "x"]]));
    dual := loop_parent / [loop_parent.x^2];
    dual_simple := SimpleModules(dual)[1];
    dual_regular := IndecProjectiveModules(dual)[1];

    Print("{\"zero_0\":", IsNthSyzygy(zero, 0));
    Print(",\"zero_1\":", IsNthSyzygy(zero, 1));
    Print(",\"simple_one_0\":", IsNthSyzygy(simples[1], 0));
    Print(",\"simple_one_1\":", IsNthSyzygy(simples[1], 1));
    Print(",\"simple_one_2\":", IsNthSyzygy(simples[1], 2));
    Print(",\"simple_two_1\":", IsNthSyzygy(simples[2], 1));
    Print(",\"projective_one_1\":", IsNthSyzygy(projectives[1], 1));
    Print(",\"projective_two_2\":", IsNthSyzygy(projectives[2], 2));
    Print(",\"dual_simple_1\":", IsNthSyzygy(dual_simple, 1));
    Print(",\"dual_simple_2\":", IsNthSyzygy(dual_simple, 2));
    Print(",\"dual_regular_1\":", IsNthSyzygy(dual_regular, 1));
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"modulehomalg_is_nth_syzygy\",\"oracle\":{");
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
