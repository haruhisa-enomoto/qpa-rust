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
    local before, tilting, after, cotilting;
    before := HasFaithfulDimension(module);
    tilting := NumberOfComplementsOfAlmostCompleteTiltingModule(module);
    after := HasFaithfulDimension(module);
    cotilting := NumberOfComplementsOfAlmostCompleteCotiltingModule(module);
    Print("{\"has_before\":", before, ",\"tilting\":"); PrintValue(tilting);
    Print(",\"has_after\":", after, ",\"cotilting\":"); PrintValue(cotilting);
    Print("}");
end;

PrintCase := function(field)
    local a2, projectives, regular, simples, zero;
    a2 := PathAlgebra(field,
        Quiver(["v1", "v2"], [["v1", "v2", "a"]]));
    projectives := IndecProjectiveModules(a2);
    regular := DirectSumOfQPAModules(projectives);
    simples := SimpleModules(a2);
    zero := ZeroModule(a2);

    Print("{\"regular\":"); PrintCheck(regular);
    Print(",\n\"projective_1\":"); PrintCheck(projectives[1]);
    Print(",\n\"simple_1\":"); PrintCheck(simples[1]);
    Print(",\n\"zero\":"); PrintCheck(zero);
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"modulehomalg_number_complements\",\"oracle\":{");
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
