if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

# IsOmegaPeriodic prints one progress line per syzygy.  Keep the source
# computation, but route JSON through the original printer while suppressing
# those diagnostic lines so the oracle remains machine-readable.
JsonPrint := Print;
MakeReadWriteGlobal("Print");
Print := function(arg...) end;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    JsonPrint("\"", escaped, "\"");
end;

PrintPeriod := function(value)
    if value = false then
        JsonPrint("false");
    else
        JsonPrint(value);
    fi;
end;

PrintCase := function(field)
    local a2, zero, simple_one, simple_two, loop_parent, dual, dual_simple;
    a2 := PathAlgebra(field,
        Quiver(["v1", "v2"], [["v1", "v2", "a"]]));
    zero := ZeroModule(a2);
    simple_one := RightModuleOverPathAlgebra(a2, [1, 0], []);
    simple_two := RightModuleOverPathAlgebra(a2, [0, 1], []);

    loop_parent := PathAlgebra(field,
        Quiver(["w"], [["w", "w", "x"]]));
    dual := loop_parent / [loop_parent.x^2];
    dual_simple := RightModuleOverPathAlgebra(dual, [1],
        [["x", [[0]] * One(field)]]);

    JsonPrint("{\"omega_zero_n0\":"); PrintPeriod(IsOmegaPeriodic(zero, 0));
    JsonPrint(",\"omega_zero_negative\":"); PrintPeriod(IsOmegaPeriodic(zero, -2));
    JsonPrint(",\"omega_zero_n3\":"); PrintPeriod(IsOmegaPeriodic(zero, 3));
    JsonPrint(",\"omega_a2_simple_one_n4\":"); PrintPeriod(IsOmegaPeriodic(simple_one, 4));
    JsonPrint(",\"omega_a2_simple_two_n4\":"); PrintPeriod(IsOmegaPeriodic(simple_two, 4));
    JsonPrint(",\"omega_dual_simple_n1\":"); PrintPeriod(IsOmegaPeriodic(dual_simple, 1));
    JsonPrint(",\"omega_dual_simple_n4\":"); PrintPeriod(IsOmegaPeriodic(dual_simple, 4));
    JsonPrint(",\"tau_zero_n0\":"); PrintPeriod(IsTauPeriodic(zero, 0));
    JsonPrint(",\"tau_zero_negative\":"); PrintPeriod(IsTauPeriodic(zero, -2));
    JsonPrint(",\"tau_zero_n3\":"); PrintPeriod(IsTauPeriodic(zero, 3));
    JsonPrint(",\"tau_a2_simple_one_n4\":"); PrintPeriod(IsTauPeriodic(simple_one, 4));
    JsonPrint(",\"tau_a2_simple_two_n4\":"); PrintPeriod(IsTauPeriodic(simple_two, 4));
    JsonPrint(",\"tau_dual_simple_n1\":"); PrintPeriod(IsTauPeriodic(dual_simple, 1));
    JsonPrint(",\"tau_dual_simple_n4\":"); PrintPeriod(IsTauPeriodic(dual_simple, 4));
    JsonPrint("}");
end;

JsonPrint("{\"schema_version\":1,\"kind\":\"modulehomalg_periodicity\",\"oracle\":{");
JsonPrint("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
JsonPrint(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
JsonPrint(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
JsonPrint(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
JsonPrint(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
JsonPrint(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
JsonPrint("},\"binary\":"); PrintCase(GF(2));
JsonPrint(",\"quinary\":"); PrintCase(GF(5));
JsonPrint(",\"rational\":"); PrintCase(Rationals);
JsonPrint("}\n");

QUIT_GAP(0);
