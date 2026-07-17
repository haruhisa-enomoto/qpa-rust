if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintCase := function(field)
    local a2, a2simple, zero, dual, dualsimple, dualprojective,
          cycle, cyclesimple;

    a2 := PathAlgebra(field,
        Quiver(["a1", "a2"], [["a1", "a2", "a"]]));
    a2simple := SimpleModules(a2)[1];
    zero := ZeroModule(a2);

    dual := PathAlgebra(field,
        Quiver(["d"], [["d", "d", "x"]]));
    dual := dual / [dual.x^2];
    dualsimple := SimpleModules(dual)[1];
    dualprojective := IndecProjectiveModules(dual)[1];

    cycle := PathAlgebra(field,
        Quiver(["c1", "c2"],
               [["c1", "c2", "p"], ["c2", "c1", "q"]]));
    cycle := cycle / [cycle.p * cycle.q, cycle.q * cycle.p];
    cyclesimple := SimpleModules(cycle)[1];

    Print("{\"a2_simple_negative\":", N_RigidModule(a2simple, -1));
    Print(",\"a2_simple_0\":", N_RigidModule(a2simple, 0));
    Print(",\"a2_simple_3\":", N_RigidModule(a2simple, 3));
    Print(",\"zero_3\":", N_RigidModule(zero, 3));
    Print(",\"dual_simple_0\":", N_RigidModule(dualsimple, 0));
    Print(",\"dual_simple_1\":", N_RigidModule(dualsimple, 1));
    Print(",\"dual_projective_3\":", N_RigidModule(dualprojective, 3));
    Print(",\"cycle_simple_1\":", N_RigidModule(cyclesimple, 1));
    Print(",\"cycle_simple_2\":", N_RigidModule(cyclesimple, 2));
    Print(",\"cycle_simple_3\":", N_RigidModule(cyclesimple, 3));
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"modulehomalg_n_rigid\",\"oracle\":{");
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
