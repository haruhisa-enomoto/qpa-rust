if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintBoolean := function(value)
    if value then Print("true"); else Print("false"); fi;
end;

CaseData := function(field)
    local a2, zero, simple_one, projective, semisimple_pair,
          loop_parent, dual, dual_simple, dual_regular;
    a2 := PathAlgebra(field,
        Quiver(["v1", "v2"], [["v1", "v2", "a"]]));
    zero := ZeroModule(a2);
    simple_one := RightModuleOverPathAlgebra(a2, [1, 0], []);
    projective := RightModuleOverPathAlgebra(a2, [1, 1],
        [["a", [[1]] * One(field)]]);
    semisimple_pair := RightModuleOverPathAlgebra(a2, [1, 1], []);

    loop_parent := PathAlgebra(field,
        Quiver(["w"], [["w", "w", "x"]]));
    dual := loop_parent / [loop_parent.x^2];
    dual_simple := RightModuleOverPathAlgebra(dual, [1],
        [["x", [[0]] * One(field)]]);
    dual_regular := RightModuleOverPathAlgebra(dual, [2],
        [["x", [[0, 1], [0, 0]] * One(field)]]);

    return rec(
        zero := zero,
        simple_one := simple_one,
        projective := projective,
        semisimple_pair := semisimple_pair,
        dual := dual,
        dual_simple := dual_simple,
        dual_regular := dual_regular
    );
end;

PrintCommonCase := function(field)
    local data;
    data := CaseData(field);
    Print("{\"exceptional_zero\":");
    PrintBoolean(IsExceptionalModule(data.zero));
    Print(",\"exceptional_simple_one\":");
    PrintBoolean(IsExceptionalModule(data.simple_one));
    Print(",\"exceptional_projective\":");
    PrintBoolean(IsExceptionalModule(data.projective));
    Print(",\"exceptional_dual_simple\":");
    PrintBoolean(IsExceptionalModule(data.dual_simple));
    Print(",\"complexity_a2_simple_n2\":", ComplexityOfModule(data.simple_one, 2));
    Print(",\"complexity_dual_simple_n2\":", ComplexityOfModule(data.dual_simple, 2));
    Print(",\"complexity_dual_simple_n3\":", ComplexityOfModule(data.dual_simple, 3));
    Print(",\"complexity_dual_simple_n10\":", ComplexityOfModule(data.dual_simple, 10));
    Print(",\"complexity_dual_regular_n2\":", ComplexityOfModule(data.dual_regular, 2));
    Print(",\"complexity_dual_algebra_n2\":", ComplexityOfAlgebra(data.dual, 2));
    Print(",\"complexity_dual_algebra_n10\":", ComplexityOfAlgebra(data.dual, 10));
    Print("}");
end;

PrintFiniteCase := function(field)
    local data;
    data := CaseData(field);
    Print("{\"common\":");
    PrintCommonCase(field);
    Print(",\"exceptional_semisimple_pair\":");
    PrintBoolean(IsExceptionalModule(data.semisimple_pair));
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"module_complexity\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"binary\":"); PrintFiniteCase(GF(2));
Print(",\"quinary\":"); PrintFiniteCase(GF(5));
Print(",\"rational\":"); PrintCommonCase(Rationals);
Print("}\n");

QUIT_GAP(0);
