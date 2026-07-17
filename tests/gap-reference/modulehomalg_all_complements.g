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

PrintSide := function(side)
    local degree, first, differential;
    if IsList(side) then
        Print("[]");
        return;
    fi;
    Print("{"); first := true;
    for degree in [-3..4] do
        if not first then Print(","); fi;
        differential := DifferentialOfComplex(side, degree);
        Print("\"", degree, "\":{");
        Print("\"source\":"); PrintIntegerList(DimensionVector(Source(differential)));
        Print(",\"range\":"); PrintIntegerList(DimensionVector(Range(differential)));
        Print("}"); first := false;
    od;
    Print("}");
end;

PrintResult := function(result)
    Print("{\"resolution\":"); PrintSide(result[1]);
    Print(",\"coresolution\":"); PrintSide(result[2]);
    Print("}");
end;

PrintCase := function(field)
    local a2, projectives, simples, zero,
          loop_parent, dual, dual_regular, dual_zero;
    a2 := PathAlgebra(field,
        Quiver(["v1", "v2"], [["v1", "v2", "a"]]));
    projectives := IndecProjectiveModules(a2);
    simples := SimpleModules(a2);
    zero := ZeroModule(a2);

    loop_parent := PathAlgebra(field,
        Quiver(["w"], [["w", "w", "x"]]));
    dual := loop_parent / [loop_parent.x^2];
    dual_regular := IndecProjectiveModules(dual)[1];
    dual_zero := ZeroModule(dual);

    Print("{\"p1_complement_p2\":");
    PrintResult(AllComplementsOfAlmostCompleteTiltingModule(projectives[1], projectives[2]));
    Print(",\n\"p2_complement_p1\":");
    PrintResult(AllComplementsOfAlmostCompleteTiltingModule(projectives[2], projectives[1]));
    Print(",\n\"p1_complement_s1\":");
    PrintResult(AllComplementsOfAlmostCompleteTiltingModule(projectives[1], simples[1]));
    Print(",\n\"zero_complement_simple\":");
    PrintResult(AllComplementsOfAlmostCompleteTiltingModule(zero, simples[1]));
    Print(",\n\"dual_zero_complement_regular\":");
    PrintResult(AllComplementsOfAlmostCompleteTiltingModule(dual_zero, dual_regular));
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"modulehomalg_all_complements\",\"oracle\":{");
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
