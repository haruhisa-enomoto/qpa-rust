if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintStringList := function(values)
    local first, value;
    Print("["); first := true;
    for value in values do
        if not first then Print(","); fi;
        PrintJsonString(value); first := false;
    od;
    Print("]");
end;

PrintMaybeString := function(value)
    if IsString(value) then
        PrintJsonString(value);
    else
        Print(value);
    fi;
end;

PrintMaybeStrings := function(values)
    local first, value;
    Print("["); first := true;
    for value in values do
        if not first then Print(","); fi;
        PrintMaybeString(value); first := false;
    od;
    Print("]");
end;

RelationStrings := function(relations)
    return List(relations, relation ->
        Concatenation(List(WalkOfPath(relation), String)));
end;

PrintCase := function(field)
    local Q, kQ, A, array, sigma, eps, inputs, valid, extensions, zeroExtensions;
    Q := Quiver(4,
        [[1,2,"a"], [1,2,"b"], [2,3,"c"], [3,4,"d"], [3,4,"e"]]);
    kQ := PathAlgebra(field, Q);
    A := kQ / [kQ.b*kQ.c, kQ.c*kQ.d];
    array := QPAStringSigmaEps(A);
    sigma := array[1]; eps := array[2];
    inputs := ["eca", "Ab", "bc", "CB", "cb", "BC", "dc", "CD",
        "(1,1)", "(4,-1)",
        "(5,1)", "aa", "aA", "ca", "1a", "-1a", "a1", "a-1"];
    valid := List(inputs, value -> IsValidString(A, ShallowCopy(value)));
    extensions := [
        QPAStringDirectLeft(A, "eca", sigma, eps),
        QPAStringInverseLeft(A, "eca", sigma, eps),
        QPAStringDirectRight(A, "eca", sigma, eps),
        QPAStringInverseRight(A, "eca", sigma, eps)
    ];
    zeroExtensions := [
        QPAStringDirectLeft(A, "(1,1)", sigma, eps),
        QPAStringInverseLeft(A, "(2,-1)", sigma, eps),
        QPAStringDirectRight(A, "(4,1)", sigma, eps),
        QPAStringInverseRight(A, "(3,-1)", sigma, eps)
    ];

    Print("{\"rho\":");
    PrintStringList(RelationStrings(QPAStringQuiverRho(Q, ["bc", "cd"])));
    Print(",\"sigma\":", sigma, ",\"epsilon\":", eps);
    Print(",\"inputs\":"); PrintStringList(inputs);
    Print(",\"valid\":", valid);
    Print(",\"extensions\":"); PrintMaybeStrings(extensions);
    Print(",\"zero_extensions\":"); PrintMaybeStrings(zeroExtensions);
    Print(",\"local_ar\":"); PrintStringList(LocalARQuiver(A, "eca"));
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"validity\",\"oracle\":{");
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
Print(",\"rational\":"); PrintCase(Rationals);
Print("}\n");

QUIT_GAP(0);
