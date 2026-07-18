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

PrintCase := function(field)
    local Q, kQ, A, endpoints, bandInputs, bands;
    Q := Quiver(4,
        [[1,2,"a"], [1,2,"b"], [2,3,"c"], [3,4,"d"], [3,4,"e"]]);
    kQ := PathAlgebra(field, Q);
    A := kQ / [kQ.b*kQ.c, kQ.c*kQ.d];
    endpoints := List(['a','A','e','E'], ch ->
        [String(SourceOfArrow(Q, ch)), String(TargetOfArrow(Q, ch))]);
    bandInputs := ["eca", "Ab", "aB", "De", "(1,1)", "A", "AbAb"];
    bands := List(bandInputs, word -> IsABand(A, word));

    Print("{\"lps\":[", computeLPSarray("ababaca"), ",",
        computeLPSarray("aaaa"), ",", computeLPSarray(""), "]");
    Print(",\"endpoints\":", endpoints);
    Print(",\"cyclic\":"); PrintStringList(CyclicPermutationOfABand("Abc"));
    Print(",\"maximum_direct\":", MaximumLengthOfDirectString(A));
    Print(",\"joints\":", NumberOfJoints(A));
    Print(",\"band_inputs\":"); PrintStringList(bandInputs);
    Print(",\"is_band\":", bands);
    Print(",\"strings_1\":"); PrintStringList(StringsLessThan(A, 1));
    Print(",\"strings_2\":"); PrintStringList(StringsLessThan(A, 2));
    Print(",\"bands_3\":"); PrintStringList(BandsLessThan(A, 3));
    Print(",\"representatives_3\":");
    PrintStringList(BandRepresentativesLessThan(A, 3));
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"stringalgebra\",\"oracle\":{");
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
