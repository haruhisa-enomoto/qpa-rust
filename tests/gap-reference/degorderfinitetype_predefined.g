if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

cases := [
    ["D4 subspace", ["D4 subspace"]],
    ["E6 subspace", ["E6 subspace"]],
    ["A3 zero", ["A3 zero"]],
    ["A1,2 zero", ["A1,2 zero"]],
    ["A2,2 comm", ["A2,2 comm"]],
    ["R nilp", ["R nilp"]],
    ["BG2", ["BG", 2]],
    ["BG5", ["BG", 5]],
    ["BG7", ["BG", 7]],
    ["BG8", ["BG", 8]],
    ["BG9", ["BG", 9]],
    ["BG11", ["BG", 11]],
    ["BG12", ["BG", 12]],
    ["BG13", ["BG", 13]],
    ["BG14", ["BG", 14]],
    ["BGo2", ["BGo", 2]],
    ["BGo13", ["BGo", 13]],
    ["BGo14", ["BGo", 14]],
    ["BG1r1", ["BG", 1, 1]],
    ["BG1r2", ["BG", 1, 2]],
    ["BG3r1", ["BG", 3, 1]],
    ["BG3r2", ["BG", 3, 2]],
    ["BG3r3", ["BG", 3, 3]],
    ["BG4r1", ["BG", 4, 1]],
    ["BG4r2", ["BG", 4, 2]],
    ["BG4r3", ["BG", 4, 3]],
    ["BG4r4", ["BG", 4, 4]],
    ["BG4r5", ["BG", 4, 5]],
    ["BG4r6", ["BG", 4, 6]],
    ["BG4r7", ["BG", 4, 7]],
    ["BG6r1", ["BG", 6, 1]],
    ["BG6r2", ["BG", 6, 2]],
    ["BG6r3", ["BG", 6, 3]],
    ["BG6r4", ["BG", 6, 4]],
    ["BG10r1", ["BG", 10, 1]],
    ["BG10r2", ["BG", 10, 2]]
];;

Print("{\"schema_version\":1,\"kind\":\"degorderfinitetype_predefined\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"cases\":[");
for i in [1..Length(cases)] do
    if i > 1 then Print(","); fi;
    Print("[\""); Print(cases[i][1]); Print("\",");
    Print(CallFuncList(PredefARQuivers, cases[i][2]));
    Print("]");
od;
Print("]}\n");

QUIT_GAP(0);
