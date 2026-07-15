if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintJsonStringList := function(values)
    local first, value;
    Print("[");
    first := true;
    for value in values do
        if not first then
            Print(",");
        fi;
        PrintJsonString(value);
        first := false;
    od;
    Print("]");
end;

PrintPairList := function(pairs)
    local first, pair;
    Print("[");
    first := true;
    for pair in pairs do
        if not first then
            Print(",");
        fi;
        Print("[");
        PrintJsonString(String(pair[1]));
        Print(",");
        PrintJsonString(String(pair[2]));
        Print("]");
        first := false;
    od;
    Print("]");
end;

PrintPoset := function(poset)
    local order, points;
    points := UnderlyingSet(poset);
    order := PartialOrderOfPoset(poset);
    Print("{\"size\":", Size(poset), ",\"underlying\":");
    PrintJsonStringList(List(points, String));
    Print(",\"minimal\":");
    PrintPairList(poset!.minimal_relations);
    Print(",\"matrix\":", List(points,
        left -> List(points, right -> order(left, right))), "}");
end;

example := Poset(["d", "b", "a", "c"],
    [["a", "b", "c"], ["b", "d"], ["c", "d"]]);;
duplicates := Poset([3, 1, 2, 1], [[1, 2], [2, 3]]);;
self_relation := Poset([2, 1], [[1, 1]]);;

Print("{\"schema_version\":1,\"kind\":\"poset\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},\"example\":");
PrintPoset(example);
Print(",\"duplicates\":");
PrintPoset(duplicates);
Print(",\"self_relation\":");
PrintPoset(self_relation);
Print("}\n");

QUIT_GAP(0);
