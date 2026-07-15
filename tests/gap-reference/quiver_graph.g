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

PrintQuiverRecord := function(quiver)
    local arrow, first;
    Print("{\"vertices\":");
    PrintJsonStringList(List(VerticesOfQuiver(quiver), String));
    Print(",\"arrows\":[");
    first := true;
    for arrow in ArrowsOfQuiver(quiver) do
        if not first then
            Print(",");
        fi;
        Print("{\"name\":");
        PrintJsonString(String(arrow));
        Print(",\"source\":");
        PrintJsonString(String(SourceOfPath(arrow)));
        Print(",\"target\":");
        PrintJsonString(String(TargetOfPath(arrow)));
        Print("}");
        first := false;
    od;
    Print("]}");
end;

line := Quiver([[0, 1, 0], [0, 0, 1], [0, 0, 0]]);;
parallel := Quiver([[0, 2, 0], [0, 0, 1], [0, 0, 0]]);;
cyclic := Quiver(3, [[1, 2, "a"], [2, 3, "b"], [3, 2, "c"]]);;
disconnected := Quiver(["p", "q", "r"], [["p", "q", "x"]]);;
components := ConnectedComponentsOfQuiver(disconnected);;
subquiver := FullSubquiver(cyclic, [cyclic.v3, cyclic.v2]);;
dynkin_a := DynkinQuiverAn(4, ["r", "l", "r"]);;
dynkin_d := DynkinQuiverDn(5, ["r", "l", "r", "l"]);;
dynkin_e := DynkinQuiverEn(6, ["r", "l", "r", "l", "d"]);;
separated := SeparatedQuiver(cyclic);;
doubled := DoubleQuiver(cyclic);;

Print("{\"schema_version\":1,\"kind\":\"quiver_graph\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},\"properties\":{");
Print("\"line\":{\"size\":", Size(line),
      ",\"connected\":", IsConnectedQuiver(line),
      ",\"unoriented_acyclic\":", IsUAcyclicQuiver(line),
      ",\"tree\":", IsTreeQuiver(line), "},");
Print("\"parallel\":{\"size\":", Size(parallel),
      ",\"connected\":", IsConnectedQuiver(parallel),
      ",\"unoriented_acyclic\":", IsUAcyclicQuiver(parallel),
      ",\"tree\":", IsTreeQuiver(parallel), "},");
Print("\"cyclic\":{\"size\":\"infinity\",\"connected\":",
      IsConnectedQuiver(cyclic),
      ",\"unoriented_acyclic\":", IsUAcyclicQuiver(cyclic),
      ",\"tree\":", IsTreeQuiver(cyclic), "}}");
Print(",\"components\":[");
PrintQuiverRecord(components[1]);
Print(",");
PrintQuiverRecord(components[2]);
Print("],\"full_subquiver\":");
PrintQuiverRecord(subquiver);
Print(",\"dynkin\":{");
Print("\"a4\":");
PrintQuiverRecord(dynkin_a);
Print(",\"d5\":");
PrintQuiverRecord(dynkin_d);
Print(",\"e6\":");
PrintQuiverRecord(dynkin_e);
Print("},\"separated\":");
PrintQuiverRecord(separated);
Print(",\"double\":");
PrintQuiverRecord(doubled);
Print("}\n");

QUIT_GAP(0);
