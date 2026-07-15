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

PathText := function(path)
    if IsQuiverVertex(path) then
        return String(path);
    elif IsZeroPath(path) then
        return "0";
    fi;
    return JoinStringsWithSeparator(List(WalkOfPathOrVertex(path), String), "*");
end;

q1 := Quiver(["u1", "u2"], [["u1", "u2", "a"]]);;
q2 := Quiver(["v1", "v2", "v3"],
             [["v1", "v2", "b"], ["v2", "v3", "c"]]);;
product := QuiverProduct(q1, q2);;
decomposition := QuiverProductDecomposition(product);;
included := IncludeInProductQuiver([q1.a, q2.b * q2.c], product);;
included_vertices := IncludeInProductQuiver([q1.u1, q2.v2], product);;

line := Quiver(3, [[1, 2, "x"], [2, 3, "y"]]);;
line_path := line.x * line.y;;
loops := Quiver(1, [[1, 1, "s"], [1, 1, "t"]]);;

Print("{\"schema_version\":1,\"kind\":\"pathalgtensor_quiver\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},\"product\":{\"vertices\":");
PrintJsonStringList(List(VerticesOfQuiver(product), String));
Print(",\"arrows\":");
PrintJsonStringList(List(ArrowsOfQuiver(product), String));
Print(",\"arrow_endpoints\":[");
for i in [1..Length(ArrowsOfQuiver(product))] do
    if i > 1 then
        Print(",");
    fi;
    PrintJsonStringList([
        String(SourceOfPath(ArrowsOfQuiver(product)[i])),
        String(TargetOfPath(ArrowsOfQuiver(product)[i]))
    ]);
od;
Print("]},\"decomposition\":{");
Print("\"length\":", Length(decomposition),
      ",\"bound\":[", IsBound(decomposition[1]), ",",
      IsBound(decomposition[2]), ",", IsBound(decomposition[3]), "]",
      ",\"factor_identity\":[", IsIdenticalObj(decomposition[1], q1), ",",
      IsIdenticalObj(decomposition[2], q2), "]},");
Print("\"walks\":{");
Print("\"vertex\":");
PrintJsonStringList(List(WalkOfPathOrVertex(q1.u1), String));
Print(",\"arrow\":");
PrintJsonStringList(List(WalkOfPathOrVertex(q1.a), String));
Print(",\"path\":");
PrintJsonStringList(List(WalkOfPathOrVertex(q2.b * q2.c), String));
Print("},\"include\":{");
Print("\"path\":");
PrintJsonString(PathText(included));
Print(",\"vertices\":");
PrintJsonString(PathText(included_vertices));
Print("},\"project\":{");
Print("\"first\":");
PrintJsonString(PathText(ProjectFromProductQuiver(1, included)));
Print(",\"second\":");
PrintJsonString(PathText(ProjectFromProductQuiver(2, included)));
Print("},\"reverse\":{");
Print("\"arrow\":");
PrintJsonString(PathText(ReversePath(line.x)));
Print(",\"line_path\":");
PrintJsonString(PathText(ReversePath(line_path)));
Print(",\"loop_path\":");
PrintJsonString(PathText(ReversePath(loops.s * loops.t)));
Print("}}\n");

QUIT_GAP(0);
