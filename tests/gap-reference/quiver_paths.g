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

PrintPathRecord := function(label, path)
    Print("{\"label\":");
    PrintJsonString(label);
    Print(",\"display\":");
    PrintJsonString(StringPrint(path));
    Print(",\"external_representation\":", ExtRepOfObj(path));
    Print(",\"length\":");
    if HasLengthOfPath(path) then
        Print(LengthOfPath(path));
    else
        Print("null");
    fi;
    Print(",\"walk\":");
    if HasWalkOfPath(path) then
        Print(List(WalkOfPath(path), ArrowIndex));
    else
        Print("null");
    fi;
    Print(",\"source\":");
    PrintJsonString(StringPrint(SourceOfPath(path)));
    Print(",\"target\":");
    PrintJsonString(StringPrint(TargetOfPath(path)));
    Print("}");
end;

Q := Quiver(3, [[1, 2, "a"], [2, 3, "b"], [3, 2, "c"]]);;
a := Q.a;;
b := Q.b;;
c := Q.c;;
ab := a * b;;
abc := ab * c;;

R := Quiver(["x", "y"], [["x", "y", "a"], ["x", "x", "b"]]);;
iterator := Iterator(R);;
iterator_first := List([1..8], index -> NextIterator(iterator));;

Qop := OppositeQuiver(Q);;
abop := OppositePath(ab);;

Print("{\"schema_version\":1,\"kind\":\"quiver_paths\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":");
PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"quiver\":{");
Print("\"vertices\":");
PrintJsonStringList(List(VerticesOfQuiver(Q), String));
Print(",\"arrows\":[");
for i in [1..NumberOfArrows(Q)] do
    if i > 1 then
        Print(",");
    fi;
    Print("{\"name\":");
    PrintJsonString(String(ArrowsOfQuiver(Q)[i]));
    Print(",\"index\":", ArrowIndex(ArrowsOfQuiver(Q)[i]));
    Print(",\"generator_index\":", GeneratorIndex(ArrowsOfQuiver(Q)[i]));
    Print(",\"source\":", VertexIndex(SourceOfPath(ArrowsOfQuiver(Q)[i])));
    Print(",\"target\":", VertexIndex(TargetOfPath(ArrowsOfQuiver(Q)[i])), "}");
od;
Print("],\"adjacency\":", AdjacencyMatrixOfQuiver(Q));
Print(",\"special_arrow_names\":", IsQuiverSA(Q));
Print(",\"acyclic\":", IsAcyclicQuiver(Q), "}");
Print(",\"paths\":[");
PrintPathRecord("zero", Zero(Q));
Print(",");
PrintPathRecord("v2", Q.v2);
Print(",");
PrintPathRecord("a", a);
Print(",");
PrintPathRecord("ab", ab);
Print(",");
PrintPathRecord("abc", abc);
Print("],\"products\":{");
Print("\"a_b\":", ExtRepOfObj(a * b));
Print(",\"b_a\":", ExtRepOfObj(b * a));
Print(",\"v1_a\":", ExtRepOfObj(Q.v1 * a));
Print(",\"v2_a\":", ExtRepOfObj(Q.v2 * a));
Print(",\"a_v2\":", ExtRepOfObj(a * Q.v2));
Print(",\"a_v1\":", ExtRepOfObj(a * Q.v1), "}");
Print(",\"ordering\":{");
Print("\"zero_lt_v1\":", Zero(Q) < Q.v1);
Print(",\"v1_lt_v2\":", Q.v1 < Q.v2);
Print(",\"v2_lt_a\":", Q.v2 < a);
Print(",\"a_lt_ab\":", a < ab, "}");
Print(",\"iterator_first\":");
PrintJsonStringList(List(iterator_first, StringPrint));
Print(",\"opposite\":{");
Print("\"vertices\":");
PrintJsonStringList(List(VerticesOfQuiver(Qop), String));
Print(",\"arrows\":");
PrintJsonStringList(List(ArrowsOfQuiver(Qop), String));
Print(",\"ab_display\":");
PrintJsonString(StringPrint(abop));
Print(",\"ab_external_representation\":", ExtRepOfObj(abop));
Print(",\"double_opposite_is_identical\":", IsIdenticalObj(OppositeQuiver(Qop), Q));
Print("}}\n");

QUIT_GAP(0);
