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

generic := StaticDictionary(["he", "she", "his", "hers"]);;
generic_overlap := StaticDictionary([[1, 2], [2], [1, 2]]);;

Q := Quiver(1, [[1, 1, "x"], [1, 1, "y"]]);;
x := Q.x;;
y := Q.y;;
finite := QuiverStaticDictionary(Q, [x*x, y*x, y*y]);;
infinite := QuiverStaticDictionary(Q, [x*x]);;
vertex_prefix := QuiverStaticDictionary(Q, [Q.v1, x]);;
overlap := QuiverStaticDictionary(Q, [x, x*y]);;

Print("{\"schema_version\":1,\"kind\":\"static_dictionary\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},\"generic\":{\"patterns\":");
PrintJsonStringList(Patterns(generic));
Print(",\"ushers\":", Search(generic, "ushers"),
      ",\"overlap\":", Search(generic_overlap, [1, 2, 1, 2]), "}");
Print(",\"quiver\":{\"patterns\":");
PrintJsonStringList(List(Patterns(finite), StringPrint));
Print(",\"search\":", Search(finite, y*x*x),
      ",\"prefix\":", PrefixSearch(finite, x*x),
      ",\"finite\":", IsFiniteDifference(finite),
      ",\"size\":", DifferenceSize(finite),
      ",\"words\":");
PrintJsonStringList(List(DifferenceWords(finite), StringPrint));
Print(",\"infinite\":", IsFiniteDifference(infinite),
      ",\"infinite_size\":\"", DifferenceSize(infinite), "\"",
      ",\"infinite_words_fail\":", DifferenceWords(infinite) = fail,
      ",\"vertex_prefix\":", PrefixSearch(vertex_prefix, x),
      ",\"overlap_finite\":", IsFiniteDifference(overlap),
      ",\"overlap_size\":");
if IsFiniteDifference(overlap) then
    Print(DifferenceSize(overlap));
else
    Print("\"infinity\"");
fi;
Print("}}\n");

QUIT_GAP(0);
