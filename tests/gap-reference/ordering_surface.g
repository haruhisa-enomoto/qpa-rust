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

OrderingFilters := function(ordering)
    return [
        IsQuiverOrdering(ordering),
        IsLexicographicOrdering(ordering),
        IsLengthOrdering(ordering),
        IsReverseOrdering(ordering),
        IsVectorOrdering(ordering),
        IsWeightOrdering(ordering),
        IsBlockOrdering(ordering),
        IsWreathOrdering(ordering),
        IsLeftVectorOrdering(ordering),
        IsRightVectorOrdering(ordering),
        IsLeftLexicographicOrdering(ordering),
        IsRightLexicographicOrdering(ordering)
    ];
end;

Q := Quiver(2, [[1, 1, "a"], [1, 2, "b"]]);;
a := Q.a;;
b := Q.b;;
default := OrderingOfQuiver(Q);;
length := LengthOrdering(Q);;
left_lex := LeftLexicographicOrdering(Q, [b, a]);;
right_lex := RightLexicographicOrdering(Q, [b, a]);;
reverse := ReverseOrdering(Q, default);;
left_vector := LeftVectorOrdering(Q, [b, a]);;
right_vector := RightVectorOrdering(Q, [b, a]);;
weight := WeightOrdering(Q, [1, 5]);;
block := BlockOrdering(Q, [[weight, GeneratorsOfQuiver(Q)]]);;
wreath := WreathOrdering(Q, [[weight, GeneratorsOfQuiver(Q)]]);;
orderings := [length, left_lex, right_lex, reverse, left_vector,
              right_vector, weight, block, wreath];;

empty := PathFromArrowList(Q, []);;
singleton_vertex := PathFromArrowList(Q, [Q.v1]);;
formal_vertices := PathFromArrowList(Q, [Q.v1, Q.v2]);;
formal_noncomposable := PathFromArrowList(Q, [b, a]);;

QC := Quiver(2, [[1, 2, "x"], [2, 1, "y"]]);;
cycle_ordering := LengthOrdering(QC);;

Print("{\"schema_version\":1,\"kind\":\"ordering_surface\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},\"category_filter_order\":");
PrintJsonStringList([
    "quiver", "lexicographic", "length", "reverse", "vector", "weight",
    "block", "wreath", "left_vector", "right_vector", "left_lexicographic",
    "right_lexicographic"
]);
Print(",\"categories\":", List(orderings, OrderingFilters));
Print(",\"names\":");
PrintJsonStringList(List(orderings, OrderingName));
Print(",\"views\":");
PrintJsonStringList(List(orderings, StringPrint));
Print(",\"attributes\":{");
Print("\"default_view\":");
PrintJsonString(StringPrint(default));
Print(",\"default_name\":");
PrintJsonString(OrderingName(default));
Print(",\"default_has_next\":", HasNextOrdering(default));
Print(",\"default_next_identity\":",
      IsIdenticalObj(NextOrdering(default), NextOrdering(default)));
Print(",\"default_lex_table\":");
PrintJsonStringList(List(LexicographicTable(default), StringPrint));
Print(",\"reverse_next_identity\":",
      IsIdenticalObj(NextOrdering(reverse), default));
Print(",\"reverse_lex_identity\":",
      IsIdenticalObj(LexicographicTable(reverse), LexicographicTable(default)));
Print("}");
Print(",\"comparison_values\":{");
Print("\"weight_a_b\":");
PrintJsonString(String(ComparisonFunction(weight)(a, b)));
Print(",\"weight_aa_b\":");
PrintJsonString(String(CompareArrowLists(weight, [a, a], [b])));
Print(",\"block_a_b\":");
PrintJsonString(String(ComparisonFunction(block)(a, b)));
Print(",\"wreath_a_b\":");
PrintJsonString(String(ComparisonFunction(wreath)(a, b)));
Print("}");
Print(",\"path_from_arrow_list\":{");
Print("\"empty_is_zero\":", IsZeroPath(empty));
Print(",\"singleton_vertex_identity\":",
      IsIdenticalObj(singleton_vertex, Q.v1));
Print(",\"formal_vertices\":{");
Print("\"view\":");
PrintJsonString(StringPrint(formal_vertices));
Print(",\"extrep\":", ExtRepOfObj(formal_vertices));
Print(",\"walk\":");
PrintJsonStringList(List(WalkOfPath(formal_vertices), StringPrint));
Print(",\"source\":");
PrintJsonString(StringPrint(SourceOfPath(formal_vertices)));
Print(",\"target\":");
PrintJsonString(StringPrint(TargetOfPath(formal_vertices)));
Print(",\"length\":", LengthOfPath(formal_vertices), "}");
Print(",\"formal_noncomposable\":{");
Print("\"view\":");
PrintJsonString(StringPrint(formal_noncomposable));
Print(",\"extrep\":", ExtRepOfObj(formal_noncomposable));
Print(",\"walk\":");
PrintJsonStringList(List(WalkOfPath(formal_noncomposable), StringPrint));
Print(",\"source\":");
PrintJsonString(StringPrint(SourceOfPath(formal_noncomposable)));
Print(",\"target\":");
PrintJsonString(StringPrint(TargetOfPath(formal_noncomposable)));
Print(",\"length\":", LengthOfPath(formal_noncomposable), "}}");
Print(",\"block_cycle_checks\":[",
      CheckBlocksForCycles(cycle_ordering, [1, 1, 1, 1]), ",",
      CheckBlocksForCycles(cycle_ordering, [1, 1, 1, 2]), ",",
      CheckBlocksForCycles(cycle_ordering, [0, 0, 1, 2]), "]}");
Print("\n");

QUIT_GAP(0);
