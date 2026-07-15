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

Sign := function(value)
    if value < 0 then
        return -1;
    elif value > 0 then
        return 1;
    else
        return 0;
    fi;
end;

ComparisonMatrix := function(ordering, paths)
    return List(paths,
        left -> List(paths,
            right -> Sign(ComparisonFunction(ordering)(left, right))));
end;

Q := Quiver(2, [[1, 1, "a"], [1, 2, "b"]]);;
a := Q.a;;
b := Q.b;;
paths := [Zero(Q), Q.v1, Q.v2, a, b, a*a, a*b, a*a*b];;
default_ordering := OrderingOfQuiver(Q);;
left_lex := LeftLexicographicOrdering(Q, [b, a]);;
right_lex := RightLexicographicOrdering(Q, [b, a]);;
reverse_ordering := ReverseOrdering(Q, default_ordering);;
left_vector := LeftVectorOrdering(Q, [b, a]);;
right_vector := RightVectorOrdering(Q, [b, a]);;
weight := WeightOrdering(Q, [1, 5]);;

Print("{\"schema_version\":1,\"kind\":\"orderings\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},\"paths\":");
PrintJsonStringList(List(paths, StringPrint));
Print(",\"matrices\":{");
Print("\"default\":", ComparisonMatrix(default_ordering, paths));
Print(",\"left_lex_ba\":", ComparisonMatrix(left_lex, paths));
Print(",\"right_lex_ba\":", ComparisonMatrix(right_lex, paths));
Print(",\"reverse_default\":", ComparisonMatrix(reverse_ordering, paths));
Print(",\"left_vector_ba\":", ComparisonMatrix(left_vector, paths));
Print(",\"right_vector_ba\":", ComparisonMatrix(right_vector, paths));
Print(",\"weight_1_5\":", ComparisonMatrix(weight, paths));
Print("}}\n");

QUIT_GAP(0);
