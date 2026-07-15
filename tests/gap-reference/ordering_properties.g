if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

OrderingProperties := function(ordering)
    return [IsWellOrdering(ordering), IsWellReversedOrdering(ordering),
            IsTotalOrdering(ordering), IsAdmissibleOrdering(ordering)];
end;

SubsetProperties := function(ordering, subset)
    return [IsWellSubset(ordering, subset),
            IsWellReversedSubset(ordering, subset),
            IsTotalSubset(ordering, subset)];
end;

QA := Quiver(2, [[1, 2, "a"]]);;
qa_default := OrderingOfQuiver(QA);;
qa_length := LengthOrdering(QA);;
qa_lex := LeftLexicographicOrdering(QA);;
qa_reverse := ReverseOrdering(QA, qa_length);;
qa_vector := LeftVectorOrdering(QA, []);;
qa_weight := WeightOrdering(QA, [0]);;
qa_block := BlockOrdering(QA, [[qa_default, GeneratorsOfQuiver(QA)]]);;

QL := Quiver(1, [[1, 1, "x"]]);;
ql_default := OrderingOfQuiver(QL);;
ql_length := LengthOrdering(QL);;
ql_vector := RightVectorOrdering(QL, []);;
ql_weight_zero := WeightOrdering(QL, [0]);;
ql_weight_one := WeightOrdering(QL, [1]);;
ql_block := BlockOrdering(QL, [[ql_default, GeneratorsOfQuiver(QL)]]);;

QS := Quiver(2, [[1, 1, "a"], [1, 2, "b"]]);;
a := QS.a;;
b := QS.b;;
v := QS.v1;;
qs_default := OrderingOfQuiver(QS);;
qs_length := LengthOrdering(QS);;
qs_lex := LeftLexicographicOrdering(QS, [b, a]);;
qs_weight_zero := WeightOrdering(QS, [0, 0]);;
qs_weight_mixed := WeightOrdering(QS, [1, 0]);;
qs_block := BlockOrdering(QS,
    [[qs_default, [a]], [qs_default, [b]]]);;
subsets := [[], [v], [a], [b], [a, b], GeneratorsOfQuiver(QS)];;

universe := [[1, 1, "a"], [1, 2, "b"],
             [2, 1, "c"], [2, 2, "d"]];;
exhaustive := [];;
for quiver_mask in [0..15] do
    specifications := Filtered([1..4],
        i -> QuoInt(quiver_mask, 2^(i-1)) mod 2 = 1);;
    QX := Quiver(2, universe{specifications});;
    arrows := ArrowsOfQuiver(QX);;
    lex := LeftLexicographicOrdering(QX);;
    reverse_lex := LeftLexicographicOrdering(QX, Reversed(arrows));;
    weights := List([1..Length(arrows)], i -> i mod 2);;
    weight := WeightOrdering(QX, weights);;
    for subset_mask in [0..2^Length(arrows)-1] do
        subset := Filtered([1..Length(arrows)],
            i -> QuoInt(subset_mask, 2^(i-1)) mod 2 = 1);;
        subset := arrows{subset};;
        Add(exhaustive, [quiver_mask, subset_mask,
            CheckQuiverSubsetForCycles(QX, subset),
            IsWellSubset(lex, subset), IsWellReversedSubset(lex, subset),
            IsWellSubset(reverse_lex, subset),
            IsWellReversedSubset(reverse_lex, subset),
            IsWellReversedSubset(weight, subset)]);
    od;
od;

Print("{\"schema_version\":1,\"kind\":\"ordering_properties\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},\"acyclic\":{");
Print("\"length\":", OrderingProperties(qa_length));
Print(",\"default\":", OrderingProperties(qa_default));
Print(",\"lex\":", OrderingProperties(qa_lex));
Print(",\"reverse_length\":", OrderingProperties(qa_reverse));
Print(",\"vector\":", OrderingProperties(qa_vector));
Print(",\"weight\":", OrderingProperties(qa_weight));
Print(",\"block\":", OrderingProperties(qa_block));
Print("},\"loop\":{");
Print("\"length\":", OrderingProperties(ql_length));
Print(",\"default\":", OrderingProperties(ql_default));
Print(",\"vector\":", OrderingProperties(ql_vector));
Print(",\"weight_zero\":", OrderingProperties(ql_weight_zero));
Print(",\"weight_one\":", OrderingProperties(ql_weight_one));
Print(",\"block\":", OrderingProperties(ql_block));
Print("},\"subsets\":{");
Print("\"cycle_checks\":", List(subsets,
    subset -> CheckQuiverSubsetForCycles(QS, subset)));
Print(",\"length\":", List(subsets,
    subset -> SubsetProperties(qs_length, subset)));
Print(",\"lex_ba\":", List(subsets,
    subset -> SubsetProperties(qs_lex, subset)));
Print(",\"weight_zero\":", List(subsets,
    subset -> SubsetProperties(qs_weight_zero, subset)));
Print(",\"weight_mixed\":", List(subsets,
    subset -> SubsetProperties(qs_weight_mixed, subset)));
Print(",\"block\":", List(subsets,
    subset -> SubsetProperties(qs_block, subset)));
Print("},\"exhaustive\":", exhaustive);
Print("}\n");

QUIT_GAP(0);
