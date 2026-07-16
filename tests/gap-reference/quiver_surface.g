if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    escaped := ReplacedString(escaped, "\n", "\\n");
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

PrintBool := function(value)
    if value then
        Print("true");
    else
        Print("false");
    fi;
end;

Q := Quiver(2, [[1, 2, "a"], [1, 1, "b"]]);;
indexed := Quiver(2, [[1, 2], [2, 1]]);;
named := Quiver(["p", "q"], [["p", "q"], ["q", "p"]]);;
adjacency := Quiver([[0, 2], [0, 0]]);;

O := RightLexicographicOrdering(Q, Reversed(GeneratorsOfQuiver(Q)));;
ordered := OrderedBy(Q, O);;
ordered_ba := ordered.b * ordered.a;;
Qop := OppositeQuiver(Q);;

next1 := NextPath(Q, 0);;
next2 := NextPath(Q, next1[2]);;
enumerator := Enumerator(Q);;
empty := Quiver([], []);;
empty_iterator := Iterator(empty);;
empty_done_before := IsDoneIterator(empty_iterator);;
empty_next := NextIterator(empty_iterator);;
empty_done_after := IsDoneIterator(empty_iterator);;

unique1 := UniqueQuiverName();;
unique2 := UniqueQuiverName();;

Print("{\"schema_version\":1,\"kind\":\"quiver_surface\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},\"constructors\":{");
Print("\"indexed_vertices\":");
PrintJsonStringList(List(VerticesOfQuiver(indexed), String));
Print(",\"indexed_arrows\":");
PrintJsonStringList(List(ArrowsOfQuiver(indexed), String));
Print(",\"named_vertices\":");
PrintJsonStringList(List(VerticesOfQuiver(named), String));
Print(",\"named_arrows\":");
PrintJsonStringList(List(ArrowsOfQuiver(named), String));
Print(",\"adjacency_vertices\":");
PrintJsonStringList(List(VerticesOfQuiver(adjacency), String));
Print(",\"adjacency_arrows\":");
PrintJsonStringList(List(ArrowsOfQuiver(adjacency), String));
Print(",\"adjacency_matrix\":", AdjacencyMatrixOfQuiver(adjacency));
Print(",\"special_named_arrows\":");
PrintBool(IsQuiverSA(Q));
Print(",\"default_ordering_cached\":");
PrintBool(IsIdenticalObj(OrderingOfQuiver(Q), OrderingOfQuiver(Q)));
Print("},\"attributes\":{");
Print("\"incoming_y\":");
PrintJsonStringList(List(IncomingArrowsOfVertex(Q.v2), String));
Print(",\"outgoing_x\":");
PrintJsonStringList(List(OutgoingArrowsOfVertex(Q.v1), String));
Print(",\"in_degree_y\":", InDegreeOfVertex(Q.v2));
Print(",\"out_degree_x\":", OutDegreeOfVertex(Q.v1));
Print(",\"neighbors_x\":");
PrintJsonStringList(List(NeighborsOfVertex(Q.v1), String));
Print(",\"piece_lookup_is_identical\":");
PrintBool(IsIdenticalObj(Q.a, ArrowsOfQuiver(Q)[1]));
Print("},\"printing\":{");
Print("\"quiver\":");
PrintJsonString(StringPrint(Q));
Print(",\"view\":");
PrintJsonString(StringView(Q));
Print(",\"power_path\":");
PrintJsonString(StringPrint(Q.b * Q.b * Q.a));
Print("},\"ordered\":{");
Print("\"ordering_is_identical\":");
PrintBool(IsIdenticalObj(OrderingOfQuiver(ordered), O));
Print(",\"new_family\":");
PrintBool(not IsIdenticalObj(Q.v1, ordered.v1));
Print(",\"vertices\":");
PrintJsonStringList(List(VerticesOfQuiver(ordered), String));
Print(",\"arrows\":");
PrintJsonStringList(List(ArrowsOfQuiver(ordered), String));
Print(",\"v1_lt_v2\":");
PrintBool(ordered.v1 < ordered.v2);
Print(",\"a_lt_b\":");
PrintBool(ordered.a < ordered.b);
Print(",\"ba\":");
PrintJsonString(StringPrint(ordered_ba));
Print("},\"opposite_name_map\":{");
Print("\"append\":");
PrintJsonString(OppositeQuiverNameMap(Q)("plain"));
Print(",\"remove_without_suffix_check\":");
PrintJsonString(OppositeQuiverNameMap(Qop)("abcd"));
Print("},\"next_path\":{");
Print("\"first\":");
PrintJsonString(StringPrint(next1[1]));
Print(",\"first_queue\":");
PrintJsonStringList(List(next1[2], StringPrint));
Print(",\"second\":");
PrintJsonString(StringPrint(next2[1]));
Print(",\"second_queue\":");
PrintJsonStringList(List(next2[2], StringPrint));
Print(",\"enumerator_5\":");
PrintJsonString(StringPrint(enumerator[5]));
Print(",\"empty_next_is_fail\":");
PrintBool(empty_next = fail);
Print(",\"empty_done_before\":");
PrintBool(empty_done_before);
Print(",\"empty_done_after\":");
PrintBool(empty_done_after);
Print("},\"indices\":{");
Print("\"vertex\":", VertexIndex(Q.v2));
Print(",\"arrow\":", ArrowIndex(Q.b));
Print(",\"generator_vertex\":", GeneratorIndex(Q.v2));
Print(",\"generator_arrow\":", GeneratorIndex(Q.b));
Print(",\"quiver_containing_path\":");
PrintBool(IsIdenticalObj(QuiverContainingPath(Q.b * Q.a), Q));
Print("},\"unique_names\":{");
Print("\"different\":");
PrintBool(unique1 <> unique2);
Print(",\"prefixes\":");
PrintBool(PositionSublist(unique1, "Quiver_") = 1 and
          PositionSublist(unique2, "Quiver_") = 1);
Print("}}\n");

QUIT_GAP(0);
