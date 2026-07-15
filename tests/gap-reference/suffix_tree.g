if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintToken := function(token)
    if IsInt(token) then
        Print(token);
    elif IsList(token) and IsEmpty(token) then
        Print("[]");
    else
        Error("unexpected suffix-tree token");
    fi;
end;

PrintTokenList := function(tokens)
    local index;
    Print("[");
    for index in [1..Length(tokens)] do
        if index > 1 then
            Print(",");
        fi;
        PrintToken(tokens[index]);
    od;
    Print("]");
end;

PrintPairs := function(pairs)
    local index;
    Print("[");
    for index in [1..Length(pairs)] do
        if index > 1 then
            Print(",");
        fi;
        Print("[", pairs[index][1], ",", pairs[index][2], "]");
    od;
    Print("]");
end;

PrintNode := function(node)
    local edge, first;
    if IsExternalSuffixTreeNode(node) then
        Print("{\"kind\":\"external\",\"label\":");
        PrintTokenList(LabelOfSuffixTreeNode(node));
        Print(",\"occurrences\":[");
        first := true;
        for edge in node!.label do
            if not first then
                Print(",");
            fi;
            Print("[", edge[1][1], ",", edge[2], "]");
            first := false;
        od;
        Print("]}");
    else
        Print("{\"kind\":\"internal\",\"label\":");
        PrintTokenList(LabelOfSuffixTreeNode(node));
        Print(",\"marked\":");
        if node!.marked then Print("true"); else Print("false"); fi;
        Print(",\"suffix\":");
        if node!.suffixLink = 0 then
            Print("null");
        else
            PrintTokenList(LabelOfSuffixTreeNode(node!.suffixLink));
        fi;
        Print(",\"edges\":[");
        first := true;
        for edge in node!.edges do
            if not first then
                Print(",");
            fi;
            Print("{\"label\":");
            PrintTokenList(LabelOfSuffixTreeEdge(edge));
            Print(",\"child\":");
            PrintNode(edge!.child);
            Print("}");
            first := false;
        od;
        Print("]}");
    fi;
end;

PrintQueries := function(tree)
    Print("{\"all_in\":");
    PrintPairs(AllPatternsInSequence(tree, [0,1,2,1,2,1,3]));
    Print(",\"containing_12\":");
    PrintPairs(SequenceInPatterns(tree, [1,2]));
    Print(",\"containing_21\":");
    PrintPairs(SequenceInPatterns(tree, [2,1]));
    Print(",\"containing_empty\":");
    PrintPairs(SequenceInPatterns(tree, []));
    Print("}");
end;

tree := CreateSuffixTree();;
indices := [
    InsertPatternIntoSuffixTree(tree, [1,2,1]),
    InsertPatternIntoSuffixTree(tree, [2,1]),
    InsertPatternIntoSuffixTree(tree, [1,2]),
    InsertPatternIntoSuffixTree(tree, [3])
];;

Print("{\"schema_version\":1,\"kind\":\"suffix_tree\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},\"indices\":");
PrintTokenList(indices);
Print(",\"patterns\":[");
PrintTokenList(tree[1]); Print(",");
PrintTokenList(tree[2]); Print(",");
PrintTokenList(tree[3]); Print(",");
PrintTokenList(tree[4]);
Print("],\"before_delete\":");
PrintQueries(tree);
Print(",\"tree_before_delete\":");
PrintNode(tree!.root);

DeletePatternFromSuffixTree(tree, 2);;
Print(",\"after_delete_2\":");
PrintQueries(tree);
DeletePatternFromSuffixTree(tree, 1);;
Print(",\"after_delete_1\":");
PrintQueries(tree);
DeletePatternFromSuffixTree(tree, 3);;
Print(",\"after_delete_3\":");
PrintQueries(tree);
Print(",\"tree_after_delete\":");
PrintNode(tree!.root);
Print("}\n");

QUIT_GAP(0);
