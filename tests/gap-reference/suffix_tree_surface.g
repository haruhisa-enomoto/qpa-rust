if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

tree := CreateSuffixTree();;
pattern := [1, [1, 2, 3, []]];;
Add(tree!.patterns, pattern);;
family := tree!.family;;
root := tree!.root;;

internal := CreateInternalSuffixTreeNode(family, pattern, 1, 1);;
external := CreateExternalSuffixTreeNode(family, pattern, 2);;
edge := CreateSuffixTreeEdge(family, internal, external, pattern, 2, -1);;

split_parent := CreateInternalSuffixTreeNode(family, pattern, 1, 1);;
split_child := CreateExternalSuffixTreeNode(family, pattern, 1);;
split_edge := CreateSuffixTreeEdge(
    family, split_parent, split_child, pattern, 1, 2);;
split_node := CreateInternalSuffixTreeNode(family, pattern, 1, 1);;
AddInternalSuffixTreeNode(split_edge, 1, split_node);;
split_remainder := split_node!.edges[1];;

added_external := CreateExternalSuffixTreeNode(family, pattern, 1);;
AddExternalSuffixTreeNode(root, added_external);;
added_edge := added_external!.parentEdge;;

other_tree := CreateSuffixTree();;
InsertPatternIntoSuffixTree(other_tree, [1, 2, 3]);;
other_edge := other_tree!.root!.edges[1];;

Print("{\"schema_version\":1,\"kind\":\"suffix_tree_surface\",",
      "\"oracle\":{\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},",
      "\"info_class\":", IsInfoClass(InfoDynamicDictionary), ",",
      "\"tree_filters\":[", IsSuffixTree(tree), ",",
          IsDefaultSuffixTreeRep(tree), "],",
      "\"root_filters\":[", IsSuffixTreeElement(root), ",",
          IsSuffixTreeNode(root), ",", IsInternalSuffixTreeNode(root), ",",
          IsExternalSuffixTreeNode(root), ",",
          IsDefaultSuffixTreeNodeRep(root), ",",
          IsDefaultInternalSuffixTreeNodeRep(root), "],",
      "\"external_filters\":[", IsSuffixTreeElement(external), ",",
          IsSuffixTreeNode(external), ",", IsInternalSuffixTreeNode(external), ",",
          IsExternalSuffixTreeNode(external), ",",
          IsDefaultSuffixTreeNodeRep(external), "],",
      "\"edge_filters\":[", IsSuffixTreeElement(edge), ",",
          IsSuffixTreeEdge(edge), ",", IsDefaultSuffixTreeEdgeRep(edge), "],",
      "\"components\":{",
          "\"dead_length\":", tree!.deadLength, ",",
          "\"pattern_count\":", Length(tree!.patterns), ",",
          "\"pattern_1\":", tree[1], "},",
      "\"manual_edge\":{",
          "\"internal_label\":", LabelOfSuffixTreeNode(internal), ",",
          "\"external_label\":", LabelOfSuffixTreeNode(external), ",",
          "\"edge_label\":", LabelOfSuffixTreeEdge(edge), ",",
          "\"root_parent_is_zero\":", ParentNode(root) = 0, ",",
          "\"child_parent_identity\":", IsIdenticalObj(ParentNode(external), internal), ",",
          "\"edge_lookup_identity\":",
              IsIdenticalObj(SuffixTreeEdgeStartingAt(internal, 2), edge), ",",
          "\"missing_lookup_is_fail\":",
              SuffixTreeEdgeStartingAt(internal, 9) = fail, ",",
          "\"self_equal\":", edge = edge, ",",
          "\"distinct_equal\":", edge = split_edge, ",",
          "\"other_tree_root_equal\":", root = other_tree!.root, ",",
          "\"other_tree_edge_equal\":", edge = other_edge, "},",
      "\"split\":{",
          "\"prefix\":", LabelOfSuffixTreeEdge(split_edge), ",",
          "\"remainder\":", LabelOfSuffixTreeEdge(split_remainder), ",",
          "\"old_child_parent_identity\":",
              IsIdenticalObj(ParentNode(split_child), split_node), "},",
      "\"add_external\":{",
          "\"edge_label\":", LabelOfSuffixTreeEdge(added_edge), ",",
          "\"parent_identity\":",
              IsIdenticalObj(ParentNode(added_external), root), "}}\n");

QUIT_GAP(0);
