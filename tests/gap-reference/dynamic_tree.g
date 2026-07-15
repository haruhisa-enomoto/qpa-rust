if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

NodeItem := function(node)
    if node = 0 then
        return false;
    fi;
    return node!.item;
end;

PrintNodeItem := function(node)
    if node = 0 then
        Print("null");
    else
        PrintJsonString(node!.item);
    fi;
end;

PrintSnapshot := function(nodes)
    local first, node;
    Print("[");
    first := true;
    for node in nodes do
        if not first then
            Print(",");
        fi;
        Print("[");
        PrintNodeItem(node);
        Print(",");
        PrintNodeItem(node!.parent);
        Print(",");
        PrintNodeItem(node!.left);
        Print(",");
        PrintNodeItem(node!.right);
        Print("]");
        first := false;
    od;
    Print("]");
end;

left := CreateDynamicTree("left");;
middle := CreateDynamicTree("middle");;
right := CreateDynamicTree("right");;
middle!.left := left;;
middle!.right := right;;
left!.parent := middle;;
right!.parent := middle;;
RightRotateDynamicTree(middle);;
right_rotation := [left, middle, right];;

left2 := CreateDynamicTree("left");;
middle2 := CreateDynamicTree("middle");;
right2 := CreateDynamicTree("right");;
middle2!.left := left2;;
middle2!.right := right2;;
left2!.parent := middle2;;
right2!.parent := middle2;;
LeftRotateDynamicTree(middle2);;
left_rotation := [left2, middle2, right2];;

splice_parent := CreateDynamicTree("parent");;
splice_child := CreateDynamicTree("child");;
splice_child!.parent := splice_parent;;
SpliceDynamicTree(splice_child);;

a := CreateDynamicTree("a");;
b := CreateDynamicTree("b");;
c := CreateDynamicTree("c");;
d := CreateDynamicTree("d");;
LinkDynamicTrees(a, b);;
LinkDynamicTrees(b, c);;
LinkDynamicTrees(d, c);;
roots_before := List([a, b, c, d], node -> RootOfDynamicTree(node)!.item);;
CutDynamicTree(b);;
roots_after_cut := List([a, b, c, d], node -> RootOfDynamicTree(node)!.item);;
LinkDynamicTrees(b, d);;
roots_after_relink := List([a, b, c, d], node -> RootOfDynamicTree(node)!.item);;

Print("{\"schema_version\":1,\"kind\":\"dynamic_tree\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},\"right_rotation\":");
PrintSnapshot(right_rotation);
Print(",\"left_rotation\":");
PrintSnapshot(left_rotation);
Print(",\"splice\":");
PrintSnapshot([splice_parent, splice_child]);
Print(",\"roots_before\":");
Print("[");
PrintJsonString(roots_before[1]);
Print(","); PrintJsonString(roots_before[2]);
Print(","); PrintJsonString(roots_before[3]);
Print(","); PrintJsonString(roots_before[4]);
Print("]");
Print(",\"roots_after_cut\":");
Print("[");
PrintJsonString(roots_after_cut[1]);
Print(","); PrintJsonString(roots_after_cut[2]);
Print(","); PrintJsonString(roots_after_cut[3]);
Print(","); PrintJsonString(roots_after_cut[4]);
Print("]");
Print(",\"roots_after_relink\":");
Print("[");
PrintJsonString(roots_after_relink[1]);
Print(","); PrintJsonString(roots_after_relink[2]);
Print(","); PrintJsonString(roots_after_relink[3]);
Print(","); PrintJsonString(roots_after_relink[4]);
Print("]}\n");

QUIT_GAP(0);
