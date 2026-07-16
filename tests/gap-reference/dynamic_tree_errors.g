if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

invalid_arity := CALL_WITH_CATCH(CreateDynamicTree, [1, 2]);;
right_leaf := CreateDynamicTree("right-leaf");;
right_without_left := CALL_WITH_CATCH(RightRotateDynamicTree, [right_leaf]);;
left_leaf := CreateDynamicTree("left-leaf");;
left_without_right := CALL_WITH_CATCH(LeftRotateDynamicTree, [left_leaf]);;
parent := CreateDynamicTree("parent");;
child := CreateDynamicTree("child");;
parent!.left := child;;
child!.parent := parent;;
SetAssertionLevel(1);;
splice_auxiliary_child := CALL_WITH_CATCH(SpliceDynamicTree, [child]);;

Print("{\"schema_version\":1,\"source\":{",
      "\"gap\":\"4.17dev@4efe877859ab0bc0ac041bf635a996fb897ee0c9\",",
      "\"qpa\":\"1.37@8d5f847258d4853700ce3a5d947db2a1441109f1\"},",
      "\"failures\":{",
      "\"invalid_arity_success\":", invalid_arity[1], ",",
      "\"right_without_left_success\":", right_without_left[1], ",",
      "\"left_without_right_success\":", left_without_right[1], ",",
      "\"splice_auxiliary_child_success\":", splice_auxiliary_child[1], "}}\n");

QUIT_GAP(0);
