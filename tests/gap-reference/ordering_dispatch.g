if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

Succeeded := result -> result[1];

Q := Quiver(2, [[1, 1, "a"], [1, 2, "b"]]);;
a := Q.a;;
b := Q.b;;
default := OrderingOfQuiver(Q);;

empty_block := CALL_WITH_CATCH(BlockOrdering, [Q, []]);;
empty_wreath := CALL_WITH_CATCH(WreathOrdering, [Q, []]);;
conflicting_vector := CALL_WITH_CATCH(LeftVectorOrdering,
    [Q, [b, a], default]);;
conflicting_block := CALL_WITH_CATCH(BlockOrdering,
    [Q, [[default, [a]], [LeftLexicographicOrdering(Q, [b, a]), [b]]]]);;
lex_well := CALL_WITH_CATCH(IsWellOrdering,
    [LeftLexicographicOrdering(Q)]);;
weight_next := WeightOrdering(Q, [0, 0], default);;
weight_total_subset := CALL_WITH_CATCH(IsTotalSubset, [weight_next, [a]]);;
wreath := WreathOrdering(Q, [[default, GeneratorsOfQuiver(Q)]]);;
wreath_subset := CALL_WITH_CATCH(IsWellSubset,
    [wreath, GeneratorsOfQuiver(Q)]);;

QT := Quiver(1, [[1, 1, "x"], [1, 1, "y"]]);;
td := OrderingOfQuiver(QT);;
two_loop_block := BlockOrdering(QT, [[td, [QT.x]], [td, [QT.y]]]);;
block_total := CALL_WITH_CATCH(IsTotalOrdering, [two_loop_block]);;
short_blocks := CALL_WITH_CATCH(CheckBlocksForCycles,
    [LengthOrdering(Q), [1]]);;

R := Quiver(1, [[1, 1, "z"]]);;
foreign_default := OrderingOfQuiver(R);;
foreign_length := CALL_WITH_CATCH(LengthOrdering, [Q, foreign_default]);;
foreign_reverse := CALL_WITH_CATCH(ReverseOrdering, [Q, foreign_default]);;
foreign_weight := CALL_WITH_CATCH(WeightOrdering,
    [Q, [1, 2], foreign_default]);;
heterogeneous_path := CALL_WITH_CATCH(PathFromArrowList,
    [Q, [Q.v1, R.v1]]);;
heterogeneous_compare := CALL_WITH_CATCH(CompareArrowLists,
    [default, [Q.a, R.z], [Q.a]]);;

Print("{\"schema_version\":1,\"kind\":\"ordering_dispatch\",",
      "\"source\":{",
      "\"gap\":\"4.17dev@4efe877859ab0bc0ac041bf635a996fb897ee0c9\",",
      "\"qpa\":\"1.37@8d5f847258d4853700ce3a5d947db2a1441109f1\"},",
      "\"success\":{",
      "\"empty_block\":", Succeeded(empty_block), ",",
      "\"empty_wreath\":", Succeeded(empty_wreath), ",",
      "\"conflicting_vector\":", Succeeded(conflicting_vector), ",",
      "\"conflicting_block\":", Succeeded(conflicting_block), ",",
      "\"lex_well\":", Succeeded(lex_well), ",",
      "\"weight_total_subset\":", Succeeded(weight_total_subset), ",",
      "\"wreath_subset\":", Succeeded(wreath_subset), ",",
      "\"block_total\":", Succeeded(block_total), ",",
      "\"short_blocks\":", Succeeded(short_blocks), ",",
      "\"foreign_length\":", Succeeded(foreign_length), ",",
      "\"foreign_reverse\":", Succeeded(foreign_reverse), ",",
      "\"foreign_weight\":", Succeeded(foreign_weight), ",",
      "\"heterogeneous_path\":", Succeeded(heterogeneous_path), ",",
      "\"heterogeneous_compare\":", Succeeded(heterogeneous_compare), "}}\n");

QUIT_GAP(0);
