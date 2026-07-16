if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintCaughtDecision := function(caught)
    Print("{\"success\":", caught[1]);
    if caught[1] and caught[2] = fail then
        Print(",\"value\":null");
    elif caught[1] then
        Print(",\"value\":", caught[2]);
    else
        Print(",\"value\":null");
    fi;
    Print("}");
end;

loop_quiver := Quiver(["v"], [["v", "v", "x"]]);;
loop_parent := PathAlgebra(Rationals, loop_quiver);;
dual_numbers := loop_parent / [loop_parent.x^2];;

bound_zero := CALL_WITH_CATCH(BongartzTest, [dual_numbers, 0]);;
bound_one := CALL_WITH_CATCH(BongartzTest, [dual_numbers, 1]);;
negative_bound := CALL_WITH_CATCH(BongartzTest, [dual_numbers, -1]);;
path_input := CALL_WITH_CATCH(BongartzTest, [loop_parent, 0]);;

wild_quiver := Quiver(["w"],
    [["w", "w", "x"], ["w", "w", "y"]]);;
wild_parent := PathAlgebra(GF(2), wild_quiver);;
wild_quotient := wild_parent /
    [wild_parent.x^2, wild_parent.x * wild_parent.y,
     wild_parent.y * wild_parent.x, wild_parent.y^2];;
wild := CALL_WITH_CATCH(BongartzTest, [wild_quotient, 20]);;

Print("{\"schema_version\":1,\"source\":{",
      "\"gap\":\"4.17dev@4efe877859ab0bc0ac041bf635a996fb897ee0c9\",",
      "\"qpa\":\"1.37@8d5f847258d4853700ce3a5d947db2a1441109f1\"},",
      "\"bound_zero\":"); PrintCaughtDecision(bound_zero);
Print(",\"bound_one\":"); PrintCaughtDecision(bound_one);
Print(",\"negative_bound\":"); PrintCaughtDecision(negative_bound);
Print(",\"path_input\":"); PrintCaughtDecision(path_input);
Print(",\"wild\":"); PrintCaughtDecision(wild);
Print("}\n");

QUIT_GAP(0);
