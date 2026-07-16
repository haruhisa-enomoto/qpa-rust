if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

cycle_parent := PathAlgebra(GF(2),
    Quiver(["c1", "c2", "c3"],
        [["c1", "c2", "x"], ["c2", "c3", "y"],
         ["c3", "c1", "z"]]));;
cycle := cycle_parent /
    [cycle_parent.x * cycle_parent.y,
     cycle_parent.y * cycle_parent.z,
     cycle_parent.z * cycle_parent.x];;
cycle_permutation := NakayamaPermutation(cycle);;
cycle_simples := SimpleModules(cycle);;
cycle_nonsimple := DirectSumOfQPAModules(cycle_simples{[1, 2]});;
cycle_nonsimple_result :=
    CALL_WITH_CATCH(cycle_permutation[1], [cycle_nonsimple]);;
cycle_zero_index_result := CALL_WITH_CATCH(cycle_permutation[2], [0]);;
cycle_negative_index_result := CALL_WITH_CATCH(cycle_permutation[2], [-1]);;
cycle_large_index_result := CALL_WITH_CATCH(cycle_permutation[2], [4]);;

foreign := PathAlgebra(GF(2), Quiver(["f1", "f2", "f3", "f4"], []));;
foreign_simple := SimpleModules(foreign)[4];;
foreign_result := CALL_WITH_CATCH(cycle_permutation[1], [foreign_simple]);;

discrete := PathAlgebra(GF(2), Quiver(["d1", "d2"], []));;
discrete_permutation := NakayamaPermutation(discrete);;
discrete_nonsimple :=
    DirectSumOfQPAModules(SimpleModules(discrete));;
discrete_nonsimple_result :=
    CALL_WITH_CATCH(discrete_permutation[1], [discrete_nonsimple]);;
discrete_zero_index_result :=
    CALL_WITH_CATCH(discrete_permutation[2], [0]);;

nonadmissible_parent := PathAlgebra(GF(2),
    Quiver(["v1", "v2"], [["v1", "v2", "a"]]));;
nonadmissible := nonadmissible_parent / [nonadmissible_parent.a];;
nonadmissible_result := CALL_WITH_CATCH(NakayamaPermutation, [nonadmissible]);;

Print("{\"schema_version\":1,\"source\":{",
      "\"gap\":\"4.17dev@4efe877859ab0bc0ac041bf635a996fb897ee0c9\",",
      "\"qpa\":\"1.37@8d5f847258d4853700ce3a5d947db2a1441109f1\"},",
      "\"cycle_nonsimple_success\":", cycle_nonsimple_result[1], ",",
      "\"cycle_zero_index_success\":", cycle_zero_index_result[1], ",",
      "\"cycle_negative_index_success\":", cycle_negative_index_result[1], ",",
      "\"cycle_large_index_success\":", cycle_large_index_result[1], ",",
      "\"foreign_simple_success\":", foreign_result[1], ",",
      "\"discrete_nonsimple_success\":", discrete_nonsimple_result[1], ",",
      "\"discrete_zero_index_success\":", discrete_zero_index_result[1], ",",
      "\"nonadmissible_success\":", nonadmissible_result[1], "}\n");

QUIT_GAP(0);
