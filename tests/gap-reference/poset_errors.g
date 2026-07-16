if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

empty := CALL_WITH_CATCH(Poset, [[], []]);;
outside := CALL_WITH_CATCH(Poset, [[1, 2], [[1, 3]]]);;
cycle := CALL_WITH_CATCH(Poset, [[1, 2], [[1, 2], [2, 1]]]);;
poset := Poset([1, 2], [[1, 2]]);;
order := PartialOrderOfPoset(poset);;
outside_query := CALL_WITH_CATCH(order, [1, 3]);;

Print("{\"schema_version\":1,\"source\":{",
      "\"gap\":\"4.17dev@4efe877859ab0bc0ac041bf635a996fb897ee0c9\",",
      "\"qpa\":\"1.37@8d5f847258d4853700ce3a5d947db2a1441109f1\"},",
      "\"empty_success\":", empty[1], ",",
      "\"outside_success\":", outside[1], ",",
      "\"cycle_success\":", cycle[1], ",",
      "\"outside_query_success\":", outside_query[1], "}\n");

QUIT_GAP(0);
