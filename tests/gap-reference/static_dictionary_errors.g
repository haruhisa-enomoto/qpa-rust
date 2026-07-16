if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

sparse := [[1], [2]];;
Unbind(sparse[1]);;
sparse_patterns := CALL_WITH_CATCH(StaticDictionary, [sparse]);;
generic := StaticDictionary([[1]]);;
empty_search := CALL_WITH_CATCH(Search, [generic, []]);;

Q := Quiver(1, [[1, 1, "x"]]);;
R := Quiver(1, [[1, 1, "y"]]);;
foreign_pattern := CALL_WITH_CATCH(QuiverStaticDictionary, [Q, [R.y]]);;
zero_pattern := CALL_WITH_CATCH(QuiverStaticDictionary, [Q, [Zero(Q)]]);;
quiver_dictionary := QuiverStaticDictionary(Q, [Q.x]);;
foreign_search := CALL_WITH_CATCH(Search, [quiver_dictionary, R.y]);;
zero_search := CALL_WITH_CATCH(Search, [quiver_dictionary, Zero(Q)]);;

Print("{\"schema_version\":1,\"source\":{",
      "\"gap\":\"4.17dev@4efe877859ab0bc0ac041bf635a996fb897ee0c9\",",
      "\"qpa\":\"1.37@8d5f847258d4853700ce3a5d947db2a1441109f1\"},",
      "\"failures\":{",
      "\"sparse_patterns_success\":", sparse_patterns[1], ",",
      "\"empty_search_success\":", empty_search[1], ",",
      "\"foreign_pattern_success\":", foreign_pattern[1], ",",
      "\"zero_pattern_success\":", zero_pattern[1], ",",
      "\"foreign_search_success\":", foreign_search[1], ",",
      "\"zero_search_success\":", zero_search[1], "}}\n");

QUIT_GAP(0);
