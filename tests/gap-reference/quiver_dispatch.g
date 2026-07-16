if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

Succeeded := result -> result[1];

Q := Quiver(2, [[1, 2, "a"]]);;
foreign := Quiver(1, []);;
foreign_ordering := RightLexicographicOrdering(
    foreign, GeneratorsOfQuiver(foreign));;
foreign_ordered := OrderedBy(Q, foreign_ordering);;

zero_vertices := CALL_WITH_CATCH(Quiver, [0, []]);;
mixed_rows := CALL_WITH_CATCH(Quiver,
    [2, [[1, 2], [2, 1, "b"]]]);;
indexed_endpoint := CALL_WITH_CATCH(Quiver, [2, [[1, 3, "a"]]]);;
named_endpoint := CALL_WITH_CATCH(Quiver,
    [["x", "y"], [["x", "z", "a"]]]);;
nonsquare := CALL_WITH_CATCH(Quiver, [[[0, 1], [0]]]);;
path_empty := CALL_WITH_CATCH(Path, [FamilyObj(Q.a), []]);;
path_foreign_product := CALL_WITH_CATCH(Path,
    [FamilyObj(Q.a), [Q.a, foreign.v1]]);;
extrep_empty := CALL_WITH_CATCH(ObjByExtRep,
    [FamilyObj(Q.a), []]);;
extrep_mixed_zero := CALL_WITH_CATCH(ObjByExtRep,
    [FamilyObj(Q.a), [0, 1, 3]]);;
next_zero := CALL_WITH_CATCH(NextPath, [Q, [Zero(Q)]]);;
opposite_zero := CALL_WITH_CATCH(OppositePath, [Zero(Q)]);;
foreign_order_compare := CALL_WITH_CATCH(\<,
    [foreign_ordered.v1, foreign_ordered.v2]);;

dynkin_a_rank := CALL_WITH_CATCH(DynkinQuiverAn, [0, []]);;
dynkin_a_few := CALL_WITH_CATCH(DynkinQuiverAn, [3, ["r"]]);;
dynkin_a_many := CALL_WITH_CATCH(DynkinQuiverAn,
    [2, ["r", "l"]]);;
dynkin_a_symbol := CALL_WITH_CATCH(DynkinQuiverAn,
    [2, ["u"]]);;
dynkin_d_rank := CALL_WITH_CATCH(DynkinQuiverDn,
    [3, ["r", "r"]]);;
dynkin_e_rank := CALL_WITH_CATCH(DynkinQuiverEn, [5, []]);;
dynkin_e_symbol := CALL_WITH_CATCH(DynkinQuiverEn,
    [6, ["r", "r", "r", "r", "r"]]);;
dynkin_unknown := CALL_WITH_CATCH(DynkinQuiver,
    ["X", 1, []]);;

Print("{\"schema_version\":1,\"kind\":\"quiver_dispatch\",",
      "\"source\":{",
      "\"gap\":\"4.17dev@4efe877859ab0bc0ac041bf635a996fb897ee0c9\",",
      "\"qpa\":\"1.37@8d5f847258d4853700ce3a5d947db2a1441109f1\"},",
      "\"constructor_success\":[",
      Succeeded(zero_vertices), ",",
      Succeeded(mixed_rows), ",",
      Succeeded(indexed_endpoint), ",",
      Succeeded(named_endpoint), ",",
      Succeeded(nonsquare), "],",
      "\"path_success\":[",
      Succeeded(path_empty), ",",
      Succeeded(path_foreign_product), ",",
      Succeeded(extrep_empty), ",",
      Succeeded(next_zero), ",",
      Succeeded(opposite_zero), "],",
      "\"path_empty_result_is_one\":",
      Succeeded(path_empty) and path_empty[2] = 1, ",",
      "\"extrep_mixed_zero\":{",
      "\"success\":", Succeeded(extrep_mixed_zero), ",",
      "\"is_zero\":", Succeeded(extrep_mixed_zero) and
          IsZeroPath(extrep_mixed_zero[2]), "},",
      "\"foreign_ordered_construction_success\":true,",
      "\"foreign_order_compare_success\":",
      Succeeded(foreign_order_compare), ",",
      "\"dynkin_success\":[",
      Succeeded(dynkin_a_rank), ",",
      Succeeded(dynkin_a_few), ",",
      Succeeded(dynkin_a_many), ",",
      Succeeded(dynkin_a_symbol), ",",
      Succeeded(dynkin_d_rank), ",",
      Succeeded(dynkin_e_rank), ",",
      Succeeded(dynkin_e_symbol), ",",
      Succeeded(dynkin_unknown), "]}\n");

QUIT_GAP(0);
