if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

loop_path := PathAlgebra(GF(2),
    Quiver(["l"], [["l", "l", "t"]]));;
loop_schurian := CALL_WITH_CATCH(IsSchurianAlgebra, [loop_path]);;
loop_semicommutative := CALL_WITH_CATCH(IsSemicommutativeAlgebra, [loop_path]);;

parent := PathAlgebra(GF(2),
    Quiver(["v1", "v2"], [["v1", "v2", "a"]]));;
nonadmissible := parent / [parent.a];;
nonadmissible_symmetric := CALL_WITH_CATCH(IsSymmetricAlgebra, [nonadmissible]);;
nonadmissible_weak := CALL_WITH_CATCH(IsWeaklySymmetricAlgebra, [nonadmissible]);;
nonadmissible_schurian := CALL_WITH_CATCH(IsSchurianAlgebra, [nonadmissible]);;
nonadmissible_semicommutative :=
    CALL_WITH_CATCH(IsSemicommutativeAlgebra, [nonadmissible]);;

Print("{\"schema_version\":1,\"source\":{",
      "\"gap\":\"4.17dev@4efe877859ab0bc0ac041bf635a996fb897ee0c9\",",
      "\"qpa\":\"1.37@8d5f847258d4853700ce3a5d947db2a1441109f1\"},",
      "\"loop_schurian_success\":", loop_schurian[1], ",",
      "\"loop_semicommutative_success\":", loop_semicommutative[1], ",",
      "\"nonadmissible\":", IsAdmissibleQuotientOfPathAlgebra(nonadmissible), ",",
      "\"nonadmissible_symmetric_success\":", nonadmissible_symmetric[1], ",",
      "\"nonadmissible_weak_success\":", nonadmissible_weak[1], ",",
      "\"nonadmissible_schurian_success\":", nonadmissible_schurian[1], ",",
      "\"nonadmissible_semicommutative_success\":",
      nonadmissible_semicommutative[1], "}\n");

QUIT_GAP(0);
