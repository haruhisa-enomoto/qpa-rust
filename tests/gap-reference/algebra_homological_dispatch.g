if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);;
parent := PathAlgebra(GF(2), quiver);;
nonadmissible := parent / [parent.a];;
selfinjective := CALL_WITH_CATCH(IsSelfinjectiveAlgebra, [nonadmissible]);;
hereditary := CALL_WITH_CATCH(IsHereditaryAlgebra, [nonadmissible]);;

Print("{\"schema_version\":1,\"source\":{",
      "\"gap\":\"4.17dev@4efe877859ab0bc0ac041bf635a996fb897ee0c9\",",
      "\"qpa\":\"1.37@8d5f847258d4853700ce3a5d947db2a1441109f1\"},",
      "\"nonadmissible\":", IsAdmissibleQuotientOfPathAlgebra(nonadmissible), ",",
      "\"selfinjective_success\":", selfinjective[1], ",",
      "\"hereditary_success\":", hereditary[1], "}\n");

QUIT_GAP(0);
