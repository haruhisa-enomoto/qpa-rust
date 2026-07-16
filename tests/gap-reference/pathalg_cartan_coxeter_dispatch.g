if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

ReturnsFail := function(operation, algebra)
    local result;
    result := CALL_WITH_CATCH(operation, [algebra]);
    return result[1] and result[2] = fail;
end;

loop_path := PathAlgebra(GF(2),
    Quiver(["l"], [["l", "l", "t"]]));;
loop_cartan_initially_known := HasCartanMatrix(loop_path);;
loop_cartan_fail := ReturnsFail(CartanMatrix, loop_path);;
loop_coxeter_fail := ReturnsFail(CoxeterMatrix, loop_path);;
loop_polynomial_fail := ReturnsFail(CoxeterPolynomial, loop_path);;

nonadmissible_parent := PathAlgebra(GF(2),
    Quiver(["v1", "v2"], [["v1", "v2", "a"]]));;
nonadmissible := nonadmissible_parent / [nonadmissible_parent.a];;
nonadmissible_cartan_initially_known := HasCartanMatrix(nonadmissible);;
nonadmissible_cartan_fail := ReturnsFail(CartanMatrix, nonadmissible);;
nonadmissible_coxeter_fail := ReturnsFail(CoxeterMatrix, nonadmissible);;
nonadmissible_polynomial_fail := ReturnsFail(CoxeterPolynomial, nonadmissible);;

cycle_parent := PathAlgebra(GF(2),
    Quiver(["c1", "c2"],
        [["c1", "c2", "x"], ["c2", "c1", "y"]]));;
singular := cycle_parent /
    [cycle_parent.x * cycle_parent.y,
     cycle_parent.y * cycle_parent.x];;
singular_cartan_initially_known := HasCartanMatrix(singular);;
singular_cartan := CartanMatrix(singular);;
singular_cartan_known_after := HasCartanMatrix(singular);;
singular_matrix_initially_known := HasCoxeterMatrix(singular);;
singular_coxeter_fail := ReturnsFail(CoxeterMatrix, singular);;
singular_matrix_known_after := HasCoxeterMatrix(singular);;
singular_polynomial_initially_known := HasCoxeterPolynomial(singular);;
singular_polynomial_fail := ReturnsFail(CoxeterPolynomial, singular);;
singular_polynomial_known_after := HasCoxeterPolynomial(singular);;

Print("{\"schema_version\":1,\"source\":{",
      "\"gap\":\"4.17dev@4efe877859ab0bc0ac041bf635a996fb897ee0c9\",",
      "\"qpa\":\"1.37@8d5f847258d4853700ce3a5d947db2a1441109f1\"},",
      "\"loop\":{",
      "\"cartan_initially_known\":", loop_cartan_initially_known, ",",
      "\"cartan_fail\":", loop_cartan_fail, ",",
      "\"cartan_known_after\":", HasCartanMatrix(loop_path), ",",
      "\"coxeter_fail\":", loop_coxeter_fail, ",",
      "\"polynomial_fail\":", loop_polynomial_fail, ",",
      "\"matrix_known_after\":", HasCoxeterMatrix(loop_path), ",",
      "\"polynomial_known_after\":", HasCoxeterPolynomial(loop_path), "},",
      "\"nonadmissible\":{",
      "\"is_admissible\":", IsAdmissibleQuotientOfPathAlgebra(nonadmissible), ",",
      "\"cartan_initially_known\":", nonadmissible_cartan_initially_known, ",",
      "\"cartan_fail\":", nonadmissible_cartan_fail, ",",
      "\"cartan_known_after\":", HasCartanMatrix(nonadmissible), ",",
      "\"coxeter_fail\":", nonadmissible_coxeter_fail, ",",
      "\"polynomial_fail\":", nonadmissible_polynomial_fail, ",",
      "\"matrix_known_after\":", HasCoxeterMatrix(nonadmissible), ",",
      "\"polynomial_known_after\":", HasCoxeterPolynomial(nonadmissible), "},",
      "\"singular\":{",
      "\"cartan_initially_known\":", singular_cartan_initially_known, ",",
      "\"cartan\":", singular_cartan, ",",
      "\"cartan_known_after\":", singular_cartan_known_after, ",",
      "\"matrix_initially_known\":", singular_matrix_initially_known, ",",
      "\"coxeter_fail\":", singular_coxeter_fail, ",",
      "\"matrix_known_after\":", singular_matrix_known_after, ",",
      "\"polynomial_initially_known\":", singular_polynomial_initially_known, ",",
      "\"polynomial_fail\":", singular_polynomial_fail, ",",
      "\"polynomial_known_after\":", singular_polynomial_known_after, "}}\n");

QUIT_GAP(0);
