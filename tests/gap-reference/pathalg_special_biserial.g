if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

Predicates := function(A)
    return [IsSpecialBiserialAlgebra(A), IsStringAlgebra(A), IsGentleAlgebra(A)];
end;

PrintFieldPredicates := function(K)
    local Q, A, q1, q2, q3, q4, L, square, cubic, linear, minimal,
          gorenstein_before, gentle, gorenstein_after;

    Q := Quiver(1, [[1,1,"a"], [1,1,"b"]]);
    A := PathAlgebra(K, Q);
    q1 := A/[A.a*A.b, A.a^2, A.b^2];
    q2 := A/[A.a*A.b-A.b*A.a, A.a^2, A.b^2];
    q3 := A/[A.a*A.b+A.b*A.a, A.a^2, A.b^2, A.b*A.a];
    q4 := A/[A.a*A.b, A.a^2, A.b^3];

    L := PathAlgebra(K, Quiver(1, [[1,1,"x"]]));
    square := L/[L.x^2];
    cubic := L/[L.x^3];
    linear := L/[L.x];

    gorenstein_before := HasIsGorensteinAlgebra(square);
    gentle := IsGentleAlgebra(square);
    gorenstein_after := HasIsGorensteinAlgebra(square) and IsGorensteinAlgebra(square);

    Print("{\"doc_q1\":", Predicates(q1));
    minimal := MinimalGeneratingSetOfIdeal(ElementsFamily(FamilyObj(q1))!.ideal);
    minimal := List(minimal, relation ->
        List(WalkOfPath(CoefficientsAndMagmaElements(relation)[1]), String));
    Print(",\"doc_q1_minimal_walks\":", minimal);
    Print(",\"doc_q2\":", Predicates(q2));
    Print(",\"doc_q3\":", Predicates(q3));
    Print(",\"doc_q4\":", Predicates(q4));
    Print(",\"square_loop\":", Predicates(square));
    Print(",\"cubic_loop\":", Predicates(cubic));
    Print(",\"linear_loop\":", Predicates(linear));
    Print(",\"gentle_implication\":[", gorenstein_before, ",", gentle,
          ",", gorenstein_after, "]}");
end;

one_loop := Quiver(1, [[1,1,"x"]]);
three_out := Quiver(4, [[1,2,"a"], [1,3,"b"], [1,4,"c"]]);
chain := Quiver(3, [[1,2,"a"], [2,3,"b"]]);
chain_path := PathAlgebra(Rationals, chain);
loop_path := PathAlgebra(Rationals, one_loop);

Print("{\"schema_version\":1,\"kind\":\"pathalg_special_biserial\",\"oracle\":{");
Print("\"gap_version\":\"", GAPInfo.Version, "\",");
Print("\"gap_revision\":\"",
    GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION, "\",");
Print("\"qpa_version\":\"", PackageInfo("qpa")[1].Version, "\",");
Print("\"qpa_revision\":\"",
    GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION, "\"},");
Print("\"quivers\":{\"one_loop\":", IsSpecialBiserialQuiver(one_loop));
Print(",\"three_out\":", IsSpecialBiserialQuiver(three_out), "},");
Print("\"paths\":{\"chain\":", Predicates(chain_path));
Print(",\"one_loop\":", Predicates(loop_path), "},");
Print("\"rational\":");
PrintFieldPredicates(Rationals);
Print(",\"f2\":");
PrintFieldPredicates(GF(2));
Print("}\n");

QUIT_GAP(0);
