if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintIntegerList := function(values)
    local first, value;
    Print("["); first := true;
    for value in values do
        if not first then Print(","); fi;
        Print(value); first := false;
    od;
    Print("]");
end;

PrintIntegerMatrix := function(matrix)
    local first, row;
    Print("["); first := true;
    for row in matrix do
        if not first then Print(","); fi;
        PrintIntegerList(row); first := false;
    od;
    Print("]");
end;

PrintScalarList := function(values)
    local first, value;
    Print("["); first := true;
    for value in values do
        if not first then Print(","); fi;
        PrintJsonString(String(value)); first := false;
    od;
    Print("]");
end;

PrintScalarMatrix := function(matrix)
    local first, row;
    Print("["); first := true;
    for row in matrix do
        if not first then Print(","); fi;
        PrintScalarList(row); first := false;
    od;
    Print("]");
end;

PrintCase := function(algebra)
    local cartan_initially_known, matrix_initially_known, polynomial_initially_known,
          cartan, coxeter, polynomial;

    cartan_initially_known := HasCartanMatrix(algebra);
    matrix_initially_known := HasCoxeterMatrix(algebra);
    polynomial_initially_known := HasCoxeterPolynomial(algebra);
    cartan := CartanMatrix(algebra);
    coxeter := CoxeterMatrix(algebra);
    polynomial := CoxeterPolynomial(algebra);

    Print("{\"cartan_initially_known\":", cartan_initially_known);
    Print(",\"matrix_initially_known\":", matrix_initially_known);
    Print(",\"polynomial_initially_known\":", polynomial_initially_known);
    Print(",\"cartan\":"); PrintIntegerMatrix(cartan);
    Print(",\"coxeter\":"); PrintScalarMatrix(coxeter);
    Print(",\"polynomial_coefficients\":");
    PrintScalarList(CoefficientsOfUnivariatePolynomial(polynomial));
    Print(",\"cartan_known_after\":", HasCartanMatrix(algebra));
    Print(",\"matrix_known_after\":", HasCoxeterMatrix(algebra));
    Print(",\"polynomial_known_after\":", HasCoxeterPolynomial(algebra));
    Print("}");
end;

PrintCases := function(field)
    local discrete, a3, loop_parent, dual_numbers, diamond_parent,
          commutative_diamond;

    discrete := PathAlgebra(field, Quiver(["d1", "d2"], []));
    a3 := PathAlgebra(field,
        Quiver(["a1", "a2", "a3"],
            [["a1", "a2", "x"], ["a2", "a3", "y"]]));
    loop_parent := PathAlgebra(field,
        Quiver(["l"], [["l", "l", "t"]]));
    dual_numbers := loop_parent / [loop_parent.t * loop_parent.t];
    diamond_parent := PathAlgebra(field,
        Quiver(["q1", "q2", "q3", "q4"],
            [["q1", "q2", "a"], ["q1", "q3", "b"],
             ["q2", "q4", "c"], ["q3", "q4", "d"]]));
    commutative_diamond := diamond_parent /
        [diamond_parent.a * diamond_parent.c
         - diamond_parent.b * diamond_parent.d];

    Print("{\"path_discrete\":"); PrintCase(discrete);
    Print(",\"path_a3\":"); PrintCase(a3);
    Print(",\"quotient_dual_numbers\":"); PrintCase(dual_numbers);
    Print(",\"quotient_commutative_diamond\":");
    PrintCase(commutative_diamond);
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"pathalg_cartan_coxeter\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"f2\":"); PrintCases(GF(2));
Print(",\"f5\":"); PrintCases(GF(5));
Print(",\"rational\":"); PrintCases(Rationals);
Print("}\n");

QUIT_GAP(0);
