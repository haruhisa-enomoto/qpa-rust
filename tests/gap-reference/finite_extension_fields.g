if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintIntegerList := function(values)
    Print("[", JoinStringsWithSeparator(List(values, String), ","), "]");
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

PrintCase := function(order)
    local field, characteristic, degree, polynomial, quiver, algebra,
          simple_one, simple_two, projective_one, mixed, repeated,
          semisimple, submodules, endomorphism_result, loop_quiver,
          loop_parent, dual_numbers, regular, regular_result;
    field := GF(order);
    characteristic := Characteristic(field);
    degree := DegreeFFE(Z(order));
    polynomial := ConwayPolynomial(characteristic, degree);

    quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
    algebra := PathAlgebra(field, quiver);
    simple_one := RightModuleOverPathAlgebra(algebra, [1, 0], []);
    simple_two := RightModuleOverPathAlgebra(algebra, [0, 1], []);
    projective_one := IndecProjectiveModules(algebra)[1];
    mixed := DirectSumOfQPAModules([simple_one, projective_one]);
    repeated := DirectSumOfQPAModules([simple_one, simple_one]);
    semisimple := RightModuleOverPathAlgebra(
        algebra, [2, 1], [["a", [[0], [0]] * One(field)]]);
    submodules := AllSubmodulesOfModule(semisimple);
    endomorphism_result := EndOfBasicModuleAsQuiverAlgebra(
        [simple_one, projective_one, simple_two]);

    loop_quiver := Quiver(["w"], [["w", "w", "x"]]);
    loop_parent := PathAlgebra(field, loop_quiver);
    dual_numbers := loop_parent / [loop_parent.x^2];
    regular := RightModuleOverPathAlgebra(dual_numbers, [2],
        [["x", [[0, 1], [0, 0]] * One(field)]]);
    regular_result := EndOfModuleAsQuiverAlgebra(regular);

    Print("{\"order\":", order);
    Print(",\"characteristic\":", characteristic);
    Print(",\"degree\":", degree);
    Print(",\"conway\":");
    PrintIntegerList(List(CoefficientsOfUnivariatePolynomial(polynomial), IntFFE));
    Print(",\"syzygy_dimension\":");
    PrintIntegerList(DimensionVector(1stSyzygy(simple_one)));
    Print(",\"projective_dimension\":", ProjDimensionOfModule(simple_one, 3));
    Print(",\"decomposition_dimensions\":");
    PrintIntegerMatrix(List(DecomposeModule(mixed), DimensionVector));
    Print(",\"multiplicities\":");
    PrintIntegerList(DecomposeModuleWithMultiplicities(repeated)[2]);
    Print(",\"simple_submodule_count\":",
          Length(AllSimpleSubmodulesOfModule(repeated)));
    Print(",\"all_submodule_levels\":");
    PrintIntegerList(List(submodules, Length));
    Print(",\"basic_endomorphism_adjacency\":");
    PrintIntegerMatrix(AdjacencyMatrixOfQuiver(QuiverOfPathAlgebra(endomorphism_result)));
    Print(",\"basic_endomorphism_dimension\":", Dimension(endomorphism_result));
    Print(",\"regular_endomorphism_dimension\":", Dimension(regular_result[1]));
    Print(",\"regular_endomorphism_adjacency\":");
    PrintIntegerMatrix(regular_result[2]);
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"finite_extension_fields\",\"cases\":[");
PrintCase(4); Print(","); PrintCase(8); Print(","); PrintCase(9);
Print("]}\n");

QUIT_GAP(0);
