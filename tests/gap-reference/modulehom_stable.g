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

PrintScalarList := function(values)
    local first, value;
    Print("["); first := true;
    for value in values do
        if not first then Print(","); fi;
        PrintJsonString(String(value)); first := false;
    od;
    Print("]");
end;

PrintMatrix := function(matrix)
    local first, row;
    Print("{\"shape\":"); PrintIntegerList(DimensionsMat(matrix));
    Print(",\"rows\":["); first := true;
    for row in matrix do
        if not first then Print(","); fi;
        PrintScalarList(row); first := false;
    od;
    Print("]}");
end;

PrintMatrixList := function(matrices)
    local first, matrix;
    Print("["); first := true;
    for matrix in matrices do
        if not first then Print(","); fi;
        PrintMatrix(matrix); first := false;
    od;
    Print("]");
end;

PrintMap := function(map)
    Print("{\"source_dimension\":");
    PrintIntegerList(DimensionVector(Source(map)));
    Print(",\"range_dimension\":");
    PrintIntegerList(DimensionVector(Range(map)));
    Print(",\"matrices\":");
    PrintMatrixList(MatricesOfPathAlgebraMatModuleHomomorphism(map));
    Print("}");
end;

PrintCoordinateList := function(vectors)
    local first, vector;
    Print("["); first := true;
    for vector in vectors do
        if not first then Print(","); fi;
        PrintScalarList(vector); first := false;
    od;
    Print("]");
end;

PrintMultiplicationTable := function(table)
    local first_row, row;
    Print("["); first_row := true;
    for row in table do
        if not first_row then Print(","); fi;
        PrintCoordinateList(row); first_row := false;
    od;
    Print("]");
end;

PrintStable := function(module, print_zero_vertex_map)
    local endomorphisms, end_basis, factoring, factor_matrices, ideal,
          natural, stable, stable_basis, images, table;
    endomorphisms := EndOverAlgebra(module);
    end_basis := BasisVectors(Basis(endomorphisms));
    factoring := HomFactoringThroughProjOverAlgebra(module, module);
    factor_matrices := List(factoring, FromHomMMToEndM);
    ideal := Ideal(endomorphisms, factor_matrices);
    natural := EndModuloProjOverAlgebra(module);
    stable := Range(natural);
    stable_basis := BasisVectors(Basis(stable));
    images := List(end_basis,
        x -> Coefficients(Basis(stable), Image(natural, x)));
    table := List(stable_basis,
        x -> List(stable_basis,
            y -> Coefficients(Basis(stable), x * y)));

    Print("{\"end_basis\":"); PrintMatrixList(end_basis);
    Print(",\"factor_matrices\":"); PrintMatrixList(factor_matrices);
    Print(",\"ideal_basis\":"); PrintMatrixList(BasisVectors(Basis(ideal)));
    Print(",\"kernel_basis\":");
    PrintMatrixList(BasisVectors(Basis(Kernel(natural))));
    Print(",\"source_cache_identity\":",
          IsIdenticalObj(Source(natural), EndOverAlgebra(module)));
    Print(",\"dimension\":", Dimension(stable));
    Print(",\"images\":"); PrintCoordinateList(images);
    Print(",\"one\":");
    if Dimension(stable) = 0 then
        Print("null");
    else
        PrintScalarList(Coefficients(Basis(stable), One(stable)));
    fi;
    Print(",\"multiplication\":"); PrintMultiplicationTable(table);
    Print(",\"roundtrip_matrices\":");
    PrintMatrixList(List(end_basis,
        x -> FromHomMMToEndM(FromEndMToHomMM(module, x))));
    if print_zero_vertex_map then
        Print(",\"zero_vertex_map\":");
        PrintMap(FromEndMToHomMM(module, end_basis[1]));
    fi;
    Print("}");
end;

PrintCase := function(field)
    local quiver, algebra, module, simple_one, simple_two, loop_quiver,
          path_algebra, loop, dual_numbers, dual_simple, dual_projective;
    quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
    algebra := PathAlgebra(field, quiver);
    module := RightModuleOverPathAlgebra(
        algebra, [2, 2], [["a", [[1,0],[0,0]] * One(field)]]);
    simple_one := RightModuleOverPathAlgebra(algebra, [1, 0], []);
    simple_two := RightModuleOverPathAlgebra(algebra, [0, 1], []);

    loop_quiver := Quiver(["w"], [["w", "w", "x"]]);
    path_algebra := PathAlgebra(field, loop_quiver);
    loop := ArrowsOfQuiver(loop_quiver)[1] * One(path_algebra);
    dual_numbers := path_algebra / [loop * loop];
    dual_simple := RightModuleOverPathAlgebra(dual_numbers, [1], []);
    dual_projective := IndecProjectiveModules(dual_numbers)[1];

    Print("{\"a2_module\":"); PrintStable(module, false);
    Print(",\"a2_simple_one\":"); PrintStable(simple_one, true);
    Print(",\"a2_simple_two\":"); PrintStable(simple_two, true);
    Print(",\"dual_simple\":"); PrintStable(dual_simple, false);
    Print(",\"dual_projective\":"); PrintStable(dual_projective, false);
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"modulehom_stable\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("GBNP")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"rational\":"); PrintCase(Rationals);
Print(",\"binary\":"); PrintCase(GF(2));
Print("}\n");

QUIT_GAP(0);
