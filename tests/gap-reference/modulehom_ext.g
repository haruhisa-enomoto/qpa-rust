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

PrintMap := function(map)
    local first, matrix;
    Print("{\"source_dimension\":");
    PrintIntegerList(DimensionVector(Source(map)));
    Print(",\"range_dimension\":");
    PrintIntegerList(DimensionVector(Range(map)));
    Print(",\"matrices\":["); first := true;
    for matrix in MatricesOfPathAlgebraMatModuleHomomorphism(map) do
        if not first then Print(","); fi;
        PrintMatrix(matrix); first := false;
    od;
    Print("]}");
end;

PrintMapList := function(maps)
    local first, map;
    Print("["); first := true;
    for map in maps do
        if not first then Print(","); fi;
        PrintMap(map); first := false;
    od;
    Print("]");
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

PrintExt := function(ext, field)
    local basis, entered, terms, combo;
    basis := ext[2];
    Print("{\"syzygy_inclusion\":"); PrintMap(ext[1]);
    Print(",\"dimension\":", Length(basis));
    Print(",\"basis\":"); PrintMapList(basis);
    if Length(basis) = 0 then
        Print(",\"basis_coordinates\":null");
        Print(",\"entered_coefficients\":null");
        Print(",\"combo_coordinates\":null}");
    else
        Print(",\"basis_coordinates\":");
        PrintCoordinateList(List(basis, ext[3]));
        entered := List([1..Length(basis)], i -> (i + 1) * One(field));
        terms := List([1..Length(basis)], i -> entered[i] * basis[i]);
        combo := Sum(terms);
        Print(",\"entered_coefficients\":"); PrintScalarList(entered);
        Print(",\"combo_coordinates\":"); PrintScalarList(ext[3](combo));
        Print("}");
    fi;
end;

PrintRigid := function(module)
    local has_before, value, has_after;
    has_before := HasIsRigidModule(module);
    value := IsRigidModule(module);
    has_after := HasIsRigidModule(module);
    Print("{\"has_before\":", has_before);
    Print(",\"value\":", value);
    Print(",\"has_after\":", has_after, "}");
end;

PrintCase := function(field)
    local quiver, algebra, simples, zero, a2_12, a2_11, a2_22, a2_02,
          loop_quiver, parent, x, quotient, quotient_simple,
          quotient_projective, quotient_sum, dual_ss, dual_sp, dual_sum,
          kronecker_quiver,
          kronecker, kronecker_simples, kronecker_ext;

    quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
    algebra := PathAlgebra(field, quiver);
    simples := SimpleModules(algebra);
    zero := ZeroModule(algebra);
    a2_12 := ExtOverAlgebra(simples[1], simples[2]);
    a2_11 := ExtOverAlgebra(simples[1], simples[1]);
    a2_22 := ExtOverAlgebra(simples[2], simples[2]);
    a2_02 := ExtOverAlgebra(zero, simples[2]);

    loop_quiver := Quiver(["w"], [["w", "w", "x"]]);
    parent := PathAlgebra(field, loop_quiver);
    x := parent.x;
    quotient := parent / [x * x];
    quotient_simple := SimpleModules(quotient)[1];
    quotient_projective := IndecProjectiveModules(quotient)[1];
    dual_ss := ExtOverAlgebra(quotient_simple, quotient_simple);
    dual_sp := ExtOverAlgebra(quotient_simple, quotient_projective);
    quotient_sum := DirectSumOfQPAModules(
        [quotient_simple, quotient_projective]);
    dual_sum := ExtOverAlgebra(quotient_simple, quotient_sum);

    kronecker_quiver := Quiver(["k1", "k2"],
        [["k1", "k2", "p"], ["k1", "k2", "q"]]);
    kronecker := PathAlgebra(field, kronecker_quiver);
    kronecker_simples := SimpleModules(kronecker);
    kronecker_ext := ExtOverAlgebra(kronecker_simples[1],
        kronecker_simples[2]);

    Print("{\"a2_s1_s2\":"); PrintExt(a2_12, field);
    Print(",\"a2_s1_s1\":"); PrintExt(a2_11, field);
    Print(",\"a2_s2_s2\":"); PrintExt(a2_22, field);
    Print(",\"a2_zero_s2\":"); PrintExt(a2_02, field);
    Print(",\"dual_simple_simple\":"); PrintExt(dual_ss, field);
    Print(",\"dual_simple_projective\":"); PrintExt(dual_sp, field);
    Print(",\"dual_simple_sum\":"); PrintExt(dual_sum, field);
    Print(",\"kronecker_s1_s2\":"); PrintExt(kronecker_ext, field);
    Print(",\"rigid\":{\"a2_s1\":"); PrintRigid(simples[1]);
    Print(",\"a2_s2\":"); PrintRigid(simples[2]);
    Print(",\"a2_zero\":"); PrintRigid(zero);
    Print(",\"dual_simple\":"); PrintRigid(quotient_simple);
    Print(",\"dual_projective\":"); PrintRigid(quotient_projective);
    Print(",\"dual_sum\":"); PrintRigid(quotient_sum);
    Print(",\"kronecker_s1\":"); PrintRigid(kronecker_simples[1]);
    Print(",\"kronecker_s2\":"); PrintRigid(kronecker_simples[2]);
    Print("}");
    Print("}");
end;

Print("{\"schema_version\":1,\"source\":{");
Print("\"gap\":\"4.17dev@4efe877859ab0bc0ac041bf635a996fb897ee0c9\",");
Print("\"qpa\":\"1.37@8d5f847258d4853700ce3a5d947db2a1441109f1\",");
Print("\"gbnp\":\"1.1.0@e0b46c1aac16f3fbabe743e3169dc6bd86bdc654\"},");
Print("\"rational\":"); PrintCase(Rationals);
Print(",\"binary\":"); PrintCase(GF(2));
Print("}\n");
