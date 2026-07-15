if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\\\", "\\\\\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\\\"");
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

PrintMapOrFalse := function(map)
    if map = false then
        Print("false");
    else
        PrintMap(map);
    fi;
end;

PrintCase := function(field)
    local quiver, algebra, one, two, equal_one, equal_two, zero_module,
          zero_iso, equal_iso, different_dimension, other_quiver, other_algebra,
          other_module, different_algebra, quiver_a2, algebra_a2, rank_one,
          conjugate_rank_one, rank_two, recursive_iso, reverse_iso,
          recursive_inverse, nonisomorphic;

    quiver := Quiver(["v"], []);
    algebra := PathAlgebra(field, quiver);
    one := RightModuleOverPathAlgebra(algebra, [1], []);
    two := RightModuleOverPathAlgebra(algebra, [2], []);
    equal_one := RightModuleOverPathAlgebra(algebra, [2], []);
    equal_two := RightModuleOverPathAlgebra(algebra, [2], []);
    zero_module := ZeroModule(algebra);
    zero_iso := IsomorphismOfModules(zero_module, zero_module);
    equal_iso := IsomorphismOfModules(equal_one, equal_two);
    different_dimension := IsomorphismOfModules(one, two);

    other_quiver := Quiver(["w1", "w2"], []);
    other_algebra := PathAlgebra(field, other_quiver);
    other_module := RightModuleOverPathAlgebra(other_algebra, [1, 0], []);
    different_algebra := IsomorphismOfModules(one, other_module);

    quiver_a2 := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
    algebra_a2 := PathAlgebra(field, quiver_a2);
    rank_one := RightModuleOverPathAlgebra(
        algebra_a2, [2, 2], [["a", [[1, 0], [0, 0]] * One(field)]]);
    conjugate_rank_one := RightModuleOverPathAlgebra(
        algebra_a2, [2, 2], [["a", [[1, 1], [1, 1]] * One(field)]]);
    rank_two := RightModuleOverPathAlgebra(
        algebra_a2, [2, 2], [["a", [[1, 0], [0, 1]] * One(field)]]);
    recursive_iso := IsomorphismOfModules(rank_one, conjugate_rank_one);
    reverse_iso := IsomorphismOfModules(conjugate_rank_one, rank_one);
    recursive_inverse := InverseOfIsomorphism(recursive_iso);
    nonisomorphic := IsomorphismOfModules(rank_one, rank_two);

    Print("{\"zero\":"); PrintMap(zero_iso);
    Print(",\"equal\":"); PrintMap(equal_iso);
    Print(",\"equal_endpoints_identical\":[",
          IsIdenticalObj(Source(equal_iso), equal_one), ",",
          IsIdenticalObj(Range(equal_iso), equal_two), "]");
    Print(",\"different_dimension\":", different_dimension = false);
    Print(",\"different_algebra\":", different_algebra = false);
    Print(",\"recursive\":"); PrintMapOrFalse(recursive_iso);
    Print(",\"recursive_isomorphism\":", IsIsomorphism(recursive_iso));
    Print(",\"recursive_triangles\":[",
          recursive_iso * recursive_inverse = IdentityMapping(rank_one), ",",
          recursive_inverse * recursive_iso = IdentityMapping(conjugate_rank_one), "]");
    Print(",\"reverse\":"); PrintMapOrFalse(reverse_iso);
    Print(",\"nonisomorphic\":", nonisomorphic = false);
    Print("}");
end;

Print("{\"schema_version\":1,\"source\":{",
      "\"gap\":\"4.17dev@4efe877859ab0bc0ac041bf635a996fb897ee0c9\",",
      "\"qpa\":\"1.37@8d5f847258d4853700ce3a5d947db2a1441109f1\",",
      "\"gbnp\":\"1.1.0@e0b46c1aac16f3fbabe743e3169dc6bd86bdc654\"},",
      "\"rational\":");
PrintCase(Rationals);
Print(",\"binary\":");
PrintCase(GF(2));
Print("}\n");
