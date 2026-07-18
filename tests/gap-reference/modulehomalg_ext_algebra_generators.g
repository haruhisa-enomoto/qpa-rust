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
        Print(Int(value)); first := false;
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

PrintGeneratorDegrees := function(degrees)
    local first, maps;
    Print("["); first := true;
    for maps in degrees do
        if not first then Print(","); fi;
        PrintMapList(maps); first := false;
    od;
    Print("]");
end;

PrintExtAlgebra := function(result)
    Print("{\"dimensions\":"); PrintIntegerList(result[1]);
    Print(",\"generator_counts\":"); PrintIntegerList(result[2]);
    Print(",\"generators\":"); PrintGeneratorDegrees(result[3]);
    Print("}");
end;

PrintCase := function(field)
    local a2_quiver, a2, a2_sum, vector_quiver, vector_algebra, vector2,
          loop_quiver, parent, x, dual, dual_simple, dual_projective,
          dual_sum, truncated, truncated_projective;

    a2_quiver := Quiver(["a1", "a2"], [["a1", "a2", "a"]]);
    a2 := PathAlgebra(field, a2_quiver);
    a2_sum := DirectSumOfQPAModules(SimpleModules(a2));

    vector_quiver := Quiver(["v"], []);
    vector_algebra := PathAlgebra(field, vector_quiver);
    vector2 := RightModuleOverPathAlgebra(vector_algebra, [2], []);

    loop_quiver := Quiver(["w"], [["w", "w", "x"]]);
    parent := PathAlgebra(field, loop_quiver);
    x := parent.x;
    dual := parent / [x^2];
    dual_simple := SimpleModules(dual)[1];
    dual_projective := IndecProjectiveModules(dual)[1];
    dual_sum := DirectSumOfQPAModules([dual_simple, dual_projective]);

    truncated := parent / [x^3];
    truncated_projective := IndecProjectiveModules(truncated)[1];

    Print("{\"a2_sum_n2\":");
    PrintExtAlgebra(ExtAlgebraGenerators(a2_sum, 2));
    Print(",\"vector2_n2\":");
    PrintExtAlgebra(ExtAlgebraGenerators(vector2, 2));
    Print(",\"dual_simple_n3\":");
    PrintExtAlgebra(ExtAlgebraGenerators(dual_simple, 3));
    Print(",\"dual_sum_n3\":");
    PrintExtAlgebra(ExtAlgebraGenerators(dual_sum, 3));
    Print(",\"truncated_projective_n2\":");
    PrintExtAlgebra(ExtAlgebraGenerators(truncated_projective, 2));
    Print("}");
end;

Print("{\"schema_version\":1,\"source\":{");
Print("\"gap\":\"4.17dev@4efe877859ab0bc0ac041bf635a996fb897ee0c9\",");
Print("\"qpa\":\"1.37@8d5f847258d4853700ce3a5d947db2a1441109f1\",");
Print("\"gbnp\":\"1.1.0@e0b46c1aac16f3fbabe743e3169dc6bd86bdc654\"},");
Print("\"rational\":"); PrintCase(Rationals);
Print(",\"binary\":"); PrintCase(GF(2));
Print("}\n");
