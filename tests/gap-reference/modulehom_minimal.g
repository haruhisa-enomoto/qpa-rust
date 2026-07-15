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

PrintMore := function(result, original)
    if IsList(result) then
        Print("{\"changed\":true,\"minimal_part\":");
        PrintMap(result[1]);
        Print(",\"removed_part\":"); PrintMap(result[2]);
        Print("}");
    else
        Print("{\"changed\":false,\"same_object\":",
              IsIdenticalObj(result, original), ",\"map\":");
        PrintMap(result);
        Print("}");
    fi;
end;

PrintModuleDimensions := function(modules)
    local first, module;
    Print("["); first := true;
    for module in modules do
        if not first then Print(","); fi;
        PrintIntegerList(DimensionVector(module)); first := false;
    od;
    Print("]");
end;

PrintMinimal := function(result, original, side)
    Print("{\"map\":"); PrintMap(result[1]);
    Print(",\"removed_dimensions\":"); PrintModuleDimensions(result[2]);
    Print(",\"map_is_original\":", IsIdenticalObj(result[1], original));
    if side = "right" then
        Print(",\"property_bound\":", HasIsRightMinimal(result[1]));
        Print(",\"property_value\":", IsRightMinimal(result[1]));
    else
        Print(",\"property_bound\":", HasIsLeftMinimal(result[1]));
        Print(",\"property_value\":", IsLeftMinimal(result[1]));
    fi;
    Print("}");
end;

PrintCase := function(field)
    local quiver, algebra, one, three, right_more_map, right_more,
          left_more_map, left_more, right_map, right_result, right_again,
          left_map, left_result, left_again, identity, right_identity,
          left_identity, right_more_unchanged, left_more_unchanged, zero_map,
          right_zero, quiver_a2, algebra_a2, projective_one, simple_two,
          direct_sum, projection, inclusion, a2_right, a2_left;

    quiver := Quiver(["v"], []);
    algebra := PathAlgebra(field, quiver);
    one := RightModuleOverPathAlgebra(algebra, [1], []);
    three := RightModuleOverPathAlgebra(algebra, [3], []);

    right_more_map := RightModuleHomOverAlgebra(
        three, one, [[[1], [0], [0]] * One(field)]);
    right_more := MoreRightMinimalVersion(right_more_map);
    left_more_map := RightModuleHomOverAlgebra(
        one, three, [[[1, 0, 0]] * One(field)]);
    left_more := MoreLeftMinimalVersion(left_more_map);

    right_map := RightModuleHomOverAlgebra(
        three, one, [[[1], [0], [0]] * One(field)]);
    right_result := RightMinimalVersion(right_map);
    right_again := RightMinimalVersion(right_map);
    left_map := RightModuleHomOverAlgebra(
        one, three, [[[1, 0, 0]] * One(field)]);
    left_result := LeftMinimalVersion(left_map);
    left_again := LeftMinimalVersion(left_map);
    identity := IdentityMapping(one);
    right_more_unchanged := MoreRightMinimalVersion(identity);
    left_more_unchanged := MoreLeftMinimalVersion(identity);
    right_identity := RightMinimalVersion(identity);
    left_identity := LeftMinimalVersion(identity);
    zero_map := ZeroMapping(one, one);
    right_zero := RightMinimalVersion(zero_map);

    quiver_a2 := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
    algebra_a2 := PathAlgebra(field, quiver_a2);
    projective_one := IndecProjectiveModules(algebra_a2)[1];
    simple_two := RightModuleOverPathAlgebra(algebra_a2, [0, 1], []);
    direct_sum := DirectSumOfQPAModules([projective_one, simple_two]);
    projection := DirectSumProjections(direct_sum)[1];
    inclusion := DirectSumInclusions(direct_sum)[1];
    a2_right := RightMinimalVersion(projection);
    a2_left := LeftMinimalVersion(inclusion);

    Print("{\"right_more\":"); PrintMore(right_more, right_more_map);
    Print(",\"left_more\":"); PrintMore(left_more, left_more_map);
    Print(",\"right_more_unchanged\":");
    PrintMore(right_more_unchanged, identity);
    Print(",\"left_more_unchanged\":");
    PrintMore(left_more_unchanged, identity);
    Print(",\"right_minimal\":"); PrintMinimal(right_result, right_map, "right");
    Print(",\"left_minimal\":"); PrintMinimal(left_result, left_map, "left");
    Print(",\"attribute_cache_identity\":[",
          IsIdenticalObj(right_result[1], right_again[1]), ",",
          ForAll([1..Length(right_result[2])],
                 i -> IsIdenticalObj(right_result[2][i], right_again[2][i])), ",",
          IsIdenticalObj(left_result[1], left_again[1]), ",",
          ForAll([1..Length(left_result[2])],
                 i -> IsIdenticalObj(left_result[2][i], left_again[2][i])), "]");
    Print(",\"right_identity\":"); PrintMinimal(right_identity, identity, "right");
    Print(",\"left_identity\":"); PrintMinimal(left_identity, identity, "left");
    Print(",\"right_zero\":"); PrintMinimal(right_zero, zero_map, "right");
    Print(",\"a2_right\":"); PrintMinimal(a2_right, projection, "right");
    Print(",\"a2_left\":"); PrintMinimal(a2_left, inclusion, "left");
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
