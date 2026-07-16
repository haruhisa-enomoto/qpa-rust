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

PrintTransposeMap := function(map)
    local has_before, value, has_after, repeated, source_attribute,
          range_attribute;
    has_before := HasTransposeOfModuleHomomorphism(map);
    value := TransposeOfModuleHomomorphism(map);
    has_after := HasTransposeOfModuleHomomorphism(map);
    repeated := TransposeOfModuleHomomorphism(map);
    source_attribute := TransposeOfModule(Range(map));
    range_attribute := TransposeOfModule(Source(map));
    Print("{\"has_before\":", has_before, ",\"value\":");
    PrintMap(value);
    Print(",\"has_after\":", has_after);
    Print(",\"repeat_identity\":", IsIdenticalObj(value, repeated));
    Print(",\"source_attribute_identity\":",
          IsIdenticalObj(Source(value), source_attribute));
    Print(",\"range_attribute_identity\":",
          IsIdenticalObj(Range(value), range_attribute), "}");
end;

PrintDTrMap := function(map)
    local has_before, value, has_after;
    has_before := HasDTr(map);
    value := DTr(map);
    has_after := HasDTr(map);
    Print("{\"has_before\":", has_before, ",\"value\":");
    PrintMap(value);
    Print(",\"has_after\":", has_after);
    Print(",\"repeat_identity\":", IsIdenticalObj(value, DTr(map)));
    Print(",\"source_attribute_identity\":",
          IsIdenticalObj(Source(value), DTr(Source(map))));
    Print(",\"range_attribute_identity\":",
          IsIdenticalObj(Range(value), DTr(Range(map))), "}");
end;

PrintTrDMap := function(map)
    local has_before, value, has_after;
    has_before := HasTrD(map);
    value := TrD(map);
    has_after := HasTrD(map);
    Print("{\"has_before\":", has_before, ",\"value\":");
    PrintMap(value);
    Print(",\"has_after\":", has_after);
    Print(",\"repeat_identity\":", IsIdenticalObj(value, TrD(map)));
    Print(",\"source_attribute_identity\":",
          IsIdenticalObj(Source(value), TrD(Source(map))));
    Print(",\"range_attribute_identity\":",
          IsIdenticalObj(Range(value), TrD(Range(map))), "}");
end;

PrintFunctorCase := function(map)
    Print("{\"entered\":"); PrintMap(map);
    Print(",\"transpose\":"); PrintTransposeMap(map);
    Print(",\"dtr\":"); PrintDTrMap(map);
    Print(",\"trd\":"); PrintTrDMap(map);
    Print("}");
end;

PrintCase := function(field)
    local a2_quiver, a2, a2_module, a2_map, kronecker_quiver,
          kronecker, kronecker_module, kronecker_map, loop_quiver,
          loop_parent, x, dual_numbers, dual_simple, dual_identity,
          bound_quiver, bound_parent, a, b, bound, bound_simple,
          bound_identity, shear, physical_zero;

    shear := [[1, 1], [0, 1]] * One(field);
    physical_zero := [[0]] * One(field);

    a2_quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
    a2 := PathAlgebra(field, a2_quiver);
    a2_module := RightModuleOverPathAlgebra(a2, [2, 0], []);
    a2_map := RightModuleHomOverAlgebra(a2_module, a2_module,
        [shear, physical_zero]);

    kronecker_quiver := Quiver(["k1", "k2"],
        [["k1", "k2", "p"], ["k1", "k2", "q"]]);
    kronecker := PathAlgebra(field, kronecker_quiver);
    kronecker_module := RightModuleOverPathAlgebra(kronecker, [2, 0], []);
    kronecker_map := RightModuleHomOverAlgebra(
        kronecker_module, kronecker_module, [shear, physical_zero]);

    loop_quiver := Quiver(["w"], [["w", "w", "x"]]);
    loop_parent := PathAlgebra(field, loop_quiver);
    x := loop_parent.x;
    dual_numbers := loop_parent / [x * x];
    dual_simple := SimpleModules(dual_numbers)[1];
    dual_identity := IdentityMapping(dual_simple);

    bound_quiver := Quiver(["r1", "r2", "r3"],
        [["r1", "r2", "a"], ["r2", "r3", "b"]]);
    bound_parent := PathAlgebra(field, bound_quiver);
    a := bound_parent.a; b := bound_parent.b;
    bound := bound_parent / [a * b];
    bound_simple := SimpleModules(bound)[1];
    bound_identity := IdentityMapping(bound_simple);

    Print("{\"a2_shear\":"); PrintFunctorCase(a2_map);
    Print(",\"kronecker_shear\":"); PrintFunctorCase(kronecker_map);
    Print(",\"dual_simple_identity\":"); PrintFunctorCase(dual_identity);
    Print(",\"bound_a3_identity\":"); PrintFunctorCase(bound_identity);
    Print("}");
end;

Print("{\"schema_version\":1,\"source\":{",
      "\"gap\":\"4.17dev@4efe877859ab0bc0ac041bf635a996fb897ee0c9\",",
      "\"qpa\":\"1.37@8d5f847258d4853700ce3a5d947db2a1441109f1\",",
      "\"gbnp\":\"1.1.0@e0b46c1aac16f3fbabe743e3169dc6bd86bdc654\"},",
      "\"rational\":");
PrintCase(Rationals);
Print(",\"binary\":"); PrintCase(GF(2));
Print("}\n");
