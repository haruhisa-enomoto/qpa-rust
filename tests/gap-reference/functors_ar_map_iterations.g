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

ComputeOne := function(map)
    return rec(
        map := map,
        d_zero := DTr(map, 0),
        d_one := DTr(map, 1),
        d_minus_one := DTr(map, -1),
        d_two := DTr(map, 2),
        t_zero := TrD(map, 0),
        t_one := TrD(map, 1),
        t_minus_one := TrD(map, -1),
        t_two := TrD(map, 2));
end;

ComputeCase := function(field)
    local a2_quiver, a2, a2_module, a2_map, kronecker_quiver,
          kronecker, kronecker_module, kronecker_map, loop_quiver,
          loop_parent, x, dual_numbers, dual_simple, dual_identity,
          shear, physical_zero;

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

    return rec(
        a2_shear := ComputeOne(a2_map),
        kronecker_shear := ComputeOne(kronecker_map),
        dual_simple_identity := ComputeOne(dual_identity));
end;

PrintOne := function(data)
    Print("{\"dtr_zero_identity\":",
          IsIdenticalObj(data.d_zero, data.map));
    Print(",\"dtr_one_identity\":",
          IsIdenticalObj(data.d_one, DTr(data.map)));
    Print(",\"dtr_minus_one_identity\":",
          IsIdenticalObj(data.d_minus_one, TrD(data.map)));
    Print(",\"dtr_two\":"); PrintMap(data.d_two);
    Print(",\"trd_zero_identity\":",
          IsIdenticalObj(data.t_zero, data.map));
    Print(",\"trd_one_identity\":",
          IsIdenticalObj(data.t_one, TrD(data.map)));
    Print(",\"trd_minus_one_identity\":",
          IsIdenticalObj(data.t_minus_one, DTr(data.map)));
    Print(",\"trd_two\":"); PrintMap(data.t_two);
    Print("}");
end;

PrintCase := function(data)
    Print("{\"a2_shear\":"); PrintOne(data.a2_shear);
    Print(",\"kronecker_shear\":"); PrintOne(data.kronecker_shear);
    Print(",\"dual_simple_identity\":");
    PrintOne(data.dual_simple_identity);
    Print("}");
end;

rational := ComputeCase(Rationals);;
binary := ComputeCase(GF(2));;

SizeScreen([100000, 100000]);
Print("{\"schema_version\":1,\"source\":{",
      "\"gap\":\"4.17dev@4efe877859ab0bc0ac041bf635a996fb897ee0c9\",",
      "\"qpa\":\"1.37@8d5f847258d4853700ce3a5d947db2a1441109f1\",",
      "\"gbnp\":\"1.1.0@e0b46c1aac16f3fbabe743e3169dc6bd86bdc654\"},",
      "\"rational\":");
PrintCase(rational);
Print(",\"binary\":"); PrintCase(binary);
Print("}\n");
