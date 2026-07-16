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

PrintNestedScalars := function(rows)
    local first, row;
    Print("["); first := true;
    for row in rows do
        if not first then Print(","); fi;
        PrintScalarList(row); first := false;
    od;
    Print("]");
end;

PrintElements := function(elements)
    local first, element, external;
    Print("["); first := true;
    for element in elements do
        if not first then Print(","); fi;
        external := ExtRepOfObj(element);
        if IsPathModuleElem(external) then
            external := ExtRepOfObj(external);
        fi;
        PrintNestedScalars(external); first := false;
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

PrintModule := function(module)
    local first, entry;
    Print("{\"dimension\":"); PrintIntegerList(DimensionVector(module));
    Print(",\"matrices\":["); first := true;
    for entry in MatricesOfPathAlgebraModule(module) do
        if not first then Print(","); fi;
        if IsMatrix(entry) then
            PrintMatrix(entry);
        else
            PrintJsonString(String(entry));
        fi;
        first := false;
    od;
    Print("]}");
end;

ComputeOne := function(module)
    local d_one;
    d_one := DTr(module, 1);
    return rec(
        module := module,
        d_zero := DTr(module, 0),
        d_one := d_one,
        d_one_generators := MinimalGeneratingSetOfModule(d_one),
        d_minus_one := DTr(module, -1),
        d_two := DTr(module, 2),
        t_zero := TrD(module, 0),
        t_one := TrD(module, 1),
        t_minus_one := TrD(module, -1),
        t_two := TrD(module, 2));
end;

ComputeCase := function(field)
    local quiver, algebra, a2_simple, kronecker_quiver, kronecker,
          kronecker_simple, loop_quiver, parent, x, dual_numbers,
          dual_simple;
    quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
    algebra := PathAlgebra(field, quiver);
    a2_simple := SimpleModules(algebra)[1];
    kronecker_quiver := Quiver(["k1", "k2"],
        [["k1", "k2", "p"], ["k1", "k2", "q"]]);
    kronecker := PathAlgebra(field, kronecker_quiver);
    kronecker_simple := SimpleModules(kronecker)[1];
    loop_quiver := Quiver(["w"], [["w", "w", "x"]]);
    parent := PathAlgebra(field, loop_quiver);
    x := parent.x;
    dual_numbers := parent / [x * x];
    dual_simple := SimpleModules(dual_numbers)[1];
    return rec(
        a2_s1 := ComputeOne(a2_simple),
        kronecker_s1 := ComputeOne(kronecker_simple),
        dual_simple := ComputeOne(dual_simple));
end;

PrintOne := function(data)
    Print("{\"dtr_zero_identity\":",
          IsIdenticalObj(data.d_zero, data.module));
    Print(",\"dtr_one_identity\":",
          IsIdenticalObj(data.d_one, DTr(data.module)));
    Print(",\"dtr_one_generators\":");
    PrintElements(data.d_one_generators);
    Print(",\"dtr_minus_one_identity\":",
          IsIdenticalObj(data.d_minus_one, TrD(data.module)));
    Print(",\"dtr_two\":"); PrintModule(data.d_two);
    Print(",\"trd_zero_identity\":",
          IsIdenticalObj(data.t_zero, data.module));
    Print(",\"trd_one_identity\":",
          IsIdenticalObj(data.t_one, TrD(data.module)));
    Print(",\"trd_minus_one_identity\":",
          IsIdenticalObj(data.t_minus_one, DTr(data.module)));
    Print(",\"trd_two\":"); PrintModule(data.t_two);
    Print("}");
end;

PrintCase := function(data)
    Print("{\"a2_s1\":"); PrintOne(data.a2_s1);
    Print(",\"kronecker_s1\":"); PrintOne(data.kronecker_s1);
    Print(",\"dual_simple\":"); PrintOne(data.dual_simple);
    Print("}");
end;

rational := ComputeCase(Rationals);;
binary := ComputeCase(GF(2));;

Print("{\"schema_version\":1,\"source\":{",
      "\"gap\":\"4.17dev@4efe877859ab0bc0ac041bf635a996fb897ee0c9\",",
      "\"qpa\":\"1.37@8d5f847258d4853700ce3a5d947db2a1441109f1\",",
      "\"gbnp\":\"1.1.0@e0b46c1aac16f3fbabe743e3169dc6bd86bdc654\"},",
      "\"rational\":");
PrintCase(rational);
Print(",\"binary\":"); PrintCase(binary);
Print("}\n");
