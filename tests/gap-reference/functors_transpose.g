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

PrintTranspose := function(module)
    local has_before, value, has_after;
    has_before := HasTransposeOfModule(module);
    value := TransposeOfModule(module);
    has_after := HasTransposeOfModule(module);
    Print("{\"has_before\":", has_before, ",\"value\":");
    PrintModule(value);
    Print(",\"has_after\":", has_after);
    Print(",\"repeat_identity\":", IsIdenticalObj(value,
          TransposeOfModule(module)));
    Print(",\"indecomposable_bound\":",
          HasIsIndecomposableModule(value));
    if HasIsIndecomposableModule(value) then
        Print(",\"indecomposable_value\":", IsIndecomposableModule(value));
    else
        Print(",\"indecomposable_value\":null");
    fi;
    Print("}");
end;

PrintDTr := function(module)
    local has_before, value, has_after;
    has_before := HasDTr(module);
    value := DTr(module);
    has_after := HasDTr(module);
    Print("{\"has_before\":", has_before, ",\"value\":");
    PrintModule(value);
    Print(",\"has_after\":", has_after);
    Print(",\"repeat_identity\":", IsIdenticalObj(value, DTr(module)), "}");
end;

PrintTrD := function(module)
    local has_before, value, has_after;
    has_before := HasTrD(module);
    value := TrD(module);
    has_after := HasTrD(module);
    Print("{\"has_before\":", has_before, ",\"value\":");
    PrintModule(value);
    Print(",\"has_after\":", has_after);
    Print(",\"repeat_identity\":", IsIdenticalObj(value, TrD(module)), "}");
end;

PrintTauRigid := function(module)
    local has_before, value, has_after;
    has_before := HasIsTauRigidModule(module);
    value := IsTauRigidModule(module);
    has_after := HasIsTauRigidModule(module);
    Print("{\"has_before\":", has_before);
    Print(",\"value\":", value);
    Print(",\"has_after\":", has_after, "}");
end;

PrintFunctorCase := function(module)
    Print("{\"transpose\":"); PrintTranspose(module);
    Print(",\"dtr\":"); PrintDTr(module);
    Print(",\"trd\":"); PrintTrD(module);
    Print(",\"tau_rigid\":"); PrintTauRigid(module);
    Print("}");
end;

PrintCase := function(field)
    local quiver, algebra, simples, projectives, mixed, kronecker_quiver,
          kronecker, kronecker_simple, loop_quiver, loop_parent, x,
          dual_numbers, dual_simple, bound_quiver, bound_parent, a, b,
          bound, bound_simple;

    quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
    algebra := PathAlgebra(field, quiver);
    simples := SimpleModules(algebra);
    projectives := IndecProjectiveModules(algebra);
    mixed := RightModuleOverPathAlgebra(algebra, [2, 1],
        [["a", [[1], [0]] * One(field)]]);

    kronecker_quiver := Quiver(["k1", "k2"],
        [["k1", "k2", "p"], ["k1", "k2", "q"]]);
    kronecker := PathAlgebra(field, kronecker_quiver);
    kronecker_simple := SimpleModules(kronecker)[1];

    loop_quiver := Quiver(["w"], [["w", "w", "x"]]);
    loop_parent := PathAlgebra(field, loop_quiver);
    x := loop_parent.x;
    dual_numbers := loop_parent / [x * x];
    dual_simple := SimpleModules(dual_numbers)[1];

    bound_quiver := Quiver(["r1", "r2", "r3"],
        [["r1", "r2", "a"], ["r2", "r3", "b"]]);
    bound_parent := PathAlgebra(field, bound_quiver);
    a := bound_parent.a; b := bound_parent.b;
    bound := bound_parent / [a * b];
    bound_simple := SimpleModules(bound)[1];

    Print("{\"a2_s1\":"); PrintFunctorCase(simples[1]);
    Print(",\"a2_s2_projective\":"); PrintFunctorCase(simples[2]);
    Print(",\"a2_mixed\":"); PrintFunctorCase(mixed);
    Print(",\"kronecker_s1\":"); PrintFunctorCase(kronecker_simple);
    Print(",\"dual_simple\":"); PrintFunctorCase(dual_simple);
    Print(",\"bound_a3_s1\":"); PrintFunctorCase(bound_simple);
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
