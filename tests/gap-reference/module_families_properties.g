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

PrintModuleList := function(modules)
    local first, module;
    Print("["); first := true;
    for module in modules do
        if not first then Print(","); fi;
        PrintModule(module); first := false;
    od;
    Print("]");
end;

PrintCase := function(field)
    local quiver, algebra, injectives, injectives_again, simples,
          simples_again, projectives, zero_module, semisimple_sum,
          injective_installed, simple_installed, injective_values,
          simple_values, semisimple_values, property_bound,
          loop_quiver, parent, x, quotient, quotient_injectives,
          quotient_injectives_again, quotient_simples,
          quotient_simples_again, quotient_projective, quotient_zero,
          quotient_injective_installed, quotient_simple_installed,
          quotient_injective_values, quotient_simple_values,
          quotient_semisimple_values, quotient_property_bound;

    quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
    algebra := PathAlgebra(field, quiver);
    injectives := IndecInjectiveModules(algebra);
    injectives_again := IndecInjectiveModules(algebra);
    simples := SimpleModules(algebra);
    simples_again := SimpleModules(algebra);
    projectives := IndecProjectiveModules(algebra);
    zero_module := ZeroModule(algebra);
    semisimple_sum := DirectSumOfQPAModules(simples);
    injective_installed := List(injectives, module -> [
        HasIsInjectiveModule(module), HasIsIndecomposableModule(module),
        HasIsProjectiveModule(module)]);
    simple_installed := List(simples, module -> [
        HasIsSimpleQPAModule(module), HasIsIndecomposableModule(module),
        HasIsSemisimpleModule(module)]);
    injective_values := List(
        Concatenation(injectives, simples, projectives, [zero_module]),
        IsInjectiveModule);
    simple_values := List(
        Concatenation(simples, projectives, [zero_module, semisimple_sum]),
        IsSimpleQPAModule);
    semisimple_values := List(
        Concatenation(simples, projectives, [zero_module, semisimple_sum]),
        IsSemisimpleModule);
    property_bound := List(
        Concatenation(injectives, simples, projectives, [zero_module,
            semisimple_sum]),
        module -> [HasIsInjectiveModule(module),
            HasIsSimpleQPAModule(module), HasIsSemisimpleModule(module)]);

    loop_quiver := Quiver(["w"], [["w", "w", "x"]]);
    parent := PathAlgebra(field, loop_quiver);
    x := parent.x;
    quotient := parent / [x * x];
    quotient_injectives := IndecInjectiveModules(quotient);
    quotient_injectives_again := IndecInjectiveModules(quotient);
    quotient_simples := SimpleModules(quotient);
    quotient_simples_again := SimpleModules(quotient);
    quotient_projective := IndecProjectiveModules(quotient)[1];
    quotient_zero := ZeroModule(quotient);
    quotient_injective_installed := List(quotient_injectives,
        module -> [HasIsInjectiveModule(module),
            HasIsIndecomposableModule(module), HasIsProjectiveModule(module)]);
    quotient_simple_installed := List(quotient_simples,
        module -> [HasIsSimpleQPAModule(module),
            HasIsIndecomposableModule(module), HasIsSemisimpleModule(module)]);
    quotient_injective_values := List(Concatenation(quotient_injectives,
        quotient_simples, [quotient_projective, quotient_zero]),
        IsInjectiveModule);
    quotient_simple_values := List(Concatenation(quotient_simples,
        [quotient_projective, quotient_zero]), IsSimpleQPAModule);
    quotient_semisimple_values := List(Concatenation(quotient_simples,
        [quotient_projective, quotient_zero]), IsSemisimpleModule);
    quotient_property_bound := List(Concatenation(quotient_injectives,
        quotient_simples, [quotient_projective, quotient_zero]),
        module -> [HasIsInjectiveModule(module),
            HasIsSimpleQPAModule(module), HasIsSemisimpleModule(module)]);

    Print("{\"path\":{");
    Print("\"injectives\":"); PrintModuleList(injectives);
    Print(",\"injective_cache_identity\":",
        List([1..Length(injectives)], i ->
            IsIdenticalObj(injectives[i], injectives_again[i])));
    Print(",\"injective_installed\":", injective_installed);
    Print(",\"simples\":"); PrintModuleList(simples);
    Print(",\"simple_cache_identity\":",
        List([1..Length(simples)], i ->
            IsIdenticalObj(simples[i], simples_again[i])));
    Print(",\"simple_installed\":", simple_installed);
    Print(",\"injective_values\":", injective_values);
    Print(",\"simple_values\":", simple_values);
    Print(",\"semisimple_values\":", semisimple_values);
    Print(",\"property_bound\":", property_bound);
    Print("},\"quotient\":{");
    Print("\"injectives\":"); PrintModuleList(quotient_injectives);
    Print(",\"injective_cache_identity\":",
        List([1..Length(quotient_injectives)], i -> IsIdenticalObj(
            quotient_injectives[i], quotient_injectives_again[i])));
    Print(",\"injective_installed\":", quotient_injective_installed);
    Print(",\"simples\":"); PrintModuleList(quotient_simples);
    Print(",\"simple_cache_identity\":",
        List([1..Length(quotient_simples)], i -> IsIdenticalObj(
            quotient_simples[i], quotient_simples_again[i])));
    Print(",\"simple_installed\":", quotient_simple_installed);
    Print(",\"injective_values\":", quotient_injective_values);
    Print(",\"simple_values\":", quotient_simple_values);
    Print(",\"semisimple_values\":", quotient_semisimple_values);
    Print(",\"property_bound\":", quotient_property_bound);
    Print("}}");
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
