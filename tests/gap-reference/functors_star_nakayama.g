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

PrintStarModule := function(module)
    local has_before, value, has_after;
    has_before := HasStarOfModule(module);
    value := StarOfModule(module);
    has_after := HasStarOfModule(module);
    Print("{\"has_before\":", has_before, ",\"value\":");
    PrintModule(value);
    Print(",\"has_after\":", has_after);
    Print(",\"repeat_identity\":",
          IsIdenticalObj(value, StarOfModule(module)), "}");
end;

PrintNakayamaModule := function(module)
    local has_before, value, has_after;
    has_before := HasNakayamaFunctorOfModule(module);
    value := NakayamaFunctorOfModule(module);
    has_after := HasNakayamaFunctorOfModule(module);
    Print("{\"has_before\":", has_before, ",\"value\":");
    PrintModule(value);
    Print(",\"has_after\":", has_after);
    Print(",\"repeat_identity\":",
          IsIdenticalObj(value, NakayamaFunctorOfModule(module)), "}");
end;

PrintOppositeNakayamaModule := function(module)
    local has_before, value, has_after;
    has_before := HasOppositeNakayamaFunctorOfModule(module);
    value := OppositeNakayamaFunctorOfModule(module);
    has_after := HasOppositeNakayamaFunctorOfModule(module);
    Print("{\"has_before\":", has_before, ",\"value\":");
    PrintModule(value);
    Print(",\"has_after\":", has_after);
    Print(",\"repeat_identity\":", IsIdenticalObj(value,
          OppositeNakayamaFunctorOfModule(module)), "}");
end;

PrintModuleFunctorCase := function(module)
    Print("{\"star\":"); PrintStarModule(module);
    Print(",\"nakayama\":"); PrintNakayamaModule(module);
    Print(",\"opposite_nakayama\":");
    PrintOppositeNakayamaModule(module);
    Print("}");
end;

PrintStarMap := function(map)
    local has_before, value, has_after;
    has_before := HasStarOfModuleHomomorphism(map);
    value := StarOfModuleHomomorphism(map);
    has_after := HasStarOfModuleHomomorphism(map);
    Print("{\"has_before\":", has_before, ",\"value\":");
    PrintMap(value);
    Print(",\"has_after\":", has_after);
    Print(",\"repeat_identity\":", IsIdenticalObj(value,
          StarOfModuleHomomorphism(map)), "}");
end;

PrintNakayamaMap := function(map)
    local has_before, value, has_after;
    has_before := HasNakayamaFunctorOfModuleHomomorphism(map);
    value := NakayamaFunctorOfModuleHomomorphism(map);
    has_after := HasNakayamaFunctorOfModuleHomomorphism(map);
    Print("{\"has_before\":", has_before, ",\"value\":");
    PrintMap(value);
    Print(",\"has_after\":", has_after);
    Print(",\"repeat_identity\":", IsIdenticalObj(value,
          NakayamaFunctorOfModuleHomomorphism(map)), "}");
end;

PrintOppositeNakayamaMap := function(map)
    local has_before, value, has_after;
    has_before := HasOppositeNakayamaFunctorOfModuleHomomorphism(map);
    value := OppositeNakayamaFunctorOfModuleHomomorphism(map);
    has_after := HasOppositeNakayamaFunctorOfModuleHomomorphism(map);
    Print("{\"has_before\":", has_before, ",\"value\":");
    PrintMap(value);
    Print(",\"has_after\":", has_after);
    Print(",\"repeat_identity\":", IsIdenticalObj(value,
          OppositeNakayamaFunctorOfModuleHomomorphism(map)), "}");
end;

PrintMapFunctorCase := function(map)
    Print("{\"star\":"); PrintStarMap(map);
    Print(",\"nakayama\":"); PrintNakayamaMap(map);
    Print(",\"opposite_nakayama\":"); PrintOppositeNakayamaMap(map);
    Print("}");
end;

PrintCase := function(field)
    local quiver, algebra, simples, projectives, zero, inclusion,
          module, module_map, loop_quiver, parent, x, quotient,
          quotient_simple, quotient_projective, quotient_radical_map;

    quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
    algebra := PathAlgebra(field, quiver);
    simples := SimpleModules(algebra);
    projectives := IndecProjectiveModules(algebra);
    zero := ZeroModule(algebra);
    inclusion := KernelInclusion(ProjectiveCover(simples[1]));
    module := RightModuleOverPathAlgebra(algebra, [2, 1],
        [["a", [[1], [0]] * One(field)]]);
    module_map := RightModuleHomOverAlgebra(module, module,
        [[[1, 0], [0, 0]] * One(field), [[1]] * One(field)]);

    loop_quiver := Quiver(["w"], [["w", "w", "x"]]);
    parent := PathAlgebra(field, loop_quiver);
    x := parent.x;
    quotient := parent / [x * x];
    quotient_simple := SimpleModules(quotient)[1];
    quotient_projective := IndecProjectiveModules(quotient)[1];
    quotient_radical_map := HomFromProjective(
        MinimalGeneratingSetOfModule(quotient_projective)[1] ^ quotient.x,
        quotient_projective);

    Print("{\"modules\":{\"a2_s2\":");
    PrintModuleFunctorCase(simples[2]);
    Print(",\"a2_p1\":"); PrintModuleFunctorCase(projectives[1]);
    Print(",\"a2_zero\":"); PrintModuleFunctorCase(zero);
    Print(",\"a2_mixed\":"); PrintModuleFunctorCase(module);
    Print(",\"dual_simple\":"); PrintModuleFunctorCase(quotient_simple);
    Print(",\"dual_projective\":");
    PrintModuleFunctorCase(quotient_projective);
    Print("},\"maps\":{\"a2_inclusion\":");
    PrintMapFunctorCase(inclusion);
    Print(",\"a2_endomorphism\":"); PrintMapFunctorCase(module_map);
    Print(",\"dual_radical\":");
    PrintMapFunctorCase(quotient_radical_map);
    Print("}}");
end;

Print("{\"schema_version\":1,\"source\":{",
      "\"gap\":\"4.17dev@4efe877859ab0bc0ac041bf635a996fb897ee0c9\",",
      "\"qpa\":\"1.37@8d5f847258d4853700ce3a5d947db2a1441109f1\",",
      "\"gbnp\":\"1.1.0@e0b46c1aac16f3fbabe743e3169dc6bd86bdc654\"},",
      "\"rational\":");
PrintCase(Rationals);
Print(",\"binary\":"); PrintCase(GF(2));
Print("}\n");
