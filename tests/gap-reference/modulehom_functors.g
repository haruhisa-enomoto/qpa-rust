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
    Print("[");
    first := true;
    for value in values do
        if not first then Print(","); fi;
        Print(value);
        first := false;
    od;
    Print("]");
end;

PrintScalarList := function(values)
    local first, value;
    Print("[");
    first := true;
    for value in values do
        if not first then Print(","); fi;
        PrintJsonString(String(value));
        first := false;
    od;
    Print("]");
end;

PrintNestedScalars := function(rows)
    local first, row;
    Print("[");
    first := true;
    for row in rows do
        if not first then Print(","); fi;
        PrintScalarList(row);
        first := false;
    od;
    Print("]");
end;

PrintMatrix := function(matrix)
    local first_row, row;
    Print("{\"shape\":");
    PrintIntegerList(DimensionsMat(matrix));
    Print(",\"rows\":[");
    first_row := true;
    for row in matrix do
        if not first_row then Print(","); fi;
        PrintScalarList(row);
        first_row := false;
    od;
    Print("]}");
end;

PrintMatrixList := function(matrices)
    local first, matrix;
    Print("[");
    first := true;
    for matrix in matrices do
        if not first then Print(","); fi;
        PrintMatrix(matrix);
        first := false;
    od;
    Print("]");
end;

PrintElements := function(elements)
    local first, element, external;
    Print("[");
    first := true;
    for element in elements do
        if not first then Print(","); fi;
        external := ExtRepOfObj(element);
        if IsPathModuleElem(external) then
            external := ExtRepOfObj(external);
        fi;
        PrintNestedScalars(external);
        first := false;
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

PrintCase := function(field)
    local quiver, algebra, module, map, minimal_generators, socle_inclusion,
          top_projection, top_module, dual_map, target, simple_one, simple_two,
          zero_module, zero_socle,
          zero_generators, zero_trace, zero_source_trace, zero_target_trace,
          nonzero_trace, legacy_method, legacy_trace, ordered_trace,
          ordered_legacy;
    quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
    algebra := PathAlgebra(field, quiver);
    module := RightModuleOverPathAlgebra(
        algebra, [2, 1], [["a", [[1],[0]] * One(field)]]);
    map := RightModuleHomOverAlgebra(
        module, module,
        [[[1,0],[0,0]] * One(field), [[1]] * One(field)]);
    minimal_generators := MinimalGeneratingSetOfModule(module);
    top_projection := TopOfModuleProjection(module);
    top_module := TopOfModule(module);
    socle_inclusion := SocleOfModuleInclusion(module);
    dual_map := DualOfModuleHomomorphism(map);
    target := RightModuleOverPathAlgebra(
        algebra, [1, 2], [["a", [[1,1]] * One(field)]]);
    simple_one := RightModuleOverPathAlgebra(algebra, [1, 0], []);
    simple_two := RightModuleOverPathAlgebra(algebra, [0, 1], []);
    zero_module := ZeroModule(algebra);
    zero_socle := SocleOfModuleInclusion(zero_module);
    zero_generators := MinimalGeneratingSetOfModule(zero_module);
    zero_trace := TraceOfModule(simple_one, target);
    zero_source_trace := TraceOfModule(zero_module, target);
    zero_target_trace := TraceOfModule(simple_two, zero_module);
    nonzero_trace := TraceOfModule(simple_two, target);
    legacy_method := First(MethodsOperation(TraceOfModule, 2),
                           method -> method.location[2] = 1735);
    legacy_trace := legacy_method.func(simple_two, target);
    ordered_trace := TraceOfModule(module, module);
    ordered_legacy := legacy_method.func(module, module);

    Print("{\"minimal_generators\":{");
    Print("\"count\":", Length(minimal_generators), ",\"elements\":");
    PrintElements(minimal_generators);
    Print(",\"zero_count\":", Length(zero_generators));
    Print(",\"top_projection_cached_identity\":",
          IsIdenticalObj(TopOfModuleProjection(module), top_projection));
    Print(",\"top_cached_identity\":",
          IsIdenticalObj(TopOfModule(module), top_module));
    Print(",\"top_projection_range_identity\":",
          IsIdenticalObj(Range(top_projection), top_module));
    Print("},\"socle\":{");
    Print("\"inclusion\":");
    PrintMap(socle_inclusion);
    Print(",\"object_dimension\":");
    PrintIntegerList(DimensionVector(SocleOfModule(module)));
    Print(",\"zero_inclusion\":");
    PrintMap(zero_socle);
    Print(",\"zero_source_identity\":",
          IsIdenticalObj(Source(zero_socle), zero_module));
    Print("},\"dual_map\":{");
    Print("\"cached_identity\":",
          IsIdenticalObj(DualOfModuleHomomorphism(map), dual_map), ",");
    Print("\"source_dimension\":");
    PrintIntegerList(DimensionVector(Source(dual_map)));
    Print(",\"range_dimension\":");
    PrintIntegerList(DimensionVector(Range(dual_map)));
    Print(",\"matrices\":");
    PrintMatrixList(MatricesOfPathAlgebraMatModuleHomomorphism(dual_map));
    Print("},\"trace\":{");
    Print("\"zero\":");
    PrintMap(zero_trace);
    Print(",\"zero_source\":");
    PrintMap(zero_source_trace);
    Print(",\"zero_target\":");
    PrintMap(zero_target_trace);
    Print(",\"nonzero\":");
    PrintMap(nonzero_trace);
    Print(",\"legacy\":");
    PrintMap(legacy_trace);
    Print(",\"legacy_equal\":", legacy_trace = nonzero_trace);
    Print(",\"ordered\":");
    PrintMap(ordered_trace);
    Print(",\"ordered_legacy\":");
    PrintMap(ordered_legacy);
    Print(",\"ordered_legacy_equal\":", ordered_legacy = ordered_trace);
    Print(",\"right_fac_equal\":",
          RightFacApproximation(simple_two, target) = nonzero_trace);
    Print("}}");
end;

Print("{\"schema_version\":1,\"kind\":\"modulehom_functors\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":");
PrintJsonString(PackageInfo("GBNP")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"rational\":");
PrintCase(Rationals);
Print(",\"binary\":");
PrintCase(GF(2));
Print("}\n");

QUIT_GAP(0);
