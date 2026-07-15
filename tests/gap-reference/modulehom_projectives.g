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

PrintMatrixList := function(matrices)
    local first, matrix;
    Print("["); first := true;
    for matrix in matrices do
        if not first then Print(","); fi;
        PrintMatrix(matrix); first := false;
    od;
    Print("]");
end;

PrintModule := function(module)
    Print("{\"dimension\":"); PrintIntegerList(DimensionVector(module));
    Print(",\"matrices\":");
    if Dimension(module) = 0 then
        Print("[]");
    else
        PrintMatrixList(MatricesOfPathAlgebraModule(module));
    fi;
    Print("}");
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

PrintMapList := function(maps)
    local first, map;
    Print("["); first := true;
    for map in maps do
        if not first then Print(","); fi;
        PrintMap(map); first := false;
    od;
    Print("]");
end;

PrintProjectiveBasisPaths := function(bases)
    local first_projective, first_vertex, first_element, projective,
          vertex_basis, element;
    Print("["); first_projective := true;
    for projective in bases do
        if not first_projective then Print(","); fi;
        Print("["); first_vertex := true;
        for vertex_basis in projective do
            if not first_vertex then Print(","); fi;
            Print("["); first_element := true;
            for element in vertex_basis do
                if not first_element then Print(","); fi;
                PrintJsonString(String(TipMonomial(element)));
                first_element := false;
            od;
            Print("]"); first_vertex := false;
        od;
        Print("]"); first_projective := false;
    od;
    Print("]");
end;

PrintCase := function(field)
    local quiver, algebra, projectives, bases, module, generators, cover,
          from_generator, projective_sum, lifting, nondirect_lifting,
          zero_module, zero_to_module, zero_lifting, factoring, simple_one,
          simple_two, layer, first, second, kernel_module, kernel_inclusion,
          inclusion_lift, inclusion_fail, loop_quiver, path_algebra, loop,
          dual_numbers, dual_projective, dual_basis, dual_simple, dual_cover;
    quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
    algebra := PathAlgebra(field, quiver);
    projectives := IndecProjectiveModules(algebra);
    bases := BasisOfProjectives(algebra);
    module := RightModuleOverPathAlgebra(
        algebra, [2, 2], [["a", [[1,0],[0,0]] * One(field)]]);
    generators := MinimalGeneratingSetOfModule(module);
    cover := ProjectiveCover(module);
    from_generator := HomFromProjective(generators[1], module);
    projective_sum := DirectSumOfQPAModules(projectives);
    lifting := LiftingMorphismFromProjective(cover, cover);
    nondirect_lifting := LiftingMorphismFromProjective(cover, from_generator);
    zero_module := ZeroModule(algebra);
    zero_to_module := ZeroMapping(zero_module, module);
    zero_lifting := LiftingMorphismFromProjective(cover, zero_to_module);
    factoring := HomFactoringThroughProjOverAlgebra(module, module);
    simple_one := RightModuleOverPathAlgebra(algebra, [1, 0], []);
    simple_two := RightModuleOverPathAlgebra(algebra, [0, 1], []);
    layer := RightModuleOverPathAlgebra(algebra, [1, 1], []);
    first := RightModuleHomOverAlgebra(
        layer, module, [[[0,1]] * One(field), [[1,0]] * One(field)]);
    second := RightModuleHomOverAlgebra(
        layer, module, [[[0,1]] * One(field), [[0,1]] * One(field)]);
    kernel_module := RightModuleOverPathAlgebra(algebra, [1, 0], []);
    kernel_inclusion := RightModuleHomOverAlgebra(
        kernel_module, module,
        [[[0,1]] * One(field), NullMat(1, 2, field)]);
    inclusion_lift := LiftingInclusionMorphisms(first, kernel_inclusion);
    inclusion_fail := LiftingInclusionMorphisms(first, second);

    loop_quiver := Quiver(["w"], [["w", "w", "x"]]);
    path_algebra := PathAlgebra(field, loop_quiver);
    loop := ArrowsOfQuiver(loop_quiver)[1] * One(path_algebra);
    dual_numbers := path_algebra / [loop * loop];
    dual_projective := IndecProjectiveModules(dual_numbers)[1];
    dual_basis := BasisOfProjectives(dual_numbers);
    dual_simple := RightModuleOverPathAlgebra(dual_numbers, [1], []);
    dual_cover := ProjectiveCover(dual_simple);

    Print("{\"a2\":{");
    Print("\"basis_paths\":"); PrintProjectiveBasisPaths(bases);
    Print(",\"projectives\":["); PrintModule(projectives[1]); Print(",");
    PrintModule(projectives[2]); Print("]");
    Print(",\"cache_identity\":[",
          IsIdenticalObj(IndecProjectiveModules(algebra)[1], projectives[1]), ",",
          IsIdenticalObj(IndecProjectiveModules(algebra)[2], projectives[2]), "]");
    Print(",\"installed_properties\":[",
          IsProjectiveModule(projectives[1]), ",",
          IsProjectiveModule(projectives[2]), ",",
          IsIndecomposableModule(projectives[1]), ",",
          IsIndecomposableModule(projectives[2]), "]");
    Print(",\"projectivity\":[",
          IsProjectiveModule(projective_sum), ",",
          IsProjectiveModule(simple_one), ",",
          IsProjectiveModule(simple_two), ",",
          IsProjectiveModule(zero_module), "]");
    Print(",\"minimal_generator_count\":", Length(generators));
    Print(",\"from_generator\":"); PrintMap(from_generator);
    Print(",\"cover\":"); PrintMap(cover);
    Print(",\"cover_cache_identity\":",
          IsIdenticalObj(ProjectiveCover(module), cover));
    Print(",\"lifting\":"); PrintMap(lifting);
    Print(",\"lifting_triangle\":", lifting * cover = cover);
    Print(",\"nondirect_lifting_is_fail\":", nondirect_lifting = fail);
    Print(",\"zero_lifting\":"); PrintMap(zero_lifting);
    Print(",\"factoring_count\":", Length(factoring));
    Print(",\"factoring\":"); PrintMapList(factoring);
    Print(",\"factoring_zero_counts\":[",
          Length(HomFactoringThroughProjOverAlgebra(zero_module, module)), ",",
          Length(HomFactoringThroughProjOverAlgebra(module, zero_module)), "]");
    Print(",\"inclusion_lift\":"); PrintMap(inclusion_lift);
    Print(",\"inclusion_triangle\":", inclusion_lift * first = kernel_inclusion);
    Print(",\"inclusion_fail\":", inclusion_fail = fail);
    Print("},\"dual_numbers\":{");
    Print("\"basis_paths\":"); PrintProjectiveBasisPaths(dual_basis);
    Print(",\"projective\":"); PrintModule(dual_projective);
    Print(",\"generator_count\":",
          Length(MinimalGeneratingSetOfModule(dual_projective)));
    Print(",\"is_projective\":", IsProjectiveModule(dual_projective));
    Print(",\"simple_is_projective\":", IsProjectiveModule(dual_simple));
    Print(",\"cover\":"); PrintMap(dual_cover);
    Print("}}");
end;

Print("{\"schema_version\":1,\"kind\":\"modulehom_projectives\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("GBNP")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"rational\":"); PrintCase(Rationals);
Print(",\"binary\":"); PrintCase(GF(2));
Print("}\n");

QUIT_GAP(0);
