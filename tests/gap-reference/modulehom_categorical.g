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

PrintMatrix := function(matrix)
    local first, row;
    Print("{\"shape\":");
    PrintIntegerList(DimensionsMat(matrix));
    Print(",\"rows\":[");
    first := true;
    for row in matrix do
        if not first then Print(","); fi;
        PrintScalarList(row);
        first := false;
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
    Print("[");
    first := true;
    for map in maps do
        if not first then Print(","); fi;
        PrintMap(map);
        first := false;
    od;
    Print("]");
end;

PrintCase := function(field)
    local quiver, algebra, module, f, t, identity, zero, layer, first,
          second, kernel_module, kernel_to_module, simple_one, simple_two,
          direct_sum, projections, inclusions, singleton, induced_kernel,
          induced_image, induced_cokernel, intersection, intersection_list,
          sum, sum_list, pullback, pushout, reject_zero, reject_one,
          reject_many;
    quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
    algebra := PathAlgebra(field, quiver);
    module := RightModuleOverPathAlgebra(
        algebra, [2, 2], [["a", [[1,1],[0,0]] * One(field)]]);
    f := RightModuleHomOverAlgebra(
        module, module,
        [[[1,0],[0,0]] * One(field), IdentityMat(2, field)]);
    t := RightModuleHomOverAlgebra(
        module, module,
        [[[2,0],[0,3]] * One(field), 2 * IdentityMat(2, field)]);
    identity := IdentityMapping(module);
    zero := ZeroModule(algebra);
    layer := RightModuleOverPathAlgebra(algebra, [1, 1], []);
    first := RightModuleHomOverAlgebra(
        layer, module, [[[0,1]] * One(field), [[1,0]] * One(field)]);
    second := RightModuleHomOverAlgebra(
        layer, module, [[[0,1]] * One(field), [[0,1]] * One(field)]);
    kernel_module := RightModuleOverPathAlgebra(algebra, [1, 0], []);
    kernel_to_module := RightModuleHomOverAlgebra(
        kernel_module, module,
        [[[0,1]] * One(field), NullMat(1, 2, field)]);
    simple_one := RightModuleOverPathAlgebra(algebra, [1, 0], []);
    simple_two := RightModuleOverPathAlgebra(algebra, [0, 1], []);

    direct_sum := DirectSumOfQPAModules([layer, kernel_module, simple_two]);
    projections := DirectSumProjections(direct_sum);
    inclusions := DirectSumInclusions(direct_sum);
    singleton := DirectSumOfQPAModules([layer]);
    induced_kernel := MorphismOnKernel(f, f, t, t);
    induced_image := MorphismOnImage(f, f, t, t);
    induced_cokernel := MorphismOnCoKernel(f, f, t, t);
    intersection := IntersectionOfSubmodules(first, second);
    intersection_list := IntersectionOfSubmodules(
        [first, second, kernel_to_module]);
    sum := SumOfSubmodules(first, second);
    sum_list := SumOfSubmodules([first, second, kernel_to_module]);
    pullback := PullBack(first, second);
    pushout := PushOut(intersection[1], intersection[2]);
    reject_zero := RejectOfModule(simple_one, simple_two);
    reject_one := RejectOfModule(module, simple_two);
    reject_many := RejectOfModule(module, simple_one);

    Print("{\"direct_sum\":{");
    Print("\"empty_is_fail\":", DirectSumOfQPAModules([]) = fail);
    Print(",\"filter\":", IsDirectSumOfModules(direct_sum));
    Print(",\"dimension\":");
    PrintIntegerList(DimensionVector(direct_sum));
    Print(",\"matrices\":");
    PrintMatrixList(MatricesOfPathAlgebraModule(direct_sum));
    Print(",\"projections\":");
    PrintMapList(projections);
    Print(",\"inclusions\":");
    PrintMapList(inclusions);
    Print(",\"projection_cache_identity\":",
          IsIdenticalObj(DirectSumProjections(direct_sum)[1], projections[1]));
    Print(",\"inclusion_cache_identity\":",
          IsIdenticalObj(DirectSumInclusions(direct_sum)[1], inclusions[1]));
    Print(",\"singleton_equal\":", singleton = layer);
    Print(",\"singleton_identical\":", IsIdenticalObj(singleton, layer));
    Print("},\"induced\":{");
    Print("\"kernel\":"); PrintMap(induced_kernel);
    Print(",\"image\":"); PrintMap(induced_image);
    Print(",\"cokernel\":"); PrintMap(induced_cokernel);
    Print(",\"noncommuting_fails\":[",
          MorphismOnKernel(f, f, identity, t) = fail, ",",
          MorphismOnImage(f, f, identity, t) = fail, ",",
          MorphismOnCoKernel(f, f, identity, t) = fail, "]");
    Print("},\"pullback\":"); PrintMapList(pullback);
    Print(",\"pushout\":"); PrintMapList(pushout);
    Print(",\"intersection\":{");
    Print("\"pair\":"); PrintMapList(intersection);
    Print(",\"list\":"); PrintMap(intersection_list);
    Print(",\"singleton_identity\":",
          IsIdenticalObj(IntersectionOfSubmodules([first]), first));
    Print("},\"sum\":{");
    Print("\"pair\":"); PrintMapList(sum);
    Print(",\"list\":"); PrintMap(sum_list);
    Print(",\"singleton_identity\":",
          IsIdenticalObj(SumOfSubmodules([first]), first));
    Print("},\"reject\":{");
    Print("\"hom_zero\":"); PrintMap(reject_zero);
    Print(",\"hom_one\":"); PrintMap(reject_one);
    Print(",\"hom_many\":"); PrintMap(reject_many);
    Print("}}");
end;

Print("{\"schema_version\":1,\"kind\":\"modulehom_categorical\",\"oracle\":{");
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
