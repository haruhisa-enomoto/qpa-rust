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

PrintDimensionVectors := function(modules)
    local first, module;
    Print("["); first := true;
    for module in modules do
        if not first then Print(","); fi;
        PrintIntegerList(DimensionVector(module)); first := false;
    od;
    Print("]");
end;

PrintCase := function(field)
    local quiver, algebra, zero_module, one, two, two_copy, common, maximal,
          quiver_a2, algebra_a2, projective_one, simple_one, simple_two,
          disjoint, mixed, mixed_common, identity, projection, inclusion,
          zero, top_projection, radical_inclusion, right_before, left_before;

    quiver := Quiver(["v"], []);
    algebra := PathAlgebra(field, quiver);
    zero_module := RightModuleOverPathAlgebra(algebra, [0], []);
    one := RightModuleOverPathAlgebra(algebra, [1], []);
    two := RightModuleOverPathAlgebra(algebra, [2], []);
    two_copy := DirectSumOfQPAModules([one, one]);

    common := CommonDirectSummand(two, one);
    maximal := MaximalCommonDirectSummand(two, one);

    quiver_a2 := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
    algebra_a2 := PathAlgebra(field, quiver_a2);
    projective_one := IndecProjectiveModules(algebra_a2)[1];
    simple_one := RightModuleOverPathAlgebra(algebra_a2, [1, 0], []);
    simple_two := RightModuleOverPathAlgebra(algebra_a2, [0, 1], []);
    disjoint := CommonDirectSummand(simple_one, simple_two);
    mixed := DirectSumOfQPAModules([simple_one, simple_two]);
    mixed_common := CommonDirectSummand(mixed, simple_one);

    identity := IdentityMapping(one);
    projection := RightModuleHomOverAlgebra(two, one,
        [[[1], [0]] * One(field)]);
    inclusion := RightModuleHomOverAlgebra(one, two,
        [[[1, 0]] * One(field)]);
    zero := ZeroMapping(one, one);
    top_projection := TopOfModuleProjection(projective_one);
    radical_inclusion := RadicalOfModuleInclusion(projective_one);
    right_before := HasIsRightMinimal(projection);
    left_before := HasIsLeftMinimal(inclusion);

    Print("{\"common\":"); PrintDimensionVectors(common);
    Print(",\"maximal\":{\"common\":"); PrintDimensionVectors(maximal[1]);
    Print(",\"source_remainder\":"); PrintIntegerList(DimensionVector(maximal[2]));
    Print(",\"range_remainder\":"); PrintIntegerList(DimensionVector(maximal[3]));
    Print("}");
    Print(",\"disjoint_is_false\":", disjoint = false);
    Print(",\"zero_common_is_false\":",
          CommonDirectSummand(zero_module, zero_module) = false);
    Print(",\"mixed_common\":"); PrintDimensionVectors(mixed_common);
    Print(",\"decisions\":{");
    Print("\"isomorphic\":", IsomorphicModules(two, two_copy));
    Print(",\"not_isomorphic\":", IsomorphicModules(simple_one, simple_two));
    Print(",\"direct_summand\":", IsDirectSummand(one, two));
    Print(",\"not_direct_summand\":", IsDirectSummand(simple_one, simple_two));
    Print(",\"additive_vector_space\":", IsInAdditiveClosure(two, one));
    Print(",\"not_additive\":", IsInAdditiveClosure(simple_one, simple_two));
    Print("}");
    Print(",\"minimality\":{");
    Print("\"right_before\":", right_before,
          ",\"right_projection\":", IsRightMinimal(projection),
          ",\"right_after\":", HasIsRightMinimal(projection));
    Print(",\"left_before\":", left_before,
          ",\"left_inclusion\":", IsLeftMinimal(inclusion),
          ",\"left_after\":", HasIsLeftMinimal(inclusion));
    Print(",\"identity\":[", IsRightMinimal(identity), ",",
          IsLeftMinimal(identity), "]");
    Print(",\"zero\":[", IsRightMinimal(zero), ",",
          IsLeftMinimal(zero), "]");
    Print(",\"projective_cover_right\":", IsRightMinimal(top_projection));
    Print(",\"radical_inclusion_left\":", IsLeftMinimal(radical_inclusion));
    Print("}}");
end;

Print("{\"schema_version\":1,\"kind\":\"modulehom_common_summand\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"rational\":"); PrintCase(Rationals);
Print(",\"binary\":"); PrintCase(GF(2));
Print("}\n");

QUIT_GAP(0);
