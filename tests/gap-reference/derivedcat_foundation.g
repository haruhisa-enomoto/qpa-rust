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
    Print("[", JoinStringsWithSeparator(List(values, String), ","), "]");
end;

PrintMaybeInteger := function(value)
    if value = false then
        PrintJsonString("false");
    else
        Print(value);
    fi;
end;

PrintMaybeIntegerList := function(values)
    local first, value;
    Print("["); first := true;
    for value in values do
        if not first then Print(","); fi;
        PrintMaybeInteger(value); first := false;
    od;
    Print("]");
end;

PrintIntegerLists := function(values)
    local first, value;
    Print("["); first := true;
    for value in values do
        if not first then Print(","); fi;
        PrintIntegerList(value); first := false;
    od;
    Print("]");
end;

PrintScalar := function(value)
    if IsFFE(value) then
        Print(IntFFE(value));
    elif IsRat(value) then
        Print("[", NumeratorRat(value), ",", DenominatorRat(value), "]");
    else
        Error("unsupported fixture scalar");
    fi;
end;

PrintMatrix := function(matrix)
    local firstRow, firstEntry, row, entry;
    Print("["); firstRow := true;
    for row in matrix do
        if not firstRow then Print(","); fi;
        Print("["); firstEntry := true;
        for entry in row do
            if not firstEntry then Print(","); fi;
            PrintScalar(entry); firstEntry := false;
        od;
        Print("]"); firstRow := false;
    od;
    Print("]");
end;

PrintMatrices := function(matrices)
    local first, matrix;
    Print("["); first := true;
    for matrix in matrices do
        if not first then Print(","); fi;
        PrintMatrix(matrix); first := false;
    od;
    Print("]");
end;

PrintMap := function(map)
    Print("["); PrintIntegerList(DimensionVector(Source(map)));
    Print(","); PrintIntegerList(DimensionVector(Range(map)));
    Print(","); PrintMatrices(MatricesOfPathAlgebraMatModuleHomomorphism(map));
    Print("]");
end;

PrintMapPair := function(pair)
    Print("["); PrintMap(pair[1]); Print(","); PrintMap(pair[2]); Print("]");
end;

PrintMapGrid := function(grid)
    local firstRow, firstMap, row, map;
    Print("["); firstRow := true;
    for row in grid do
        if not firstRow then Print(","); fi;
        Print("["); firstMap := true;
        for map in row do
            if not firstMap then Print(","); fi;
            PrintMap(map); firstMap := false;
        od;
        Print("]"); firstRow := false;
    od;
    Print("]");
end;

PrintCase := function(field)
    local algebra, cat, projectives, injectives, simples, nonprojective,
          noninjective, projectiveSum, injectiveSum, projectiveComplex,
          injectiveComplex, nonprojectiveComplex, noninjectiveComplex,
          identity, components, complement, allRemoved, projectiveCachedBefore,
          injectiveCachedBefore;
    algebra := PathAlgebra(field,
        Quiver(["1", "2"], [["1", "2", "a"]]));
    cat := CatOfRightAlgebraModules(algebra);
    projectives := IndecProjectiveModules(algebra);
    injectives := IndecInjectiveModules(algebra);
    simples := SimpleModules(algebra);
    nonprojective := First(simples, module -> not IsProjectiveModule(module));
    noninjective := First(simples, module -> not IsInjectiveModule(module));
    projectiveSum := DirectSumOfQPAModules(
        [projectives[1], projectives[2], projectives[1]]);
    injectiveSum := DirectSumOfQPAModules(
        [injectives[2], injectives[1], injectives[2]]);
    projectiveComplex := StalkComplex(cat, projectiveSum, 1);
    injectiveComplex := StalkComplex(cat, injectiveSum, -1);
    nonprojectiveComplex := StalkComplex(cat, nonprojective, 0);
    noninjectiveComplex := StalkComplex(cat, noninjective, 0);
    projectiveCachedBefore := HasIsProjectiveComplex(projectiveComplex);
    injectiveCachedBefore := HasIsInjectiveComplex(injectiveComplex);
    identity := IdentityMapping(projectiveSum);
    components := FindAllMapComponents(identity);
    complement := DirectSumMinusSummands(projectiveSum, [2]);
    allRemoved := DirectSumMinusSummands(projectiveSum, [1, 2, 3]);

    Print("{\"projective_complex\":", IsProjectiveComplex(projectiveComplex));
    Print(",\"injective_complex\":", IsInjectiveComplex(injectiveComplex));
    Print(",\"nonprojective_complex\":", IsProjectiveComplex(nonprojectiveComplex));
    Print(",\"noninjective_complex\":", IsInjectiveComplex(noninjectiveComplex));
    Print(",\"cache\":[", projectiveCachedBefore, ",",
          HasIsProjectiveComplex(projectiveComplex), ",", injectiveCachedBefore,
          ",", HasIsInjectiveComplex(injectiveComplex), "]");
    Print(",\"compare_projective\":[");
    PrintMaybeInteger(CompareWithIndecProjective(projectives[1])); Print(",");
    PrintMaybeInteger(CompareWithIndecProjective(projectiveSum)); Print(",");
    PrintMaybeInteger(CompareWithIndecProjective(nonprojective)); Print("]");
    Print(",\"compare_injective\":[");
    PrintMaybeInteger(CompareWithIndecInjective(injectives[2])); Print(",");
    PrintMaybeInteger(CompareWithIndecInjective(injectiveSum)); Print(",");
    PrintMaybeInteger(CompareWithIndecInjective(noninjective)); Print("]");
    Print(",\"projective_descriptions\":[");
    PrintMaybeIntegerList(DescriptionOfProjectiveModule(projectives[1])); Print(",");
    PrintMaybeIntegerList(DescriptionOfProjectiveModule(projectiveSum)); Print(",");
    PrintMaybeIntegerList(DescriptionOfProjectiveModule(ZeroModule(algebra))); Print("]");
    Print(",\"degree_descriptions\":[");
    PrintIntegerList(DescriptionOfProjComplexInDegree(projectiveComplex, 0)); Print(",");
    PrintIntegerList(DescriptionOfProjComplexInDegree(projectiveComplex, 1)); Print(",");
    PrintIntegerList(DescriptionOfProjComplexInDegree(projectiveComplex, 2)); Print(",");
    PrintIntegerList(DescriptionOfInjComplexInDegree(injectiveComplex, -1)); Print("]");
    Print(",\"finite_descriptions\":[");
    PrintIntegerLists(DescriptionOfFiniteProjComplex(projectiveComplex)); Print(",");
    PrintIntegerLists(DescriptionOfFiniteInjComplex(injectiveComplex)); Print("]");
    Print(",\"components\":"); PrintMapGrid(components);
    Print(",\"dual_matrices\":"); PrintMatrices(MatricesOfDualMap(complement[1]));
    Print(",\"complement\":"); PrintMapPair(complement);
    Print(",\"all_removed\":"); PrintMapPair(allRemoved);
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"derivedcat_foundation\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"binary\":"); PrintCase(GF(2));
Print(",\"rational\":"); PrintCase(Rationals);
Print("}\n");

QUIT_GAP(0);
