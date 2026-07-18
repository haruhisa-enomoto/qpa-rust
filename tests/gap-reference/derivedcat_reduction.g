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

PrintIso := function(value)
    if IsEmpty(value) then
        Print("[]");
    else
        Print("["); PrintMap(value[1]); Print(",", value[2], ",",
              value[3], ",", value[4], "]");
    fi;
end;

PrintCase := function(field)
    local algebra, projective, source, range, one, zero, differential,
          noIso, iso, sourceReduced, rangeReduced, reduced, reducedIso,
          sourceZero, rangeZero, reducedZero, tMap, reducedTMap,
          previous, reducedPrevious, previousT, reducedPreviousT;
    algebra := PathAlgebra(field, Quiver(["v"], []));
    projective := IndecProjectiveModules(algebra)[1];
    source := DirectSumOfQPAModules([projective, projective]);
    range := DirectSumOfQPAModules([projective, projective]);
    one := IdentityMapping(projective);
    zero := ZeroMapping(projective, projective);
    differential := MultiplyListsOfMaps(
        DirectSumProjections(source),
        [[one, one], [one, zero]],
        DirectSumInclusions(range));
    noIso := CheckForIsomorphism(ZeroMapping(source, range));
    iso := CheckForIsomorphism(differential);
    sourceReduced := DirectSumMinusSummands(source, [iso[3]]);
    rangeReduced := DirectSumMinusSummands(range, [iso[4]]);
    reduced := ReduceDifferential(differential, iso,
        [Source(sourceReduced[1]), Source(rangeReduced[1])]);
    reducedIso := CheckForIsomorphism(reduced);
    sourceZero := DirectSumMinusSummands(Source(reduced), [reducedIso[3]]);
    rangeZero := DirectSumMinusSummands(Range(reduced), [reducedIso[4]]);
    reducedZero := ReduceDifferential(reduced, reducedIso,
        [Source(sourceZero[1]), Source(rangeZero[1])]);

    tMap := DirectSumProjections(source)[1];
    reducedTMap := ReduceTMap(tMap, differential, iso, sourceReduced[1]);
    previous := DirectSumProjections(range)[1] + DirectSumProjections(range)[2];
    reducedPrevious := ReducePreviousDifferential(previous, rangeReduced[1]);
    previousT := DirectSumProjections(source)[1] + DirectSumProjections(source)[2];
    reducedPreviousT := ReducePreviousTMap(previousT, sourceReduced[1]);

    Print("{\"differential\":"); PrintMap(differential);
    Print(",\"no_iso\":"); PrintIso(noIso);
    Print(",\"iso\":"); PrintIso(iso);
    Print(",\"reduced\":"); PrintMap(reduced);
    Print(",\"reduced_iso\":"); PrintIso(reducedIso);
    Print(",\"reduced_zero\":"); PrintMap(reducedZero);
    Print(",\"reduced_t\":"); PrintMap(reducedTMap);
    Print(",\"reduced_previous\":"); PrintMap(reducedPrevious);
    Print(",\"reduced_previous_t\":"); PrintMap(reducedPreviousT);
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"derivedcat_reduction\",\"oracle\":{");
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
