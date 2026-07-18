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

PrintComplex := function(complex, degrees)
    local first, degree;
    Print("{\"upper\":"); PrintJsonString(String(UpperBound(complex)));
    Print(",\"lower\":"); PrintJsonString(String(LowerBound(complex)));
    Print(",\"length\":"); PrintJsonString(String(LengthOfComplex(complex)));
    Print(",\"projective\":", IsProjectiveComplex(complex));
    Print(",\"objects\":["); first := true;
    for degree in degrees do
        if not first then Print(","); fi;
        PrintIntegerList(DimensionVector(ObjectOfComplex(complex, degree)));
        first := false;
    od;
    Print("],\"differentials\":["); first := true;
    for degree in degrees do
        if not first then Print(","); fi;
        PrintMap(DifferentialOfComplex(complex, degree)); first := false;
    od;
    Print("]}");
end;

PrintCase := function(field)
    local a2, cat, simple, projective, projectiveComplex, already,
          stalk, stalkResolution, zeroComplex, zeroResolution, loopParent, dual, dualCat,
          dualSimple, dualZeroComplex, dualZeroResolution, degrees;
    a2 := PathAlgebra(field,
        Quiver(["v1", "v2"], [["v1", "v2", "a"]]));
    cat := CatOfRightAlgebraModules(a2);
    simple := RightModuleOverPathAlgebra(a2, [1, 0], []);
    projective := IndecProjectiveModules(a2)[1];
    projectiveComplex := StalkComplex(cat, projective, 0);
    already := ProjectiveResolutionOfComplex(projectiveComplex);
    stalk := StalkComplex(cat, simple, 2);
    stalkResolution := ProjectiveResolutionOfComplex(stalk);
    zeroComplex := FiniteComplex(cat, 1, [ZeroMapping(simple, simple)]);
    zeroResolution := ProjectiveResolutionOfComplex(zeroComplex);

    loopParent := PathAlgebra(field,
        Quiver(["w"], [["w", "w", "x"]]));
    dual := loopParent / [loopParent.x^2];
    dualCat := CatOfRightAlgebraModules(dual);
    dualSimple := RightModuleOverPathAlgebra(dual, [1],
        [["x", [[0]] * One(field)]]);
    dualZeroComplex := FiniteComplex(dualCat, 1,
        [ZeroMapping(dualSimple, dualSimple)]);
    dualZeroResolution := ProjectiveResolutionOfComplex(dualZeroComplex);
    degrees := [-2 .. 6];

    Print("{\"already_identical\":", IsIdenticalObj(already, projectiveComplex));
    Print(",\"already\":"); PrintComplex(already, degrees);
    Print(",\"stalk\":"); PrintComplex(stalkResolution, degrees);
    Print(",\"zero_map\":"); PrintComplex(zeroResolution, degrees);
    Print(",\"dual_zero_map\":"); PrintComplex(dualZeroResolution, degrees);
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"derivedcat_projective_resolution\",\"oracle\":{");
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
