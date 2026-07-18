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

PrintMap := function(map)
    local first, matrix;
    Print("["); PrintIntegerList(DimensionVector(Source(map)));
    Print(","); PrintIntegerList(DimensionVector(Range(map)));
    Print(",["); first := true;
    for matrix in MatricesOfPathAlgebraMatModuleHomomorphism(map) do
        if not first then Print(","); fi;
        PrintMatrix(matrix); first := false;
    od;
    Print("]]");
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

PrintComplexSample := function(complex, degrees)
    local degree, first;
    Print("{\"objects\":["); first := true;
    for degree in degrees do
        if not first then Print(","); fi;
        PrintIntegerList(DimensionVector(ObjectOfComplex(complex, degree)));
        first := false;
    od;
    Print("],\"differentials\":");
    PrintMapList(List(degrees, degree -> DifferentialOfComplex(complex, degree)));
    Print("}");
end;

PrintComplexProperties := function(complex)
    local value, finiteness;
    Print("{\"upper\":"); PrintJsonString(String(UpperBound(complex)));
    Print(",\"lower\":"); PrintJsonString(String(LowerBound(complex)));
    Print(",\"finite\":"); finiteness := IsFiniteComplex(complex);
    if finiteness = fail then PrintJsonString("fail"); else Print(finiteness); fi;
    Print(",\"length\":"); PrintJsonString(String(LengthOfComplex(complex)));
    Print(",\"exact\":");
    if finiteness = fail then
        PrintJsonString("property_error");
    else
        value := IsExactSequence(complex); Print(value);
    fi;
    Print(",\"short_exact\":");
    if finiteness = fail then
        PrintJsonString("property_error");
    else
        value := IsShortExactSequence(complex); Print(value);
    fi;
    Print(",\"positive_repeating\":", IsPositiveRepeating(complex));
    Print(",\"negative_repeating\":", IsNegativeRepeating(complex));
    Print(",\"positive_repeat_degrees\":");
    value := PositiveRepeatDegrees(complex);
    if value = fail then PrintJsonString("fail"); else PrintIntegerList(value); fi;
    Print(",\"negative_repeat_degrees\":");
    value := NegativeRepeatDegrees(complex);
    if value = fail then PrintJsonString("fail"); else PrintIntegerList(value); fi;
    Print(",\"highest_known\":"); PrintJsonString(String(HighestKnownDegree(complex)));
    Print(",\"lowest_known\":"); PrintJsonString(String(LowestKnownDegree(complex)));
    Print("}");
end;

PrintChainMapSample := function(map, degrees)
    Print("{\"morphisms\":");
    PrintMapList(List(degrees, degree -> MorphismOfChainMap(map, degree)));
    Print(",\"highest_known\":"); PrintJsonString(String(HighestKnownDegree(map)));
    Print(",\"lowest_known\":"); PrintJsonString(String(LowestKnownDegree(map)));
    Print("}");
end;

PrintHomologyDimensions := function(complex, degrees)
    local first, degree;
    Print("["); first := true;
    for degree in degrees do
        if not first then Print(","); fi;
        PrintIntegerList(DimensionVector(HomologyOfComplex(complex, degree)));
        first := false;
    od;
    Print("]");
end;

PrintCase := function(field)
    local algebra, cat, zero, one, two, f, g, zeroEndo, ses, stalk,
          periodic, calls, positioned, shiftOne, shiftMinusOne, brutalBelow,
          brutalAbove, brutal, goodBelow, goodAbove, good, syzygy, cosyzygy,
          both, yoneda, identityChain, zeroChain, cone, comparison, degrees,
          coneDegrees, cutAbove, cutBelow;
    algebra := PathAlgebra(field, Quiver(["v"], []));
    cat := CatOfRightAlgebraModules(algebra);
    zero := ZeroModule(algebra);
    one := RightModuleOverPathAlgebra(algebra, [1], []);
    two := RightModuleOverPathAlgebra(algebra, [2], []);
    f := RightModuleHomOverAlgebra(one, two, [[[1, 0]] * One(field)]);
    g := RightModuleHomOverAlgebra(two, one, [[[0], [1]] * One(field)]);
    zeroEndo := ZeroMapping(one, one);
    ses := ShortExactSequence(cat, f, g);
    stalk := StalkComplex(cat, one, 2);
    periodic := Complex(cat, 0, [zeroEndo],
                        ["repeat", [zeroEndo]], ["repeat", [zeroEndo]]);
    calls := 0;
    positioned := Complex(cat, 0, [zeroEndo],
        ["pos", function(complex, degree)
                    calls := calls + 1; return zeroEndo;
                end, false],
        ["pos", function(complex, degree)
                    calls := calls + 1; return zeroEndo;
                end, false]);
    DifferentialOfComplex(positioned, 5);
    DifferentialOfComplex(positioned, -5);

    shiftOne := Shift(ses, 1);
    shiftMinusOne := Shift(ses, -1);
    brutalBelow := BrutalTruncationBelow(ses, 0);
    brutalAbove := BrutalTruncationAbove(ses, 0);
    brutal := BrutalTruncation(ses, 0, 0);
    goodBelow := GoodTruncationBelow(ses, 0);
    goodAbove := GoodTruncationAbove(ses, 0);
    good := GoodTruncation(ses, 0, 0);
    syzygy := SyzygyTruncation(ses, 0);
    cosyzygy := CosyzygyTruncation(ses, 0);
    both := SyzygyCosyzygyTruncation(ses, 0, 0);
    yoneda := YonedaProduct(ses, ses);
    identityChain := FiniteChainMap(ses, ses, -1,
        List([-1 .. 1], degree -> IdentityMapping(ObjectOfComplex(ses, degree))));
    zeroChain := ZeroChainMap(ses, ses);
    cone := MappingCone(identityChain);
    comparison := ComparisonLiftingToProjectiveResolution(IdentityMapping(one));
    cutAbove := CutComplexAbove(ProjectiveResolution(one));
    cutBelow := CutComplexBelow(InjectiveResolution(one));
    degrees := [-2 .. 3];
    coneDegrees := [-3 .. 3];

    Print("{\"category_equal\":", cat = CatOfRightAlgebraModules(algebra));
    Print(",\"zero\":"); PrintComplexSample(ZeroComplex(cat), [-1,0,1]);
    Print(",\"ses\":"); PrintComplexSample(ses, degrees);
    Print(",\"ses_homology\":"); PrintHomologyDimensions(ses, [-2 .. 2]);
    Print(",\"stalk\":"); PrintComplexSample(stalk, [0 .. 4]);
    Print(",\"periodic\":"); PrintComplexSample(periodic, [-2 .. 2]);
    Print(",\"position_calls\":", calls);
    Print(",\"position_known\":[");
    PrintJsonString(String(LowestKnownDegree(positioned))); Print(",");
    PrintJsonString(String(HighestKnownDegree(positioned))); Print("]");
    Print(",\"shift_one\":"); PrintComplexSample(shiftOne, degrees);
    Print(",\"shift_minus_one\":"); PrintComplexSample(shiftMinusOne, degrees);
    Print(",\"brutal_below\":"); PrintComplexSample(brutalBelow, degrees);
    Print(",\"brutal_above\":"); PrintComplexSample(brutalAbove, degrees);
    Print(",\"brutal\":"); PrintComplexSample(brutal, degrees);
    Print(",\"good_below\":"); PrintComplexSample(goodBelow, degrees);
    Print(",\"good_above\":"); PrintComplexSample(goodAbove, degrees);
    Print(",\"good\":"); PrintComplexSample(good, degrees);
    Print(",\"syzygy\":"); PrintComplexSample(syzygy, degrees);
    Print(",\"cosyzygy\":"); PrintComplexSample(cosyzygy, degrees);
    Print(",\"syzygy_cosyzygy\":"); PrintComplexSample(both, degrees);
    Print(",\"yoneda\":"); PrintComplexSample(yoneda, [-3 .. 3]);
    Print(",\"identity_chain\":"); PrintChainMapSample(identityChain, degrees);
    Print(",\"zero_chain\":"); PrintChainMapSample(zeroChain, degrees);
    Print(",\"cone\":"); PrintComplexSample(cone[1], coneDegrees);
    Print(",\"cone_inclusion\":"); PrintChainMapSample(cone[2], coneDegrees);
    Print(",\"cone_projection\":"); PrintChainMapSample(cone[3], coneDegrees);
    Print(",\"comparison\":"); PrintChainMapSample(comparison, [-2 .. 2]);
    Print(",\"cut_above\":"); PrintComplexSample(cutAbove, [-2 .. 2]);
    Print(",\"cut_below\":"); PrintComplexSample(cutBelow, [-2 .. 2]);
    Print(",\"properties\":{\"zero\":"); PrintComplexProperties(ZeroComplex(cat));
    Print(",\"ses\":"); PrintComplexProperties(ses);
    Print(",\"stalk\":"); PrintComplexProperties(stalk);
    Print(",\"periodic\":"); PrintComplexProperties(periodic);
    Print(",\"good\":"); PrintComplexProperties(good);
    Print(",\"yoneda\":"); PrintComplexProperties(yoneda);
    Print("}");
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"complex\",\"oracle\":{");
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
