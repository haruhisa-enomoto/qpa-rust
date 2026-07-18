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

PrintComplex := function(complex, degrees, property)
    local first, degree;
    Print("{\"upper\":"); PrintJsonString(String(UpperBound(complex)));
    Print(",\"lower\":"); PrintJsonString(String(LowerBound(complex)));
    Print(",\"length\":"); PrintJsonString(String(LengthOfComplex(complex)));
    Print(",\"property\":", property(complex));
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
    local a2, cat, projectives, sum, identity, starDecomp, finite,
          finiteInjective, stalk, stalkInjective, loopParent, dual, dualCat,
          dualProjective, generator, xmap, starIndec, compute, infinite,
          infiniteInjective, semisimple, semisimpleCat, semisimpleProjective,
          tau, degrees;
    a2 := PathAlgebra(field,
        Quiver(["v1", "v2"], [["v1", "v2", "a"]]));
    cat := CatOfRightAlgebraModules(a2);
    projectives := IndecProjectiveModules(a2);
    sum := DirectSumOfQPAModules(projectives);
    identity := IdentityMapping(sum);
    starDecomp := StarOfMapBetweenProjectives(identity, [1, 2], [1, 2]);
    finite := FiniteComplex(cat, 1, [identity]);
    finiteInjective := ProjectiveToInjectiveComplex(finite);
    stalk := StalkComplex(cat, projectives[1], 2);
    stalkInjective := ProjectiveToInjectiveFiniteComplex(stalk);

    loopParent := PathAlgebra(field,
        Quiver(["w"], [["w", "w", "x"]]));
    dual := loopParent / [loopParent.x^2];
    dualCat := CatOfRightAlgebraModules(dual);
    dualProjective := IndecProjectiveModules(dual)[1];
    generator := MinimalGeneratingSetOfModule(dualProjective)[1];
    xmap := HomFromProjective(generator^dual.x, dualProjective);
    starIndec := StarOfMapBetweenProjectives(xmap, [1], [1]);
    compute := function(C, i)
        return xmap;
    end;
    infinite := Complex(dualCat, 0, [xmap], ["pos", compute, true], "zero");
    infiniteInjective := ProjectiveToInjectiveComplex(infinite);

    semisimple := PathAlgebra(field, Quiver(["s"], []));
    semisimpleCat := CatOfRightAlgebraModules(semisimple);
    semisimpleProjective := IndecProjectiveModules(semisimple)[1];
    tau := TauOfComplex(StalkComplex(semisimpleCat, semisimpleProjective, 0));
    degrees := [-2 .. 5];

    Print("{\"star_decomposable\":"); PrintMap(starDecomp);
    Print(",\"star_indecomposable\":"); PrintMap(starIndec);
    Print(",\"finite\":"); PrintComplex(finiteInjective, degrees, IsInjectiveComplex);
    Print(",\"stalk\":"); PrintComplex(stalkInjective, degrees, IsInjectiveComplex);
    Print(",\"infinite\":"); PrintComplex(infiniteInjective, degrees, IsInjectiveComplex);
    Print(",\"tau\":"); PrintComplex(tau, degrees, IsProjectiveComplex);
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"derivedcat_tau\",\"oracle\":{");
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
