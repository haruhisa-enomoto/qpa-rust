if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

SizeScreen([100000, 100000]);

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintWalks := function(paths)
    local firstPath, firstArrow, path, arrow;
    Print("[");
    firstPath := true;
    for path in paths do
        if not firstPath then Print(","); fi;
        firstPath := false;
        Print("[");
        firstArrow := true;
        for arrow in WalkOfPath(path) do
            if not firstArrow then Print(","); fi;
            firstArrow := false;
            PrintJsonString(String(arrow));
        od;
        Print("]");
    od;
    Print("]");
end;

PermutationImages := function(permutation, size)
    return List([1..size], point -> point^permutation);
end;

PrintGenus := function(genus)
    Print("[", NumeratorRat(genus), ",", DenominatorRat(genus), "]");
end;

PrintMap := function(map)
    Print("{\"size\":", Size(map));
    Print(",\"ordering\":", PermutationImages(
        OrderingOfCombinatorialMap(map), Size(map)));
    Print(",\"pairing\":", PermutationImages(
        PairingOfCombinatorialMap(map), Size(map)));
    Print(",\"marked\":", MarkedEdgesOfCombinatorialMap(map));
    Print(",\"faces\":", FacesOfCombinatorialMap(map));
    Print(",\"genus\":"); PrintGenus(GenusOfCombinatorialMap(map));
    Print("}");
end;

PrintCase := function(algebra)
    local paths;
    paths := MaximalPathsOfGentleAlgebra(algebra);
    Print("{\"maximal_paths\":"); PrintWalks(paths);
    Print(",\"map\":"); PrintMap(CombinatorialMapOfGentleAlgebra(algebra));
    Print("}");
end;

PrintFieldCases := function(field)
    local Q, A, chain, discrete, loop, square, branch;

    Q := Quiver(3, [[1,2,"a"], [2,3,"b"]]);
    A := PathAlgebra(field, Q);
    chain := A;

    discrete := PathAlgebra(field, Quiver(2, []));

    Q := Quiver(1, [[1,1,"x"]]);
    A := PathAlgebra(field, Q);
    loop := A;
    square := A/[A.x^2];

    Q := Quiver(4, [[1,2,"a"], [2,3,"b"], [2,4,"c"]]);
    A := PathAlgebra(field, Q);
    branch := A/[A.a*A.b];

    Print("{\"chain\":"); PrintCase(chain);
    Print(",\"discrete\":"); PrintCase(discrete);
    Print(",\"square_loop\":{\"maximal_paths\":");
    PrintWalks(MaximalPathsOfGentleAlgebra(square)); Print("}");
    Print(",\"branch\":"); PrintCase(branch);
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"combinatorialmap_gentle\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"rational\":"); PrintFieldCases(Rationals);
Print(",\"f2\":"); PrintFieldCases(GF(2));
Print("}\n");

QUIT_GAP(0);
