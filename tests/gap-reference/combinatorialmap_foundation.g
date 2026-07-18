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
    Print(",\"marked_boundaries\":",
        MarkedBoundariesOfCombinatorialMap(map));
    Print("}");
end;

torus := CombinatorialMap(4, (1,2,3,4), (1,3)(2,4), [1,3]);
dual := DualOfCombinatorialMap(torus);
removed := RemoveEdgeOfCombinatorialMap(torus, [1,3]);
searchmap := CombinatorialMap(
    8, (1,2,3)(4,5)(6,7), (1,4)(2,6)(3,8)(5,7), [1,5,8]);
fixed := CombinatorialMap(3, (), (1,2), [0,4]);

Print("{\"schema_version\":1,\"kind\":\"combinatorialmap_foundation\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"torus\":"); PrintMap(torus);
Print(",\"dual\":"); PrintMap(dual);
Print(",\"dual_identity\":", IsIdenticalObj(dual, DualOfCombinatorialMap(torus)));
Print(",\"removed\":"); PrintMap(removed);
Print(",\"search_map\":"); PrintMap(searchmap);
Print(",\"searches\":{");
Print("\"depth_1\":", DepthSearchCombinatorialMap(searchmap, 1));
Print(",\"width_1\":", WidthSearchCombinatorialMap(searchmap, 1));
Print(",\"depth_9\":", DepthSearchCombinatorialMap(searchmap, 9));
Print("}");
Print(",\"boundary_curves\":{");
Print("\"torus\":", BoundaryCurvesOfCombinatorialMap(torus));
Print(",\"search_map\":", BoundaryCurvesOfCombinatorialMap(searchmap));
Print("}");
Print(",\"winding\":{");
Print("\"equal\":", WindingNumber(torus, [1,1]));
Print(",\"marked\":", WindingNumber(torus, [1,2]));
Print(",\"rotated\":", WindingNumber(torus, [2,4]));
Print("}");
Print(",\"fixed_pairing\":"); PrintMap(fixed);
Print("}\n");

QUIT_GAP(0);
