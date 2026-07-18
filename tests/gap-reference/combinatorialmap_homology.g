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
    Print("}");
end;

PrintTopology := function(map)
    local alpha, cut, cutmap, bound1, bound2, index, beta, joined, basis,
          winding;

    alpha := NonSeperatingCurve(map);
    cut := CutNonSepCombinatorialMap(map, alpha, [1..Size(map)]);
    cutmap := cut[1];
    bound1 := cut[2][1];
    bound2 := cut[2][2];
    index := cut[3];
    beta := JoinCurveCombinatorialMap(cutmap, bound1, bound2, index);
    joined := CutJoinCurveCombinatorialMap(
        cutmap, bound1, bound2, beta, index);
    basis := HomologyBasisOfCombinatorialMap(map);
    winding := List(basis, pair -> [
        WindingNumber(map, pair[1]), WindingNumber(map, pair[2])]);

    Print("{\"map\":"); PrintMap(map);
    Print(",\"nonseparating\":", alpha);
    Print(",\"nonseparating_identity\":",
        IsIdenticalObj(alpha, NonSeperatingCurve(map)));
    Print(",\"cut\":{\"map\":"); PrintMap(cutmap);
    Print(",\"boundaries\":[", bound1, ",", bound2, "]");
    Print(",\"index\":", index, "}");
    Print(",\"join_curve\":", beta);
    Print(",\"cut_join\":{\"map\":"); PrintMap(joined[1]);
    Print(",\"index\":", joined[2], "}");
    Print(",\"homology_basis\":", basis);
    Print(",\"homology_identity\":",
        IsIdenticalObj(basis, HomologyBasisOfCombinatorialMap(map)));
    Print(",\"basis_winding\":", winding, "}");
end;

torus := CombinatorialMap(
    4, (1,3,2,4), (1,2)(3,4), [1]);
genus_two := CombinatorialMap(
    8, (1,3,5,7,2,4,6,8), (1,2)(3,4)(5,6)(7,8), [1,5]);

Print("{\"schema_version\":1,\"kind\":\"combinatorialmap_homology\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"torus\":"); PrintTopology(torus);
Print(",\"genus_two\":"); PrintTopology(genus_two);
Print("}\n");

QUIT_GAP(0);
