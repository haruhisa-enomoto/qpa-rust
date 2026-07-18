if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    escaped := ReplacedString(escaped, "\n", "\\n");
    Print("\"", escaped, "\"");
end;

d4desc := [
    [0], [1,0], [1,0], [1,0],
    [2,3,4,1], [5,2], [5,3], [5,4],
    [6,7,8,5], [9,6], [9,7], [9,8]
];;
d4 := ARQuiverNumerical(12, 4, d4desc);;
orbitdesc := ["orbits", [
    [0,[4,0]],
    [1,[3,2]],
    [3,[4,0],[2,1]],
    [1,[3,1],[3,3]]
]];;
orbit := ARQuiverNumerical(9, 2, orbitdesc);;
a3 := ARQuiverNumerical(5, 3, [[0],[1,0],[4,0],[2,1],[3,4]]);;

indec := 5;;
decomposable := [0,0,0,1,0,0,0,1,0,0,0,0];;
modules_number := ModulesOfDimVect(d4, indec);;
modules_vector := ModulesOfDimVect(d4, DimensionVector(d4, indec));;
predecessors := DegOrderPredecessors(d4, decomposable);;
direct_predecessors := DegOrderDirectPredecessors(d4, decomposable);;
predecessors_with_direct := DegOrderPredecessorsWithDirect(d4, decomposable);;
successors := DegOrderSuccessors(d4, indec);;
direct_successors := DegOrderDirectSuccessors(d4, indec);;
successors_with_direct := DegOrderSuccessorsWithDirect(d4, indec);;
different_dimension_leq := DegOrderLEQ(d4, 1, 2);;
conjecture := TestCodimConjecture(a3, false, false);;
conjecture_counts := List([1..5], i ->
    List(DegOrderPredecessorsWithDirect(a3, i), Length));;

Print("{\"schema_version\":1,\"kind\":\"degorderfinitetype\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"mesh\":{");
Print("\"number_of_indecomposables\":", NumberOfIndecomposables(d4));
Print(",\"number_of_projectives\":", NumberOfProjectives(d4));
Print(",\"dim_hom_matrix\":", d4!.DimHomMat);
Print(",\"simples\":", d4!.Simples);
Print(",\"display\":"); PrintJsonString(StringPrint(d4));
Print(",\"dimension_vectors\":", List([1..12], i -> DimensionVector(d4, i)));
Print(",\"dimension_vector_of_vector\":", DimensionVector(d4, decomposable));
Print(",\"modules_number\":", modules_number);
Print(",\"modules_vector\":", modules_vector);
Print(",\"predecessors\":", predecessors);
Print(",\"direct_predecessors\":", direct_predecessors);
Print(",\"predecessors_with_direct\":", predecessors_with_direct);
Print(",\"successors\":", successors);
Print(",\"direct_successors\":", direct_successors);
Print(",\"successors_with_direct\":", successors_with_direct);
Print(",\"dim_hom_dispatch\":[", DimHom(d4, indec, indec), ",",
    DimHom(d4, indec, decomposable), ",",
    DimHom(d4, decomposable, indec), ",",
    DimHom(d4, decomposable, decomposable), "]");
Print(",\"dim_end\":[", DimEnd(d4, indec), ",", DimEnd(d4, decomposable), "]");
Print(",\"orbit_dimension\":[", OrbitDim(d4, indec), ",",
    OrbitDim(d4, decomposable), "]");
Print(",\"orbit_codimension\":", OrbitCodim(d4, indec, decomposable));
Print(",\"degeneration_leq\":[", DegOrderLEQ(d4, indec, decomposable), ",",
    DegOrderLEQ(d4, decomposable, indec), "]");
Print(",\"degeneration_leq_nc\":[", DegOrderLEQNC(d4, indec, decomposable), ",",
    DegOrderLEQNC(d4, decomposable, indec), "]");
Print(",\"different_dimension_leq\":", different_dimension_leq);
Print("},\"orbit\":{");
Print("\"dim_hom_matrix\":", orbit!.DimHomMat);
Print(",\"simples\":", orbit!.Simples);
Print("},\"conjecture\":{");
Print("\"counterexamples\":", conjecture);
Print(",\"predecessor_direct_counts\":", conjecture_counts);
Print("}}\n");

QUIT_GAP(0);
