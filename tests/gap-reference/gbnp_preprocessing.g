if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintRationalNP := function(polynomial)
    local i;
    Print("[");
    for i in [1 .. Length(polynomial[1])] do
        if i > 1 then Print(","); fi;
        Print("{\"monomial\":", polynomial[1][i]);
        Print(",\"coefficient\":[", NumeratorRat(polynomial[2][i]), ",",
              DenominatorRat(polynomial[2][i]), "]}");
    od;
    Print("]");
end;

PrintRationalNPList := function(polynomials)
    local i;
    Print("[");
    for i in [1 .. Length(polynomials)] do
        if i > 1 then Print(","); fi;
        PrintRationalNP(polynomials[i]);
    od;
    Print("]");
end;

PrintPrimeNP := function(polynomial)
    local i;
    Print("[");
    for i in [1 .. Length(polynomial[1])] do
        if i > 1 then Print(","); fi;
        Print("{\"monomial\":", polynomial[1][i]);
        Print(",\"coefficient\":", IntFFE(polynomial[2][i]), "}");
    od;
    Print("]");
end;

PrintPrimeNPList := function(polynomials)
    local i;
    Print("[");
    for i in [1 .. Length(polynomials)] do
        if i > 1 then Print(","); fi;
        PrintPrimeNP(polynomials[i]);
    od;
    Print("]");
end;

rational_input := [
    [[[1, 2], [1]], [1, 1]],
    [[[1, 2], [2]], [1, 1]]
];;
rational_reduced := GBNP.ReducePol(rational_input);;
duplicate_reduced := GBNP.ReducePol([rational_input[1], rational_input[1]]);;

GBNP.SetOption("CancellativeMonoid", true);;
GBNP.SetOption("Alphabetsize", 2);;
cancellation_input := [[[1, 2, 1], [1, 1, 1]], [1, -1]];;
cancellation_reduced := GBNP.ReduceCancellation(StructuralCopy(cancellation_input));;
GBNP.ClearOptions();;

F5 := GF(5);;
f5_input := [
    [[[1, 2], [1]], [One(F5), One(F5)]],
    [[[1, 2], [2]], [One(F5), One(F5)]]
];;
f5_reduced := GBNP.ReducePol(f5_input);;

F2 := GF(2);;
f2_input := [
    [[[1, 2], [1]], [One(F2), One(F2)]],
    [[[1, 2], [2]], [One(F2), One(F2)]]
];;
f2_reduced := GBNP.ReducePol(f2_input);;

patterns := [[1, 2], [1], [2, 1]];;
lookup_tree := GBNP.CreateOccurTreePTSLR(patterns, 0, true);;
lookup_word := [2, 1, 2];;
first_lookup := GBNP.OccurInLstT(lookup_word, lookup_tree);;
all_lookup := GBNP.LookUpOccurTreeAllLstPTSLR(lookup_word, lookup_tree, true);;

mutable_tree := GBNP.CreateOccurTreePTSLR([[1], [2]], 0, true);;
GBNP.AddMonToTreePTSLR([1, 2], 2, mutable_tree, true);;
after_insert := GBNP.LookUpOccurTreeAllLstPTSLR([1, 2], mutable_tree, true);;
GBNP.RemoveMonFromTreePTSLR([1, 2], 2, mutable_tree, true);;
after_remove := GBNP.LookUpOccurTreeAllLstPTSLR([1, 2], mutable_tree, true);;

Print("{\"schema_version\":1,\"kind\":\"gbnp_preprocessing\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"gbnp_version\":");
PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"occurrence_tree\":{");
Print("\"first_lookup\":", first_lookup);
Print(",\"all_lookup\":", all_lookup);
Print(",\"after_insert\":", after_insert);
Print(",\"after_remove\":", after_remove, "}");
Print(",\"rational\":{");
Print("\"reduced\":");
PrintRationalNPList(rational_reduced);
Print(",\"duplicate_reduced\":");
PrintRationalNPList(duplicate_reduced);
Print(",\"cancellation\":");
PrintRationalNP(cancellation_reduced);
Print("},\"f5\":{");
Print("\"reduced\":");
PrintPrimeNPList(f5_reduced);
Print("},\"f2\":{");
Print("\"reduced\":");
PrintPrimeNPList(f2_reduced);
Print("}}\n");

QUIT_GAP(0);
