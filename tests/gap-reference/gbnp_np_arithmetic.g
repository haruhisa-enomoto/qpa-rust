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
        if i > 1 then
            Print(",");
        fi;
        Print("{\"monomial\":", polynomial[1][i]);
        Print(",\"coefficient\":[", NumeratorRat(polynomial[2][i]), ",",
              DenominatorRat(polynomial[2][i]), "]}");
    od;
    Print("]");
end;

PrintPrimeNP := function(polynomial)
    local i;
    Print("[");
    for i in [1 .. Length(polynomial[1])] do
        if i > 1 then
            Print(",");
        fi;
        Print("{\"monomial\":", polynomial[1][i]);
        Print(",\"coefficient\":", IntFFE(polynomial[2][i]), "}");
    od;
    Print("]");
end;

p1 := [[[1, 1, 2], []], [1, -3]];;
p2 := [[[1, 2, 2], []], [1, -4]];;
dirty := [[[1, 1, 2], [], [1, 1, 2], []], [1, -3, -2, 3]];;
monic_input := [[[1, 1, 2], []], [2, -1]];;
module_polynomial := [[[-2, 1], [-3, 2]], [1, -1]];;
normal_input := [[[1, 1, 1, 2], [2, 1], []], [1, -1, 3]];;

F5 := GF(5);;
f5_cancelled := CleanNP([[[1], [1]], [2 * One(F5), 3 * One(F5)]]);;
f5_monic := MkMonicNP([[[2, 1], [1]], [2 * One(F5), 4 * One(F5)]]);;
f5_reducer := [[[2, 1], [1, 2]], [One(F5), 2 * One(F5)]];;
f5_normal := StrongNormalFormNP([[[2, 1]], [One(F5)]], [f5_reducer]);;

F2 := GF(2);;
f2_cancelled := CleanNP([[[1], [1]], [One(F2), One(F2)]]);;
f2_clean := CleanNP([[[2, 1], [], [2, 1]],
                     [One(F2), One(F2), One(F2)]]);;
f2_reducer := [[[2, 1], [1, 2]], [One(F2), One(F2)]];;
f2_normal := StrongNormalFormNP([[[2, 1]], [One(F2)]], [f2_reducer]);;

Print("{\"schema_version\":1,\"kind\":\"gbnp_np_arithmetic\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"gbnp_version\":");
PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"rational\":{");
Print("\"clean\":");
PrintRationalNP(CleanNP(dirty));
Print(",\"add_p1_minus_p2\":");
PrintRationalNP(AddNP(p1, p2, 1, -1));
Print(",\"bimul_right_2\":");
PrintRationalNP(BimulNP([], p1, [2]));
Print(",\"bimul_left_1\":");
PrintRationalNP(BimulNP([1], p2, []));
Print(",\"multiply\":");
PrintRationalNP(MulNP(p1, p2));
Print(",\"commutator\":");
PrintRationalNP(AddNP(MulNP(p1, p2), MulNP(p2, p1), 1, -1));
Print(",\"monic\":");
PrintRationalNP(MkMonicNP(monic_input));
Print(",\"leading_monomial\":", LMonNP(p1));
Print(",\"leading_term\":");
PrintRationalNP(LTermNP(p1));
Print(",\"generator_counts\":{");
Print("\"algebra\":", NumAlgGensNP(module_polynomial));
Print(",\"module\":", NumModGensNP(module_polynomial), "}");
Print(",\"occurrences\":{");
Print("\"occur_ab_in_bab\":", GBNP.Occur([1, 2], [2, 1, 2]));
Print(",\"left_ab_in_bab\":", GBNP.LeftOccur([1, 2], [2, 1, 2]));
Print(",\"right_ab_in_bab\":", GBNP.RightOccur([1, 2], [2, 1, 2]));
Print(",\"empty_in_bab\":", GBNP.Occur([], [2, 1, 2]));
Print(",\"first_in_list\":", GBNP.OccurInLst([2, 1, 2], [[1, 2], [2]]), "}");
Print(",\"strong_normal_form\":");
PrintRationalNP(StrongNormalFormNP(normal_input, [p1, p2]));
Print("},\"f5\":{");
Print("\"cancelled\":");
PrintPrimeNP(f5_cancelled);
Print(",\"monic\":");
PrintPrimeNP(f5_monic);
Print(",\"strong_normal_form\":");
PrintPrimeNP(f5_normal);
Print("},\"f2\":{");
Print("\"cancelled\":");
PrintPrimeNP(f2_cancelled);
Print(",\"clean\":");
PrintPrimeNP(f2_clean);
Print(",\"strong_normal_form\":");
PrintPrimeNP(f2_normal);
Print("}}\n");

QUIT_GAP(0);
