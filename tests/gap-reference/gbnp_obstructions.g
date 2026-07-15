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

leading := [[1, 2, 1], [2, 1, 2]];;
self_first := GBNP.SelfObs(1, leading);;
self_second := GBNP.SelfObs(2, leading);;
left_second := GBNP.LeftObs(2, leading, self_second);;
right_second := GBNP.RightObs(2, leading, self_second);;

rational_basis := [
    [[[1, 2, 1], [1]], [1, -1]],
    [[[2, 1, 2], [2]], [1, 2]]
];;
rational_spoly := GBNP.Spoly(left_second[1], rational_basis);;

F5 := GF(5);;
f5_basis := [
    [[[1, 2, 1], [1]], [One(F5), -One(F5)]],
    [[[2, 1, 2], [2]], [One(F5), 2 * One(F5)]]
];;
f5_spoly := GBNP.Spoly(left_second[1], f5_basis);;

F2 := GF(2);;
f2_basis := [
    [[[1, 2, 1], [1]], [One(F2), One(F2)]],
    [[[2, 1, 2], [2]], [One(F2), Zero(F2)]]
];;
f2_spoly := GBNP.Spoly(left_second[1], f2_basis);;

Print("{\"schema_version\":1,\"kind\":\"gbnp_obstructions\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"gbnp_version\":");
PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"leading_monomials\":", leading);
Print(",\"self_first\":", self_first);
Print(",\"self_second\":", self_second);
Print(",\"left_second\":", left_second);
Print(",\"right_second\":", right_second);
Print(",\"rational_spolynomial\":");
PrintRationalNP(rational_spoly);
Print(",\"f5_spolynomial\":");
PrintPrimeNP(f5_spoly);
Print(",\"f2_spolynomial\":");
PrintPrimeNP(f2_spoly);
Print("}\n");

QUIT_GAP(0);
