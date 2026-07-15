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

p1 := [[[1, 1, 2], []], [1, -1]];;
p2 := [[[1, 2, 2], []], [1, -1]];;
rational_input := [p1, p2];;
rational_phase1 := GBNP.ReducePol(rational_input);;
funcs := ShallowCopy(GBNP.SGrobnerLoopRec);;
funcs.pg := GBNP.CalculatePG(rational_phase1);;
rational_initial_todo := GBNP.AllObs(rational_phase1, funcs);;
rational_full := SGrobner(rational_input);;
rational_one := SGrobner(rational_input, 1);;

F5 := GF(5);;
f5_input := [
    [[[1, 1, 2], []], [One(F5), -One(F5)]],
    [[[1, 2, 2], []], [One(F5), -One(F5)]]
];;
f5_full := SGrobner(f5_input);;

F2 := GF(2);;
f2_input := [
    [[[1, 1, 2], []], [One(F2), One(F2)]],
    [[[1, 2, 2], []], [One(F2), One(F2)]]
];;
f2_full := SGrobner(f2_input);;

Print("{\"schema_version\":1,\"kind\":\"gbnp_sgrobner\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"gbnp_version\":");
PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"rational\":{");
Print("\"phase1\":");
PrintRationalNPList(rational_phase1);
Print(",\"initial_todo\":");
PrintRationalNPList(rational_initial_todo);
Print(",\"full\":");
PrintRationalNPList(rational_full);
Print(",\"one_iteration\":{\"basis\":");
PrintRationalNPList(rational_one.G);
Print(",\"todo\":");
PrintRationalNPList(rational_one.todo);
Print(",\"iterations\":", rational_one.iterations);
if rational_one.completed then
    Print(",\"completed\":true}");
else
    Print(",\"completed\":false}");
fi;
Print(",\"empty\":");
PrintRationalNPList(SGrobner([]));
Print(",\"identity\":");
PrintRationalNPList(SGrobner([[[[]], [1]]]));
Print("},\"f5\":{\"full\":");
PrintPrimeNPList(f5_full);
Print("},\"f2\":{\"full\":");
PrintPrimeNPList(f2_full);
Print("}}\n");

QUIT_GAP(0);
