if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

ComputeStructure := function(M, n)
    local P, dimlistP, revdimlistP, top, K, degrees, origin, i,
          target, dimV, products, zero, i_first, j, f, r, first, s,
          g, gpowers, u, h, temp, i_times_j, ir_first, v, A;

    P := List([0..n], i -> HomOverAlgebraWithBasisFunction(M, TrD(M, i)));
    if Dimension(TrD(M, n + 1)) > 0 then
        Error("Must increase the second argument.\n");
    fi;
    dimlistP := List(P, p -> Length(p[1]));
    revdimlistP := Reversed(dimlistP);
    top := n + 2 - PositionNonZero(revdimlistP);
    P := P{[1..top]};
    dimlistP := dimlistP{[1..top]};
    K := LeftActingDomain(M);
    degrees := [];
    origin := 1;
    for i in [1..top] do
        target := origin + dimlistP[i] - 1;
        Add(degrees, [origin..target]);
        origin := target + 1;
    od;
    dimV := target;
    products := EmptySCTable(dimV, Zero(K));
    zero := Flat(List([1..dimV], u -> [Zero(K), u]));

    for i in [1..top] do
        if i = 1 then i_first := 0;
        else i_first := Maximum(degrees[i - 1]); fi;
        for j in degrees[i] do
            f := P[i][1][j - i_first];
            for r in [1..top] do
                if r = 1 then first := 0;
                else first := Maximum(degrees[r - 1]); fi;
                for s in degrees[r] do
                    if i + r - 1 > top then
                        SetEntrySCTable(products, j, s, zero);
                    else
                        g := P[r][1][s - first];
                        gpowers := g;
                        for u in [1..i - 1] do
                            gpowers := DualOfModuleHomomorphism(gpowers);
                            gpowers := TransposeOfModuleHomomorphism(gpowers);
                        od;
                        h := f * gpowers;
                        temp := P[i + r - 1][2](h);
                        i_times_j := [];
                        if i + r = 2 then ir_first := 0;
                        else ir_first := Maximum(degrees[i + r - 2]); fi;
                        for v in [1..dimV] do
                            if v in degrees[i + r - 1] then
                                Append(i_times_j, [temp[v - ir_first], v]);
                            else
                                Append(i_times_j, [Zero(K), v]);
                            fi;
                        od;
                        SetEntrySCTable(products, j, s, i_times_j);
                    fi;
                od;
            od;
        od;
    od;
    A := AlgebraByStructureConstants(K, products);
    return rec(algebra := A, hom_dimensions := dimlistP, degrees := degrees);
end;

PrintIntegerList := function(values)
    Print("[", JoinStringsWithSeparator(List(values, String), ","), "]");
end;

PrintRationalRow := function(values)
    local first, i, value;
    Print("["); first := true;
    for i in [1..Length(values)] do
        value := values[i];
        if not IsZero(value) then
            if not first then Print(","); fi;
            Print("[", i, ",[", NumeratorRat(value), ",", DenominatorRat(value), "]]");
            first := false;
        fi;
    od;
    Print("]");
end;

PrintPrimeRow := function(values)
    local first, i, value;
    Print("["); first := true;
    for i in [1..Length(values)] do
        value := values[i];
        if not IsZero(value) then
            if not first then Print(","); fi;
            Print("[", i, ",", IntFFE(value), "]");
            first := false;
        fi;
    od;
    Print("]");
end;

PrintStructure := function(data, row_printer)
    local A, basis, first, i, j, coordinates;
    A := data.algebra;
    basis := BasisVectors(Basis(A));
    Print("{\"dimension\":", Dimension(A), ",\"hom_dimensions\":");
    PrintIntegerList(data.hom_dimensions);
    Print(",\"degrees\":["); first := true;
    for i in data.degrees do
        if not first then Print(","); fi;
        PrintIntegerList(i); first := false;
    od;
    Print("],\"nonzero_products\":["); first := true;
    for i in [1..Length(basis)] do
        for j in [1..Length(basis)] do
            coordinates := Coefficients(Basis(A), basis[i] * basis[j]);
            if not IsZero(coordinates) then
                if not first then Print(","); fi;
                Print("{\"left\":", i, ",\"right\":", j, ",\"coordinates\":");
                row_printer(coordinates);
                Print("}"); first := false;
            fi;
        od;
    od;
    Print("]}");
end;

PrintFinalSummary := function(result)
    local algebra, quiver, first, arrow;
    algebra := result[1];
    quiver := QuiverOfPathAlgebra(algebra);
    Print("{\"dimension\":", Dimension(algebra));
    Print(",\"vertices\":", NumberOfVertices(quiver), ",\"arrows\":[");
    first := true;
    for arrow in ArrowsOfQuiver(quiver) do
        if not first then Print(","); fi;
        Print("[", VertexIndex(SourceOfPath(arrow)), ",",
              VertexIndex(TargetOfPath(arrow)), "]");
        first := false;
    od;
    Print("],\"relations\":");
    PrintJsonString(String(RelationsOfAlgebra(algebra)));
    Print("}");
end;

BuildCase := function(field)
    local quiver, algebra, simples, projectives, generator;
    quiver := Quiver(2, [[1, 2, "a"]]);
    algebra := PathAlgebra(field, quiver);
    simples := SimpleModules(algebra);
    projectives := IndecProjectiveModules(algebra);
    generator := DirectSumOfQPAModules([simples[1], simples[2], projectives[1]]);
    return rec(module := generator, structure := ComputeStructure(generator, 1));
end;

rational := BuildCase(Rationals);;
binary := BuildCase(GF(2));;
prime_five := BuildCase(GF(5));;
binary_final := PreprojectiveAlgebra(binary.module, 1);;
prime_five_final := PreprojectiveAlgebra(prime_five.module, 1);;

Print("{\"schema_version\":1,\"kind\":\"pathalg_module_preprojective_structure\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":"); PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":"); PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":"); PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"rational\":"); PrintStructure(rational.structure, PrintRationalRow);
Print(",\"f2\":"); PrintStructure(binary.structure, PrintPrimeRow);
Print(",\"f5\":"); PrintStructure(prime_five.structure, PrintPrimeRow);
Print(",\"final\":{\"f2\":"); PrintFinalSummary(binary_final);
Print(",\"f5\":"); PrintFinalSummary(prime_five_final);
Print("}}\n");

QUIT_GAP(0);
