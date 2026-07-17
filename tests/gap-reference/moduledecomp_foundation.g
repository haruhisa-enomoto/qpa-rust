if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

CoefficientInteger := function(coefficient)
    if IsFFE(coefficient) then
        return IntFFE(coefficient);
    fi;
    return Int(coefficient);
end;

PrintRows := function(rows)
    local first, row;
    Print("["); first := true;
    for row in rows do
        if not first then Print(","); fi;
        Print(List(row, CoefficientInteger)); first := false;
    od;
    Print("]");
end;

ComplementRows := function(field, rows)
    local space, complement;
    space := VectorSpace(field, rows * One(field));
    complement := ComplementInFullRowSpace(space);
    return List(BasisVectors(Basis(complement)),
                row -> List(row, CoefficientInteger));
end;

ShiftedUpperTriangularTwoByTwo := function(field)
    local table;
    # b1 = e11 + e12, b2 = e22 + e12, and b3 = e12.
    table := EmptySCTable(3, Zero(field));
    SetEntrySCTable(table, 1, 1, [One(field), 1]);
    SetEntrySCTable(table, 1, 2, [2 * One(field), 3]);
    SetEntrySCTable(table, 1, 3, [One(field), 3]);
    SetEntrySCTable(table, 2, 2, [One(field), 2]);
    SetEntrySCTable(table, 3, 2, [One(field), 3]);
    return AlgebraByStructureConstants(field, table);
end;

PrintLiftingCase := function(field)
    local algebra, basis, ideal, map, quotient_basis, quotient_ids,
          lifted, primitive, first, pair;
    algebra := ShiftedUpperTriangularTwoByTwo(field);
    basis := Basis(algebra);
    ideal := Ideal(algebra, [basis[3]]);
    map := NaturalHomomorphismByIdeal(algebra, ideal);
    quotient_basis := Basis(Range(map));
    quotient_ids := [ImageElm(map, basis[2]), ImageElm(map, basis[1])];
    lifted := LiftIdempotentsForDecomposition(map, quotient_ids);
    primitive := IdempotentsForDecomposition(algebra);
    first := LiftIdempotent(map, quotient_ids[1]);
    pair := LiftTwoOrthogonalIdempotents(map, first, quotient_ids[2]);

    Print("{\"quotient_ids\":");
    PrintRows(List(quotient_ids, x -> Coefficients(quotient_basis, x)));
    Print(",\"lifted_for_decomposition\":");
    PrintRows(List(lifted, x -> Coefficients(basis, x)));
    Print(",\"idempotents_for_decomposition\":");
    PrintRows(List(primitive, x -> Coefficients(basis, x)));
    Print(",\"lift_idempotent\":");
    PrintRows([Coefficients(basis, first)]);
    Print(",\"lift_two\":");
    PrintRows(List(pair, x -> Coefficients(basis, x)));
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"moduledecomp_foundation\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"complement\":{\"rational\":");
PrintRows(ComplementRows(Rationals, [[1, 2, 0, 1], [0, 0, 1, 3]]));
Print(",\"f2\":");
PrintRows(ComplementRows(GF(2), [[1, 1, 0, 1], [0, 0, 1, 1]]));
Print(",\"full_f5\":");
PrintRows(ComplementRows(GF(5), IdentityMat(3)));
Print("},\"f5\":"); PrintLiftingCase(GF(5));
Print("}\n");

QUIT_GAP(0);
