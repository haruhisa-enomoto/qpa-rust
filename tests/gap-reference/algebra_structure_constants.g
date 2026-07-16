if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

IdentityCoordinates := function(field)
    local table, algebra, identity;
    table := EmptySCTable(2, Zero(field));
    SetEntrySCTable(table, 1, 1, [One(field), 1]);
    SetEntrySCTable(table, 1, 2, [One(field), 2]);
    algebra := AlgebraByStructureConstants(field, table);
    identity := FindMultiplicativeIdentity(algebra);
    return Coefficients(Basis(algebra), identity);
end;

NoIdentity := function(field)
    local table, algebra;
    table := EmptySCTable(1, Zero(field));
    algebra := AlgebraByStructureConstants(field, table);
    return FindMultiplicativeIdentity(algebra) = fail;
end;

DualNumbers := function(field)
    local table;
    table := EmptySCTable(2, Zero(field));
    SetEntrySCTable(table, 1, 1, [One(field), 1]);
    SetEntrySCTable(table, 1, 2, [One(field), 2]);
    SetEntrySCTable(table, 2, 1, [One(field), 2]);
    return AlgebraByStructureConstants(field, table);
end;

ShiftedDualNumbers := function(field)
    local table;
    # b1 = 1 + epsilon and b2 = epsilon, so b1^2 = b1 + b2.
    table := EmptySCTable(2, Zero(field));
    SetEntrySCTable(table, 1, 1, [One(field), 1, One(field), 2]);
    SetEntrySCTable(table, 1, 2, [One(field), 2]);
    SetEntrySCTable(table, 2, 1, [One(field), 2]);
    return AlgebraByStructureConstants(field, table);
end;

UpperTriangularTwoByTwo := function(field)
    local table;
    table := EmptySCTable(3, Zero(field));
    SetEntrySCTable(table, 1, 1, [One(field), 1]);
    SetEntrySCTable(table, 1, 2, [One(field), 2]);
    SetEntrySCTable(table, 2, 3, [One(field), 2]);
    SetEntrySCTable(table, 3, 3, [One(field), 3]);
    return AlgebraByStructureConstants(field, table);
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

ProductOfTwoFields := function(field)
    local table;
    table := EmptySCTable(2, Zero(field));
    SetEntrySCTable(table, 1, 1, [One(field), 1]);
    SetEntrySCTable(table, 2, 2, [One(field), 2]);
    return AlgebraByStructureConstants(field, table);
end;

IntCoordinates := function(basis, element)
    return List(Coefficients(basis, element), IntFFE);
end;

PrintIntRows := function(rows)
    local row;
    Print("[");
    for row in [1..Length(rows)] do
        if row > 1 then
            Print(",");
        fi;
        Print(List(rows[row], IntFFE));
    od;
    Print("]");
end;

PrintPowerSubalgebra := function(algebra)
    local subalgebra;
    subalgebra := PowerSubalgebra(algebra, 1);
    PrintIntRows(List(BasisVectors(Basis(subalgebra)),
        element -> Coefficients(Basis(algebra), element)));
end;

PrintIdempotentLifting := function(field)
    local algebra, basis, ideal, quotient_map, quotient_idempotent,
          lifting, lift, complete;
    algebra := ShiftedDualNumbers(field);
    basis := Basis(algebra);
    ideal := Ideal(algebra, [basis[2]]);
    quotient_map := NaturalHomomorphismByIdeal(algebra, ideal);
    quotient_idempotent := One(Range(quotient_map));
    lifting := LiftingIdempotent(quotient_map, quotient_idempotent);
    lift := LiftIdempotent(quotient_map, quotient_idempotent);
    complete := LiftingCompleteSetOfOrthogonalIdempotents(
        quotient_map, [quotient_idempotent]);
    Print("{\"quotient_idempotent\":");
    Print(List(Coefficients(Basis(Range(quotient_map)), quotient_idempotent), IntFFE));
    Print(",\"lifting\":");
    Print(IntCoordinates(basis, lifting));
    Print(",\"lift\":");
    Print(IntCoordinates(basis, lift));
    Print(",\"complete\":");
    PrintIntRows(List(complete, element -> Coefficients(basis, element)));
    Print("}");
end;

PrintOrthogonalIdempotentLifting := function(field)
    local algebra, basis, ideal, quotient_map, idempotents, pair, complete;
    algebra := ShiftedUpperTriangularTwoByTwo(field);
    basis := Basis(algebra);
    ideal := Ideal(algebra, [basis[3]]);
    quotient_map := NaturalHomomorphismByIdeal(algebra, ideal);
    # Reverse the quotient-basis order so g*v is nonzero in the second lift.
    idempotents := [ImageElm(quotient_map, basis[2]),
                    ImageElm(quotient_map, basis[1])];
    pair := LiftTwoOrthogonalIdempotents(
        quotient_map, basis[2], idempotents[2]);
    complete := LiftingCompleteSetOfOrthogonalIdempotents(
        quotient_map, idempotents);
    Print("{\"idempotents\":");
    PrintIntRows(List(idempotents,
        element -> Coefficients(Basis(Range(quotient_map)), element)));
    Print(",\"pair\":");
    PrintIntRows(List(pair, element -> Coefficients(basis, element)));
    Print(",\"complete\":");
    PrintIntRows(List(complete, element -> Coefficients(basis, element)));
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"algebra_structure_constants\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":");
PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"rational\":{\"identity\":");
Print(IdentityCoordinates(Rationals));
Print(",\"no_identity\":", NoIdentity(Rationals));
Print("},\"f2\":{\"identity\":[");
Print(JoinStringsWithSeparator(List(IdentityCoordinates(GF(2)), x -> String(IntFFE(x))), ","));
Print("],\"no_identity\":", NoIdentity(GF(2)));
Print("},\"power_subalgebra\":{\"dual_numbers_f5\":");
PrintPowerSubalgebra(DualNumbers(GF(5)));
Print(",\"upper_triangular_f2\":");
PrintPowerSubalgebra(UpperTriangularTwoByTwo(GF(2)));
Print(",\"product_f2\":");
PrintPowerSubalgebra(ProductOfTwoFields(GF(2)));
Print("},\"idempotent_lifting\":{\"shifted_dual_numbers_f5\":");
PrintIdempotentLifting(GF(5));
Print(",\"shifted_upper_triangular_f5\":");
PrintOrthogonalIdempotentLifting(GF(5));
Print("}}\n");

QUIT_GAP(0);
