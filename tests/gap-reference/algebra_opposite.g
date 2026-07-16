if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintObjectString := function(value)
    local result, stream;
    result := "";
    stream := OutputTextString(result, true);
    SetPrintFormattingStatus(stream, false);
    PrintTo(stream, value);
    CloseStream(stream);
    return result;
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

DualNumbers := function(field)
    local table;
    table := EmptySCTable(2, Zero(field));
    SetEntrySCTable(table, 1, 1, [One(field), 1]);
    SetEntrySCTable(table, 1, 2, [One(field), 2]);
    SetEntrySCTable(table, 2, 1, [One(field), 2]);
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

PrintOpposite := function(field)
    local algebra, basis, opposite, opposite2, opposite_basis,
          canonical_basis, family, x, y, z, recreated;
    algebra := UpperTriangularTwoByTwo(field);
    basis := Basis(algebra);
    IsCommutative(algebra);
    opposite := OppositeAlgebra(algebra);
    opposite2 := OppositeAlgebra(algebra);
    opposite_basis := Basis(opposite);
    canonical_basis := CanonicalBasis(opposite);
    family := ElementsFamily(FamilyObj(opposite));
    x := BasisVectors(opposite_basis)[1];
    y := BasisVectors(opposite_basis)[2];
    z := BasisVectors(opposite_basis)[3];
    recreated := ObjByExtRep(family, basis[2]);

    Print("{\"is_opposite\":", IsOppositeAlgebra(opposite));
    Print(",\"is_algebra\":", IsAlgebra(opposite));
    Print(",\"underlying_with_one\":", IsAlgebraWithOne(algebra));
    Print(",\"opposite_with_one\":", IsAlgebraWithOne(opposite));
    Print(",\"underlying_identical\":",
        IsIdenticalObj(UnderlyingAlgebra(opposite), algebra));
    Print(",\"fresh\":", not IsIdenticalObj(opposite, opposite2));
    Print(",\"dimension\":", Dimension(opposite));
    Print(",\"finite_dimensional\":", IsFiniteDimensional(opposite));
    Print(",\"basis\":");
    PrintIntRows(List(BasisVectors(opposite_basis),
        element -> Coefficients(basis, element![1])));
    Print(",\"canonical_basis\":");
    PrintIntRows(List(BasisVectors(canonical_basis),
        element -> Coefficients(basis, element![1])));
    Print(",\"coefficients\":", List(Coefficients(opposite_basis, y), IntFFE));
    Print(",\"member\":", y in opposite);
    Print(",\"products\":");
    PrintIntRows(List([x * y, y * x, y * z, z * y],
        element -> Coefficients(basis, element![1])));
    Print(",\"sum\":", IntCoordinates(basis, (y + z)![1]));
    Print(",\"left_scalar\":", IntCoordinates(basis, (2 * y)![1]));
    Print(",\"right_scalar\":", IntCoordinates(basis, (y * 2)![1]));
    Print(",\"negative\":", IntCoordinates(basis, (-y)![1]));
    Print(",\"zero\":", IntCoordinates(basis, Zero(y)![1]));
    Print(",\"equal\":", recreated = y);
    Print(",\"less\":[", x < y, ",", y < x, "]");
    Print(",\"underlying_less\":[", basis[1] < basis[2], ",",
        basis[2] < basis[1], "]");
    Print(",\"external_equal\":", ExtRepOfObj(y) = ExtRepOfObj(basis[2]));
    Print(",\"recreated\":", IntCoordinates(basis, recreated![1]));
    Print(",\"element_print\":"); PrintJsonString(PrintObjectString(y));
    Print(",\"sum_print\":"); PrintJsonString(PrintObjectString(y + z));
    Print(",\"scalar_print\":"); PrintJsonString(PrintObjectString(2 * y));
    Print(",\"negative_print\":"); PrintJsonString(PrintObjectString(-y));
    Print(",\"zero_print\":"); PrintJsonString(PrintObjectString(Zero(y)));
    Print("}");
end;

PrintRationalForms := function()
    local algebra, opposite, basis, x, y;
    algebra := UpperTriangularTwoByTwo(Rationals);
    IsCommutative(algebra);
    opposite := OppositeAlgebra(algebra);
    basis := BasisVectors(Basis(opposite));
    x := basis[1];
    y := basis[2];
    Print("{\"sum\":"); PrintJsonString(PrintObjectString(x + y));
    Print(",\"double\":"); PrintJsonString(PrintObjectString(2 * y));
    Print(",\"half\":"); PrintJsonString(PrintObjectString((1/2) * y));
    Print(",\"negative\":"); PrintJsonString(PrintObjectString(-y));
    Print(",\"mixed\":"); PrintJsonString(PrintObjectString(x - 2 * y));
    Print(",\"zero\":"); PrintJsonString(PrintObjectString(Zero(x)));
    Print("}");
end;

PrintCommutative := function(field)
    local algebra, opposite;
    algebra := DualNumbers(field);
    IsCommutative(algebra);
    opposite := OppositeAlgebra(algebra);
    Print("{\"identical\":", IsIdenticalObj(algebra, opposite));
    Print(",\"has_is_opposite\":", HasIsOppositeAlgebra(opposite));
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"algebra_opposite\",\"oracle\":{");
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
Print("},\"upper_triangular_f5\":");
PrintOpposite(GF(5));
Print(",\"rational_print\":");
PrintRationalForms();
Print(",\"dual_numbers_f5\":");
PrintCommutative(GF(5));
Print("}\n");

QUIT_GAP(0);
