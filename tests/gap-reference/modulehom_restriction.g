if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintIntegerList := function(values)
    Print("[", JoinStringsWithSeparator(List(values, String), ","), "]");
end;

PrintPrimeRow := function(row)
    PrintIntegerList(List(row, IntFFE));
end;

PrintRationalRow := function(row)
    local first, entry;
    Print("["); first := true;
    for entry in row do
        if not first then Print(","); fi;
        Print("[", NumeratorRat(entry), ",", DenominatorRat(entry), "]");
        first := false;
    od;
    Print("]");
end;

PrintMatrix := function(matrix, row_printer)
    local first, row;
    Print("["); first := true;
    for row in matrix do
        if not first then Print(","); fi;
        row_printer(row); first := false;
    od;
    Print("]");
end;

PrintModule := function(module, row_printer)
    local first, matrix;
    Print("{\"dimension\":"); PrintIntegerList(DimensionVector(module));
    Print(",\"matrices\":["); first := true;
    for matrix in MatricesOfPathAlgebraModule(module) do
        if not first then Print(","); fi;
        PrintMatrix(matrix, row_printer); first := false;
    od;
    Print("]}");
end;

PrintMap := function(map, row_printer)
    local first, matrix;
    Print("{\"source\":"); PrintModule(Source(map), row_printer);
    Print(",\"range\":"); PrintModule(Range(map), row_printer);
    Print(",\"matrices\":["); first := true;
    for matrix in MatricesOfPathAlgebraMatModuleHomomorphism(map) do
        if not first then Print(","); fi;
        PrintMatrix(matrix, row_printer); first := false;
    od;
    Print("]}");
end;

PrintCase := function(field, row_printer)
    local loop_quiver, loop_parent, dual, dual_generators, zero_arrow_map,
          regular, radical_map, restricted_regular, restricted_radical_map,
          a2_quiver, a2, field_algebra, field_generator, projective,
          identity_map, field_inclusion, restricted_projective,
          restricted_identity;

    loop_quiver := Quiver(["w"], [["w", "w", "x"]]);
    loop_parent := PathAlgebra(field, loop_quiver);
    dual := loop_parent / [loop_parent.x^2];
    dual_generators := GeneratorsOfAlgebraWithOne(dual);
    zero_arrow_map := AlgebraHomomorphismByImages(dual, dual,
        dual_generators, [dual_generators[1], Zero(dual)]);
    zero_arrow_map!.generators := dual_generators;
    zero_arrow_map!.genimages := [dual_generators[1], Zero(dual)];
    regular := RightModuleOverPathAlgebra(dual, [2],
        [["x", [[0, 1], [0, 0]] * One(field)]]);
    radical_map := RightModuleHomOverAlgebra(regular, regular,
        [[[0, 1], [0, 0]] * One(field)]);
    restricted_regular :=
        RestrictionViaAlgebraHomomorphism(zero_arrow_map, regular);
    restricted_radical_map :=
        RestrictionViaAlgebraHomomorphismMap(zero_arrow_map, radical_map);

    a2_quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
    a2 := PathAlgebra(field, a2_quiver);
    field_algebra := PathAlgebra(field, Quiver(["u"], []));
    field_generator := GeneratorsOfAlgebraWithOne(field_algebra);
    field_inclusion := AlgebraHomomorphismByImages(field_algebra, a2,
        field_generator, [One(a2)]);
    field_inclusion!.generators := field_generator;
    field_inclusion!.genimages := [One(a2)];
    projective := RightModuleOverPathAlgebra(a2, [1, 1],
        [["a", [[1]] * One(field)]]);
    identity_map := IdentityMapping(projective);
    restricted_projective :=
        RestrictionViaAlgebraHomomorphism(field_inclusion, projective);
    restricted_identity :=
        RestrictionViaAlgebraHomomorphismMap(field_inclusion, identity_map);

    Print("{\"zero_arrow_module\":");
    PrintModule(restricted_regular, row_printer);
    Print(",\"zero_arrow_map\":");
    PrintMap(restricted_radical_map, row_printer);
    Print(",\"field_module\":");
    PrintModule(restricted_projective, row_printer);
    Print(",\"field_identity_map\":");
    PrintMap(restricted_identity, row_printer);
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"modulehom_restriction\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"binary\":"); PrintCase(GF(2), PrintPrimeRow);
Print(",\"quinary\":"); PrintCase(GF(5), PrintPrimeRow);
Print(",\"rational\":"); PrintCase(Rationals, PrintRationalRow);
Print("}\n");

QUIT_GAP(0);
