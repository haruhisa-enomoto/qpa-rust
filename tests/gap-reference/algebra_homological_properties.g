if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintGlobalDimension := function(value)
    if value = infinity then
        PrintJsonString("infinity");
    else
        Print(value);
    fi;
end;

ProductOfTwoFields := function(field)
    local table;
    table := EmptySCTable(2, Zero(field));
    SetEntrySCTable(table, 1, 1, [One(field), 1]);
    SetEntrySCTable(table, 2, 2, [One(field), 2]);
    return AlgebraByStructureConstants(field, table);
end;

StructureDualNumbers := function(field)
    local table;
    table := EmptySCTable(2, Zero(field));
    SetEntrySCTable(table, 1, 1, [One(field), 1]);
    SetEntrySCTable(table, 1, 2, [One(field), 2]);
    SetEntrySCTable(table, 2, 1, [One(field), 2]);
    return AlgebraByStructureConstants(field, table);
end;

PrintCases := function(field)
    local semisimple_structure, dual_structure, discrete_quiver, discrete,
          a2_quiver, a2, a2_quotient, loop_quiver, loop_path, dual_numbers, a3_quiver,
          dual_numbers_reverse, a3_parent, radical_square_zero_a3,
          initial_global, selfinjective,
          after_selfinjective, semisimple, after_semisimple,
          global_dimension, after_global, hereditary, after_hereditary,
          bounded_dimension;

    semisimple_structure := ProductOfTwoFields(field);
    dual_structure := StructureDualNumbers(field);

    discrete_quiver := Quiver(["d"], []);
    discrete := PathAlgebra(field, discrete_quiver);
    a2_quiver := Quiver(["a1", "a2"], [["a1", "a2", "a"]]);
    a2 := PathAlgebra(field, a2_quiver);
    a2_quotient := a2 / [];
    loop_quiver := Quiver(["l"], [["l", "l", "t"]]);
    loop_path := PathAlgebra(field, loop_quiver);
    dual_numbers := loop_path / [loop_path.t * loop_path.t];
    dual_numbers_reverse := loop_path / [loop_path.t * loop_path.t];
    a3_quiver := Quiver(["r1", "r2", "r3"],
        [["r1", "r2", "x"], ["r2", "r3", "y"]]);
    a3_parent := PathAlgebra(field, a3_quiver);
    radical_square_zero_a3 := a3_parent / [a3_parent.x * a3_parent.y];

    Print("{\"structure_semisimple\":{");
    Print("\"known_before\":", HasIsSemisimpleAlgebra(semisimple_structure));
    semisimple := IsSemisimpleAlgebra(semisimple_structure);
    Print(",\"value\":", semisimple);
    Print(",\"known_after\":", HasIsSemisimpleAlgebra(semisimple_structure));
    Print(",\"global_known_before\":", HasGlobalDimension(semisimple_structure));
    global_dimension := GlobalDimension(semisimple_structure);
    Print(",\"global_dimension\":"); PrintGlobalDimension(global_dimension);
    Print(",\"global_known_after\":", HasGlobalDimension(semisimple_structure));
    Print(",\"finite_global_dimension\":",
          IsFiniteGlobalDimensionAlgebra(semisimple_structure));
    Print("},\"structure_dual_numbers\":{");
    Print("\"known_before\":", HasIsSemisimpleAlgebra(dual_structure));
    semisimple := IsSemisimpleAlgebra(dual_structure);
    Print(",\"value\":", semisimple);
    Print(",\"known_after\":", HasIsSemisimpleAlgebra(dual_structure));
    Print(",\"global_known\":", HasGlobalDimension(dual_structure));

    Print("},\"path_discrete\":{");
    Print("\"selfinjective\":", IsSelfinjectiveAlgebra(discrete));
    Print(",\"semisimple\":", IsSemisimpleAlgebra(discrete));
    Print(",\"hereditary\":", IsHereditaryAlgebra(discrete));
    Print(",\"finite_global_dimension\":",
          IsFiniteGlobalDimensionAlgebra(discrete));
    Print(",\"global_dimension\":"); PrintGlobalDimension(GlobalDimension(discrete));

    Print("},\"path_a2\":{");
    Print("\"selfinjective\":", IsSelfinjectiveAlgebra(a2));
    Print(",\"semisimple\":", IsSemisimpleAlgebra(a2));
    Print(",\"hereditary\":", IsHereditaryAlgebra(a2));
    Print(",\"finite_global_dimension\":", IsFiniteGlobalDimensionAlgebra(a2));
    Print(",\"global_dimension\":"); PrintGlobalDimension(GlobalDimension(a2));

    Print("},\"path_loop\":{");
    Print("\"selfinjective\":", IsSelfinjectiveAlgebra(loop_path));
    Print(",\"hereditary\":", IsHereditaryAlgebra(loop_path));
    Print(",\"finite_global_dimension\":",
          IsFiniteGlobalDimensionAlgebra(loop_path));
    Print(",\"global_dimension\":"); PrintGlobalDimension(GlobalDimension(loop_path));

    semisimple := IsSemisimpleAlgebra(a2_quotient);
    hereditary := IsHereditaryAlgebra(a2_quotient);
    Print("},\"quotient_a2_zero_ideal\":{");
    Print("\"semisimple\":", semisimple);
    Print(",\"hereditary\":", hereditary);
    Print(",\"finite_global_dimension\":",
          IsFiniteGlobalDimensionAlgebra(a2_quotient));
    Print(",\"global_dimension\":");
    PrintGlobalDimension(GlobalDimension(a2_quotient));

    initial_global := HasGlobalDimension(dual_numbers);
    selfinjective := IsSelfinjectiveAlgebra(dual_numbers);
    after_selfinjective := HasGlobalDimension(dual_numbers);
    semisimple := IsSemisimpleAlgebra(dual_numbers);
    after_semisimple := HasGlobalDimension(dual_numbers);
    global_dimension := GlobalDimension(dual_numbers);
    after_global := HasGlobalDimension(dual_numbers);
    hereditary := IsHereditaryAlgebra(dual_numbers);
    Print("},\"quotient_dual_numbers\":{");
    Print("\"global_known_initially\":", initial_global);
    Print(",\"selfinjective\":", selfinjective);
    Print(",\"global_known_after_selfinjective\":", after_selfinjective);
    Print(",\"semisimple\":", semisimple);
    Print(",\"global_known_after_semisimple\":", after_semisimple);
    Print(",\"global_dimension\":"); PrintGlobalDimension(global_dimension);
    Print(",\"global_known_after_global\":", after_global);
    Print(",\"hereditary\":", hereditary);

    semisimple := IsSemisimpleAlgebra(dual_numbers_reverse);
    after_semisimple := HasGlobalDimension(dual_numbers_reverse);
    selfinjective := IsSelfinjectiveAlgebra(dual_numbers_reverse);
    after_selfinjective := HasGlobalDimension(dual_numbers_reverse);
    Print("},\"quotient_dual_numbers_semisimple_first\":{");
    Print("\"semisimple\":", semisimple);
    Print(",\"global_known_after_semisimple\":", after_semisimple);
    Print(",\"selfinjective\":", selfinjective);
    Print(",\"global_known_after_selfinjective\":", after_selfinjective);
    Print(",\"global_dimension\":");
    PrintGlobalDimension(GlobalDimension(dual_numbers_reverse));

    selfinjective := IsSelfinjectiveAlgebra(radical_square_zero_a3);
    semisimple := IsSemisimpleAlgebra(radical_square_zero_a3);
    hereditary := IsHereditaryAlgebra(radical_square_zero_a3);
    after_hereditary := HasGlobalDimension(radical_square_zero_a3);
    bounded_dimension := GlobalDimensionOfAlgebra(radical_square_zero_a3, 2);
    after_global := HasGlobalDimension(radical_square_zero_a3);
    Print("},\"quotient_radical_square_zero_a3\":{");
    Print("\"selfinjective\":", selfinjective);
    Print(",\"semisimple\":", semisimple);
    Print(",\"hereditary\":", hereditary);
    Print(",\"global_known_after_hereditary\":", after_hereditary);
    Print(",\"bounded_global_dimension\":");
    PrintGlobalDimension(bounded_dimension);
    Print(",\"global_known_after_bounded\":", after_global);
    Print("}}");
end;

Print("{\"schema_version\":1,\"kind\":\"algebra_homological_properties\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"f2\":"); PrintCases(GF(2));
Print(",\"rational\":"); PrintCases(Rationals);
Print("}\n");

QUIT_GAP(0);
