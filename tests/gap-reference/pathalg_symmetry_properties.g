if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintProperties := function(algebra)
    Print("{\"symmetric\":", IsSymmetricAlgebra(algebra));
    Print(",\"weakly_symmetric\":", IsWeaklySymmetricAlgebra(algebra));
    Print(",\"schurian\":", IsSchurianAlgebra(algebra));
    Print(",\"semicommutative\":", IsSemicommutativeAlgebra(algebra));
    Print("}");
end;

PrintCases := function(field)
    local discrete, a2, loop_parent, dual_numbers, diamond_parent,
          diamond_commutative, diamond_one_path_zero, cycle_parent,
          radical_square_zero_cycle;

    discrete := PathAlgebra(field, Quiver(["d"], []));
    a2 := PathAlgebra(field,
        Quiver(["a1", "a2"], [["a1", "a2", "a"]]));
    loop_parent := PathAlgebra(field,
        Quiver(["l"], [["l", "l", "t"]]));
    dual_numbers := loop_parent / [loop_parent.t * loop_parent.t];

    diamond_parent := PathAlgebra(field,
        Quiver(["q1", "q2", "q3", "q4"],
            [["q1", "q2", "a"], ["q1", "q3", "b"],
             ["q2", "q4", "c"], ["q3", "q4", "d"]]));
    diamond_commutative := diamond_parent /
        [diamond_parent.a * diamond_parent.c - diamond_parent.b * diamond_parent.d];
    diamond_one_path_zero := diamond_parent /
        [diamond_parent.a * diamond_parent.c];

    cycle_parent := PathAlgebra(field,
        Quiver(["c1", "c2", "c3"],
            [["c1", "c2", "x"], ["c2", "c3", "y"],
             ["c3", "c1", "z"]]));
    radical_square_zero_cycle := cycle_parent /
        [cycle_parent.x * cycle_parent.y,
         cycle_parent.y * cycle_parent.z,
         cycle_parent.z * cycle_parent.x];

    Print("{\"path_discrete\":"); PrintProperties(discrete);
    Print(",\"path_a2\":"); PrintProperties(a2);
    Print(",\"path_loop\":{\"symmetric\":", IsSymmetricAlgebra(loop_parent));
    Print(",\"weakly_symmetric\":", IsWeaklySymmetricAlgebra(loop_parent));
    Print("}");
    Print(",\"path_diamond\":"); PrintProperties(diamond_parent);
    Print(",\"quotient_dual_numbers\":"); PrintProperties(dual_numbers);
    Print(",\"quotient_radical_square_zero_cycle\":");
    PrintProperties(radical_square_zero_cycle);
    Print(",\"quotient_diamond_commutative\":");
    PrintProperties(diamond_commutative);
    Print(",\"quotient_diamond_one_path_zero\":");
    PrintProperties(diamond_one_path_zero);
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"pathalg_symmetry_properties\",\"oracle\":{");
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
