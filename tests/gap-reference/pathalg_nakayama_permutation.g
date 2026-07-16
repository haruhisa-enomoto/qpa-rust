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
    local first, value;
    Print("["); first := true;
    for value in values do
        if not first then Print(","); fi;
        Print(value); first := false;
    od;
    Print("]");
end;

PrintIntegerMatrix := function(matrix)
    local first, row;
    Print("["); first := true;
    for row in matrix do
        if not first then Print(","); fi;
        PrintIntegerList(row); first := false;
    od;
    Print("]");
end;

PrintBooleanList := function(values)
    local first, value;
    Print("["); first := true;
    for value in values do
        if not first then Print(","); fi;
        Print(value); first := false;
    od;
    Print("]");
end;

PrintCase := function(algebra)
    local initially_known, permutation, simples, indices, module_images,
          module_identical;

    initially_known := HasNakayamaPermutation(algebra);
    permutation := NakayamaPermutation(algebra);
    Print("{\"initially_known\":", initially_known);
    if permutation = false then
        Print(",\"value\":false");
    else
        simples := SimpleModules(algebra);
        indices := List([1..Length(simples)], permutation[2]);
        module_images := List(simples,
            simple -> DimensionVector(permutation[1](simple)));
        module_identical := List(simples,
            simple -> IsIdenticalObj(permutation[1](simple), simple));
        Print(",\"index_images\":"); PrintIntegerList(indices);
        Print(",\"module_images\":"); PrintIntegerMatrix(module_images);
        Print(",\"module_identical\":"); PrintBooleanList(module_identical);
    fi;
    Print(",\"known_after\":", HasNakayamaPermutation(algebra));
    Print("}");
end;

PrintCases := function(field)
    local discrete, a2, loop_parent, dual_numbers, a3_parent,
          nonselfinjective_quotient, cycle_parent, radical_square_zero_cycle;

    discrete := PathAlgebra(field, Quiver(["d1", "d2", "d3"], []));
    a2 := PathAlgebra(field,
        Quiver(["a1", "a2"], [["a1", "a2", "a"]]));
    loop_parent := PathAlgebra(field,
        Quiver(["l"], [["l", "l", "t"]]));
    dual_numbers := loop_parent / [loop_parent.t * loop_parent.t];
    a3_parent := PathAlgebra(field,
        Quiver(["q1", "q2", "q3"],
            [["q1", "q2", "u"], ["q2", "q3", "v"]]));
    nonselfinjective_quotient := a3_parent / [a3_parent.u * a3_parent.v];
    cycle_parent := PathAlgebra(field,
        Quiver(["c1", "c2", "c3"],
            [["c1", "c2", "x"], ["c2", "c3", "y"],
             ["c3", "c1", "z"]]));
    radical_square_zero_cycle := cycle_parent /
        [cycle_parent.x * cycle_parent.y,
         cycle_parent.y * cycle_parent.z,
         cycle_parent.z * cycle_parent.x];

    Print("{\"path_discrete\":"); PrintCase(discrete);
    Print(",\"path_a2\":"); PrintCase(a2);
    Print(",\"quotient_dual_numbers\":"); PrintCase(dual_numbers);
    Print(",\"quotient_nonselfinjective\":");
    PrintCase(nonselfinjective_quotient);
    Print(",\"quotient_radical_square_zero_cycle\":");
    PrintCase(radical_square_zero_cycle);
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"pathalg_nakayama_permutation\",\"oracle\":{");
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
Print(",\"f5\":"); PrintCases(GF(5));
Print(",\"rational\":"); PrintCases(Rationals);
Print("}\n");

QUIT_GAP(0);
