if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintResult := function(value)
    if value = false then Print("false"); else Print(value); fi;
end;

RadicalSquareZeroA3 := function(field, prefix)
    local algebra;
    algebra := PathAlgebra(field,
        Quiver([Concatenation(prefix, "1"), Concatenation(prefix, "2"),
                Concatenation(prefix, "3")],
               [[Concatenation(prefix, "1"), Concatenation(prefix, "2"),
                 Concatenation(prefix, "x")],
                [Concatenation(prefix, "2"), Concatenation(prefix, "3"),
                 Concatenation(prefix, "y")]]));
    return algebra / [algebra.(Concatenation(prefix, "x")) *
                      algebra.(Concatenation(prefix, "y"))];
end;

DualNumbers := function(field, prefix)
    local algebra, arrow;
    arrow := Concatenation(prefix, "x");
    algebra := PathAlgebra(field,
        Quiver([Concatenation(prefix, "v")],
               [[Concatenation(prefix, "v"), Concatenation(prefix, "v"), arrow]]));
    return algebra / [algebra.(arrow)^2];
end;

PrintCase := function(field)
    local discrete, a2, radical_square, known_finite, dual, known_infinite,
          setup;
    discrete := PathAlgebra(field, Quiver(["d"], []));
    a2 := PathAlgebra(field,
        Quiver(["v1", "v2"], [["v1", "v2", "a"]]));
    radical_square := RadicalSquareZeroA3(field, "r");
    known_finite := RadicalSquareZeroA3(field, "k");
    dual := DualNumbers(field, "u");
    known_infinite := DualNumbers(field, "i");

    Print("{\"discrete_bound_0\":");
    PrintResult(GorensteinDimensionOfAlgebra(discrete, 0));
    Print(",\"a2_bound_0\":");
    PrintResult(GorensteinDimensionOfAlgebra(a2, 0));
    Print(",\"a2_attribute\":", GorensteinDimension(a2));
    Print(",\"a2_cached_negative\":");
    PrintResult(GorensteinDimensionOfAlgebra(a2, -1));

    Print(",\"radical_initially_bound\":", HasGorensteinDimension(radical_square));
    Print(",\"radical_bound_0\":");
    PrintResult(GorensteinDimensionOfAlgebra(radical_square, 0));
    Print(",\"radical_bound_after_0\":", HasGorensteinDimension(radical_square));
    Print(",\"radical_bound_2\":");
    PrintResult(GorensteinDimensionOfAlgebra(radical_square, 2));
    Print(",\"radical_attribute\":", GorensteinDimension(radical_square));
    Print(",\"radical_cached_bound_0\":");
    PrintResult(GorensteinDimensionOfAlgebra(radical_square, 0));

    setup := GlobalDimensionOfAlgebra(known_finite, 2);
    Print(",\"known_finite_setup\":"); PrintResult(setup);
    Print(",\"known_finite_negative\":");
    PrintResult(GorensteinDimensionOfAlgebra(known_finite, -1));

    Print(",\"dual_bound_0\":");
    PrintResult(GorensteinDimensionOfAlgebra(dual, 0));
    Print(",\"dual_attribute\":", GorensteinDimension(dual));
    setup := GlobalDimensionOfAlgebra(known_infinite, 0);
    Print(",\"known_infinite_setup\":\"", setup, "\"");
    Print(",\"known_infinite_bound_0\":");
    PrintResult(GorensteinDimensionOfAlgebra(known_infinite, 0));
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"modulehomalg_gorenstein_dimension\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"binary\":"); PrintCase(GF(2));
Print(",\"quinary\":"); PrintCase(GF(5));
Print(",\"rational\":"); PrintCase(Rationals);
Print("}\n");

QUIT_GAP(0);
