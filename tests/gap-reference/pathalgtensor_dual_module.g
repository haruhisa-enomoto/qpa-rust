if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintJsonStringList := function(values)
    local first, value;
    Print("[");
    first := true;
    for value in values do
        if not first then Print(","); fi;
        PrintJsonString(value);
        first := false;
    od;
    Print("]");
end;

PrintIntegerList := function(values)
    local first, value;
    Print("[");
    first := true;
    for value in values do
        if not first then Print(","); fi;
        Print(value);
        first := false;
    od;
    Print("]");
end;

PrintMatrix := function(matrix)
    local first_row, first_entry, row, entry;
    Print("{\"shape\":");
    PrintIntegerList(DimensionsMat(matrix));
    Print(",\"rows\":[");
    first_row := true;
    for row in matrix do
        if not first_row then Print(","); fi;
        Print("[");
        first_entry := true;
        for entry in row do
            if not first_entry then Print(","); fi;
            PrintJsonString(String(entry));
            first_entry := false;
        od;
        Print("]");
        first_row := false;
    od;
    Print("]}");
end;

PrintMatrixList := function(matrices)
    local first, matrix;
    Print("[");
    first := true;
    for matrix in matrices do
        if not first then Print(","); fi;
        PrintMatrix(matrix);
        first := false;
    od;
    Print("]");
end;

PrintModule := function(module)
    local algebra, arrows;
    algebra := RightActingAlgebra(module);
    arrows := ArrowsOfQuiver(QuiverOfPathAlgebra(algebra));
    Print("{\"dimension_vector\":");
    PrintIntegerList(DimensionVector(module));
    Print(",\"dimension\":", Dimension(module));
    Print(",\"arrows\":");
    PrintJsonStringList(List(arrows, String));
    Print(",\"matrices\":");
    PrintMatrixList(MatricesOfPathAlgebraModule(module));
    Print("}");
end;

qloop := Quiver(["v"], [["v", "v", "x"]]);;
parent := PathAlgebra(Rationals, qloop);;
rational_algebra := parent / [parent.x * parent.x];;
rational_cached_before :=
    HasDualOfAlgebraAsModuleOverEnvelopingAlgebra(rational_algebra);;
rational_dual :=
    DualOfAlgebraAsModuleOverEnvelopingAlgebra(rational_algebra);;
rational_again :=
    DualOfAlgebraAsModuleOverEnvelopingAlgebra(rational_algebra);;
rational_regular := AlgebraAsModuleOverEnvelopingAlgebra(rational_algebra);;
rational_generic_dual := DualOfModule(rational_regular);;
rational_generic_back := DualOfModule(rational_generic_dual);;
rational_enveloping := EnvelopingAlgebra(rational_algebra);;

qpath := Quiver(["u1", "u2"], [["u1", "u2", "a"]]);;
binary_algebra := PathAlgebra(GF(2), qpath);;
binary_dual := DualOfAlgebraAsModuleOverEnvelopingAlgebra(binary_algebra);;
binary_again := DualOfAlgebraAsModuleOverEnvelopingAlgebra(binary_algebra);;
binary_enveloping := EnvelopingAlgebra(binary_algebra);;

infinite_algebra := PathAlgebra(Rationals, qloop);;
infinite_result :=
    DualOfAlgebraAsModuleOverEnvelopingAlgebra(infinite_algebra);;

Print("{\"schema_version\":1,\"kind\":\"pathalgtensor_dual_module\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":");
PrintJsonString(PackageInfo("GBNP")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"rational_quotient\":{");
Print("\"attribute_identity\":[", rational_cached_before, ",",
      HasDualOfAlgebraAsModuleOverEnvelopingAlgebra(rational_algebra), ",",
      IsIdenticalObj(rational_dual, rational_again), ",",
      IsIdenticalObj(RightActingAlgebra(rational_dual), rational_enveloping), "]");
Print(",\"generic_identity\":[",
      IsIdenticalObj(DualOfModule(rational_regular), rational_generic_dual), ",",
      IsIdenticalObj(rational_generic_back, rational_regular), "]");
Print(",\"generic_dual\":");
PrintModule(rational_generic_dual);
Print(",\"bimodule_dual\":");
PrintModule(rational_dual);
Print("},\"binary_path\":{");
Print("\"attribute_identity\":[",
      IsIdenticalObj(binary_dual, binary_again), ",",
      IsIdenticalObj(RightActingAlgebra(binary_dual), binary_enveloping), "]");
Print(",\"bimodule_dual\":");
PrintModule(binary_dual);
Print("},\"infinite_path\":{");
Print("\"is_fail\":", infinite_result = fail, ",");
Print("\"attribute_bound\":",
      HasDualOfAlgebraAsModuleOverEnvelopingAlgebra(infinite_algebra), ",");
Print("\"regular_bound\":",
      HasAlgebraAsModuleOverEnvelopingAlgebra(infinite_algebra), ",");
Print("\"enveloping_bound\":", HasEnvelopingAlgebra(infinite_algebra));
Print("}}\n");

QUIT_GAP(0);
