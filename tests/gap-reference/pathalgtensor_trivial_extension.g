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

PrintRelations := function(algebra)
    local first_relation, first_term, relation, data, i;
    Print("[");
    first_relation := true;
    for relation in RelatorsOfFpAlgebra(algebra) do
        if not first_relation then Print(","); fi;
        Print("[");
        data := CoefficientsAndMagmaElements(relation);
        first_term := true;
        for i in [1,3..Length(data)-1] do
            if not first_term then Print(","); fi;
            Print("{\"walk\":");
            PrintJsonStringList(List(WalkOfPath(data[i]), String));
            Print(",\"coefficient\":");
            PrintJsonString(String(data[i+1]));
            Print("}");
            first_term := false;
        od;
        Print("]");
        first_relation := false;
    od;
    Print("]");
end;

PrintAlgebra := function(algebra)
    local quiver, first, arrow;
    quiver := QuiverOfPathAlgebra(algebra);
    Print("{\"level\":", TrivialExtensionOfQuiverAlgebraLevel(algebra));
    Print(",\"dimension\":", Dimension(algebra));
    Print(",\"arrows\":[");
    first := true;
    for arrow in ArrowsOfQuiver(quiver) do
        if not first then Print(","); fi;
        Print("[");
        PrintJsonString(String(SourceOfPath(arrow)));
        Print(",");
        PrintJsonString(String(TargetOfPath(arrow)));
        Print(",");
        PrintJsonString(String(arrow));
        Print("]");
        first := false;
    od;
    Print("],\"relations\":");
    PrintRelations(algebra);
    Print("}");
end;

PrintTopData := function(field)
    local quiver, algebra, module, radical, top;
    quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
    algebra := PathAlgebra(field, quiver);
    module := RightModuleOverPathAlgebra(
        algebra, [2, 2], [["a", [[1,1],[0,0]] * One(field)]]);
    radical := RadicalOfModuleInclusion(module);
    top := TopOfModuleProjection(module);
    Print("{\"radical_dimension\":");
    PrintIntegerList(DimensionVector(Source(radical)));
    Print(",\"radical_map\":");
    PrintMatrixList(MatricesOfPathAlgebraMatModuleHomomorphism(radical));
    Print(",\"top_dimension\":");
    PrintIntegerList(DimensionVector(Range(top)));
    Print(",\"top_arrow_matrices\":");
    PrintMatrixList(MatricesOfPathAlgebraModule(Range(top)));
    Print(",\"top_map\":");
    PrintMatrixList(MatricesOfPathAlgebraMatModuleHomomorphism(top));
    Print(",\"top_lifts\":");
    PrintJsonStringList(List(
        BasisVectors(Basis(Range(top))),
        x -> String(ExtRepOfObj(ExtRepOfObj(
            PreImagesRepresentativeNC(top, x))))));
    Print("}");
end;

qa2 := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);;
rational_a2 := PathAlgebra(Rationals, qa2);;
rational_ta2 := TrivialExtensionOfQuiverAlgebra(rational_a2);;
rational_projection :=
    TrivialExtensionOfQuiverAlgebraProjection(rational_a2);;
rational_gens := GeneratorsOfAlgebraWithOne(rational_a2);;
rational_identity := AlgebraHomomorphismByImages(
    rational_a2, rational_a2, rational_gens, rational_gens);;
rational_identity!.generators := rational_gens;;
rational_identity!.genimages := rational_gens;;
rational_opposite_identity :=
    OppositeAlgebraHomomorphism(rational_identity);;
rational_enveloping_identity :=
    EnvelopingAlgebraHomomorphism(rational_identity);;

binary_a2 := PathAlgebra(GF(2), qa2);;
binary_ta2 := TrivialExtensionOfQuiverAlgebra(binary_a2);;

qloop := Quiver(["v"], [["v", "v", "x"]]);;
rational_parent := PathAlgebra(Rationals, qloop);;
rational_dual_numbers := rational_parent / [rational_parent.x^2];;
rational_tdual := TrivialExtensionOfQuiverAlgebra(rational_dual_numbers);;
binary_parent := PathAlgebra(GF(2), qloop);;
binary_dual_numbers := binary_parent / [binary_parent.x^2];;
binary_tdual := TrivialExtensionOfQuiverAlgebra(binary_dual_numbers);;

infinite := PathAlgebra(Rationals, qloop);;
infinite_extension := TrivialExtensionOfQuiverAlgebra(infinite);;
infinite_projection := TrivialExtensionOfQuiverAlgebraProjection(infinite);;

Print("{\"schema_version\":1,\"kind\":\"pathalgtensor_trivial_extension\",\"oracle\":{");
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
Print("},\"rational_top\":");
PrintTopData(Rationals);
Print(",\"binary_top\":");
PrintTopData(GF(2));
Print(",\"rational_a2\":");
PrintAlgebra(rational_ta2);
Print(",\"binary_a2\":");
PrintAlgebra(binary_ta2);
Print(",\"rational_dual_numbers\":");
PrintAlgebra(rational_tdual);
Print(",\"binary_dual_numbers\":");
PrintAlgebra(binary_tdual);
Print(",\"homomorphisms\":{");
Print("\"trivial_extension_identity\":",
      IsIdenticalObj(rational_ta2,
          TrivialExtensionOfQuiverAlgebra(rational_a2)), ",");
Print("\"projection_zero_pattern\":",
      List(rational_projection!.genimages, IsZero));
Print(",\"opposite_identity\":",
      rational_opposite_identity!.generators =
      rational_opposite_identity!.genimages);
Print(",\"enveloping_identity\":",
      rational_enveloping_identity!.generators =
      rational_enveloping_identity!.genimages);
Print("},\"infinite\":{");
Print("\"level\":", TrivialExtensionOfQuiverAlgebraLevel(infinite), ",");
Print("\"extension_is_fail\":", infinite_extension = fail, ",");
Print("\"projection_is_fail\":", infinite_projection = fail);
Print("}}");
Print("\n");

QUIT_GAP(0);
