if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

SizeScreen([100000, 100000]);

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintPrimeCoefficient := function(coefficient)
    Print(IntFFE(coefficient));
end;

PrintRationalCoefficient := function(coefficient)
    Print("[", NumeratorRat(coefficient), ",", DenominatorRat(coefficient), "]");
end;

PrintElement := function(element, coefficient_printer)
    local data, first, i;
    data := CoefficientsAndMagmaElements(element);
    Print("["); first := true;
    for i in [1, 3 .. Length(data) - 1] do
        if not first then Print(","); fi;
        Print("{\"path\":", ExtRepOfObj(data[i]), ",\"coefficient\":");
        coefficient_printer(data[i + 1]);
        Print("}");
        first := false;
    od;
    Print("]");
end;

PrintElementList := function(elements, coefficient_printer)
    local element, first;
    Print("["); first := true;
    for element in elements do
        if not first then Print(","); fi;
        PrintElement(element, coefficient_printer);
        first := false;
    od;
    Print("]");
end;

PrintElementMatrix := function(matrix, coefficient_printer)
    local first, row;
    Print("["); first := true;
    for row in matrix do
        if not first then Print(","); fi;
        PrintElementList(row, coefficient_printer);
        first := false;
    od;
    Print("]");
end;

PrintElementMatrices := function(matrices, coefficient_printer)
    local first, matrix;
    Print("["); first := true;
    for matrix in matrices do
        if not first then Print(","); fi;
        PrintElementMatrix(matrix, coefficient_printer);
        first := false;
    od;
    Print("]");
end;

PrintModuleVector := function(vector, coefficient_printer)
    PrintElementList(ExtRepOfObj(vector)![1], coefficient_printer);
end;

PrintModuleVectorList := function(vectors, coefficient_printer)
    local first, vector;
    Print("["); first := true;
    for vector in vectors do
        if not first then Print(","); fi;
        PrintModuleVector(vector, coefficient_printer);
        first := false;
    od;
    Print("]");
end;

PrintModuleVectorLists := function(vector_lists, coefficient_printer)
    local first, vectors;
    Print("["); first := true;
    for vectors in vector_lists do
        if not first then Print(","); fi;
        PrintModuleVectorList(vectors, coefficient_printer);
        first := false;
    od;
    Print("]");
end;

PrintProjectiveVertices := function(projectives)
    local first, projective, vertices, vector, raw;
    Print("["); first := true;
    for projective in projectives do
        if not first then Print(","); fi;
        vertices := [];
        for vector in GeneratorsOfAlgebraModule(projective) do
            raw := ExtRepOfObj(vector);
            Add(vertices,
                ExtRepOfObj(LeadingMonomial(raw![1][LeadingPosition(raw)])));
        od;
        Print(vertices);
        first := false;
    od;
    Print("]");
end;

PrintResolution := function(resolution, coefficient_printer)
    local f_lists;
    f_lists := ProjectivesFList(resolution);
    Print("{\"maps\":");
    PrintElementMatrices(Maps(resolution), coefficient_printer);
    Print(",\"f\":");
    PrintModuleVectorLists(f_lists[1], coefficient_printer);
    Print(",\"f_prime\":");
    PrintModuleVectorLists(f_lists[2], coefficient_printer);
    Print(",\"r_projective_vertices\":");
    PrintProjectiveVertices(RProjectives(resolution));
    Print(",\"has_projectives\":", HasProjectives(resolution));
    Print("}");
end;

BuildCase := function(field)
    local quiver, path_algebra, relation, algebra, signed_module,
          signed_term_one, signed_term_three;
    quiver := Quiver(["v1", "v2", "v3"],
        [["v1", "v2", "a"], ["v2", "v3", "b"]]);
    path_algebra := PathAlgebra(field, quiver);
    relation := path_algebra.a * path_algebra.b;
    algebra := path_algebra / [relation];
    signed_module := RightModuleOverPathAlgebra(
        algebra, [2, 1, 0], [["a", [[1], [1]] * One(field)]]);
    signed_term_one := ProjectiveResolutionOfPathAlgebraModule(signed_module, 1);
    signed_term_three := ProjectiveResolutionOfPathAlgebraModule(signed_module, 3);
    return [signed_term_one, signed_term_three];
end;

PrintCase := function(case, coefficient_printer)
    Print("{\"signed_term_one\":");
    PrintResolution(case[1], coefficient_printer);
    Print(",\"signed_term_three\":");
    PrintResolution(case[2], coefficient_printer);
    Print("}");
end;

binary_case := BuildCase(GF(2));
quinary_case := BuildCase(GF(5));
rational_case := BuildCase(Rationals);

Print("{\"schema_version\":1,\"kind\":\"moduleprojres_module_resolution\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"binary\":"); PrintCase(binary_case, PrintPrimeCoefficient);
Print(",\"quinary\":"); PrintCase(quinary_case, PrintPrimeCoefficient);
Print(",\"rational\":"); PrintCase(rational_case, PrintRationalCoefficient);
Print("}\n");

QUIT_GAP(0);
