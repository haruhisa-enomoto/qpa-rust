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
    local first, projective, vertices, vector;
    Print("["); first := true;
    for projective in projectives do
        if not first then Print(","); fi;
        vertices := [];
        for vector in GeneratorsOfAlgebraModule(projective) do
            Add(vertices, ExtRepOfObj(LeadingMonomial(ExtRepOfObj(vector)![1][LeadingPosition(ExtRepOfObj(vector))])));
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
    Print(",\"projective_vertices\":");
    PrintProjectiveVertices(Projectives(resolution));
    Print("}");
end;

PrintCase := function(field, coefficient_printer)
    local q, path_algebra, x, ideal, gb, algebra, resolution, short_resolution,
          q_overlap, overlap_path_algebra, ox, oy, overlap_ideal, overlap_gb,
          overlap_algebra, overlap_resolution;
    q := Quiver(1, [[1, 1, "x"]]);
    path_algebra := PathAlgebra(field, q);
    x := path_algebra.x;
    ideal := Ideal(path_algebra, [x * x]);
    gb := GroebnerBasis(ideal, [x * x]);
    algebra := path_algebra / ideal;
    resolution := ProjectiveResolutionFpPathAlgebraModule(
        algebra, ideal, [[algebra.x]], 4);
    short_resolution := ProjectiveResolutionFpPathAlgebraModule(
        algebra, ideal, [[algebra.x]], 1);

    q_overlap := Quiver(1, [[1, 1, "x"], [1, 1, "y"]]);
    overlap_path_algebra := PathAlgebra(field, q_overlap);
    ox := overlap_path_algebra.x;
    oy := overlap_path_algebra.y;
    overlap_ideal := Ideal(overlap_path_algebra,
        [ox * ox, ox * oy, oy * oy]);
    overlap_gb := GroebnerBasis(overlap_ideal,
        [ox * ox, ox * oy, oy * oy]);
    overlap_algebra := overlap_path_algebra / overlap_ideal;
    overlap_resolution := ProjectiveResolutionFpPathAlgebraModule(
        overlap_algebra, overlap_ideal,
        [[overlap_algebra.x], [overlap_algebra.y]], 3);

    Print("{\"term_four\":");
    PrintResolution(resolution, coefficient_printer);
    Print(",\"term_one\":");
    PrintResolution(short_resolution, coefficient_printer);
    Print(",\"overlap\":");
    PrintResolution(overlap_resolution, coefficient_printer);
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"moduleprojres_resolution\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"binary\":"); PrintCase(GF(2), PrintPrimeCoefficient);
Print(",\"quinary\":"); PrintCase(GF(5), PrintPrimeCoefficient);
Print(",\"rational\":"); PrintCase(Rationals, PrintRationalCoefficient);
Print("}\n");

QUIT_GAP(0);
