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

PrintRawVectorList := function(vectors, coefficient_printer)
    local first, vector;
    Print("["); first := true;
    for vector in vectors do
        if not first then Print(","); fi;
        PrintElementList(vector![1], coefficient_printer);
        first := false;
    od;
    Print("]");
end;

PrintProjectiveVertices := function(projective)
    local first, raw, vector;
    Print("["); first := true;
    for vector in GeneratorsOfAlgebraModule(projective) do
        if not first then Print(","); fi;
        raw := ExtRepOfObj(vector);
        Print(ExtRepOfObj(LeadingMonomial(raw![1][LeadingPosition(raw)])));
        first := false;
    od;
    Print("]");
end;

PrintPresentation := function(presentation, coefficient_printer)
    local entered;
    Print("{\"result_length\":", Length(presentation));
    Print(",\"projective_vertices\":");
    PrintProjectiveVertices(presentation[1]);
    Print(",\"entered_submodule_generators\":");
    if IsEmpty(presentation[2]) then
        Print("[]");
    else
        entered := GeneratorsOfAlgebraModule(presentation[2]);
        PrintModuleVectorList(entered, coefficient_printer);
    fi;
    Print(",\"f0\":");
    PrintModuleVectorList(presentation[3], coefficient_printer);
    Print(",\"f1\":");
    PrintModuleVectorList(presentation[4], coefficient_printer);
    Print(",\"matrix\":");
    PrintElementMatrix(presentation[5], coefficient_printer);
    Print(",\"reduced_basis\":");
    if Length(presentation) = 6 then
        PrintRawVectorList(presentation[6], coefficient_printer);
    else
        Print("[]");
    fi;
    Print("}");
end;

PrintCase := function(field, coefficient_printer)
    local quiver, path_algebra, relation, bound, signed_module,
          signed_presentation, relation_projective, relation_presentation,
          sink_projective, sink_presentation;

    quiver := Quiver(["v1", "v2", "v3"],
        [["v1", "v2", "a"], ["v2", "v3", "b"]]);
    path_algebra := PathAlgebra(field, quiver);
    relation := path_algebra.a * path_algebra.b;
    bound := path_algebra / [relation];
    signed_module := RightModuleOverPathAlgebra(
        bound, [2, 1, 0], [["a", [[1], [1]] * One(field)]]);
    signed_presentation := ProjectivePathAlgebraPresentation(signed_module);
    relation_projective := IndecProjectiveModules(bound)[1];
    relation_presentation := ProjectivePathAlgebraPresentation(relation_projective);
    sink_projective := IndecProjectiveModules(bound)[3];
    sink_presentation := ProjectivePathAlgebraPresentation(sink_projective);

    Print("{\"relation_projective\":");
    PrintPresentation(relation_presentation, coefficient_printer);
    Print(",\"signed_kernel\":");
    PrintPresentation(signed_presentation, coefficient_printer);
    Print(",\"sink_projective\":");
    PrintPresentation(sink_presentation, coefficient_printer);
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"projpathalgmodule_presentation\",\"oracle\":{");
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
