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

PrintVector := function(vector, coefficient_printer)
    local component, first, raw;
    if IsPathAlgebraVector(vector) then raw := vector;
    else raw := ExtRepOfObj(vector); fi;
    Print("["); first := true;
    for component in raw![1] do
        if not first then Print(","); fi;
        PrintElement(component, coefficient_printer);
        first := false;
    od;
    Print("]");
end;

PrintVectorList := function(vectors, coefficient_printer)
    local first, vector;
    Print("["); first := true;
    for vector in vectors do
        if not first then Print(","); fi;
        PrintVector(vector, coefficient_printer);
        first := false;
    od;
    Print("]");
end;

PrintCase := function(field, coefficient_printer)
    local q, algebra, zero, map, module, explicit, generator1, generator2,
          relation, probe_a, probe_c, probe_b, two;
    q := Quiver(["v", "w", "x"],
        [["v", "w", "a"], ["w", "x", "b"], ["v", "x", "c"]]);
    algebra := PathAlgebra(field, q);
    zero := Zero(algebra);
    map := [[algebra.a, -algebra.w], [algebra.c, -algebra.b]];
    module := VertexProjectivePresentation(algebra, map);
    explicit := VertexProjectivePresentation(algebra, map, [algebra.v, algebra.w]);
    generator1 := GeneratorsOfAlgebraModule(module)[1];
    generator2 := GeneratorsOfAlgebraModule(module)[2];
    relation := Vectorize(module, map[1]);
    probe_a := Vectorize(module, [algebra.a, zero]);
    probe_c := Vectorize(module, [algebra.c, zero]);
    probe_b := Vectorize(module, [zero, algebra.b]);
    two := One(field) + One(field);

    Print("{\"dimension\":", Dimension(module));
    Print(",\"properties\":[", IsFpPathAlgebraModule(module), ",",
        IsPathAlgebraModule(module), ",", IsWholeFamily(module), "]");
    Print(",\"generators\":");
    PrintVectorList(GeneratorsOfAlgebraModule(module), coefficient_printer);
    Print(",\"basis\":");
    PrintVectorList(BasisVectors(Basis(module)), coefficient_printer);
    Print(",\"zero\":"); PrintVector(Zero(module), coefficient_printer);
    Print(",\"positions\":[", LeadingPosition(ExtRepOfObj(Zero(module))), ",",
        LeadingPosition(ExtRepOfObj(relation)), ",",
        LeadingPosition(ExtRepOfObj(generator1)), ",",
        LeadingPosition(ExtRepOfObj(generator2)), "]");
    Print(",\"normal_forms\":[");
    PrintVector(relation, coefficient_printer); Print(",");
    PrintVector(probe_a, coefficient_printer); Print(",");
    PrintVector(probe_c, coefficient_printer); Print(",");
    PrintVector(probe_b, coefficient_printer); Print("]");
    Print(",\"actions\":[");
    PrintVector(generator1^algebra.a, coefficient_printer); Print(",");
    PrintVector(generator2^algebra.b, coefficient_printer); Print("]");
    Print(",\"arithmetic\":[");
    PrintVector(generator1 + generator2, coefficient_printer); Print(",");
    PrintVector(-probe_b, coefficient_printer); Print(",");
    PrintVector(two * probe_b, coefficient_printer); Print(",");
    PrintVector(probe_b + (-probe_b), coefficient_printer); Print("]");
    Print(",\"explicit\":{\"dimension\":", Dimension(explicit));
    Print(",\"basis\":");
    PrintVectorList(BasisVectors(Basis(explicit)), coefficient_printer);
    Print("}}");
end;

Print("{\"schema_version\":1,\"kind\":\"projpathalgmodule_vertex_presentation\",\"oracle\":{");
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
