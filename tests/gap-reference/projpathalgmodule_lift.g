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
    local component, first;
    Print("["); first := true;
    for component in ExtRepOfObj(vector)![1] do
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

ProjectiveVertexPositions := function(projective)
    local generator, positions, raw;
    positions := [];
    for generator in GeneratorsOfAlgebraModule(projective) do
        raw := ExtRepOfObj(generator);
        Add(positions, VertexPosition(raw![1][LeadingPosition(raw)]));
    od;
    return positions;
end;

PrintLift := function(lifted, coefficient_printer)
    local parent, generators;
    if HasParent(lifted) then
        parent := Parent(lifted);
        generators := GeneratorsOfAlgebraModule(lifted);
        Print("{\"is_submodule\":true");
    else
        parent := lifted;
        generators := [];
        Print("{\"is_submodule\":false");
    fi;
    Print(",\"parent_vertices\":", ProjectiveVertexPositions(parent));
    Print(",\"entered_generators\":");
    PrintVectorList(generators, coefficient_printer);
    Print("}");
end;

PrintCase := function(field, coefficient_printer)
    local q, path_algebra, quotient, quotient_generators,
          source_projective, source_lift, sink_projective, sink_lift;
    q := Quiver(["v", "w", "x"],
        [["v", "w", "a"], ["w", "x", "b"]]);
    path_algebra := PathAlgebra(field, q);
    quotient := path_algebra / [path_algebra.a * path_algebra.b];
    quotient_generators := GeneratorsOfAlgebra(quotient);

    source_projective := RightProjectiveModule(
        quotient, [quotient_generators[2], quotient_generators[3]]);
    source_lift := LiftPathAlgebraModule(source_projective);
    sink_projective := RightProjectiveModule(quotient, [quotient_generators[4]]);
    sink_lift := LiftPathAlgebraModule(sink_projective);

    Print("{\"source_vertices\":[1,2],\"source_lift\":");
    PrintLift(source_lift, coefficient_printer);
    Print(",\"sink_vertices\":[3],\"sink_lift\":");
    PrintLift(sink_lift, coefficient_printer);
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"projpathalgmodule_lift\",\"oracle\":{");
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
