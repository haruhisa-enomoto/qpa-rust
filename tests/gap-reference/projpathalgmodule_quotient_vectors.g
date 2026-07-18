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
    data := CoefficientsAndMagmaElements(element![1]);
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
    local q, path_algebra, quotient, generators, projective, zero,
          vector, raw, other, other_raw, scalar, nc, canonical;
    q := Quiver(["v", "w", "x"],
        [["v", "w", "a"], ["w", "x", "b"]]);
    path_algebra := PathAlgebra(field, q);
    quotient := path_algebra / [path_algebra.a * path_algebra.b];
    generators := GeneratorsOfAlgebra(quotient);
    projective := RightProjectiveModule(quotient, [generators[2], generators[3]]);
    zero := Zero(quotient);
    vector := Vectorize(projective, [generators[5], generators[3]]);
    raw := ExtRepOfObj(vector);
    other := Vectorize(projective, [generators[2], zero]);
    other_raw := ExtRepOfObj(other);
    scalar := One(field) + One(field);
    nc := PathAlgebraVectorNC(FamilyObj(raw), raw![1], 2, false);
    canonical := CanonicalBasis(projective);

    Print("{\"dimension\":", Dimension(projective));
    Print(",\"generators\":");
    PrintVectorList(GeneratorsOfAlgebraModule(projective), coefficient_printer);
    Print(",\"canonical_basis\":");
    PrintVectorList(BasisVectors(canonical), coefficient_printer);
    Print(",\"vector\":"); PrintVector(vector, coefficient_printer);
    Print(",\"other\":"); PrintVector(other, coefficient_printer);
    Print(",\"sum\":"); PrintVector(vector + other, coefficient_printer);
    Print(",\"difference\":"); PrintVector(raw - other_raw, coefficient_printer);
    Print(",\"negated\":"); PrintVector(-vector, coefficient_printer);
    Print(",\"scaled_left\":"); PrintVector(scalar * raw, coefficient_printer);
    Print(",\"scaled_right\":"); PrintVector(raw * scalar, coefficient_printer);
    Print(",\"action_b\":"); PrintVector(vector^generators[6], coefficient_printer);
    Print(",\"leading_term\":"); PrintVector(LeadingTerm(raw), coefficient_printer);
    Print(",\"leading_component\":");
    PrintElement(LeadingComponent(raw), coefficient_printer);
    Print(",\"leading_position\":", LeadingPosition(raw));
    Print(",\"leading_coefficient\":"); coefficient_printer(LeadingCoefficient(raw));
    Print(",\"target_vertex\":", ExtRepOfObj(TargetVertex(raw)));
    Print(",\"source_vertex\":", ExtRepOfObj(SourceVertex(raw)));
    Print(",\"left_divisibility\":[", IsLeftDivisible(raw, other_raw), ",",
        IsLeftDivisible(other_raw, raw), "]");
    Print(",\"uniformity\":[", IsLeftUniform(raw), ",", IsRightUniform(raw), ",",
        IsUniform(raw), "]");
    Print(",\"normal_form\":", IsNormalForm(raw));
    Print(",\"zero_positions\":[", LeadingPosition(ExtRepOfObj(Zero(projective))), ",",
        LeadingPosition(ExtRepOfObj(Vectorize(projective, [zero, zero]))), "]");
    Print(",\"nc\":{\"vector\":"); PrintVector(nc, coefficient_printer);
    Print(",\"leading_position\":", LeadingPosition(nc), "}");
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"projpathalgmodule_quotient_vectors\",\"oracle\":{");
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
