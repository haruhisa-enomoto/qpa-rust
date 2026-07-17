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
    for component in vector![1] do
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

PrintEndpoint := function(path)
    if IsZeroPath(path) then
        Print("[0]");
    else
        Print(ExtRepOfObj(path));
    fi;
end;

PrintCase := function(field, coefficient_printer)
    local q, algebra, zero, projective,
          p1, p2, p1raw, p2raw, reducible, divisor, exact, nonuniform,
          reducibleraw, divisorraw, exactraw, nonuniformraw, vertexraw, zero_checked;
    q := Quiver(["v", "w", "x"],
        [["v", "w", "a"], ["v", "w", "b"], ["w", "x", "c"]]);
    algebra := PathAlgebra(field, q);
    zero := Zero(algebra);
    projective := RightProjectiveModule(algebra, [algebra.v, algebra.v, algebra.w]);

    p1 := Vectorize(projective, [algebra.b * algebra.c, algebra.a * algebra.c, algebra.c]);
    p2 := Vectorize(projective, [algebra.a, algebra.b, algebra.w]);
    p1raw := ExtRepOfObj(p1);
    p2raw := ExtRepOfObj(p2);

    reducible := Vectorize(projective, [algebra.b * algebra.c + algebra.a, zero, zero]);
    divisor := Vectorize(projective, [algebra.b, zero, zero]);
    exact := Vectorize(projective, [algebra.b * algebra.c, zero, zero]);
    nonuniform := Vectorize(projective, [algebra.b, algebra.v, zero]);
    reducibleraw := ExtRepOfObj(reducible);
    divisorraw := ExtRepOfObj(divisor);
    exactraw := ExtRepOfObj(exact);
    nonuniformraw := ExtRepOfObj(nonuniform);
    vertexraw := ExtRepOfObj(Vectorize(projective, [algebra.v, zero, zero]));
    zero_checked := ExtRepOfObj(Vectorize(projective, [zero, zero, zero]));

    Print("{\"dimension\":", Dimension(projective));
    Print(",\"generators\":");
    PrintVectorList(List(GeneratorsOfAlgebraModule(projective), ExtRepOfObj), coefficient_printer);
    Print(",\"p1\":"); PrintVector(p1raw, coefficient_printer);
    Print(",\"p2\":"); PrintVector(p2raw, coefficient_printer);
    Print(",\"sum\":"); PrintVector(p1raw + p2raw, coefficient_printer);
    Print(",\"difference\":"); PrintVector(p1raw - p2raw, coefficient_printer);
    Print(",\"negated\":"); PrintVector(-p2raw, coefficient_printer);
    Print(",\"action_c\":"); PrintVector(p2raw^algebra.c, coefficient_printer);
    Print(",\"leading\":{\"position\":", LeadingPosition(p1raw));
    Print(",\"coefficient\":"); coefficient_printer(LeadingCoefficient(p1raw));
    Print(",\"component\":"); PrintElement(LeadingComponent(p1raw), coefficient_printer);
    Print(",\"term\":"); PrintVector(LeadingTerm(p1raw), coefficient_printer);
    Print(",\"target\":"); PrintEndpoint(TargetVertex(p1raw));
    Print(",\"source\":"); PrintEndpoint(SourceVertex(p1raw));
    Print("}");
    Print(",\"uniformity\":{\"p1\":[", IsLeftUniform(p1raw), ",",
        IsRightUniform(p1raw), ",", IsUniform(p1raw), "],");
    Print("\"divisor\":[", IsLeftUniform(divisorraw), ",",
        IsRightUniform(divisorraw), ",", IsUniform(divisorraw), "]}");
    Print(",\"ordering\":{\"p1_lt_p2\":", p1raw < p2raw,
        ",\"p2_lt_p1\":", p2raw < p1raw, "}");
    Print(",\"zero_positions\":{\"module_zero\":",
        LeadingPosition(ExtRepOfObj(Zero(projective))),
        ",\"checked_zero\":", LeadingPosition(zero_checked), "}");
    Print(",\"division\":{\"divisible\":", IsLeftDivisible(exactraw, divisorraw),
        ",\"vertex_divisible\":", IsLeftDivisible(vertexraw, vertexraw),
        ",\"factor\":"); PrintElement(LeftDivision(exactraw, divisorraw), coefficient_printer);
    Print(",\"vertex_factor\":"); PrintElement(LeftDivision(vertexraw, vertexraw), coefficient_printer);
    Print(",\"wrong_position\":", IsLeftDivisible(p1raw, p2raw), "}");
    Print(",\"tip_reduce_list\":");
    PrintVectorList(TipReduce([reducible, divisor]), coefficient_printer);
    Print(",\"tip_reduce_exact_list\":");
    PrintVectorList(TipReduce([exact, divisor]), coefficient_printer);
    Print(",\"tip_reduce_element\":");
    PrintVector(TipReduce([divisor], reducibleraw), coefficient_printer);
    Print(",\"tip_reduce_nonuniform_exact\":");
    PrintVector(TipReduce([nonuniform], nonuniformraw), coefficient_printer);
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"projpathalgmodule_vectors\",\"oracle\":{");
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
