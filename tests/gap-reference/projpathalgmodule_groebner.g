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

PrintModuleVector := function(vector, coefficient_printer)
    PrintVector(ExtRepOfObj(vector), coefficient_printer);
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

PrintCase := function(field, coefficient_printer)
    local q, algebra, zero, projective, g1, g2, g3, submodule, zero_submodule,
          gb, probe1, probe2, probe3;
    q := Quiver(["v", "w", "x"],
        [["v", "w", "a"], ["v", "w", "b"], ["w", "x", "c"]]);
    algebra := PathAlgebra(field, q);
    zero := Zero(algebra);
    projective := RightProjectiveModule(algebra, [algebra.v, algebra.v, algebra.w]);

    g1 := Vectorize(projective, [algebra.a + algebra.a * algebra.c, algebra.b, zero]);
    g2 := Vectorize(projective, [algebra.b * algebra.c + algebra.a * algebra.c, zero, algebra.c]);
    g3 := Vectorize(projective, [algebra.a + algebra.a * algebra.c, zero, zero]);
    submodule := SubAlgebraModule(projective, [g1, g2, g3]);
    zero_submodule := SubAlgebraModule(projective, [Zero(projective)]);
    gb := RightGroebnerBasisOfModule(submodule);

    probe1 := Vectorize(projective, [algebra.a * algebra.c + algebra.a, zero, zero]);
    probe2 := Vectorize(projective,
        [algebra.b * algebra.c + algebra.a * algebra.c + algebra.a, zero, algebra.c]);
    probe3 := Vectorize(projective, [algebra.b * algebra.c, algebra.a * algebra.c, zero]);

    Print("{\"projective_dimension\":", Dimension(projective));
    Print(",\"projective_basis\":");
    PrintModuleVectorList(BasisVectors(CanonicalBasis(projective)), coefficient_printer);
    Print(",\"whole_uniform_generators\":");
    PrintModuleVectorList(UniformGeneratorsOfModule(projective), coefficient_printer);
    Print(",\"entered_generators\":");
    PrintModuleVectorList(GeneratorsOfAlgebraModule(submodule), coefficient_printer);
    Print(",\"uniform_generators\":");
    PrintModuleVectorList(UniformGeneratorsOfModule(submodule), coefficient_printer);
    Print(",\"zero_uniform_generators\":");
    PrintModuleVectorList(UniformGeneratorsOfModule(zero_submodule), coefficient_printer);
    Print(",\"groebner_basis\":");
    PrintModuleVectorList(BasisVectors(gb), coefficient_printer);
    Print(",\"reduced_probes\":[");
    PrintModuleVector(CompletelyReduce(gb, probe1), coefficient_printer); Print(",");
    PrintModuleVector(CompletelyReduce(gb, probe2), coefficient_printer); Print(",");
    PrintModuleVector(CompletelyReduce(gb, probe3), coefficient_printer); Print("]");
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"projpathalgmodule_groebner\",\"oracle\":{");
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
