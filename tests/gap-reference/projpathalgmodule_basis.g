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

PrintCoefficientList := function(coefficients, coefficient_printer)
    local coefficient, first;
    if coefficients = fail then
        Print("null");
        return;
    fi;
    Print("["); first := true;
    for coefficient in coefficients do
        if not first then Print(","); fi;
        coefficient_printer(coefficient);
        first := false;
    od;
    Print("]");
end;

PrintCase := function(field, coefficient_printer)
    local q, algebra, zero, one, projective, u, v, w, raw,
          space, entered, new_basis, unchecked, checked, domain_basis,
          dependent, outside, probe, nc, g1, g2, submodule, sub_basis,
          based_submodule, based_basis;
    q := Quiver(["v", "w", "x"],
        [["v", "w", "a"], ["w", "x", "b"]]);
    algebra := PathAlgebra(field, q);
    zero := Zero(algebra);
    one := One(field);
    projective := RightProjectiveModule(algebra, [algebra.v, algebra.v]);

    u := ExtRepOfObj(Vectorize(projective, [algebra.a, zero]));
    v := ExtRepOfObj(Vectorize(projective, [2 * one * algebra.a + algebra.v, zero]));
    w := ExtRepOfObj(Vectorize(projective, [zero, algebra.a * algebra.b]));
    raw := ExtRepOfObj(Vectorize(projective, [zero, algebra.v]));
    entered := [Zero(raw), v, u, w];
    space := LeftModuleByGenerators(field, entered);

    new_basis := NewBasis(space, entered);
    unchecked := BasisNC(space, entered);
    checked := Basis(space, entered);
    domain_basis := BasisOfDomain(space);
    dependent := Basis(space, [u, u]);
    outside := LeftModuleByGenerators(field, [u]);
    probe := v + w;
    nc := PathAlgebraVectorNC(FamilyObj(raw), raw![1], 2, false);

    g1 := Vectorize(projective, [algebra.a + algebra.a * algebra.b, algebra.a]);
    g2 := Vectorize(projective, [algebra.a * algebra.b, zero]);
    submodule := SubAlgebraModule(projective, [g1, g2]);
    sub_basis := Basis(submodule);
    based_submodule := SubAlgebraModule(projective, [g1, g2], "basis");
    based_basis := Basis(based_submodule);

    Print("{\"new_basis\":");
    PrintVectorList(BasisVectors(new_basis), coefficient_printer);
    Print(",\"unchecked_basis\":");
    if unchecked = fail then Print("null");
    else PrintVectorList(BasisVectors(unchecked), coefficient_printer); fi;
    Print(",\"checked_basis\":");
    if checked = fail then Print("null");
    else PrintVectorList(BasisVectors(checked), coefficient_printer); fi;
    Print(",\"domain_basis\":");
    PrintVectorList(BasisVectors(domain_basis), coefficient_printer);
    Print(",\"probe_coefficients\":");
    PrintCoefficientList(Coefficients(domain_basis, probe), coefficient_printer);
    Print(",\"outside_coefficients\":");
    PrintCoefficientList(Coefficients(domain_basis, raw), coefficient_printer);
    Print(",\"dependent_failed\":", dependent = fail);
    Print(",\"membership\":[", u in space, ",", raw in space, ",",
        Zero(raw) in space, ",", raw in outside, "]");
    Print(",\"nc\":{\"vector\":"); PrintVector(nc, coefficient_printer);
    Print(",\"leading_position\":", LeadingPosition(nc),
        ",\"is_normal_form\":", IsNormalForm(nc), "}");
    Print(",\"submodule_basis\":");
    PrintVectorList(List(BasisVectors(sub_basis), ExtRepOfObj), coefficient_printer);
    Print(",\"submodule_dimension\":", Dimension(submodule));
    Print(",\"based_basis\":");
    PrintVectorList(List(BasisVectors(based_basis), ExtRepOfObj), coefficient_printer);
    Print(",\"based_dimension\":", Dimension(based_submodule));
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"projpathalgmodule_basis\",\"oracle\":{");
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
