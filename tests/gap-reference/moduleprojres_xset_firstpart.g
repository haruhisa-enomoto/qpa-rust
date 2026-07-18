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

PrintOSet := function(oset, coefficient_printer)
    local entry, first;
    Print("["); first := true;
    for entry in oset do
        if not first then Print(","); fi;
        Print("{\"right_tip\":", ExtRepOfObj(entry[1]));
        Print(",\"generator_position\":", entry[2]);
        Print(",\"q_prime\":"); PrintElement(entry[3], coefficient_printer);
        Print(",\"relation\":"); PrintElement(entry[4], coefficient_printer);
        Print("}");
        first := false;
    od;
    Print("]");
end;

PrintNSet := function(nset, coefficient_printer)
    local entry, first;
    Print("["); first := true;
    for entry in nset do
        if not first then Print(","); fi;
        Print("{\"right_tip\":", ExtRepOfObj(entry[1]));
        Print(",\"z\":"); PrintElement(entry[2], coefficient_printer);
        Print(",\"relation\":"); PrintElement(entry[3], coefficient_printer);
        Print("}");
        first := false;
    od;
    Print("]");
end;

PrintCase := function(field, coefficient_printer)
    local q, algebra, vertex, x, y, ideal, gb, projective, generator,
          fx, fy, mixed, prime_only, mixed_first, prime_first, fn, xset,
          vertex_xset, rightgb;
    q := Quiver(1, [[1, 1, "x"], [1, 1, "y"]]);
    algebra := PathAlgebra(field, q);
    vertex := algebra.v1;
    x := algebra.x;
    y := algebra.y;
    ideal := Ideal(algebra, [x * x, x * y, y * y]);
    gb := GroebnerBasis(ideal, [x * x, x * y, y * y]);

    projective := RightProjectiveModule(algebra, [vertex]);
    generator := GeneratorsOfAlgebraModule(projective)[1];
    fx := Vectorize(projective, [x]);
    fy := Vectorize(projective, [y]);
    mixed := Vectorize(projective, [y * x + x * y]);
    prime_only := Vectorize(projective, [y * x]);
    mixed_first := FirstPart([fx], [fy], ExtRepOfObj(mixed));
    prime_first := FirstPart([fx], [fy], ExtRepOfObj(prime_only));

    fn := Vectorize(projective, [-y]);
    xset := XSetOfPathAlgebraVector([generator], [], ideal, ExtRepOfObj(fn));
    vertex_xset := XSetOfPathAlgebraVector(
        [generator], [], ideal, ExtRepOfObj(generator));
    rightgb := RightGroebnerBasis(ideal);

    Print("{\"first_part_mixed\":");
    PrintElementList(mixed_first, coefficient_printer);
    Print(",\"first_part_prime_only\":");
    PrintElementList(prime_first, coefficient_printer);
    Print(",\"two_sided_relations\":");
    PrintElementList(GroebnerBasisOfIdeal(ideal)!.relations, coefficient_printer);
    Print(",\"right_relations\":");
    PrintElementList(rightgb!.relations, coefficient_printer);
    Print(",\"o_set\":"); PrintOSet(xset[1], coefficient_printer);
    Print(",\"n_set\":"); PrintNSet(xset[2], coefficient_printer);
    Print(",\"vertex_o_set\":"); PrintOSet(vertex_xset[1], coefficient_printer);
    Print(",\"vertex_n_set\":"); PrintNSet(vertex_xset[2], coefficient_printer);
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"moduleprojres_xset_firstpart\",\"oracle\":{");
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
