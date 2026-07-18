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

PrintCase := function(field, coefficient_printer)
    local q, algebra, module, generator, raw, source_action, cancelled,
          target_action, vectorized, canonical, right_gb, two;
    q := Quiver(["v", "w"], [["v", "w", "a"]]);
    algebra := PathAlgebra(field, q);
    module := RightProjectiveModule(algebra, [algebra.a]);
    generator := GeneratorsOfAlgebraModule(module)[1];
    raw := ExtRepOfObj(generator);
    source_action := generator^algebra.v;
    cancelled := generator + (-generator);
    two := One(field) + One(field);
    target_action := CALL_WITH_CATCH(\^, [generator, algebra.w]);
    vectorized := CALL_WITH_CATCH(Vectorize, [module, [algebra.a]]);
    canonical := CALL_WITH_CATCH(CanonicalBasis, [module]);
    right_gb := CALL_WITH_CATCH(RightGroebnerBasisOfModule, [module]);

    Print("{\"generators\":["); PrintVector(generator, coefficient_printer); Print("]");
    Print(",\"module_zero\":"); PrintVector(Zero(module), coefficient_printer);
    Print(",\"source_action\":"); PrintVector(source_action, coefficient_printer);
    Print(",\"negative\":"); PrintVector(-generator, coefficient_printer);
    Print(",\"scaled_two\":"); PrintVector(two * generator, coefficient_printer);
    Print(",\"cancelled\":"); PrintVector(cancelled, coefficient_printer);
    Print(",\"positions\":[", LeadingPosition(ExtRepOfObj(Zero(module))), ",",
        LeadingPosition(ExtRepOfObj(cancelled)), ",", LeadingPosition(raw), "]");
    Print(",\"properties\":[", IsVertexProjectiveModule(module), ",",
        IsNormalForm(raw), "]");
    Print(",\"failure_success\":[", target_action[1], ",", vectorized[1], ",",
        canonical[1], ",", right_gb[1], "]");
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"projpathalgmodule_uniform_projective\",\"oracle\":{");
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
