if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintScalarList := function(values)
    local first, value;
    Print("["); first := true;
    for value in values do
        if not first then Print(","); fi;
        PrintJsonString(String(value)); first := false;
    od;
    Print("]");
end;

PrintAnnihilator := function(algebra, elements)
    local basis, first, element;
    basis := Basis(algebra);
    Print("["); first := true;
    for element in elements do
        if not first then Print(","); fi;
        PrintScalarList(Coefficients(basis, element)); first := false;
    od;
    Print("]");
end;

PrintCase := function(field)
    local quiver, algebra, projectives, faithful, simples, zero_module,
          arrow_zero, zero_annihilator, simple1_annihilator,
          simple2_annihilator, arrow_zero_annihilator,
          faithful_annihilator, loop_quiver, parent, x, quotient,
          quotient_projective, quotient_simple, quotient_zero,
          quotient_two_zero, quotient_projective_annihilator,
          quotient_simple_annihilator, quotient_zero_annihilator,
          quotient_two_zero_annihilator, bound_quiver, bound_parent,
          bound_a, bound_b, bound_algebra, bound_zero, bound_simple,
          bound_zero_annihilator, bound_simple_annihilator;

    quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
    algebra := PathAlgebra(field, quiver);
    projectives := IndecProjectiveModules(algebra);
    faithful := DirectSumOfQPAModules(projectives);
    simples := SimpleModules(algebra);
    zero_module := ZeroModule(algebra);
    arrow_zero := RightModuleOverPathAlgebra(algebra, [1, 1], []);
    zero_annihilator := AnnihilatorOfModule(zero_module);
    simple1_annihilator := AnnihilatorOfModule(simples[1]);
    simple2_annihilator := AnnihilatorOfModule(simples[2]);
    arrow_zero_annihilator := AnnihilatorOfModule(arrow_zero);
    faithful_annihilator := AnnihilatorOfModule(faithful);

    loop_quiver := Quiver(["w"], [["w", "w", "x"]]);
    parent := PathAlgebra(field, loop_quiver);
    x := parent.x;
    quotient := parent / [x * x];
    quotient_projective := IndecProjectiveModules(quotient)[1];
    quotient_simple := SimpleModules(quotient)[1];
    quotient_zero := ZeroModule(quotient);
    quotient_two_zero := RightModuleOverPathAlgebra(quotient, [2], []);
    quotient_projective_annihilator :=
        AnnihilatorOfModule(quotient_projective);
    quotient_simple_annihilator := AnnihilatorOfModule(quotient_simple);
    quotient_zero_annihilator := AnnihilatorOfModule(quotient_zero);
    quotient_two_zero_annihilator := AnnihilatorOfModule(quotient_two_zero);

    bound_quiver := Quiver(["u1", "u2", "u3"],
        [["u1", "u2", "c"], ["u2", "u3", "d"]]);
    bound_parent := PathAlgebra(field, bound_quiver);
    bound_a := bound_parent.c;
    bound_b := bound_parent.d;
    bound_algebra := bound_parent / [bound_a * bound_b];
    bound_zero := ZeroModule(bound_algebra);
    bound_simple := SimpleModules(bound_algebra)[3];
    bound_zero_annihilator := AnnihilatorOfModule(bound_zero);
    bound_simple_annihilator := AnnihilatorOfModule(bound_simple);

    Print("{\"path\":{");
    Print("\"basis_dimension\":", Dimension(algebra));
    Print(",\"zero\":"); PrintAnnihilator(algebra, zero_annihilator);
    Print(",\"simple1\":"); PrintAnnihilator(algebra, simple1_annihilator);
    Print(",\"simple2\":"); PrintAnnihilator(algebra, simple2_annihilator);
    Print(",\"arrow_zero\":");
    PrintAnnihilator(algebra, arrow_zero_annihilator);
    Print(",\"faithful\":"); PrintAnnihilator(algebra, faithful_annihilator);
    Print(",\"cache_bound\":[",
        HasAnnihilatorOfModule(zero_module), ",",
        HasAnnihilatorOfModule(simples[1]), ",",
        HasAnnihilatorOfModule(simples[2]), ",",
        HasAnnihilatorOfModule(arrow_zero), ",",
        HasAnnihilatorOfModule(faithful), "]");
    Print("},\"quotient\":{");
    Print("\"basis_dimension\":", Dimension(quotient));
    Print(",\"projective\":");
    PrintAnnihilator(quotient, quotient_projective_annihilator);
    Print(",\"simple\":");
    PrintAnnihilator(quotient, quotient_simple_annihilator);
    Print(",\"zero\":");
    PrintAnnihilator(quotient, quotient_zero_annihilator);
    Print(",\"two_zero\":");
    PrintAnnihilator(quotient, quotient_two_zero_annihilator);
    Print(",\"cache_bound\":[",
        HasAnnihilatorOfModule(quotient_projective), ",",
        HasAnnihilatorOfModule(quotient_simple), ",",
        HasAnnihilatorOfModule(quotient_zero), ",",
        HasAnnihilatorOfModule(quotient_two_zero), "]},");
    Print("\"bound_quotient\":{");
    Print("\"basis_dimension\":", Dimension(bound_algebra));
    Print(",\"zero\":");
    PrintAnnihilator(bound_algebra, bound_zero_annihilator);
    Print(",\"simple3\":");
    PrintAnnihilator(bound_algebra, bound_simple_annihilator);
    Print(",\"cache_bound\":[",
        HasAnnihilatorOfModule(bound_zero), ",",
        HasAnnihilatorOfModule(bound_simple), "]}}");
end;

Print("{\"schema_version\":1,\"source\":{",
      "\"gap\":\"4.17dev@4efe877859ab0bc0ac041bf635a996fb897ee0c9\",",
      "\"qpa\":\"1.37@8d5f847258d4853700ce3a5d947db2a1441109f1\",",
      "\"gbnp\":\"1.1.0@e0b46c1aac16f3fbabe743e3169dc6bd86bdc654\"},",
      "\"rational\":");
PrintCase(Rationals);
Print(",\"binary\":");
PrintCase(GF(2));
Print("}\n");
