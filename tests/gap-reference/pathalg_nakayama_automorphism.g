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

PrintScalarMatrix := function(matrix)
    local first, row;
    Print("["); first := true;
    for row in matrix do
        if not first then Print(","); fi;
        PrintScalarList(row); first := false;
    od;
    Print("]");
end;

PrintOrder := function(value)
    if value = false then
        Print("false");
    elif value = infinity then
        PrintJsonString("infinity");
    else
        Print(value);
    fi;
end;

PrintPathCase := function(algebra)
    local auto_initially_known, order_initially_known, auto, basis, order;
    auto_initially_known := HasNakayamaAutomorphism(algebra);
    order_initially_known := HasOrderOfNakayamaAutomorphism(algebra);
    auto := NakayamaAutomorphism(algebra);
    order := OrderOfNakayamaAutomorphism(algebra);
    Print("{\"auto_initially_known\":", auto_initially_known);
    Print(",\"order_initially_known\":", order_initially_known);
    if auto = false then
        Print(",\"automorphism\":false");
    else
        basis := CanonicalBasis(algebra);
        Print(",\"automorphism\":");
        PrintScalarMatrix(List(basis,
            beta -> Coefficients(basis, auto(beta))));
    fi;
    Print(",\"order\":"); PrintOrder(order);
    Print(",\"auto_known_after\":", HasNakayamaAutomorphism(algebra));
    Print(",\"order_known_after\":",
        HasOrderOfNakayamaAutomorphism(algebra));
    Print("}");
end;

PrintQuotientCase := function(algebra, trigger)
    local auto_initially_known, form_initially_known,
          functional_initially_known, order_initially_known,
          auto, form, functional, basis, order,
          auto_known_after_trigger, form_known_after_trigger,
          functional_known_after_trigger;
    auto_initially_known := HasNakayamaAutomorphism(algebra);
    form_initially_known := HasFrobeniusForm(algebra);
    functional_initially_known := HasFrobeniusLinearFunctional(algebra);
    order_initially_known := HasOrderOfNakayamaAutomorphism(algebra);

    if trigger = "form" then
        form := FrobeniusForm(algebra);
        auto := NakayamaAutomorphism(algebra);
    elif trigger = "functional" then
        functional := FrobeniusLinearFunctional(algebra);
        auto := NakayamaAutomorphism(algebra);
    else
        auto := NakayamaAutomorphism(algebra);
    fi;
    auto_known_after_trigger := HasNakayamaAutomorphism(algebra);
    form_known_after_trigger := HasFrobeniusForm(algebra);
    functional_known_after_trigger := HasFrobeniusLinearFunctional(algebra);

    if not IsBound(form) then form := FrobeniusForm(algebra); fi;
    if not IsBound(functional) then
        functional := FrobeniusLinearFunctional(algebra);
    fi;

    Print("{\"auto_initially_known\":", auto_initially_known);
    Print(",\"form_initially_known\":", form_initially_known);
    Print(",\"functional_initially_known\":",
        functional_initially_known);
    Print(",\"order_initially_known\":", order_initially_known);
    Print(",\"auto_known_after_trigger\":", auto_known_after_trigger);
    Print(",\"form_known_after_trigger\":", form_known_after_trigger);
    Print(",\"functional_known_after_trigger\":",
        functional_known_after_trigger);
    if auto = false then
        Print(",\"automorphism\":false");
        Print(",\"form\":false,\"functional\":false");
    else
        basis := CanonicalBasis(algebra);
        Print(",\"automorphism\":");
        PrintScalarMatrix(List(basis,
            beta -> Coefficients(basis, auto(beta))));
        Print(",\"form\":");
        PrintScalarMatrix(List(basis,
            left -> List(basis, right -> form(left, right))));
        Print(",\"functional\":");
        PrintScalarList(List(basis, functional));
        order := OrderOfNakayamaAutomorphism(algebra);
        Print(",\"order\":"); PrintOrder(order);
    fi;
    Print(",\"auto_known_after\":", HasNakayamaAutomorphism(algebra));
    Print(",\"form_known_after\":", HasFrobeniusForm(algebra));
    Print(",\"functional_known_after\":",
        HasFrobeniusLinearFunctional(algebra));
    Print(",\"order_known_after\":",
        HasOrderOfNakayamaAutomorphism(algebra));
    Print("}");
end;

PrintCases := function(field)
    local discrete, a2, loop_parent, dual_numbers, a3_parent,
          nonselfinjective, cycle_parent, cycle, quantum_parent,
          quantum, scalar;
    discrete := PathAlgebra(field, Quiver(["d1", "d2"], []));
    a2 := PathAlgebra(field,
        Quiver(["a1", "a2"], [["a1", "a2", "a"]]));

    loop_parent := PathAlgebra(field,
        Quiver(["l"], [["l", "l", "t"]]));
    dual_numbers := loop_parent / [loop_parent.t * loop_parent.t];

    a3_parent := PathAlgebra(field,
        Quiver(["q1", "q2", "q3"],
            [["q1", "q2", "u"], ["q2", "q3", "v"]]));
    nonselfinjective := a3_parent /
        [a3_parent.u * a3_parent.v];

    cycle_parent := PathAlgebra(field,
        Quiver(["c1", "c2", "c3"],
            [["c1", "c2", "x"], ["c2", "c3", "y"],
             ["c3", "c1", "z"]]));
    cycle := cycle_parent /
        [cycle_parent.x * cycle_parent.y,
         cycle_parent.y * cycle_parent.z,
         cycle_parent.z * cycle_parent.x];

    quantum_parent := PathAlgebra(field,
        Quiver(["w"], [["w", "w", "x"], ["w", "w", "y"]]));
    if Characteristic(field) = 2 then
        scalar := One(field);
    else
        scalar := 2 * One(field);
    fi;
    quantum := quantum_parent /
        [quantum_parent.x * quantum_parent.x,
         quantum_parent.y * quantum_parent.y,
         quantum_parent.x * quantum_parent.y
            - scalar * quantum_parent.y * quantum_parent.x];

    Print("{\"path_discrete\":"); PrintPathCase(discrete);
    Print(",\"path_a2\":"); PrintPathCase(a2);
    Print(",\"quotient_dual_numbers\":");
    PrintQuotientCase(dual_numbers, "form");
    Print(",\"quotient_nonselfinjective\":");
    PrintQuotientCase(nonselfinjective, "automorphism");
    Print(",\"quotient_radical_square_zero_cycle\":");
    PrintQuotientCase(cycle, "automorphism");
    Print(",\"quotient_quantum_complete_intersection\":");
    PrintQuotientCase(quantum, "functional");
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"pathalg_nakayama_automorphism\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"f2\":"); PrintCases(GF(2));
Print(",\"f5\":"); PrintCases(GF(5));
Print(",\"rational\":"); PrintCases(Rationals);
Print("}\n");

QUIT_GAP(0);
