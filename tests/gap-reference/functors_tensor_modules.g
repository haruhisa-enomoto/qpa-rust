if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintIntegerList := function(values)
    local first, value;
    Print("["); first := true;
    for value in values do
        if not first then Print(","); fi;
        Print(value); first := false;
    od;
    Print("]");
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

PrintVectorList := function(vectors)
    local first, vector;
    Print("["); first := true;
    for vector in vectors do
        if not first then Print(","); fi;
        PrintScalarList(vector); first := false;
    od;
    Print("]");
end;

PrintSpace := function(space)
    Print("{\"dimension\":", Dimension(space), ",\"basis\":");
    PrintVectorList(BasisVectors(Basis(space)));
    Print("}");
end;

PrintCaughtVector := function(caught)
    Print("{\"success\":", caught[1], ",\"value\":");
    if caught[1] then
        PrintScalarList(caught[2]);
    else
        Print("null");
    fi;
    Print("}");
end;

PrintCaughtGrid := function(grid)
    local first_row, first, row, caught;
    Print("["); first_row := true;
    for row in grid do
        if not first_row then Print(","); fi;
        Print("["); first := true;
        for caught in row do
            if not first then Print(","); fi;
            PrintCaughtVector(caught); first := false;
        od;
        Print("]"); first_row := false;
    od;
    Print("]");
end;

PrintCase := function(field)
    local a2_quiver, a2, simples, projectives, zero_right,
          zero_tensor, zero_element, projective_right,
          projective_tensor, projective_element, path_main,
          wrong_algebras, loop_quiver, loop_parent, cubic,
          cubic_module, cubic_right, cubic_tensor, cubic_left_basis,
          cubic_right_basis, cubic_grid, bound_quiver, bound_parent,
          a, b, bound, bound_simple, bound_projective, bound_right,
          bound_tensor, bound_element;

    a2_quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
    a2 := PathAlgebra(field, a2_quiver);
    simples := SimpleModules(a2);
    projectives := IndecProjectiveModules(a2);

    zero_right := DualOfModule(simples[2]);
    zero_tensor := TensorProductOfModules(simples[1], zero_right);
    zero_element := zero_tensor[2](
        BasisVectors(Basis(simples[1]))[1],
        BasisVectors(Basis(zero_right))[1]);

    projective_right := DualOfModule(simples[1]);
    projective_tensor := TensorProductOfModules(
        projectives[1], projective_right);
    projective_element := CALL_WITH_CATCH(projective_tensor[2], [
        BasisVectors(Basis(projectives[1]))[1],
        BasisVectors(Basis(projective_right))[1]]);

    path_main := CALL_WITH_CATCH(TensorProductOfModules,
        [simples[1], DualOfModule(projectives[1])]);
    wrong_algebras := CALL_WITH_CATCH(TensorProductOfModules,
        [simples[1], simples[1]]);

    loop_quiver := Quiver(["w"], [["w", "w", "x"]]);
    loop_parent := PathAlgebra(field, loop_quiver);
    cubic := loop_parent / [loop_parent.x^3];
    cubic_module := RightModuleOverPathAlgebra(cubic, [2],
        [["x", [[0, 1], [0, 0]] * One(field)]]);
    cubic_right := DualOfModule(IndecProjectiveModules(cubic)[1]);
    cubic_tensor := TensorProductOfModules(cubic_module, cubic_right);
    cubic_left_basis := BasisVectors(Basis(cubic_module));
    cubic_right_basis := BasisVectors(Basis(cubic_right));
    cubic_grid := List(cubic_left_basis, m -> List(cubic_right_basis,
        n -> CALL_WITH_CATCH(cubic_tensor[2], [m, n])));

    bound_quiver := Quiver(["r1", "r2", "r3"],
        [["r1", "r2", "a"], ["r2", "r3", "b"]]);
    bound_parent := PathAlgebra(field, bound_quiver);
    a := bound_parent.a; b := bound_parent.b;
    bound := bound_parent / [a * b];
    bound_simple := SimpleModules(bound)[1];
    bound_projective := IndecProjectiveModules(bound)[1];
    bound_right := DualOfModule(bound_projective);
    bound_tensor := TensorProductOfModules(bound_simple, bound_right);
    bound_element := CALL_WITH_CATCH(bound_tensor[2], [
        BasisVectors(Basis(bound_simple))[1],
        BasisVectors(Basis(bound_right))[1]]);

    Print("{\"easy_zero\":{\"space\":"); PrintSpace(zero_tensor[1]);
    Print(",\"element\":"); PrintScalarList(zero_element); Print("}");
    Print(",\"projective_branch\":{\"space\":");
    PrintSpace(projective_tensor[1]);
    Print(",\"element\":"); PrintCaughtVector(projective_element); Print("}");
    Print(",\"path_main_constructor_success\":", path_main[1]);
    Print(",\"wrong_algebra_constructor_success\":", wrong_algebras[1]);
    Print(",\"cubic_quotient\":{\"space\":"); PrintSpace(cubic_tensor[1]);
    Print(",\"elements\":"); PrintCaughtGrid(cubic_grid); Print("}");
    Print(",\"bound_capture_bug\":{\"space\":"); PrintSpace(bound_tensor[1]);
    Print(",\"element\":"); PrintCaughtVector(bound_element); Print("}}");
end;

Print("{\"schema_version\":1,\"source\":{",
      "\"gap\":\"4.17dev@4efe877859ab0bc0ac041bf635a996fb897ee0c9\",",
      "\"qpa\":\"1.37@8d5f847258d4853700ce3a5d947db2a1441109f1\",",
      "\"gbnp\":\"1.1.0@e0b46c1aac16f3fbabe743e3169dc6bd86bdc654\"},",
      "\"rational\":");
PrintCase(Rationals);
Print(",\"binary\":"); PrintCase(GF(2));
Print("}\n");
