if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

SizeScreen([100000, 100000]);

DataOf := function(module)
    return [DimensionVector(module), MatricesOfPathAlgebraModule(module)];
end;

RunCase := function(field)
    local algebra, regular, converted, custom_space, custom_regular,
          custom_converted, zero_space, zero_regular, zero_converted,
          matrix_module, reconverted,
          loop_parent, quotient, quotient_regular, quotient_converted,
          arrowless, arrowless_regular, arrowless_converted;
    algebra := PathAlgebra(field,
        Quiver(["v1", "v2"], [["v1", "v2", "a"]]));
    regular := RightAlgebraModule(algebra, \*, algebra);
    converted := RightAlgebraModuleToPathAlgebraMatModule(regular);
    custom_space := VectorSpace(field,
        [algebra.v1 + algebra.a, algebra.v2, algebra.a], "basis");
    custom_regular := RightAlgebraModule(algebra, \*, custom_space);
    custom_converted :=
        RightAlgebraModuleToPathAlgebraMatModule(custom_regular);
    zero_space := Subspace(algebra, []);
    zero_regular := RightAlgebraModule(algebra, \*, zero_space);
    zero_converted :=
        RightAlgebraModuleToPathAlgebraMatModule(zero_regular);

    matrix_module := RightModuleOverPathAlgebra(algebra, [2, 2],
        [["a", [[One(field), One(field)],
                 [Zero(field), One(field)]]]]);
    reconverted := RightAlgebraModuleToPathAlgebraMatModule(matrix_module);

    loop_parent := PathAlgebra(field,
        Quiver(["w"], [["w", "w", "t"]]));
    quotient := loop_parent / [loop_parent.t^2];
    quotient_regular := RightAlgebraModule(quotient, \*, quotient);
    quotient_converted :=
        RightAlgebraModuleToPathAlgebraMatModule(quotient_regular);

    arrowless := PathAlgebra(field, Quiver(["s"], []));
    arrowless_regular := RightAlgebraModule(arrowless, \*, arrowless);
    arrowless_converted :=
        RightAlgebraModuleToPathAlgebraMatModule(arrowless_regular);

    return [
        DataOf(converted),
        DataOf(custom_converted),
        DataOf(zero_converted),
        DataOf(reconverted),
        DataOf(quotient_converted),
        DataOf(arrowless_converted)
    ];
end;

Print("QQ=", RunCase(Rationals), "\n");
Print("GF2=", RunCase(GF(2)), "\n");

QUIT_GAP(0);
