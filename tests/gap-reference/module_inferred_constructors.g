if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

SizeScreen([100000, 100000]);

PrintModuleData := function(name, module)
    Print(name, "=", DimensionVector(module), "|",
        MatricesOfPathAlgebraModule(module), "\n");
end;

chain := PathAlgebra(Rationals,
    Quiver(["v1", "v2", "v3"],
        [["v1", "v2", "a"], ["v2", "v3", "b"]]));;
a_matrix := [[1, 2]];;
b_matrix := [[3], [4]];;
arrows := ArrowsOfQuiver(QuiverOfPathAlgebra(chain));;
named := RightModuleOverPathAlgebra(chain,
    [["b", b_matrix], ["a", a_matrix]]);;
positional := RightModuleOverPathAlgebra(chain, [a_matrix, b_matrix]);;
indexed := RightModuleOverPathAlgebra(chain,
    [[arrows[2], b_matrix], [arrows[1], a_matrix]]);;
unchecked := RightModuleOverPathAlgebraNC(chain, [a_matrix, b_matrix]);;

PrintModuleData("NAMED", named);
PrintModuleData("POSITIONAL", positional);
PrintModuleData("INDEXED", indexed);
PrintModuleData("UNCHECKED", unchecked);

zero_target := PathAlgebra(Rationals,
    Quiver(["z1", "z2"], [["z1", "z2", "z"]]));;
zero_source_module := RightModuleOverPathAlgebra(zero_target,
    [["z", [0, 2]]]);;
zero_target_module := RightModuleOverPathAlgebra(zero_target,
    [["z", [2, 0]]]);;
zero_module := RightModuleOverPathAlgebra(zero_target,
    [["z", [0, 0]]]);;
PrintModuleData("ZERO_SOURCE", zero_source_module);
PrintModuleData("ZERO_TARGET", zero_target_module);
PrintModuleData("ZERO_MODULE", zero_module);

with_isolated := PathAlgebra(Rationals,
    Quiver(["i1", "i2", "i3"], [["i1", "i2", "c"]]));;
isolated := RightModuleOverPathAlgebra(with_isolated, [[[5]]]);;
PrintModuleData("ISOLATED_CHECKED", isolated);

prime := PathAlgebra(GF(2),
    Quiver(["p1", "p2"], [["p1", "p2", "d"]]));;
prime_module := RightModuleOverPathAlgebra(prime,
    [[[One(GF(2))]]]);;
PrintModuleData("GF2", prime_module);

quotient_parent := PathAlgebra(Rationals,
    Quiver(["q1", "q2", "q3"],
        [["q1", "q2", "x"], ["q2", "q3", "y"]]));;
quotient := quotient_parent / [quotient_parent.x * quotient_parent.y];;
relation_ok := RightModuleOverPathAlgebra(quotient, [[[1]], [[0]]]);;
relation_fail := RightModuleOverPathAlgebra(quotient, [[[1]], [[1]]]);;
PrintModuleData("RELATION_OK", relation_ok);
Print("RELATION_FAIL=", relation_fail = fail, "\n");

QUIT_GAP(0);
