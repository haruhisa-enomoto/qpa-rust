if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

SizeScreen([100000, 100000]);

Print("FULL_MATRIX_QQ=", IsInFullMatrixRing([[1, -2]], Rationals), "\n");
Print("FULL_MATRIX_FOREIGN=", IsInFullMatrixRing([[Z(5)]], GF(7)), "\n");

a2 := PathAlgebra(Rationals,
    Quiver(["v1", "v2"], [["v1", "v2", "a"]]));;
module := RightModuleOverPathAlgebra(a2, [2, 1],
    [["a", [[1], [2]]]]);;
basis := BasisVectors(Basis(module));;
v := 2 * basis[1] - 3 * basis[2] + 5 * basis[3];;
w := basis[1] + basis[2] - basis[3];;
element_family := FamilyObj(ExtRepOfObj(basis[1]));;

Print("BASIS=", basis, "\n");
Print("V=", v, "\n");
Print("EXTREP=", ExtRepOfObj(v), "\n");
Print("SUM=", v + w, "\n");
Print("DIFFERENCE=", v - w, "\n");
Print("LEFT_SCALE=", 2 * w, "\n");
Print("RIGHT_SCALE=", w * 2, "\n");
Print("NEGATIVE=", -v, "\n");
Print("ZERO=", Zero(v), "\n");
Print("EQUAL_SELF=", v = v, "\n");
Print("EQUAL_OTHER=", v = w, "\n");
Print("ACTION_VERTEX_1=", v^a2.v1, "\n");
Print("ACTION_VERTEX_2=", v^a2.v2, "\n");
Print("ACTION_ARROW=", v^a2.a, "\n");
Print("ACTION_SUM=", v^(a2.v1 + 2*a2.a), "\n");
Print("REBUILT=", PathModuleElem(element_family, [[2, -3], [5]]), "\n");

loop_parent := PathAlgebra(Rationals,
    Quiver(["d"], [["d", "d", "x"]]));;
dual := loop_parent / [loop_parent.x^2];;
dual_module := RightModuleOverPathAlgebra(dual, [2],
    [["x", [[0, 1], [0, 0]]]]);;
dual_basis := BasisVectors(Basis(dual_module));;
Print("QUOTIENT_ACTION=", dual_basis[1]^dual.x, "\n");

zero_vertex_module := RightModuleOverPathAlgebra(a2, [1, 0], []);;
zero_vertex_basis := BasisVectors(Basis(zero_vertex_module));;
zero_vertex_raw := ExtRepOfObj(zero_vertex_basis[1]);;
pathological := PathModuleElem(FamilyObj(zero_vertex_raw), [[0], []]);;
Print("ZERO_DIMENSION_EXTERNAL=", ExtRepOfObj(pathological), "\n");
Print("ZERO_ROW_SUM_LEFT=", pathological + zero_vertex_raw, "\n");
Print("ZERO_ROW_SUM_RIGHT=", zero_vertex_raw + pathological, "\n");
Print("ZERO_ROW_DIFFERENCE_LEFT=", pathological - zero_vertex_raw, "\n");
Print("ZERO_ROW_DIFFERENCE_RIGHT=", zero_vertex_raw - pathological, "\n");

prime_a2 := PathAlgebra(GF(2),
    Quiver(["p1", "p2"], [["p1", "p2", "b"]]));;
prime_module := RightModuleOverPathAlgebra(prime_a2, [2, 1],
    [["b", [[One(GF(2))], [One(GF(2))]]]]);;
prime_basis := BasisVectors(Basis(prime_module));;
prime_v := prime_basis[1] + prime_basis[3];;
Print("GF2_ACTION=", prime_v^prime_a2.b, "\n");

QUIT_GAP(0);
