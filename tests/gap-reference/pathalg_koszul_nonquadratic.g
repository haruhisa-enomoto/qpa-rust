if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

Q := Quiver(1, [[1, 1, "x"]]);;
KQ := PathAlgebra(Rationals, Q);;
x := ElementOfPathAlgebra(KQ, Q.x);;

result := QuadraticPerpOfPathAlgebraIdeal([x]);;
if result <> fail then
    Error("QPA unexpectedly accepted a nonquadratic ideal");
fi;
Print("NONQUADRATIC_FAIL=true\n");

QUIT_GAP(0);
