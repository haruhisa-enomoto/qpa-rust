if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

Q := Quiver(1, []);;
KQ := PathAlgebra(Rationals, Q);;

QuadraticPerpOfPathAlgebraIdeal([Zero(KQ)]);

QUIT_GAP(0);
