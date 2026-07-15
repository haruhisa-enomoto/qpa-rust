if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

Q := Quiver(1, []);;
A := PathAlgebra(Rationals, Q);;
TipMonomialandCoefficientOfVector(A, []);

QUIT_GAP(0);
