if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

Q := Quiver(1, [[1, 1, "x"]]);;
A := PathAlgebra(Rationals, Q);;
B := A/[A.x];;
AssociatedGradedAlgebra(B);

QUIT_GAP(0);
