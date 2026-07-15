if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

Q := Quiver(1, [[1, 1, "x"]]);;
A := PathAlgebra(Rationals, Q);;
B := A/[A.x];;
if IsAdmissibleQuotientOfPathAlgebra(B) <> false then
    Error("expected a cached false admissibility property");
fi;
AssociatedGradedAlgebra(B);

QUIT_GAP(0);
