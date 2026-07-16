if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

Q := Quiver(1, [[1, 1, "x"]]);;
A := PathAlgebra(Rationals, Q);;

ProjectivePresentation(A, [[A.x], [A.x, A.x]]);

QUIT_GAP(0);
