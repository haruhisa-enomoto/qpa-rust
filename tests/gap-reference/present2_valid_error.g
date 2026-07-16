if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

Q := Quiver(3, [[1, 2, "a"], [1, 3, "b"], [2, 3, "c"]]);;
A := PathAlgebra(Rationals, Q);;

ProjectivePresentation(A, [[A.a+A.b, Zero(A)], [A.c, A.v2]]);

QUIT_GAP(0);
