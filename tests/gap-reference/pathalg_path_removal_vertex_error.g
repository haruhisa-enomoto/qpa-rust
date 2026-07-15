if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

Q := Quiver(1, []);;
A := PathAlgebra(Rationals, Q);;
PathRemoval(A.v1, [1]);

QUIT_GAP(0);
