if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

loop_algebra := PathAlgebra(Rationals, Quiver(1, [[1,1,"x"]]));
non_gentle := loop_algebra/[loop_algebra.x^3];
gentle := PathAlgebra(Rationals, Quiver(2, [[1,2,"a"]]));
AreDerivedEquivalent(non_gentle, gentle);

QUIT_GAP(0);
