if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

loop_algebra := PathAlgebra(Rationals, Quiver(1, [[1,1,"x"]]));
dual_numbers := loop_algebra/[loop_algebra.x^2];
TitsUnitFormOfAlgebra(dual_numbers);

QUIT_GAP(0);
