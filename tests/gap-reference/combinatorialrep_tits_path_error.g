if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

loop_algebra := PathAlgebra(Rationals, Quiver(3, [[1,1,"loop_arrow"]]));
TitsUnitFormOfAlgebra(loop_algebra);

QUIT_GAP(0);
