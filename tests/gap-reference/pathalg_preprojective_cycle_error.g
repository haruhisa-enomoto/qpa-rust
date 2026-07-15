if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

cyclic := PathAlgebra(Rationals,
    Quiver(2, [[1, 2, "a"], [2, 1, "b"]]));;
PreprojectiveAlgebra(cyclic);

QUIT_GAP(0);
