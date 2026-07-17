if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

algebra := PathAlgebra(GF(2), Quiver(["v1", "v2"], [["v1", "v2", "a"]]));;
module := RightModuleOverPathAlgebra(algebra, [0, 0], []);;
MakeBasicListOfModules([module]);

QUIT_GAP(0);
