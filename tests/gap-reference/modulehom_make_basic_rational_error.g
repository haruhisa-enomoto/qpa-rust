if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

algebra := PathAlgebra(Rationals, Quiver(["v"], []));;
module := RightModuleOverPathAlgebra(algebra, [1], []);;
MakeBasicListOfModules([module]);

QUIT_GAP(0);
