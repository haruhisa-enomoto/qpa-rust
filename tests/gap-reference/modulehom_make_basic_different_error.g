if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

field := GF(2);;
first_algebra := PathAlgebra(field, Quiver(["u"], []));;
second_algebra := PathAlgebra(field, Quiver(["v1", "v2"], []));;
first := RightModuleOverPathAlgebra(first_algebra, [1], []);;
second := RightModuleOverPathAlgebra(second_algebra, [1, 0], []);;
MakeBasicListOfModules([first, second]);

QUIT_GAP(0);
