if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

field := GF(2);;
algebra := PathAlgebra(field, Quiver(["v"], [["v", "v", "x"]]));;
module := RightModuleOverPathAlgebra(algebra, [1], [["x", [[0]] * One(field)]]);;
MakeBasicListOfModules([module]);

QUIT_GAP(0);
