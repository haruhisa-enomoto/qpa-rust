if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
algebra := PathAlgebra(Rationals, quiver);
module := RightModuleOverPathAlgebra(algebra, [0, 0], []);
EndModuloProjOverAlgebra(module);

QUIT_GAP(0);
