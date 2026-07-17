if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

first_algebra := PathAlgebra(GF(2), Quiver(["x"], []));;
second_algebra := PathAlgebra(GF(2), Quiver(["y"], []));;
first := RightModuleOverPathAlgebra(first_algebra, [1], []);;
second := RightModuleOverPathAlgebra(second_algebra, [1], []);;
EndOfBasicModuleAsQuiverAlgebra([first, second]);

