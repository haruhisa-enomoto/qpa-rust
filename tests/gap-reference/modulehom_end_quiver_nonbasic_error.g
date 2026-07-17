if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

algebra := PathAlgebra(GF(2), Quiver(["v"], []));;
module := RightModuleOverPathAlgebra(algebra, [2], []);;
EndOfModuleAsQuiverAlgebra(module);

QUIT_GAP(0);
