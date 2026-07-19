if LoadPackage("qpa") <> true then Error("failed to load QPA"); fi;

algebra := PathAlgebra(GF(2), Quiver(1, []));;
module := RightModuleOverPathAlgebra(algebra, [1], []);;
EndOfModuleAsQuiverAlgebra(module);
