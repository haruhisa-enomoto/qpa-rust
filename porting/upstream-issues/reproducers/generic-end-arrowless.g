if LoadPackage("qpa") <> true then Error("failed to load QPA"); fi;

algebra := PathAlgebra(Rationals, Quiver(1, []));;
module := RightModuleOverPathAlgebra(algebra, [2], []);;
End(algebra, module);
