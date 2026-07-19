if LoadPackage("qpa") <> true then Error("failed to load QPA"); fi;

algebra := PathAlgebra(Rationals, Quiver(3, [[1, 2, "a"]]));;
RightModuleOverPathAlgebraNC(algebra, [[[1]]]);
