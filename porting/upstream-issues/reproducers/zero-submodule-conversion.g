if LoadPackage("qpa") <> true then Error("failed to load QPA"); fi;

algebra := PathAlgebra(Rationals, Quiver(2, [[1, 2, "a"]]));;
module := RightModuleOverPathAlgebra(algebra, [1, 1], [["a", [[1]]]]);;
submodule := SubAlgebraModule(module, [Zero(module)]);;
SubmoduleAsModule(submodule);
