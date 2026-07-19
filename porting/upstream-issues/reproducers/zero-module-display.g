if LoadPackage("qpa") <> true then Error("failed to load QPA"); fi;

algebra := PathAlgebra(Rationals, Quiver(2, [[1, 2, "a"]]));;
Display(ZeroModule(algebra));
