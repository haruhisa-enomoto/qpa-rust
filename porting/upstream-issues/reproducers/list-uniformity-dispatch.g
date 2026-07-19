if LoadPackage("qpa") <> true then Error("failed to load QPA"); fi;

quiver := Quiver(2, [[1, 2, "a"]]);;
algebra := PathAlgebra(Rationals, quiver);;
IsLeftUniform([algebra.a], quiver.v1);
