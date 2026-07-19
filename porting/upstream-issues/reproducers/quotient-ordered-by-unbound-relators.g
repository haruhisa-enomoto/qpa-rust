if LoadPackage("qpa") <> true then Error("failed to load QPA"); fi;

quiver := Quiver(1, [[1, 1, "a"]]);;
algebra := PathAlgebra(Rationals, quiver);;
quotient := algebra / [algebra.a^2];;
OrderedBy(quotient, OrderingOfAlgebra(quotient));
