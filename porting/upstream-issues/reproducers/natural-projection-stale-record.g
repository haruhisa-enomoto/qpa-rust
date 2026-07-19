if LoadPackage("qpa") <> true then Error("failed to load QPA"); fi;

quiver := Quiver(1, [[1, 1, "x"]]);;
algebra := PathAlgebra(Rationals, quiver, HighLevelGroebnerBasis);;
ideal := Ideal(algebra, [algebra.x^2]);;
projection := NaturalHomomorphismByIdeal(algebra, ideal);;
ImagesRepresentative(projection, algebra.x);
