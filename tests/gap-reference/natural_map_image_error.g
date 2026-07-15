if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

Q := Quiver(1, [[1, 1, "x"]]);;
KQ := PathAlgebra(Rationals, Q, HighLevelGroebnerBasis);;
x := ElementOfPathAlgebra(KQ, Q.x);;
I := Ideal(KQ, [x^2]);;
projection := NaturalHomomorphismByIdeal(KQ, I);;

# QPA 1.37 reads obsolete internal fields from GAP's mapping record here.
ImagesRepresentative(projection, x);

QUIT_GAP(0);
