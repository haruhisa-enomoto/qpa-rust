if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

Q := Quiver(1, [[1,1,"x"]]);
A := PathAlgebra(Rationals, Q);
square := A/[A.x^2];
CombinatorialMapOfGentleAlgebra(square);
