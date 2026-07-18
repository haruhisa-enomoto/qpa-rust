if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

Q := Quiver(1, [[1,1,"x"]]);
A := PathAlgebra(Rationals, Q);
cubic := A/[A.x^3];
MaximalPathsOfGentleAlgebra(cubic);
