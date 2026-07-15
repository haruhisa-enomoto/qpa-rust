if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
algebra := PathAlgebra(Rationals, quiver);
v1 := One(algebra) * VerticesOfQuiver(quiver)[1];
FromMatrixToHomomorphismOfProjectives(algebra, [[v1]], [2], [1]);
