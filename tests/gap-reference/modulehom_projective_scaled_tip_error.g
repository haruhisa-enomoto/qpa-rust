if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
algebra := PathAlgebra(Rationals, quiver);
a := 2 * (One(algebra) * ArrowsOfQuiver(quiver)[1]);
ElementIn_vA_AsElementInIndecProj(algebra, a);
