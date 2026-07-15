if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

quiver := Quiver(["v"], []);
algebra := PathAlgebra(Rationals, quiver);
ElementIn_vA_AsElementInIndecProj(algebra, Zero(algebra));
