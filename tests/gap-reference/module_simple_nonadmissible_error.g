if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

quiver := Quiver(["v"], [["v", "v", "x"]]);
parent := PathAlgebra(Rationals, quiver);
algebra := parent / [parent.x];
SimpleModules(algebra);
