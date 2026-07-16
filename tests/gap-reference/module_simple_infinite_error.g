if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

quiver := Quiver(["v"], [["v", "v", "x"]]);
algebra := PathAlgebra(Rationals, quiver);
SimpleModules(algebra);
