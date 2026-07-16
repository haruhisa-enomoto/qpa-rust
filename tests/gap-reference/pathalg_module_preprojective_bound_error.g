if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

quiver := Quiver(2, [[1, 2, "a"]]);
algebra := PathAlgebra(GF(2), quiver);
module := SimpleModules(algebra)[2];
PreprojectiveAlgebra(module, 0);
