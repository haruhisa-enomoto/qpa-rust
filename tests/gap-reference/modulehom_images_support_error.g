if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
algebra := PathAlgebra(Rationals, quiver);
module := RightModuleOverPathAlgebra(algebra, [2, 1], [["a", [[1], [0]]]]);
basis := BasisVectors(Basis(module));
images := ShallowCopy(basis);
images[1] := basis[3];
HomomorphismFromImages(module, module, images);
