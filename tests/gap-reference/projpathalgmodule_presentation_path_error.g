if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
path_algebra := PathAlgebra(GF(2), quiver);
module := RightModuleOverPathAlgebra(path_algebra, [1, 0], []);
ProjectivePathAlgebraPresentation(module);

QUIT_GAP(0);
