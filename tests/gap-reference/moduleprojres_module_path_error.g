if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
path_algebra := PathAlgebra(GF(2), quiver);
module := RightModuleOverPathAlgebra(path_algebra, [1, 0], []);
ProjectiveResolutionOfPathAlgebraModule(module, 1);

QUIT_GAP(0);
