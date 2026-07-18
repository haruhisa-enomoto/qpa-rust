if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

quiver := Quiver(["v1", "v2", "v3"],
    [["v1", "v2", "a"], ["v2", "v3", "b"]]);
path_algebra := PathAlgebra(GF(2), quiver);
algebra := path_algebra / [path_algebra.a * path_algebra.b];
projective := IndecProjectiveModules(algebra)[1];
ProjectiveResolutionOfPathAlgebraModule(projective, 2);

QUIT_GAP(0);
