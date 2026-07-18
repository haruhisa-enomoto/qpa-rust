if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

a2 := PathAlgebra(Rationals,
    Quiver(["v1", "v2"], [["v1", "v2", "a"]]));;
Display(ZeroModule(a2));

QUIT_GAP(0);
