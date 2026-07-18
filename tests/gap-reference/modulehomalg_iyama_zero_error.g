if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

algebra := PathAlgebra(GF(2), Quiver(["v"], []));;
zero := RightModuleOverPathAlgebra(algebra, [0], []);;
IyamaGenerator(zero);

QUIT_GAP(0);
