if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

algebra := PathAlgebra(GF(2), Quiver(["x"], [["x", "x", "a"]]));;
module := RightModuleOverPathAlgebra(algebra, [1], []);;
ADRAlgebraOfBasicModule([module]);

