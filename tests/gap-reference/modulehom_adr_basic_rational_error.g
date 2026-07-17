if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

algebra := PathAlgebra(Rationals, Quiver(["x"], []));;
module := RightModuleOverPathAlgebra(algebra, [1], []);;
ADRAlgebraOfBasicModule([module]);

