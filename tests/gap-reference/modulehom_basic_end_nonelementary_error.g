if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

algebra := PathAlgebra(GF(2), Quiver(["x"], []));;
simple := RightModuleOverPathAlgebra(algebra, [1], []);;
repeated := DirectSumOfQPAModules([simple, simple]);;
EndOfBasicModuleAsQuiverAlgebra([repeated]);

