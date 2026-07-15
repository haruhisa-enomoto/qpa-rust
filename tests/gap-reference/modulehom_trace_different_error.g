if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

firstAlgebra := PathAlgebra(Rationals, Quiver(["v"], []));
secondAlgebra := PathAlgebra(Rationals, Quiver(["w"], []));
firstModule := RightModuleOverPathAlgebra(firstAlgebra, [1], []);
secondModule := RightModuleOverPathAlgebra(secondAlgebra, [1], []);
TraceOfModule(firstModule, secondModule);

QUIT_GAP(0);
