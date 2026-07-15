if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

quiver := Quiver(["v"], []);
algebra := PathAlgebra(Rationals, quiver);
one := RightModuleOverPathAlgebra(algebra, [1], []);
InverseOfIsomorphism(ZeroMapping(one, one));
