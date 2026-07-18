if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

cycle := PathAlgebra(Rationals, Quiver(2, [[1,2,"a"],[2,1,"b"]]));
radical_square_zero := cycle/[cycle.a*cycle.b,cycle.b*cycle.a];
TitsUnitFormOfAlgebra(radical_square_zero);

QUIT_GAP(0);
