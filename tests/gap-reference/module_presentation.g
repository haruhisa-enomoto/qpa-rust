if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

SizeScreen([100000, 100000]);

a2 := PathAlgebra(Rationals,
    Quiver(["v1", "v2"], [["v1", "v2", "a"]]));;
source := RightModuleOverPathAlgebra(a2, [1, 0], []);;
range := RightModuleOverPathAlgebra(a2, [0, 1], []);;
zero := ZeroModule(a2);;
map := ZeroMapping(source, range);;
loop_parent := PathAlgebra(Rationals,
    Quiver(["w"], [["w", "w", "x"]]));;
dual := loop_parent / [loop_parent.x^2];;
dual_simple := SimpleModules(dual)[1];;
prime_a1 := PathAlgebra(GF(5), Quiver(["p"], []));;
prime_module := RightModuleOverPathAlgebra(prime_a1, [1], []);;
wide := RightModuleOverPathAlgebra(a2, [2, 2],
    [["a", [[1, -2], [0, 3]]]]);;

Print("MODULE_VIEW="); View(source); Print("\n");
Print("ZERO_VIEW="); View(zero); Print("\n");
Print("MODULE_PRINT="); Print(source); Print("\n");
Print("QUOTIENT_MODULE_PRINT="); Print(dual_simple); Print("\n");
Print("PRIME_MODULE_PRINT="); Print(prime_module); Print("\n");
Print("MAP_VIEW="); View(map); Print("\n");
Print("MAP_PRINT="); Print(map); Print("\n");
Print("MODULE_DISPLAY_BEGIN\n"); Display(source); Print("MODULE_DISPLAY_END\n");
Print("MAP_DISPLAY_BEGIN\n"); Display(map); Print("MAP_DISPLAY_END\n");
Print("WIDE_MODULE_DISPLAY_BEGIN\n");
Display(wide);
Print("WIDE_MODULE_DISPLAY_END\n");

QUIT_GAP(0);
