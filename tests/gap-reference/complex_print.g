if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

algebra := PathAlgebra(Rationals, Quiver(["v"], []));;
cat := CatOfRightAlgebraModules(algebra);;
one := RightModuleOverPathAlgebra(algebra, [1], []);;
two := RightModuleOverPathAlgebra(algebra, [2], []);;
f := RightModuleHomOverAlgebra(one, two, [[[1, 0]]]);;
g := RightModuleHomOverAlgebra(two, one, [[[0], [1]]]);;
ses := ShortExactSequence(cat, f, g);;
zeroEndo := ZeroMapping(one, one);;
periodic := Complex(cat, 0, [zeroEndo],
                    ["repeat", [zeroEndo]], ["repeat", [zeroEndo]]);;
identityChain := FiniteChainMap(ses, ses, -1,
    List([-1 .. 1], degree -> IdentityMapping(ObjectOfComplex(ses, degree))));;

Print(cat, "\n");
Print(ZeroComplex(cat), "\n");
Print(ses, "\n");
Print(periodic, "\n");
Print(identityChain, "\n");

QUIT_GAP(0);
