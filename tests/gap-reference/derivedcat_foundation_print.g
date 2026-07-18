if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

algebra := PathAlgebra(Rationals,
    Quiver(["1", "2"], [["1", "2", "a"]]));
cat := CatOfRightAlgebraModules(algebra);
projectives := IndecProjectiveModules(algebra);
injectives := IndecInjectiveModules(algebra);
projectiveSum := DirectSumOfQPAModules(
    [projectives[1], projectives[2], projectives[1]]);
injectiveSum := DirectSumOfQPAModules(
    [injectives[2], injectives[1], injectives[2]]);
projectiveComplex := StalkComplex(cat, projectiveSum, 1);
injectiveComplex := StalkComplex(cat, injectiveSum, -1);
IsProjectiveComplex(projectiveComplex);
IsInjectiveComplex(injectiveComplex);
Print(projectiveComplex, "\n");
Print(injectiveComplex, "\n");

QUIT_GAP(0);
