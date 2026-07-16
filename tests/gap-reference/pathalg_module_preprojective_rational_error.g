if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

quiver := Quiver(2, [[1, 2, "a"]]);
algebra := PathAlgebra(Rationals, quiver);
simples := SimpleModules(algebra);
projectives := IndecProjectiveModules(algebra);
module := DirectSumOfQPAModules([simples[1], simples[2], projectives[1]]);
PreprojectiveAlgebra(module, 1);
