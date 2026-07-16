if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

field := GF(2);;
e12 := [[Zero(field), One(field)], [Zero(field), Zero(field)]];;
e21 := [[Zero(field), Zero(field)], [One(field), Zero(field)]];;
algebra := AlgebraWithOne(field, [e12, e21]);;
AlgebraAsQuiverAlgebra(algebra);

QUIT_GAP(0);
