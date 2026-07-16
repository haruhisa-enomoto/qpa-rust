if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

field := GF(2);;
space := FullRowSpace(field, 2);;
basis := Basis(space);;
identity := LeftModuleHomomorphismByMatrix(
    basis, [[1, 0], [0, 1]], basis);;
nilpotent := LeftModuleHomomorphismByMatrix(
    basis, [[0, 1], [0, 0]], basis);;
algebra := AlgebraWithOneByGenerators(
    field, [identity, nilpotent]);;

RadicalOfAlgebra(algebra);
QUIT_GAP(0);
