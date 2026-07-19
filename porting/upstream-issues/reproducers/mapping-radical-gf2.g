if LoadPackage("qpa") <> true then Error("failed to load QPA"); fi;

space := FullRowSpace(GF(2), 2);;
basis := Basis(space);;
identity := LeftModuleHomomorphismByMatrix(basis, [[1, 0], [0, 1]], basis);;
nilpotent := LeftModuleHomomorphismByMatrix(basis, [[0, 1], [0, 0]], basis);;
algebra := AlgebraWithOneByGenerators(GF(2), [identity, nilpotent]);;
RadicalOfAlgebra(algebra);
