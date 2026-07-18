if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

SizeScreen([100000, 100000]);

MatrixOf := function(mapping)
    return mapping!.matrix;
end;

NonzeroPositions := function(list)
    return Filtered([1..Length(list)], position -> not IsZero(list[position]));
end;

SupportOfMatrix := function(matrix)
    return NonzeroPositions(Concatenation(matrix));
end;

SupportsOfMappings := function(space)
    return List(GeneratorsOfLeftModule(space),
        mapping -> SupportOfMatrix(MatrixOf(mapping)));
end;

RunCase := function(field)
    local algebra, module, homspace, endalgebra, basis, generator, submodule,
          subhom, subend, to_submodule, from_submodule;
    algebra := PathAlgebra(field,
        Quiver(["v1", "v2"], [["v1", "v2", "a"]]));
    module := RightModuleOverPathAlgebra(algebra, [2, 2],
        [["a", [[One(field), Zero(field)],
                 [Zero(field), One(field)]]]]);
    homspace := Hom(algebra, module, module);
    endalgebra := End(algebra, module);
    basis := BasisVectors(Basis(module));
    generator := basis[1] + basis[4];
    submodule := SubAlgebraModule(module, [generator]);
    subhom := Hom(algebra, submodule, submodule);
    subend := End(algebra, submodule);
    to_submodule := Hom(algebra, module, submodule);
    from_submodule := Hom(algebra, submodule, module);

    return [
        Dimension(homspace),
        SupportsOfMappings(homspace),
        List(GeneratorsOfLeftModule(homspace),
            IsAlgebraModuleHomomorphism),
        Dimension(endalgebra),
        SupportsOfMappings(endalgebra),
        List(BasisVectors(CanonicalBasis(submodule)),
            element -> NonzeroPositions(
                Concatenation(ExtRepOfObj(element)![1]))),
        Dimension(subhom),
        SupportsOfMappings(subhom),
        Dimension(subend),
        SupportsOfMappings(subend),
        SupportsOfMappings(to_submodule),
        SupportsOfMappings(from_submodule)
    ];
end;

Print("QQ=", RunCase(Rationals), "\n");
Print("GF2=", RunCase(GF(2)), "\n");

QUIT_GAP(0);
