if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

SizeScreen([100000, 100000]);

RunCase := function(field)
    local algebra, module, basis, generator, submodule, standalone, projection;
    algebra := PathAlgebra(field,
        Quiver(["v1", "v2"], [["v1", "v2", "a"]]));
    module := RightModuleOverPathAlgebra(algebra, [3, 3],
        [["a", [[One(field), Zero(field), Zero(field)],
                 [Zero(field), One(field), Zero(field)],
                 [Zero(field), Zero(field), One(field)]]]]);
    basis := BasisVectors(Basis(module));
    generator := basis[1] + basis[5];
    submodule := SubAlgebraModule(module, [generator, Zero(module)]);
    standalone := SubmoduleAsModule(submodule);
    projection := NaturalHomomorphismBySubAlgebraModule(module, submodule);

    return [
        List(GeneratorsOfAlgebraModule(submodule), ExtRepOfObj),
        List(BasisVectors(CanonicalBasis(submodule)), ExtRepOfObj),
        DimensionVector(standalone),
        MatricesOfPathAlgebraModule(standalone),
        DimensionVector(Range(projection)),
        MatricesOfPathAlgebraModule(Range(projection)),
        List(BasisVectors(Basis(module)),
            x -> ExtRepOfObj(Image(projection, x))),
        IsSurjective(projection)
    ];
end;

Print("QQ=", RunCase(Rationals), "\n");
Print("GF2=", RunCase(GF(2)), "\n");

QUIT_GAP(0);
