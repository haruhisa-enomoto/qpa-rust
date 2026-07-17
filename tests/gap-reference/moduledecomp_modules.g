if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintIntegerList := function(values)
    local first, value;
    Print("["); first := true;
    for value in values do
        if not first then Print(","); fi;
        Print(Int(value)); first := false;
    od;
    Print("]");
end;

PrintMatrix := function(matrix)
    local first, row;
    Print("["); first := true;
    for row in matrix do
        if not first then Print(","); fi;
        PrintIntegerList(row); first := false;
    od;
    Print("]");
end;

PrintMatrixList := function(matrices)
    local first, matrix;
    Print("["); first := true;
    for matrix in matrices do
        if not first then Print(","); fi;
        PrintMatrix(matrix); first := false;
    od;
    Print("]");
end;

PrintModules := function(modules)
    local first, module;
    Print("["); first := true;
    for module in modules do
        if not first then Print(","); fi;
        Print("{\"dimension\":"); PrintIntegerList(DimensionVector(module));
        Print(",\"matrices\":");
        PrintMatrixList(MatricesOfPathAlgebraModule(module));
        Print("}"); first := false;
    od;
    Print("]");
end;

PrintInclusions := function(inclusions)
    local first, inclusion;
    Print("["); first := true;
    for inclusion in inclusions do
        if not first then Print(","); fi;
        Print("{\"source_dimension\":");
        PrintIntegerList(DimensionVector(Source(inclusion)));
        Print(",\"matrices\":");
        PrintMatrixList(MatricesOfPathAlgebraMatModuleHomomorphism(inclusion));
        Print("}"); first := false;
    od;
    Print("]");
end;

PrintHomomorphisms := function(homomorphisms)
    local first, homomorphism;
    Print("["); first := true;
    for homomorphism in homomorphisms do
        if not first then Print(","); fi;
        PrintMatrixList(
            MatricesOfPathAlgebraMatModuleHomomorphism(homomorphism));
        first := false;
    od;
    Print("]");
end;

PrintCase := function(field)
    local quiver, algebra, module, decomposition, inclusions, simple_one,
          simple_two, direct_sum, direct_decomposition, endomorphisms,
          endomorphism_basis, transposed_algebra, transposed_basis, idempotents,
          mixed_multiplicities, vector_quiver, vector_algebra, repeated_module,
          repeated_multiplicities, zero_module, block_idempotents,
          block_decomposition, repeated_block_idempotents,
          repeated_block_decomposition, basic_mixed, basic_repeated, basic_zero;
    quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
    algebra := PathAlgebra(field, quiver);
    module := RightModuleOverPathAlgebra(
        algebra, [2, 1], [["a", [[1], [0]] * One(field)]]);
    endomorphisms := EndOverAlgebra(module);
    endomorphism_basis := BasisVectors(Basis(endomorphisms));
    transposed_algebra := AlgebraWithOne(field,
        List(endomorphism_basis, TransposedMat));
    transposed_basis := BasisVectors(Basis(transposed_algebra));
    idempotents := IdempotentsForDecomposition(transposed_algebra);
    decomposition := DecomposeModule(module);
    inclusions := DecomposeModuleWithInclusions(module);

    simple_one := RightModuleOverPathAlgebra(algebra, [1, 0], []);
    simple_two := RightModuleOverPathAlgebra(algebra, [0, 1], []);
    direct_sum := DirectSumOfQPAModules([simple_one, simple_two]);
    direct_decomposition := DecomposeModule(direct_sum);
    mixed_multiplicities := DecomposeModuleWithMultiplicities(module);

    vector_quiver := Quiver(["w"], []);
    vector_algebra := PathAlgebra(field, vector_quiver);
    repeated_module := RightModuleOverPathAlgebra(vector_algebra, [2], []);
    repeated_multiplicities :=
        DecomposeModuleWithMultiplicities(repeated_module);
    zero_module := RightModuleOverPathAlgebra(vector_algebra, [0], []);
    block_idempotents := BlockSplittingIdempotents(module);
    block_decomposition := BlockDecompositionOfModule(module);
    repeated_block_idempotents := BlockSplittingIdempotents(repeated_module);
    repeated_block_decomposition := BlockDecompositionOfModule(repeated_module);
    basic_mixed := BasicVersionOfModule(module);
    basic_repeated := BasicVersionOfModule(repeated_module);
    basic_zero := BasicVersionOfModule(zero_module);

    Print("{\"algorithm\":{");
    Print("\"endomorphism_basis\":"); PrintMatrixList(endomorphism_basis);
    Print(",\"transposed_basis\":"); PrintMatrixList(transposed_basis);
    Print(",\"idempotents\":"); PrintMatrixList(idempotents);
    Print("}");
    Print(",\"decomposition\":"); PrintModules(decomposition);
    Print(",\"inclusions\":"); PrintInclusions(inclusions);
    Print(",\"cache\":{");
    Print("\"decomposition_first\":",
          IsIdenticalObj(DecomposeModule(module)[1], decomposition[1]));
    Print(",\"inclusion_first\":",
          IsIdenticalObj(DecomposeModuleWithInclusions(module)[1], inclusions[1]));
    Print("}");
    Print(",\"direct_sum\":"); PrintModules(direct_decomposition);
    Print(",\"multiplicities\":{");
    Print("\"mixed\":{\"basic\":"); PrintModules(mixed_multiplicities[1]);
    Print(",\"values\":"); PrintIntegerList(mixed_multiplicities[2]);
    Print("}");
    Print(",\"repeated\":{\"basic\":");
    PrintModules(repeated_multiplicities[1]);
    Print(",\"values\":"); PrintIntegerList(repeated_multiplicities[2]);
    Print(",\"cache_first\":",
          IsIdenticalObj(DecomposeModuleWithMultiplicities(repeated_module)[1][1],
                         repeated_multiplicities[1][1]));
    Print("}");
    Print(",\"zero_is_fail\":",
          DecomposeModuleWithMultiplicities(zero_module) = fail);
    Print("}");
    Print(",\"blocks\":{");
    Print("\"mixed_idempotents\":"); PrintHomomorphisms(block_idempotents);
    Print(",\"mixed_modules\":"); PrintModules(block_decomposition);
    Print(",\"repeated_idempotents\":");
    PrintHomomorphisms(repeated_block_idempotents);
    Print(",\"repeated_modules\":"); PrintModules(repeated_block_decomposition);
    Print(",\"zero_idempotents_empty\":",
          BlockSplittingIdempotents(zero_module) = []);
    Print("}");
    Print(",\"basic\":{");
    Print("\"mixed\":"); PrintModules([basic_mixed]);
    Print(",\"repeated\":"); PrintModules([basic_repeated]);
    Print(",\"zero_identical\":", IsIdenticalObj(basic_zero, zero_module));
    Print("}");
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"moduledecomp_modules\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("GBNP")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"binary\":"); PrintCase(GF(2));
Print(",\"quinary\":"); PrintCase(GF(5));
Print("}\n");

QUIT_GAP(0);
