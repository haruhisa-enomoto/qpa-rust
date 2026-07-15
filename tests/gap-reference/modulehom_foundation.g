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
    Print("[");
    first := true;
    for value in values do
        if not first then Print(","); fi;
        Print(value);
        first := false;
    od;
    Print("]");
end;

PrintMatrix := function(matrix)
    local first_row, first_entry, row, entry;
    Print("{\"shape\":");
    PrintIntegerList(DimensionsMat(matrix));
    Print(",\"rows\":[");
    first_row := true;
    for row in matrix do
        if not first_row then Print(","); fi;
        Print("[");
        first_entry := true;
        for entry in row do
            if not first_entry then Print(","); fi;
            PrintJsonString(String(entry));
            first_entry := false;
        od;
        Print("]");
        first_row := false;
    od;
    Print("]}");
end;

PrintMatrixList := function(matrices)
    local first, matrix;
    Print("[");
    first := true;
    for matrix in matrices do
        if not first then Print(","); fi;
        PrintMatrix(matrix);
        first := false;
    od;
    Print("]");
end;

PrintHomBasis := function(homomorphisms)
    local first, homomorphism;
    Print("[");
    first := true;
    for homomorphism in homomorphisms do
        if not first then Print(","); fi;
        PrintMatrixList(MatricesOfPathAlgebraMatModuleHomomorphism(homomorphism));
        first := false;
    od;
    Print("]");
end;

PrintScalarList := function(values)
    local first, value;
    Print("[");
    first := true;
    for value in values do
        if not first then Print(","); fi;
        PrintJsonString(String(value));
        first := false;
    od;
    Print("]");
end;

PrintCase := function(field)
    local quiver, algebra, module, map, kernel, factorization, projection,
          inclusion, cokernel, zero, identity, top_map, hom_basis,
          hom_with_basis, free_quiver, free_algebra, free_source, free_range,
          free_basis, disjoint_source, disjoint_range, disjoint_basis,
          disjoint_with_basis, endomorphism_algebra, endomorphism_basis;
    quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
    algebra := PathAlgebra(field, quiver);
    module := RightModuleOverPathAlgebra(
        algebra, [2, 2], [["a", [[1,1],[0,0]] * One(field)]]);
    map := RightModuleHomOverAlgebra(
        module, module,
        [[[1,0],[0,0]] * One(field), IdentityMat(2, field)]);
    kernel := KernelInclusion(map);
    factorization := ImageProjectionInclusion(map);
    projection := factorization[1];
    inclusion := factorization[2];
    cokernel := CoKernelProjection(map);
    zero := Zero(map);
    identity := IdentityMapping(module);
    top_map := TopOfModule(map);
    hom_basis := HomOverAlgebra(module, module);
    hom_with_basis := HomOverAlgebraWithBasisFunction(module, module);
    free_quiver := Quiver(["w"], []);
    free_algebra := PathAlgebra(field, free_quiver);
    free_source := RightModuleOverPathAlgebra(free_algebra, [1], []);
    free_range := RightModuleOverPathAlgebra(free_algebra, [2], []);
    free_basis := HomOverAlgebra(free_source, free_range);
    disjoint_source := RightModuleOverPathAlgebra(algebra, [1, 0], []);
    disjoint_range := RightModuleOverPathAlgebra(algebra, [0, 1], []);
    disjoint_basis := HomOverAlgebra(disjoint_source, disjoint_range);
    disjoint_with_basis :=
        HomOverAlgebraWithBasisFunction(disjoint_source, disjoint_range);
    endomorphism_algebra := EndOverAlgebra(module);
    endomorphism_basis := BasisVectors(Basis(endomorphism_algebra));

    Print("{\"predicates\":{");
    Print("\"injective\":", IsInjective(map), ",");
    Print("\"surjective\":", IsSurjective(map), ",");
    Print("\"isomorphism\":", IsIsomorphism(map), ",");
    Print("\"zero\":", IsZero(map), ",");
    Print("\"identity_isomorphism\":", IsIsomorphism(identity));
    Print("},\"kernel\":{");
    Print("\"dimension\":");
    PrintIntegerList(DimensionVector(Source(kernel)));
    Print(",\"matrices\":");
    PrintMatrixList(MatricesOfPathAlgebraMatModuleHomomorphism(kernel));
    Print(",\"cached_identity\":", IsIdenticalObj(KernelInclusion(map), kernel));
    Print(",\"of_what_identity\":", IsIdenticalObj(KernelOfWhat(kernel), map));
    Print(",\"composite_zero\":", IsZero(kernel * map));
    Print("},\"image\":{");
    Print("\"dimension\":");
    PrintIntegerList(DimensionVector(Range(projection)));
    Print(",\"projection\":");
    PrintMatrixList(MatricesOfPathAlgebraMatModuleHomomorphism(projection));
    Print(",\"inclusion\":");
    PrintMatrixList(MatricesOfPathAlgebraMatModuleHomomorphism(inclusion));
    Print(",\"projection_cached_identity\":",
          IsIdenticalObj(ImageProjection(map), projection));
    Print(",\"inclusion_cached_identity\":",
          IsIdenticalObj(ImageInclusion(map), inclusion));
    Print(",\"projection_of_what_identity\":",
          IsIdenticalObj(ImageOfWhat(projection), map));
    Print(",\"inclusion_of_what_identity\":",
          IsIdenticalObj(ImageOfWhat(inclusion), map));
    Print(",\"factorization_equal\":", projection * inclusion = map);
    Print(",\"images_set_length\":", Length(ImagesSet(map, Source(map))));
    Print("},\"top_map\":{");
    Print("\"source_dimension\":");
    PrintIntegerList(DimensionVector(Source(top_map)));
    Print(",\"range_dimension\":");
    PrintIntegerList(DimensionVector(Range(top_map)));
    Print(",\"matrices\":");
    PrintMatrixList(MatricesOfPathAlgebraMatModuleHomomorphism(top_map));
    Print("},\"cokernel\":{");
    Print("\"dimension\":");
    PrintIntegerList(DimensionVector(Range(cokernel)));
    Print(",\"matrices\":");
    PrintMatrixList(MatricesOfPathAlgebraMatModuleHomomorphism(cokernel));
    Print(",\"cached_identity\":", IsIdenticalObj(CoKernelProjection(map), cokernel));
    Print(",\"of_what_identity\":", IsIdenticalObj(CoKernelOfWhat(cokernel), map));
    Print(",\"composite_zero\":", IsZero(map * cokernel));
    Print("},\"arithmetic\":{");
    Print("\"add_zero_equal\":", map + zero = map, ",");
    Print("\"inverse_sum_zero\":", IsZero(map + (-map)), ",");
    Print("\"left_identity_equal\":", identity * map = map, ",");
    Print("\"right_identity_equal\":", map * identity = map);
    Print("},\"hom_space\":{");
    Print("\"dimension\":", Length(hom_basis), ",");
    Print("\"basis\":");
    PrintHomBasis(hom_basis);
    Print(",\"map_coefficients\":");
    PrintScalarList(hom_with_basis[2](map));
    Print(",\"free_dimension\":", Length(free_basis), ",");
    Print("\"free_basis\":");
    PrintHomBasis(free_basis);
    Print(",\"disjoint_dimension\":", Length(disjoint_basis), ",");
    Print("\"disjoint_zero_coefficients\":");
    PrintScalarList(disjoint_with_basis[2](
        ZeroMapping(disjoint_source, disjoint_range)));
    Print("},\"endomorphism_algebra\":{");
    Print("\"dimension\":", Dimension(endomorphism_algebra), ",");
    Print("\"cached_identity\":",
          IsIdenticalObj(EndOverAlgebra(module), endomorphism_algebra), ",");
    Print("\"basis\":");
    PrintMatrixList(endomorphism_basis);
    Print(",\"one\":");
    PrintMatrix(One(endomorphism_algebra));
    Print(",\"one_coefficients\":");
    PrintScalarList(Coefficients(Basis(endomorphism_algebra),
                                 One(endomorphism_algebra)));
    Print("}}");
end;

Print("{\"schema_version\":1,\"kind\":\"modulehom_foundation\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":");
PrintJsonString(PackageInfo("GBNP")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"rational\":");
PrintCase(Rationals);
Print(",\"binary\":");
PrintCase(GF(2));
Print("}\n");

QUIT_GAP(0);
