if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\\\", "\\\\\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\\\"");
    Print("\"", escaped, "\"");
end;

PrintIntegerList := function(values)
    local first, value;
    Print("["); first := true;
    for value in values do
        if not first then Print(","); fi;
        Print(value); first := false;
    od;
    Print("]");
end;

PrintScalarList := function(values)
    local first, value;
    Print("["); first := true;
    for value in values do
        if not first then Print(","); fi;
        PrintJsonString(String(value)); first := false;
    od;
    Print("]");
end;

PrintMatrix := function(matrix)
    local first, row;
    Print("{\"shape\":"); PrintIntegerList(DimensionsMat(matrix));
    Print(",\"rows\":["); first := true;
    for row in matrix do
        if not first then Print(","); fi;
        PrintScalarList(row); first := false;
    od;
    Print("]}");
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

PrintMap := function(map)
    Print("{\"source_dimension\":");
    PrintIntegerList(DimensionVector(Source(map)));
    Print(",\"range_dimension\":");
    PrintIntegerList(DimensionVector(Range(map)));
    Print(",\"matrices\":");
    PrintMatrixList(MatricesOfPathAlgebraMatModuleHomomorphism(map));
    Print("}");
end;

PrintCase := function(field)
    local quiver, algebra, two, automorphism, inverse, zero_module,
          zero_inverse, quiver_a2, algebra_a2, representation, source_map,
          basis, images, reconstructed, reconstructed_extra, p1, p2, matrix,
          i1, i2, product, product_extra, MakeMap;

    quiver := Quiver(["v"], []);
    algebra := PathAlgebra(field, quiver);
    two := RightModuleOverPathAlgebra(algebra, [2], []);
    automorphism := RightModuleHomOverAlgebra(
        two, two, [[[1, 1], [0, 1]] * One(field)]);
    inverse := InverseOfIsomorphism(automorphism);
    zero_module := ZeroModule(algebra);
    zero_inverse := InverseOfIsomorphism(ZeroMapping(zero_module, zero_module));

    quiver_a2 := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
    algebra_a2 := PathAlgebra(field, quiver_a2);
    representation := RightModuleOverPathAlgebra(
        algebra_a2, [2, 1], [["a", [[1], [0]] * One(field)]]);
    source_map := RightModuleHomOverAlgebra(
        representation, representation,
        [[[1, 1], [0, 1]] * One(field), [[1]] * One(field)]);
    basis := BasisVectors(Basis(representation));
    images := List(basis, x -> ImageElm(source_map, x));
    reconstructed := HomomorphismFromImages(representation, representation, images);
    reconstructed_extra := HomomorphismFromImages(
        representation, representation, Concatenation(images, [Zero(basis[1])]));

    MakeMap := function(matrix)
        return RightModuleHomOverAlgebra(two, two, [matrix * One(field)]);
    end;
    p1 := MakeMap([[1, 1], [0, 1]]);
    p2 := MakeMap([[1, 0], [1, 1]]);
    matrix := [
        [MakeMap([[1, 0], [0, 0]]), MakeMap([[1, 0], [0, 1]])],
        [MakeMap([[0, 1], [1, 0]]), MakeMap([[1, 1], [1, 0]])]
    ];
    i1 := MakeMap([[1, 0], [1, 1]]);
    i2 := MakeMap([[0, 1], [1, 0]]);
    product := MultiplyListsOfMaps([p1, p2], matrix, [i1, i2]);
    product_extra := MultiplyListsOfMaps(
        [p1, p2], Concatenation(matrix, [[p1, p2]]), [i1, i2]);

    Print("{\"inverse\":"); PrintMap(inverse);
    Print(",\"inverse_triangles\":[",
          automorphism * inverse = IdentityMapping(two), ",",
          inverse * automorphism = IdentityMapping(two), "]");
    Print(",\"zero_inverse\":"); PrintMap(zero_inverse);
    Print(",\"from_images\":"); PrintMap(reconstructed);
    Print(",\"from_images_equal\":", reconstructed = source_map);
    Print(",\"extra_image_ignored\":", reconstructed_extra = source_map);
    Print(",\"list_product\":"); PrintMap(product);
    Print(",\"extra_row_ignored\":", product_extra = product);
    Print("}");
end;

Print("{\"schema_version\":1,\"source\":{",
      "\"gap\":\"4.17dev@4efe877859ab0bc0ac041bf635a996fb897ee0c9\",",
      "\"qpa\":\"1.37@8d5f847258d4853700ce3a5d947db2a1441109f1\",",
      "\"gbnp\":\"1.1.0@e0b46c1aac16f3fbabe743e3169dc6bd86bdc654\"},",
      "\"rational\":");
PrintCase(Rationals);
Print(",\"binary\":");
PrintCase(GF(2));
Print("}\n");
