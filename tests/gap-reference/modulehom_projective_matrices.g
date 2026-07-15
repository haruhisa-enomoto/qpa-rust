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

PrintPathElement := function(element)
    local terms, first, i;
    terms := ExtRepOfObj(element)[2];
    Print("["); first := true;
    for i in [1, 3 .. Length(terms) - 1] do
        if not first then Print(","); fi;
        Print("{\"walk\":"); PrintIntegerList(terms[i]);
        Print(",\"coefficient\":"); PrintJsonString(String(terms[i + 1]));
        Print("}"); first := false;
    od;
    Print("]");
end;

PrintPathElementMatrix := function(matrix)
    local firstrow, firstentry, row, entry;
    Print("["); firstrow := true;
    for row in matrix do
        if not firstrow then Print(","); fi;
        Print("["); firstentry := true;
        for entry in row do
            if not firstentry then Print(","); fi;
            PrintPathElement(entry); firstentry := false;
        od;
        Print("]"); firstrow := false;
    od;
    Print("]");
end;

PrintModuleElement := function(module, element)
    local support;
    support := List(SupportModuleElement(element), VertexPosition);
    Print("{\"dimension\":"); PrintIntegerList(DimensionVector(module));
    Print(",\"support\":"); PrintIntegerList(support);
    Print(",\"coordinates\":");
    PrintScalarList(Flat(ExtRepOfObj(ExtRepOfObj(element))));
    Print("}");
end;

PrintCase := function(field)
    local q, algebra, vertices, arrows, v1, v2, a, zero, matrix, map,
          back, uniform, projective_element, recovered, projectives, sum,
          zero_module, zero_to_sum, sum_to_zero, zero_to_zero, single,
          q3, parent, qa, qb, quotient, qvertices, qv1, qv2, qv3, qzero,
          quotient_matrix, quotient_map, quotient_back, quotient_uniform,
          quotient_projective_element, quotient_recovered;

    q := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
    algebra := PathAlgebra(field, q);
    vertices := VerticesOfQuiver(q);
    arrows := ArrowsOfQuiver(q);
    v1 := One(algebra) * vertices[1];
    v2 := One(algebra) * vertices[2];
    a := One(algebra) * arrows[1];
    zero := Zero(algebra);
    matrix := [[v1, a], [zero, v2]];
    map := FromMatrixToHomomorphismOfProjectives(algebra, matrix, [1, 2], [1, 2]);
    back := MatrixOfHomomorphismBetweenProjectives(map);
    uniform := v1 + a;
    projective_element := ElementIn_vA_AsElementInIndecProj(algebra, uniform);
    recovered := ElementInIndecProjective(algebra, projective_element, 1);
    projectives := IndecProjectiveModules(algebra);
    sum := DirectSumOfQPAModules(projectives);
    zero_module := ZeroModule(algebra);
    zero_to_sum := MatrixOfHomomorphismBetweenProjectives(ZeroMapping(zero_module, sum));
    sum_to_zero := MatrixOfHomomorphismBetweenProjectives(ZeroMapping(sum, zero_module));
    zero_to_zero := MatrixOfHomomorphismBetweenProjectives(
        ZeroMapping(zero_module, zero_module));
    single := FromMatrixToHomomorphismOfProjectives(algebra, [[a]], [2], [1]);

    q3 := Quiver(["u1", "u2", "u3"],
                  [["u1", "u2", "x"], ["u2", "u3", "y"]]);
    parent := PathAlgebra(field, q3);
    qa := parent.x;
    qb := parent.y;
    quotient := parent / [qa * qb];
    qvertices := VerticesOfQuiver(q3);
    qv1 := One(quotient) * qvertices[1];
    qv2 := One(quotient) * qvertices[2];
    qv3 := One(quotient) * qvertices[3];
    qa := One(quotient) * ArrowsOfQuiver(q3)[1];
    qb := One(quotient) * ArrowsOfQuiver(q3)[2];
    qzero := Zero(quotient);
    quotient_matrix := [[qa, qzero], [qv2, qb]];
    quotient_map := FromMatrixToHomomorphismOfProjectives(
        quotient, quotient_matrix, [2, 3], [1, 2]);
    quotient_back := MatrixOfHomomorphismBetweenProjectives(quotient_map);
    quotient_uniform := qv1 + qa;
    quotient_projective_element := ElementIn_vA_AsElementInIndecProj(
        quotient, quotient_uniform);
    quotient_recovered := ElementInIndecProjective(
        quotient, quotient_projective_element, 1);

    Print("{\"path\":{");
    Print("\"map\":"); PrintMap(map);
    Print(",\"matrix\":"); PrintPathElementMatrix(back);
    Print(",\"roundtrip_equal\":", back = matrix);
    Print(",\"underlying\":"); PrintMatrix(UnderlyingLinearMap(map));
    Print(",\"uniform_projective_element\":");
    PrintModuleElement(projectives[1], projective_element);
    Print(",\"uniform_recovered\":"); PrintPathElement(recovered);
    Print(",\"single\":"); PrintMap(single);
    Print(",\"single_underlying\":"); PrintMatrix(UnderlyingLinearMap(single));
    Print(",\"zero_to_sum\":"); PrintPathElementMatrix(zero_to_sum);
    Print(",\"sum_to_zero\":"); PrintPathElementMatrix(sum_to_zero);
    Print(",\"zero_to_zero\":"); PrintPathElementMatrix(zero_to_zero);
    Print("},\"quotient\":{");
    Print("\"map\":"); PrintMap(quotient_map);
    Print(",\"matrix_parent_representatives\":");
    PrintPathElementMatrix(quotient_back);
    Print(",\"underlying\":"); PrintMatrix(UnderlyingLinearMap(quotient_map));
    Print(",\"uniform_projective_element\":");
    PrintModuleElement(IndecProjectiveModules(quotient)[1], quotient_projective_element);
    Print(",\"uniform_recovered_parent\":"); PrintPathElement(quotient_recovered);
    Print("}}");
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
