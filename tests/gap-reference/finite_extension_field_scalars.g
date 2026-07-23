if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintIntegerList := function(values)
    Print("[", JoinStringsWithSeparator(List(values, String), ","), "]");
end;

PrintIntegerMatrix := function(matrix)
    local first, row;
    Print("["); first := true;
    for row in matrix do
        if not first then Print(","); fi;
        PrintIntegerList(row); first := false;
    od;
    Print("]");
end;

PrintScalarCase := function(order)
    local field, characteristic, degree, generator, polynomial, vector_space,
          basis, coordinates, subdegrees, first, subdegree, exponent;
    field := GF(order);
    characteristic := Characteristic(field);
    generator := Z(order);
    degree := DegreeFFE(generator);
    polynomial := ConwayPolynomial(characteristic, degree);
    vector_space := AsVectorSpace(GF(characteristic), field);
    basis := Basis(vector_space,
        List([0..degree - 1], power -> generator^power));
    coordinates := element -> List(Coefficients(basis, element), IntFFE);
    subdegrees := Filtered(DivisorsInt(degree),
        subdegree -> 1 < subdegree and subdegree < degree);

    Print("{\"order\":", order);
    Print(",\"characteristic\":", characteristic);
    Print(",\"degree\":", degree);
    Print(",\"conway\":");
    PrintIntegerList(List(CoefficientsOfUnivariatePolynomial(polynomial), IntFFE));
    Print(",\"generator_order\":", Order(generator));
    Print(",\"powers\":");
    PrintIntegerMatrix(List([0..order - 2], power -> coordinates(generator^power)));
    Print(",\"inverse_generator\":");
    PrintIntegerList(coordinates(Inverse(generator)));
    Print(",\"frobenius_generator\":");
    PrintIntegerList(coordinates(generator^characteristic));
    Print(",\"proper_extension_subfields\":["); first := true;
    for subdegree in subdegrees do
        if not first then Print(","); fi;
        exponent := (order - 1) / (characteristic^subdegree - 1);
        Print("{\"degree\":", subdegree, ",\"generator\":");
        PrintIntegerList(coordinates(generator^exponent));
        Print("}"); first := false;
    od;
    Print("]}");
end;

Print("{\"schema_version\":1,\"kind\":\"finite_extension_field_scalars\",\"cases\":[");
PrintScalarCase(4); Print(","); PrintScalarCase(8); Print(",");
PrintScalarCase(9); Print(","); PrintScalarCase(16); Print(",");
PrintScalarCase(25); Print("]}\n");

QUIT_GAP(0);
