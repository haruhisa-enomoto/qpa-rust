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

PrintCase := function(field)
    local quiver, algebra, layer, simple_one, simple_two, zero_module,
          projective_one, inclusion, projection, right_inverse, left_inverse,
          radical_inclusion, top_projection, zero_inclusion, zero_projection,
          identity, identity_right, identity_left, property_inclusion,
          property_projection, property_values, property_right, property_left;
    quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
    algebra := PathAlgebra(field, quiver);
    layer := RightModuleOverPathAlgebra(algebra, [1, 1], []);
    simple_one := RightModuleOverPathAlgebra(algebra, [1, 0], []);
    simple_two := RightModuleOverPathAlgebra(algebra, [0, 1], []);
    zero_module := ZeroModule(algebra);
    projective_one := IndecProjectiveModules(algebra)[1];
    inclusion := RightModuleHomOverAlgebra(
        simple_one, layer, [[[1]] * One(field), NullMat(1,1,field)]);
    projection := RightModuleHomOverAlgebra(
        layer, simple_one, [[[1]] * One(field), NullMat(1,1,field)]);
    right_inverse := RightInverseOfHomomorphism(inclusion);
    left_inverse := LeftInverseOfHomomorphism(projection);
    radical_inclusion := RightModuleHomOverAlgebra(
        simple_two, projective_one,
        [NullMat(1,1,field), [[1]] * One(field)]);
    top_projection := RightModuleHomOverAlgebra(
        projective_one, simple_one,
        [[[1]] * One(field), NullMat(1,1,field)]);
    zero_inclusion := ZeroMapping(zero_module, layer);
    zero_projection := ZeroMapping(layer, zero_module);
    identity := IdentityMapping(layer);
    identity_right := RightInverseOfHomomorphism(identity);
    identity_left := LeftInverseOfHomomorphism(identity);
    property_inclusion := RightModuleHomOverAlgebra(
        simple_one, layer, [[[1]] * One(field), NullMat(1,1,field)]);
    property_projection := RightModuleHomOverAlgebra(
        layer, simple_one, [[[1]] * One(field), NullMat(1,1,field)]);
    property_values := [IsSplitMonomorphism(property_inclusion),
                        IsSplitEpimorphism(property_projection)];
    property_right := RightInverseOfHomomorphism(property_inclusion);
    property_left := LeftInverseOfHomomorphism(property_projection);

    Print("{\"right_inverse\":"); PrintMap(right_inverse);
    Print(",\"left_inverse\":"); PrintMap(left_inverse);
    Print(",\"triangles\":[",
          inclusion * right_inverse = IdentityMapping(simple_one), ",",
          left_inverse * projection = IdentityMapping(simple_one), "]");
    Print(",\"cache_identity\":[",
          IsIdenticalObj(RightInverseOfHomomorphism(inclusion), right_inverse), ",",
          IsIdenticalObj(LeftInverseOfHomomorphism(projection), left_inverse), ",",
          IsIdenticalObj(LeftInverseOfHomomorphism(right_inverse), inclusion), ",",
          IsIdenticalObj(RightInverseOfHomomorphism(left_inverse), projection), "]");
    Print(",\"split_properties\":[",
          IsSplitMonomorphism(inclusion), ",",
          IsSplitEpimorphism(right_inverse), ",",
          IsSplitEpimorphism(projection), ",",
          IsSplitMonomorphism(left_inverse), "]");
    Print(",\"nonsplit\":[",
          RightInverseOfHomomorphism(radical_inclusion) = false, ",",
          IsSplitMonomorphism(radical_inclusion), ",",
          LeftInverseOfHomomorphism(top_projection) = false, ",",
          IsSplitEpimorphism(top_projection), "]");
    Print(",\"zero_space_bug\":[",
          RightInverseOfHomomorphism(zero_inclusion) = false, ",",
          IsSplitMonomorphism(zero_inclusion), ",",
          LeftInverseOfHomomorphism(zero_projection) = false, ",",
          IsSplitEpimorphism(zero_projection), "]");
    Print(",\"identity_right\":"); PrintMap(identity_right);
    Print(",\"identity_left\":"); PrintMap(identity_left);
    Print(",\"identity_cache_same\":",
          IsIdenticalObj(identity_right, identity_left));
    Print(",\"property_values\":[", property_values[1], ",",
          property_values[2], "]");
    Print(",\"property_right\":"); PrintMap(property_right);
    Print(",\"property_left\":"); PrintMap(property_left);
    Print(",\"property_triangles\":[",
          property_inclusion * property_right = IdentityMapping(simple_one), ",",
          property_left * property_projection = IdentityMapping(simple_one), "]");
    Print("}");
end;

Print("{\"schema_version\":1,\"kind\":\"modulehom_split\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("GBNP")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"rational\":"); PrintCase(Rationals);
Print(",\"binary\":"); PrintCase(GF(2));
Print("}\n");

QUIT_GAP(0);
