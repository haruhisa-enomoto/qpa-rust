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

PrintModule := function(module)
    local first, entry;
    Print("{\"dimension\":"); PrintIntegerList(DimensionVector(module));
    Print(",\"matrices\":["); first := true;
    for entry in MatricesOfPathAlgebraModule(module) do
        if not first then Print(","); fi;
        if IsMatrix(entry) then
            PrintMatrix(entry);
        else
            PrintJsonString(String(entry));
        fi;
        first := false;
    od;
    Print("]}");
end;

PrintMap := function(map)
    local first, matrix;
    Print("{\"source_dimension\":");
    PrintIntegerList(DimensionVector(Source(map)));
    Print(",\"range_dimension\":");
    PrintIntegerList(DimensionVector(Range(map)));
    Print(",\"matrices\":["); first := true;
    for matrix in MatricesOfPathAlgebraMatModuleHomomorphism(map) do
        if not first then Print(","); fi;
        PrintMatrix(matrix); first := false;
    od;
    Print("]}");
end;

PrintCase := function(field)
    local quiver, algebra, zero, zero_again, simple1, simple2,
          syzygy1, syzygy2, syzygy0, envelope1, envelope2, envelope0,
          loop_quiver, parent, x, quotient, quotient_zero, quotient_simple,
          quotient_syzygy, quotient_envelope;

    quiver := Quiver(["v1", "v2"], [["v1", "v2", "a"]]);
    algebra := PathAlgebra(field, quiver);
    zero := ZeroModule(algebra);
    zero_again := ZeroModule(algebra);
    simple1 := RightModuleOverPathAlgebra(algebra, [1, 0], []);
    simple2 := RightModuleOverPathAlgebra(algebra, [0, 1], []);
    syzygy1 := 1stSyzygy(simple1);
    syzygy2 := 1stSyzygy(simple2);
    syzygy0 := 1stSyzygy(zero);
    envelope1 := InjectiveEnvelope(simple1);
    envelope2 := InjectiveEnvelope(simple2);
    envelope0 := InjectiveEnvelope(zero);

    loop_quiver := Quiver(["w"], [["w", "w", "x"]]);
    parent := PathAlgebra(field, loop_quiver);
    x := parent.x;
    quotient := parent / [x * x];
    quotient_zero := ZeroModule(quotient);
    quotient_simple := RightModuleOverPathAlgebra(quotient, [1], []);
    quotient_syzygy := 1stSyzygy(quotient_simple);
    quotient_envelope := InjectiveEnvelope(quotient_simple);

    Print("{\"path\":{");
    Print("\"zero\":"); PrintModule(zero);
    Print(",\"zero_cache_identity\":", IsIdenticalObj(zero, zero_again));
    Print(",\"syzygy_simple1\":"); PrintModule(syzygy1);
    Print(",\"syzygy_simple2\":"); PrintModule(syzygy2);
    Print(",\"syzygy_zero\":"); PrintModule(syzygy0);
    Print(",\"syzygy_cache_identity\":[",
          IsIdenticalObj(syzygy1, 1stSyzygy(simple1)), ",",
          IsIdenticalObj(syzygy2, 1stSyzygy(simple2)), ",",
          IsIdenticalObj(syzygy0, zero), "]");
    Print(",\"injective_simple1\":"); PrintMap(envelope1);
    Print(",\"injective_simple2\":"); PrintMap(envelope2);
    Print(",\"injective_zero\":"); PrintMap(envelope0);
    Print(",\"injective_cache_identity\":[",
          IsIdenticalObj(envelope1, InjectiveEnvelope(simple1)), ",",
          IsIdenticalObj(envelope2, InjectiveEnvelope(simple2)), ",",
          IsIdenticalObj(envelope0, InjectiveEnvelope(zero)), "]");
    Print(",\"injective_source_identity\":[",
          IsIdenticalObj(Source(envelope1), simple1), ",",
          IsIdenticalObj(Source(envelope2), simple2), ",",
          IsIdenticalObj(Source(envelope0), zero), "]");
    Print("},\"quotient\":{");
    Print("\"zero\":"); PrintModule(quotient_zero);
    Print(",\"syzygy_simple\":"); PrintModule(quotient_syzygy);
    Print(",\"injective_simple\":"); PrintMap(quotient_envelope);
    Print(",\"cache_identity\":[",
          IsIdenticalObj(quotient_zero, ZeroModule(quotient)), ",",
          IsIdenticalObj(quotient_syzygy, 1stSyzygy(quotient_simple)), ",",
          IsIdenticalObj(quotient_envelope, InjectiveEnvelope(quotient_simple)),
          "]}}");
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
