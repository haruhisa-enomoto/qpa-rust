if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

SizeScreen([100000, 100000]);

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintIntegerList := function(values)
    Print("[", JoinStringsWithSeparator(List(values, String), ","), "]");
end;

PrintRow := function(row)
    PrintIntegerList(List(row, Int));
end;

PrintMatrix := function(matrix)
    local first, row;
    Print("["); first := true;
    for row in matrix do
        if not first then Print(","); fi;
        PrintRow(row); first := false;
    od;
    Print("]");
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

PrintSequence := function(sequence)
    if sequence = fail then
        Print("null");
        return;
    fi;
    Print("{\"monomorphism\":"); PrintMap(sequence[1]);
    Print(",\"epimorphism\":"); PrintMap(sequence[2]);
    Print("}");
end;

PrintCase := function(field)
    local a2, injectives, simples, cotilting;
    a2 := PathAlgebra(field,
        Quiver(["v1", "v2"], [["v1", "v2", "a"]]));
    injectives := IndecInjectiveModules(a2);
    simples := SimpleModules(a2);
    cotilting := DirectSumOfQPAModules([injectives[2], simples[2]]);
    if CotiltingModule(cotilting, 1) = false then
        Error("failed to install the A2 cotilting attribute");
    fi;

    Print("{\"simple_one\":");
    PrintSequence(AlmostSplitSequenceInPerpT(cotilting, simples[1]));
    Print(",\"projective_simple_two\":");
    PrintSequence(AlmostSplitSequenceInPerpT(cotilting, simples[2]));
    Print("}");
end;

Print("{\"schema_version\":1,\"source\":{");
Print("\"gap\":\"4.17dev@4efe877859ab0bc0ac041bf635a996fb897ee0c9\",");
Print("\"qpa\":\"1.37@8d5f847258d4853700ce3a5d947db2a1441109f1\",");
Print("\"gbnp\":\"1.1.0@e0b46c1aac16f3fbabe743e3169dc6bd86bdc654\"},");
Print("\"rational\":"); PrintCase(Rationals);
Print(",\"binary\":"); PrintCase(GF(2));
Print(",\"quinary\":"); PrintCase(GF(5));
Print("}\n");

QUIT_GAP(0);
