if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintStringList := function(values)
    local first, value;
    Print("["); first := true;
    for value in values do
        if not first then Print(","); fi;
        PrintJsonString(String(value)); first := false;
    od;
    Print("]");
end;

PrintCase := function(field)
    local zero, one, values;
    zero := Zero(field); one := One(field);
    values := [zero, one + one, zero, -one, zero];
    if Characteristic(field) = 2 then
        values := [zero, one, zero, one, zero];
    fi;
    Print("{\"positions\":", PositionsNonzero(values),
          ",\"null_list\":");
    PrintStringList(NullList(4, field));
    Print("}");
end;

Print("{\"schema_version\":1,\"source\":{",
      "\"gap\":\"4.17dev@4efe877859ab0bc0ac041bf635a996fb897ee0c9\",",
      "\"qpa\":\"1.37@8d5f847258d4853700ce3a5d947db2a1441109f1\"},",
      "\"rational\":");
PrintCase(Rationals);
Print(",\"binary\":"); PrintCase(GF(2));
Print(",\"generic_matrices\":{\"positions\":",
      PositionsNonzero([[[0]], [[1]], [[0]], [[-2]]]), "}}\n");

QUIT_GAP(0);
