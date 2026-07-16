if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

OneVertexAlgebra := function(field)
    return PathAlgebra(field, Quiver(["v"], []));
end;

IsolatedAlgebra := function(field, count)
    return PathAlgebra(field, Quiver(
        List([1..count], i -> Concatenation("v", String(i))), []));
end;

PrintCaughtModule := function(caught)
    Print("{\"success\":", caught[1]);
    if caught[1] and caught[2] <> fail then
        Print(",\"is_fail\":false,\"dimension_vector\":",
              DimensionVector(caught[2]));
    elif caught[1] then
        Print(",\"is_fail\":true,\"dimension_vector\":null");
    else
        Print(",\"is_fail\":false,\"dimension_vector\":null");
    fi;
    Print("}");
end;

Reset(GlobalMersenneTwister, 1);;
binary := CALL_WITH_CATCH(RandomModule, [OneVertexAlgebra(GF(2)), 1, 1]);;
Reset(GlobalMersenneTwister, 1);;
ternary := CALL_WITH_CATCH(RandomModule, [OneVertexAlgebra(GF(3)), 1, 1]);;
Reset(GlobalMersenneTwister, 1);;
rational := CALL_WITH_CATCH(RandomModule, [OneVertexAlgebra(Rationals), 1, 1]);;
Reset(GlobalMersenneTwister, 1);;
isolated := CALL_WITH_CATCH(RandomModule, [IsolatedAlgebra(GF(2), 100), 1, 1]);;

Print("{\"schema_version\":1,\"source\":{",
      "\"gap\":\"4.17dev@4efe877859ab0bc0ac041bf635a996fb897ee0c9\",",
      "\"qpa\":\"1.37@8d5f847258d4853700ce3a5d947db2a1441109f1\"},",
      "\"binary\":"); PrintCaughtModule(binary);
Print(",\"ternary\":"); PrintCaughtModule(ternary);
Print(",\"rational\":"); PrintCaughtModule(rational);
Print(",\"isolated\":"); PrintCaughtModule(isolated);
Print("}\n");

QUIT_GAP(0);
