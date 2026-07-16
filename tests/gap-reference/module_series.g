if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintCase := function(field)
    local quiver, algebra, projectives, path_simple, zero_module, semisimple,
          loop_quiver, parent, x, quotient, quotient_projective,
          quotient_simple, quotient_zero;

    quiver := Quiver(["v1", "v2", "v3"],
        [["v1", "v2", "a"], ["v2", "v3", "b"]]);
    algebra := PathAlgebra(field, quiver);
    projectives := IndecProjectiveModules(algebra);
    path_simple := SimpleModules(algebra)[1];
    zero_module := ZeroModule(algebra);
    semisimple := RightModuleOverPathAlgebra(algebra, [1, 1, 1], []);

    loop_quiver := Quiver(["w"], [["w", "w", "x"]]);
    parent := PathAlgebra(field, loop_quiver);
    x := parent.x;
    quotient := parent / [x * x * x];
    quotient_projective := IndecProjectiveModules(quotient)[1];
    quotient_simple := SimpleModules(quotient)[1];
    quotient_zero := ZeroModule(quotient);

    Print("{\"path\":{");
    Print("\"dimensions\":", List(Concatenation(projectives,
        [zero_module, semisimple]), Dimension));
    Print(",\"partial_order\":[",
        DimensionVectorPartialOrder(projectives[3], projectives[1]), ",",
        DimensionVectorPartialOrder(projectives[1], projectives[2]), ",",
        DimensionVectorPartialOrder(projectives[1], projectives[1]), ",",
        DimensionVectorPartialOrder(projectives[2], projectives[1]), "]");
    Print(",\"radical_series\":[", RadicalSeries(projectives[1]), ",",
        RadicalSeries(projectives[2]), ",", RadicalSeries(projectives[3]),
        ",", RadicalSeries(semisimple), "]");
    Print(",\"socle_series\":[", SocleSeries(projectives[1]), ",",
        SocleSeries(projectives[2]), ",", SocleSeries(projectives[3]),
        ",", SocleSeries(semisimple), "]");
    Print(",\"zero_radical_series\":", RadicalSeries(zero_module));
    Print(",\"zero_socle_series\":", SocleSeries(zero_module));
    Print("},\"quotient\":{");
    Print("\"dimensions\":[", Dimension(quotient_projective), ",",
        Dimension(quotient_simple), ",", Dimension(quotient_zero), "]");
    Print(",\"partial_order\":[",
        DimensionVectorPartialOrder(quotient_simple, quotient_projective), ",",
        DimensionVectorPartialOrder(quotient_projective, quotient_simple), "]");
    Print(",\"cross_length_partial_order\":[",
        DimensionVectorPartialOrder(quotient_simple, path_simple), ",",
        DimensionVectorPartialOrder(path_simple, quotient_simple), ",",
        DimensionVectorPartialOrder(quotient_simple, projectives[1]), ",",
        DimensionVectorPartialOrder(projectives[1], quotient_simple), ",",
        DimensionVectorPartialOrder(quotient_simple, projectives[2]), ",",
        DimensionVectorPartialOrder(projectives[2], quotient_simple), "]");
    Print(",\"radical_projective\":", RadicalSeries(quotient_projective));
    Print(",\"socle_projective\":", SocleSeries(quotient_projective));
    Print(",\"radical_simple\":", RadicalSeries(quotient_simple));
    Print(",\"socle_simple\":", SocleSeries(quotient_simple));
    Print(",\"zero_radical_series\":", RadicalSeries(quotient_zero));
    Print(",\"zero_socle_series\":", SocleSeries(quotient_zero));
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
