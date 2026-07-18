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

A3 := function(field)
    return PathAlgebra(field, Quiver(3, [[1,2,"a"],[2,3,"b"]]));
end;

RadicalSquareZeroA3 := function(field)
    local algebra;
    algebra := A3(field);
    return algebra/[algebra.a*algebra.b];
end;

PrintEulerData := function(algebra)
    local cartan, form;
    cartan := CartanMatrix(algebra);
    form := EulerBilinearFormOfAlgebra(algebra);
    Print("{\"inverse_cartan\":", Inverse(cartan));
    Print(",\"value\":", form([1,2,-1],[3,0,2]), "}");
end;

Print("{\"schema_version\":1,\"kind\":\"combinatorialrep_euler\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"forms\":{\"path_a3\":"); PrintEulerData(A3(Rationals));
Print(",\"radical_square_zero_a3\":");
PrintEulerData(RadicalSquareZeroA3(Rationals));
Print(",\"path_a3_f2\":"); PrintEulerData(A3(GF(2)));
Print("}}\n");

QUIT_GAP(0);
