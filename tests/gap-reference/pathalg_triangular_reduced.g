if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

Qone := Quiver(1, []);;
Aone := PathAlgebra(Rationals, Qone);;
Qarrow := Quiver(2, [[1, 2, "u"]]);;
Aarrow := PathAlgebra(Rationals, Qarrow);;
Qcycle := Quiver(2, [[1, 2, "a"], [2, 1, "b"]]);;
T := TruncatedPathAlgebra(Rationals, Qcycle, 2);;
P := PathAlgebra(Rationals, Qcycle);;
F := P/[P.b];;

Aone2 := PathAlgebra(GF(2), Qone);;
Aarrow2 := PathAlgebra(GF(2), Qarrow);;
T2 := TruncatedPathAlgebra(GF(2), Qcycle, 2);;
P2 := PathAlgebra(GF(2), Qcycle);;
F2 := P2/[P2.b];;

Print("{\"schema_version\":1,\"kind\":\"pathalg_triangular_reduced\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},\"rational\":{");
Print("\"one_vertex_path\":", IsTriangularReduced(Aone));
Print(",\"arrow_path\":", IsTriangularReduced(Aarrow));
Print(",\"two_cycle_truncated\":", IsTriangularReduced(T));
Print(",\"killed_reverse_arrow\":", IsTriangularReduced(F));
Print("},\"f2\":{");
Print("\"one_vertex_path\":", IsTriangularReduced(Aone2));
Print(",\"arrow_path\":", IsTriangularReduced(Aarrow2));
Print(",\"two_cycle_truncated\":", IsTriangularReduced(T2));
Print(",\"killed_reverse_arrow\":", IsTriangularReduced(F2));
Print("}}");
Print("\n");

QUIT_GAP(0);
