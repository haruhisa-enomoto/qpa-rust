if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

q := Quiver(["v", "w"], [["v", "w", "a"]]);;
algebra := PathAlgebra(GF(2), q);;
projective := RightProjectiveModule(algebra, [algebra.v]);;
generator := Vectorize(projective, [algebra.a]);;
submodule := SubAlgebraModule(projective, [generator]);;
submodule_result := CALL_WITH_CATCH(IsVertexProjectiveModule, [submodule]);;

Print("{\"schema_version\":1,\"kind\":\"projpathalgmodule_vertex_projective_error\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"whole_projective\":", IsVertexProjectiveModule(projective));
Print(",\"submodule_success\":", submodule_result[1], "}\n");

QUIT_GAP(0);
