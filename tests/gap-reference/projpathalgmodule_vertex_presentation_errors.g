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

q := Quiver(["v", "w", "x"],
    [["v", "w", "a"], ["w", "x", "b"], ["v", "x", "c"]]);
algebra := PathAlgebra(GF(2), q);
zero := Zero(algebra);
map := [[algebra.a, algebra.w], [algebra.c, algebra.b]];

success_inferred := CALL_WITH_CATCH(VertexProjectivePresentation, [algebra, map]);
success_explicit := CALL_WITH_CATCH(VertexProjectivePresentation,
    [algebra, map, [algebra.v, algebra.w]]);
empty_inferred := CALL_WITH_CATCH(VertexProjectivePresentation, [algebra, []]);
ragged_inferred := CALL_WITH_CATCH(VertexProjectivePresentation,
    [algebra, [[algebra.a, algebra.w], [algebra.c]]]);
nonuniform_inferred := CALL_WITH_CATCH(VertexProjectivePresentation,
    [algebra, [[algebra.a + algebra.b, algebra.w]]]);
mismatch_inferred := CALL_WITH_CATCH(VertexProjectivePresentation,
    [algebra, [[algebra.a, algebra.w], [algebra.b, zero]]]);
zero_column_inferred := CALL_WITH_CATCH(VertexProjectivePresentation,
    [algebra, [[zero, algebra.w]]]);
empty_explicit := CALL_WITH_CATCH(VertexProjectivePresentation,
    [algebra, [], [algebra.v, algebra.w]]);
length_explicit := CALL_WITH_CATCH(VertexProjectivePresentation,
    [algebra, [[algebra.a]], [algebra.v, algebra.w]]);
bad_vertices_explicit := CALL_WITH_CATCH(VertexProjectivePresentation,
    [algebra, map, [algebra.a, algebra.w]]);
mismatch_explicit := CALL_WITH_CATCH(VertexProjectivePresentation,
    [algebra, map, [algebra.w, algebra.w]]);

Print("{\"schema_version\":1,\"kind\":\"projpathalgmodule_vertex_presentation_errors\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"success\":[", success_inferred[1], ",", success_explicit[1], "]");
Print(",\"failure_success\":[", empty_inferred[1], ",", ragged_inferred[1], ",",
    nonuniform_inferred[1], ",", mismatch_inferred[1], ",",
    zero_column_inferred[1], ",", empty_explicit[1], ",", length_explicit[1], ",",
    bad_vertices_explicit[1], ",", mismatch_explicit[1], "]}");
Print("\n");

QUIT_GAP(0);
