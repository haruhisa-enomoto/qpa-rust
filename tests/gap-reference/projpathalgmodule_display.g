if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

SizeScreen([100000, 100000]);

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    escaped := ReplacedString(escaped, "\n", "\\n");
    Print("\"", escaped, "\"");
end;

PrintObjectString := function(value)
    local result, stream;
    result := "";
    stream := OutputTextString(result, true);
    SetPrintFormattingStatus(stream, false);
    PrintTo(stream, value);
    CloseStream(stream);
    return result;
end;

q := Quiver(["v", "w"], [["v", "w", "a"]]);
algebra := PathAlgebra(GF(2), q);
projective := RightProjectiveModule(algebra, [algebra.v]);
vector := ExtRepOfObj(Vectorize(projective, [algebra.v + algebra.a]));
submodule := SubAlgebraModule(projective, [Vectorize(projective, [algebra.a])]);
gb := RightGroebnerBasisOfModule(submodule);
fp_module := VertexProjectivePresentation(algebra, [[algebra.a]]);
fp_vector := ExtRepOfObj(GeneratorsOfAlgebraModule(fp_module)[1]);
bare := Objectify(NewType(FamilyObj(projective),
    IsPathAlgebraModuleGroebnerBasis and
    IsPathAlgebraModuleGroebnerBasisDefaultRep),
    rec(gbasisElems := Immutable([[]]), staticDictionaries := Immutable([])));
SetUnderlyingModule(bare, projective);
SetIsLeftPathAlgebraModuleGroebnerBasis(bare, true);
two_argument := CALL_WITH_CATCH(RightGroebnerBasisOfModule,
    [algebra, [algebra.a]]);

Print("{\"schema_version\":1,\"kind\":\"projpathalgmodule_display\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"info_level\":", InfoLevel(InfoPathAlgebraModule));
Print(",\"two_argument_success\":", two_argument[1]);
Print(",\"left_properties\":[", IsLeftPathAlgebraModuleGroebnerBasis(bare), ",",
    IsPathAlgebraModuleGroebnerBasis(bare), "]");
Print(",\"vector_print\":"); PrintJsonString(PrintObjectString(vector));
Print(",\"fp_vector_print\":"); PrintJsonString(PrintObjectString(fp_vector));
Print(",\"gb_view\":\""); View(gb); Print("\"");
Print(",\"bare_gb_view\":\""); View(bare); Print("\"");
Print("}\n");

QUIT_GAP(0);
