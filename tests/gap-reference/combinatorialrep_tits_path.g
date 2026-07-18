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

a3 := PathAlgebra(Rationals, Quiver(3, [[1,2,"a"],[2,3,"b"]]));
kronecker_f2 := PathAlgebra(GF(2), Quiver(3, [[1,2,"x"],[1,2,"y"]]));
opposite_orientation := PathAlgebra(
    Rationals, Quiver(3, [[2,1,"c"],[3,2,"d"]]));

Print("{\"schema_version\":1,\"kind\":\"combinatorialrep_tits_path\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"matrices\":{");
Print("\"a3\":", SymmetricMatrixOfUnitForm(TitsUnitFormOfAlgebra(a3)));
Print(",\"kronecker_f2\":",
    SymmetricMatrixOfUnitForm(TitsUnitFormOfAlgebra(kronecker_f2)));
Print(",\"opposite_orientation\":",
    SymmetricMatrixOfUnitForm(TitsUnitFormOfAlgebra(opposite_orientation)));
Print("}}\n");

QUIT_GAP(0);
