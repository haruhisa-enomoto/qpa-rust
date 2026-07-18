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

a2 := UnitForm([[2,-1],[-1,2]]);
a3 := UnitForm([[2,-1,0],[-1,2,-1],[0,-1,2]]);
affine_a1 := UnitForm([[2,-2],[-2,2]]);
indefinite := UnitForm([[2,-3],[-3,2]]);
nonnegative_failure := UnitForm([[2,-2,-1],[-2,2,0],[-1,0,2]]);
positive_off_diagonal := UnitForm([[2,1],[1,2]]);

Print("{\"schema_version\":1,\"kind\":\"combinatorialrep_unit_form\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"a2\":{");
Print("\"matrix\":", SymmetricMatrixOfUnitForm(a2));
Print(",\"bilinear\":", BilinearFormOfUnitForm(a2)([1,2],[-1,3]));
Print(",\"quadratic\":", QuadraticFormOfUnitForm(a2)([1,2]));
Print(",\"reflection\":", ReflectionByBilinearForm(a2,2,[1,2]));
Print(",\"weakly_nonnegative\":", IsWeaklyNonnegativeUnitForm(a2));
Print(",\"weakly_positive\":", IsWeaklyPositiveUnitForm(a2));
Print(",\"positive_roots\":", PositiveRootsOfUnitForm(a2));
Print("},\"a3\":{");
Print("\"weakly_nonnegative\":", IsWeaklyNonnegativeUnitForm(a3));
Print(",\"weakly_positive\":", IsWeaklyPositiveUnitForm(a3));
Print(",\"positive_roots\":", PositiveRootsOfUnitForm(a3));
Print("},\"affine_a1\":{");
Print("\"weakly_nonnegative\":", IsWeaklyNonnegativeUnitForm(affine_a1));
Print(",\"weakly_positive\":", IsWeaklyPositiveUnitForm(affine_a1));
Print("},\"indefinite\":{");
Print("\"weakly_nonnegative\":", IsWeaklyNonnegativeUnitForm(indefinite));
Print(",\"weakly_positive\":", IsWeaklyPositiveUnitForm(indefinite));
Print("},\"nonnegative_failure\":{");
Print("\"weakly_nonnegative\":", IsWeaklyNonnegativeUnitForm(nonnegative_failure));
Print("},\"positive_off_diagonal\":{");
Print("\"weakly_positive\":", IsWeaklyPositiveUnitForm(positive_off_diagonal));
Print(",\"positive_roots\":", PositiveRootsOfUnitForm(positive_off_diagonal));
Print("}}\n");

QUIT_GAP(0);
