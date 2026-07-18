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

RadicalSquareZeroA3 := function(field)
    local algebra;
    algebra := PathAlgebra(field, Quiver(3, [[1,2,"a"],[2,3,"b"]]));
    return algebra/[algebra.a*algebra.b];
end;

MixedEndpointRelations := function(field)
    local algebra;
    algebra := PathAlgebra(field, Quiver(6, [
        [1,2,"a1"],[2,3,"a2"],[4,5,"b1"],[5,6,"b2"]]));
    return algebra/[algebra.a1*algebra.a2 + algebra.b1*algebra.b2];
end;

Print("{\"schema_version\":1,\"kind\":\"combinatorialrep_tits_quotient\",\"oracle\":{");
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
Print("\"source_nakayama\":", SymmetricMatrixOfUnitForm(
    TitsUnitFormOfAlgebra(NakayamaAlgebra(Rationals,[3,2,1]))));
Print(",\"radical_square_zero_a3\":", SymmetricMatrixOfUnitForm(
    TitsUnitFormOfAlgebra(RadicalSquareZeroA3(Rationals))));
Print(",\"radical_square_zero_a3_f2\":", SymmetricMatrixOfUnitForm(
    TitsUnitFormOfAlgebra(RadicalSquareZeroA3(GF(2)))));
Print(",\"mixed_endpoint_relations\":", SymmetricMatrixOfUnitForm(
    TitsUnitFormOfAlgebra(MixedEndpointRelations(Rationals))));
Print("}}\n");

QUIT_GAP(0);
