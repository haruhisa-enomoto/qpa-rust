if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

IdentityCoordinates := function(field)
    local table, algebra, identity;
    table := EmptySCTable(2, Zero(field));
    SetEntrySCTable(table, 1, 1, [One(field), 1]);
    SetEntrySCTable(table, 1, 2, [One(field), 2]);
    algebra := AlgebraByStructureConstants(field, table);
    identity := FindMultiplicativeIdentity(algebra);
    return Coefficients(Basis(algebra), identity);
end;

NoIdentity := function(field)
    local table, algebra;
    table := EmptySCTable(1, Zero(field));
    algebra := AlgebraByStructureConstants(field, table);
    return FindMultiplicativeIdentity(algebra) = fail;
end;

Print("{\"schema_version\":1,\"kind\":\"algebra_structure_constants\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":");
PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"rational\":{\"identity\":");
Print(IdentityCoordinates(Rationals));
Print(",\"no_identity\":", NoIdentity(Rationals));
Print("},\"f2\":{\"identity\":[");
Print(JoinStringsWithSeparator(List(IdentityCoordinates(GF(2)), x -> String(IntFFE(x))), ","));
Print("],\"no_identity\":", NoIdentity(GF(2)));
Print("}}\n");

QUIT_GAP(0);
