if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

primes := [2, 3, 5, 7, 11, 13, 17, 19, 31];;
Print("{\"schema_version\":1,\"kind\":\"pathalg_prime_field_logs\",\"oracle\":{");
Print("\"gap_version\":\"", GAPInfo.Version, "\",");
Print("\"gap_revision\":\"",
    GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION, "\",");
Print("\"qpa_version\":\"", PackageInfo("qpa")[1].Version, "\",");
Print("\"qpa_revision\":\"",
    GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION, "\"},");
Print("\"logs\":[");
for position in [1..Length(primes)] do
    p := primes[position];
    K := GF(p);
    root := PrimitiveRoot(K);
    logs := [];
    for residue in [1..p - 1] do
        Add(logs, LogFFE(residue * One(K), root));
    od;
    if position > 1 then
        Print(",");
    fi;
    Print("[", p, ",", logs, "]");
od;
Print("]}\n");

QUIT_GAP(0);
