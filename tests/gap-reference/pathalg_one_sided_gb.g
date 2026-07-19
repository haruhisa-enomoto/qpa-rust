if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

quiver := Quiver(1, [[1, 1, "a"]]);;
parent := PathAlgebra(GF(2), quiver);;
ideal := Ideal(parent, [parent.a]);;
left := CALL_WITH_CATCH(GroebnerBasisOfLeftIdeal, [ideal]);;
right := CALL_WITH_CATCH(GroebnerBasisOfRightIdeal, [ideal]);;

Print("{\"schema_version\":1,\"kind\":\"pathalg_one_sided_gb\",\"oracle\":{");
Print("\"gap_version\":\"", GAPInfo.Version, "\"");
Print(",\"gap_revision\":\"",
    GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION, "\"");
Print(",\"qpa_version\":\"", PackageInfo("qpa")[1].Version, "\"");
Print(",\"qpa_revision\":\"",
    GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION, "\"");
Print(",\"gbnp_version\":\"", PackageInfo("gbnp")[1].Version, "\"");
Print(",\"gbnp_revision\":\"",
    GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION, "\"}");
Print(",\"left_success\":", left[1]);
Print(",\"right_success\":", right[1], "}\n");

QUIT_GAP(0);
