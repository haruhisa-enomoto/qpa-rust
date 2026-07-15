if LoadPackage("qpa") <> true then
    Error("failed to load the pinned QPA package");
fi;

qpa_info := PackageInfo("qpa")[1];
gbnp_info := PackageInfo("gbnp")[1];

if qpa_info.Version <> "1.37" then
    Error(Concatenation("unexpected QPA version: ", qpa_info.Version));
fi;

if gbnp_info.Version <> "1.1.0" then
    Error(Concatenation("unexpected GBNP version: ", gbnp_info.Version));
fi;

if GAPInfo.Version <> "4.17dev" then
    Error(Concatenation("unexpected GAP version: ", GAPInfo.Version));
fi;

if PositionSublist(qpa_info.InstallationPath, "qpa-rust-gap-root") = fail then
    Error(Concatenation("QPA was not loaded through the isolated oracle root: ",
                        qpa_info.InstallationPath));
fi;

if PositionSublist(gbnp_info.InstallationPath, "qpa-rust-gap-root") = fail then
    Error(Concatenation("GBNP was not loaded through the isolated oracle root: ",
                        gbnp_info.InstallationPath));
fi;

Print("{\"schema_version\":1,\"kind\":\"oracle_metadata\",",
      "\"gap_version\":\"", GAPInfo.Version, "\",",
      "\"gap_revision\":\"", GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION, "\",",
      "\"qpa_version\":\"", qpa_info.Version, "\",",
      "\"qpa_revision\":\"", GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION, "\",",
      "\"gbnp_version\":\"", gbnp_info.Version, "\",",
      "\"gbnp_revision\":\"", GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION, "\"}\n");

QUIT_GAP(0);
