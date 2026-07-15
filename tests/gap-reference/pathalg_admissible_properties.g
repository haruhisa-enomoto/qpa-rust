if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

Qchain := Quiver(3, [[1, 2, "a"], [2, 3, "b"]]);;
Qloop := Quiver(1, [[1, 1, "x"]]);;

Cparent := PathAlgebra(Rationals, Qchain);;
C := Cparent/[Cparent.a*Cparent.b];;
L := TruncatedPathAlgebra(Rationals, Qloop, 3);;
Nparent := PathAlgebra(Rationals, Qloop);;
N := Nparent/[Nparent.x];;
nonadmissible_cached_before := HasIsAdmissibleQuotientOfPathAlgebra(N);;
nonadmissible := IsAdmissibleQuotientOfPathAlgebra(N);;
nonadmissible_cached_after := HasIsAdmissibleQuotientOfPathAlgebra(N);;

Cparent2 := PathAlgebra(GF(2), Qchain);;
C2 := Cparent2/[Cparent2.a*Cparent2.b];;
L2 := TruncatedPathAlgebra(GF(2), Qloop, 3);;

Print("{\"schema_version\":1,\"kind\":\"pathalg_admissible_properties\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},\"rational\":{");
Print("\"chain\":{\"admissible\":", IsAdmissibleQuotientOfPathAlgebra(C));
Print(",\"cached\":", HasIsAdmissibleQuotientOfPathAlgebra(C));
Print(",\"semisimple\":", IsSemisimpleAlgebra(C));
Print(",\"elementary\":", IsElementaryAlgebra(C));
Print(",\"basic\":", IsBasicAlgebra(C), "}");
Print(",\"truncated_loop\":{\"admissible\":", IsAdmissibleQuotientOfPathAlgebra(L));
Print(",\"semisimple\":", IsSemisimpleAlgebra(L));
Print(",\"elementary\":", IsElementaryAlgebra(L));
Print(",\"basic\":", IsBasicAlgebra(L), "}");
Print(",\"nonadmissible\":{\"value\":", nonadmissible);
Print(",\"cached_before\":", nonadmissible_cached_before);
Print(",\"cached_after\":", nonadmissible_cached_after, "}");
Print("},\"f2\":{");
Print("\"chain\":{\"admissible\":", IsAdmissibleQuotientOfPathAlgebra(C2));
Print(",\"cached\":", HasIsAdmissibleQuotientOfPathAlgebra(C2));
Print(",\"semisimple\":", IsSemisimpleAlgebra(C2));
Print(",\"elementary\":", IsElementaryAlgebra(C2));
Print(",\"basic\":", IsBasicAlgebra(C2), "}");
Print(",\"truncated_loop\":{\"admissible\":", IsAdmissibleQuotientOfPathAlgebra(L2));
Print(",\"semisimple\":", IsSemisimpleAlgebra(L2));
Print(",\"elementary\":", IsElementaryAlgebra(L2));
Print(",\"basic\":", IsBasicAlgebra(L2), "}");
Print("}}");
Print("\n");

QUIT_GAP(0);
