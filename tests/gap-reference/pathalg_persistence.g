if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    escaped := ReplacedString(escaped, "\n", "\\n");
    escaped := ReplacedString(escaped, "\r", "\\r");
    escaped := ReplacedString(escaped, "\t", "\\t");
    Print("\"", escaped, "\"");
end;

directory := DirectoryTemporary();;
pathfile := Filename(directory, "qpa-rust-pathalg-save-path.txt");;
qqfile := Filename(directory, "qpa-rust-pathalg-save-qq.txt");;
f2file := Filename(directory, "qpa-rust-pathalg-save-f2.txt");;
f5file := Filename(directory, "qpa-rust-pathalg-save-f5.txt");;

P := PathAlgebra(Rationals,
    DynkinQuiver("A", 11, ["r","r","r","r","r","r","r","r","r","r"]));;
path_saved := SaveAlgebra(P, pathfile, "delete");;
path_text := StringFile(pathfile);;
Pread := ReadAlgebra(pathfile);;

Q := Quiver(3, [[1, 2, "a"], [2, 3, "b"], [1, 3, "c"]]);;
A := PathAlgebra(Rationals, Q);;
QQ := A/[A.a*A.b - (1/2)*A.c];;
SaveAlgebra(QQ, qqfile, "delete old file");;
qq_text := StringFile(qqfile);;
QQread := ReadAlgebra(qqfile);;

A2 := PathAlgebra(GF(2), Q);;
F2 := A2/[A2.a*A2.b + A2.c];;
SaveAlgebra(F2, f2file, "delete");;
f2_text := StringFile(f2file);;
F2read := ReadAlgebra(f2file);;

A5 := PathAlgebra(GF(5), Q);;
F5 := A5/[A5.a*A5.b + 2*A5.c];;
SaveAlgebra(F5, f5file, "delete");;
f5_text := StringFile(f5file);;
F5read := ReadAlgebra(f5file);;

Print("{\"schema_version\":1,\"kind\":\"pathalg_persistence\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},\"path\":{\"saved\":", path_saved);
Print(",\"text\":");
PrintJsonString(path_text);
Print(",\"read_is_path\":", IsPathAlgebra(Pread));
Print(",\"read_vertices\":", NumberOfVertices(QuiverOfPathAlgebra(Pread)));
Print(",\"read_arrows\":", NumberOfArrows(QuiverOfPathAlgebra(Pread)), "}");
Print(",\"rational\":{\"text\":");
PrintJsonString(qq_text);
Print(",\"read_is_quotient\":", IsQuotientOfPathAlgebra(QQread));
Print(",\"read_dimension\":", Dimension(QQread), "}");
Print(",\"f2\":{\"text\":");
PrintJsonString(f2_text);
Print(",\"read_is_quotient\":", IsQuotientOfPathAlgebra(F2read));
Print(",\"read_dimension\":", Dimension(F2read), "}");
Print(",\"f5\":{\"text\":");
PrintJsonString(f5_text);
Print(",\"read_is_quotient\":", IsQuotientOfPathAlgebra(F5read));
Print(",\"read_dimension\":", Dimension(F5read), "}}");
Print("\n");

RemoveFile(pathfile);;
RemoveFile(qqfile);;
RemoveFile(f2file);;
RemoveFile(f5file);;
QUIT_GAP(0);
