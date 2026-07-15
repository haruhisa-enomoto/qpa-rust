if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

file := Filename(DirectoryTemporary(), "qpa-rust-pathalg-save-query-new.txt");;
if IsExistingFile(file) then
    RemoveFile(file);
fi;
A := PathAlgebra(Rationals, Quiver(1, []));;
SaveAlgebra(A, file, "query");

QUIT_GAP(0);
