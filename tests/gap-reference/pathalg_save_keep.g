if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

file := Filename(DirectoryTemporary(), "qpa-rust-pathalg-save-keep.txt");;
A := PathAlgebra(Rationals, Quiver(1, []));;
SaveAlgebra(A, file, "delete");;
result := SaveAlgebra(A, file, "keep existing");;
Print("KEEP_RESULT=", result, "\n");
RemoveFile(file);;

QUIT_GAP(0);
