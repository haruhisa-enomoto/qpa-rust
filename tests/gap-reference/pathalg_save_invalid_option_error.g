if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

file := Filename(DirectoryTemporary(), "qpa-rust-pathalg-save-invalid.txt");;
A := PathAlgebra(Rationals, Quiver(1, []));;
SaveAlgebra(A, file, "invalid");

QUIT_GAP(0);
