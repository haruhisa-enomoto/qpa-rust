if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

file := Filename(DirectoryTemporary(), "qpa-rust-pathalg-read-prefix-empty.txt");;
PrintTo(file,
    "IsQuotientOfPathAlgebra\n",
    "Vertices accepted by prefix\n",
    "v1\n",
    "Arrows accepted by prefix\n",
    "Field accepted by prefix\n",
    "Rationals\n",
    "Relations accepted by prefix\n");;
A := ReadAlgebra(file);;
Print("PREFIX_EMPTY_IS_QUOTIENT=", IsQuotientOfPathAlgebra(A), "\n");
Print("PREFIX_EMPTY_RELATIONS=", Length(RelationsOfAlgebra(A)), "\n");
RemoveFile(file);;

QUIT_GAP(0);
