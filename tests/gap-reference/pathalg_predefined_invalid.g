if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

result := NakayamaAlgebra(Rationals, [1, 1]);;
Print("INVALID_RETURN=", result, "\n");

QUIT_GAP(0);
