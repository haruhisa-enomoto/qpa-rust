if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

NakayamaAlgebra(Rationals, [0]);

QUIT_GAP(0);
