if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

NakayamaAlgebra(Rationals, []);

QUIT_GAP(0);
