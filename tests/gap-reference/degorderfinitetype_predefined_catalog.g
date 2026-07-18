if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

result := PredefARQuivers("what");;
if result <> false then
    Error("PredefARQuivers(\"what\") did not return false");
fi;

QUIT_GAP(0);
