if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

MakeBasicListOfModules([]);

QUIT_GAP(0);
