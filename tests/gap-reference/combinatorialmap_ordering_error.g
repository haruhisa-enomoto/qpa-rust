if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

CombinatorialMap(3, (1,4), (1,2), []);
QUIT_GAP(0);
