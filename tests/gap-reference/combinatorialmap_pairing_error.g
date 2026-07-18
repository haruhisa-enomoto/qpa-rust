if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

CombinatorialMap(3, (), (), []);
QUIT_GAP(0);
