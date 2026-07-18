if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

a2 := UnitForm([[2,-1],[-1,2]]);;
a3 := UnitForm([[2,-1,0],[-1,2,-1],[0,-1,2]]);;
View(a2); Print("\n");
View(a3); Print("\n");

QUIT_GAP(0);
