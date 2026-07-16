if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

Q := Quiver(1, [[1, 1, "x"]]);;

Print("GENERIC_0="); ViewObj(StaticDictionary([])); Print("\n");
Print("GENERIC_1="); ViewObj(StaticDictionary([[1]])); Print("\n");
Print("GENERIC_2="); ViewObj(StaticDictionary([[1], [2]])); Print("\n");
Print("QUIVER_0="); ViewObj(QuiverStaticDictionary(Q, [])); Print("\n");
Print("QUIVER_1="); ViewObj(QuiverStaticDictionary(Q, [Q.x])); Print("\n");
Print("QUIVER_2="); ViewObj(QuiverStaticDictionary(Q, [Q.v1, Q.x])); Print("\n");

QUIT_GAP(0);
