if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

Q := Quiver(1, [[1, 1, "x"]]);;
dict := QuiverStaticDictionary(Q, [Q.x]);;
SetInfoLevel(InfoStaticDictionary, 1);;
IsFiniteDifference(dict);;

QUIT_GAP(0);
