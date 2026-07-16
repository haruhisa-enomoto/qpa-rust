if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

example := Poset(["d", "b", "a", "c"],
    [["a", "b", "c"], ["b", "d"], ["c", "d"]]);;

ViewObj(example);
Print("\n");
Display(example);

QUIT_GAP(0);
