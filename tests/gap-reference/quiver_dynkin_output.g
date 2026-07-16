if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

Probe := function(label, quiver)
    local result;
    Print("BEGIN_", label, "\n");
    result := IsDynkinQuiver(quiver);
    Print("RESULT_", label, "=", result, "\n");
end;

Probe("A4", Quiver(4,
    [[1, 2, "a"], [2, 3, "b"], [3, 4, "c"]]));
Probe("D5", Quiver(5,
    [[1, 3, "a"], [2, 3, "b"], [3, 4, "c"], [4, 5, "d"]]));
Probe("E8", Quiver(8,
    [[1, 2, "a1"],
     [1, 3, "b1"], [3, 4, "b2"],
     [1, 5, "c1"], [5, 6, "c2"], [6, 7, "c3"], [7, 8, "c4"]]));
Probe("DISCONNECTED", Quiver(2, []));
Probe("CYCLE", Quiver(2, [[1, 2, "a"], [2, 1, "b"]]));
Probe("FOUR_ARMS", Quiver(5,
    [[1, 2, "a"], [1, 3, "b"], [1, 4, "c"], [1, 5, "d"]]));
Probe("LONG_STAR", Quiver(9,
    [[1, 2, "a1"],
     [1, 3, "b1"], [3, 4, "b2"],
     [1, 5, "c1"], [5, 6, "c2"], [6, 7, "c3"],
     [7, 8, "c4"], [8, 9, "c5"]]));

QUIT_GAP(0);
