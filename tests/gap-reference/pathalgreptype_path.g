if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

Probe := function(label, algebra)
    local result;
    Print("BEGIN_", label, "\n");
    result := IsFiniteTypeAlgebra(algebra);
    Print("RESULT_", label, "=", result, "\n");
end;

Probe("SEMISIMPLE", PathAlgebra(Rationals, Quiver(3, [])));

F2 := GF(2);;
Probe("A3", PathAlgebra(F2,
    Quiver(3, [[1, 2, "a"], [2, 3, "b"]])));

Probe("DISJOINT_DYNKIN", PathAlgebra(Rationals,
    Quiver(6, [[1, 2, "a"], [3, 4, "b"], [3, 5, "c"], [3, 6, "d"]])));

Probe("NONDYNKIN_TREE", PathAlgebra(F2,
    Quiver(5, [[1, 2, "a"], [1, 3, "b"], [1, 4, "c"], [1, 5, "d"]])));

Probe("SHORT_CIRCUIT", PathAlgebra(Rationals,
    Quiver(4, [[1, 2, "a"], [3, 4, "b"], [4, 3, "c"]])));

QUIT_GAP(0);
