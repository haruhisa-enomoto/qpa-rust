if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

Probe := function(label, configuration)
    local result;
    Print("BEGIN_", label, "\n");
    result := BrauerConfigurationAlgebra(Rationals, configuration);
    Print("FAIL_", label, "=", result = fail, "\n");
end;

Probe("OUTER_SHAPE", []);
Probe("NO_VERTICES", [[], [], []]);
Probe("VERTEX_SHAPE", [[[["v"]]], [], []]);
Probe("MULTIPLICITY", [
    [["v", 0]], [["e", "v"]], [["e"]]
]);
Probe("DUPLICATE_VERTEX", [
    [["v", 1], ["v", 1]], [["e", "v"]], [["e"], ["e"]]
]);
Probe("EMPTY_POLYGON", [
    [["v", 1]], [["e"]], [["e"]]
]);
Probe("POLYGON_VERTEX_NAME", [
    [["v", 1]], [["v", "v"]], [["v"]]
]);
Probe("UNKNOWN_VERTEX", [
    [["v", 1]], [["e", "x"]], [["e"]]
]);
Probe("DUPLICATE_POLYGON", [
    [["v", 1]], [["e", "v"], ["e", "v"]], [["e", "e"]]
]);
Probe("ORIENTATION_COUNT", [
    [["v", 1]], [["e", "v"]], []
]);
Probe("EMPTY_ORIENTATION", [
    [["v", 1]], [["e", "v"]], [[]]
]);
Probe("UNKNOWN_POLYGON", [
    [["v", 1]], [["e", "v"]], [["x"]]
]);

partial_mutation := [
    [["a", 1], ["b", 1]],
    [["z", "a"], ["x", "b"]],
    [["z"], ["z", "x"]]
];;
Probe("CONTAINMENT", partial_mutation);
Print("MUTATED_ORIENTATIONS=", partial_mutation[3], "\n");

QUIT_GAP(0);
