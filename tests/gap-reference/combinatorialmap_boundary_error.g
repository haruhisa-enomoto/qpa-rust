if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

map := CombinatorialMap(3, (), (1,2), [0,4]);
BoundaryCurvesOfCombinatorialMap(map);
