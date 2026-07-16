if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

Q := Quiver(1, [[1, 1, "x"]]);;
KQ := PathAlgebra(Rationals, Q, HighLevelGroebnerBasis);;

HighLevelGroebnerBasis([One(KQ)], KQ);

QUIT_GAP(0);
