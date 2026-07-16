if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

Q := Quiver(1, [[1, 1, "x"]]);;
KQ := PathAlgebra(Rationals, Q, HighLevelGroebnerBasis);;
LQ := PathAlgebra(Rationals, Q, HighLevelGroebnerBasis);;
x := ElementOfPathAlgebra(LQ, Q.x);;

HighLevelGroebnerBasis([x], KQ);

QUIT_GAP(0);
