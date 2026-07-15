if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

left := PathAlgebra(Rationals, Quiver(1, []));;
right := PathAlgebra(Rationals, Quiver(1, []));;
TensorProductOfAlgebras(left, right);

QUIT_GAP(0);
