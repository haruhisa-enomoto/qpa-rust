if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

isolated := PathAlgebra(Rationals, Quiver(2, []));;
PreprojectiveAlgebra(isolated);

QUIT_GAP(0);
