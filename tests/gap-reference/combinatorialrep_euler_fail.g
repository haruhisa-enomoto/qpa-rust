if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

loop_algebra := PathAlgebra(Rationals, Quiver(1, [[1,1,"x"]]));
dual_numbers := loop_algebra/[loop_algebra.x^2];
if EulerBilinearFormOfAlgebra(dual_numbers) <> fail then
    Error("expected non-unimodular Cartan matrix to return fail");
fi;

QUIT_GAP(0);
