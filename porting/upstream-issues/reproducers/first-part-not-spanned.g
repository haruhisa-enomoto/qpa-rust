# WARNING: this demonstrates nontermination in pinned QPA. Run only with an
# external timeout.
if LoadPackage("qpa") <> true then Error("failed to load QPA"); fi;

quiver := Quiver(1, [[1, 1, "x"], [1, 1, "y"]]);;
algebra := PathAlgebra(Rationals, quiver);;
projective := RightProjectiveModule(algebra, [algebra.v1]);;
fx := Vectorize(projective, [algebra.x]);;
fy := Vectorize(projective, [algebra.y]);;
FirstPart([fx], [], ExtRepOfObj(fy));
