# WARNING: this demonstrates nontermination in pinned QPA. Run only with an
# external timeout.
if LoadPackage("qpa") <> true then Error("failed to load QPA"); fi;

map := CombinatorialMap(4, (1,2)(3,4), (1,2)(3,4), [1,3]);;
WindingNumber(map, [1,3]);
