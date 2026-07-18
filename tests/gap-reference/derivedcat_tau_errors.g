if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

algebra := PathAlgebra(Rationals,
    Quiver(["1", "2"], [["1", "2", "a"]]));
cat := CatOfRightAlgebraModules(algebra);
projective := IndecProjectiveModules(algebra)[1];
identity := IdentityMapping(projective);
finite := FiniteComplex(cat, 1, [identity]);

loopParent := PathAlgebra(Rationals,
    Quiver(["w"], [["w", "w", "x"]]));
dual := loopParent / [loopParent.x^2];
dualCat := CatOfRightAlgebraModules(dual);
dualProjective := IndecProjectiveModules(dual)[1];
generator := MinimalGeneratingSetOfModule(dualProjective)[1];
xmap := HomFromProjective(generator^dual.x, dualProjective);
repeating := Complex(dualCat, 0, [xmap], ["repeat", [xmap]], "zero");

success := rec();
success.dimension_vector_as_vertices := CALL_WITH_CATCH(
    ProjectiveToInjectiveFiniteComplex, [finite])[1];
success.repeating_highest_known := CALL_WITH_CATCH(
    ProjectiveToInjectiveComplex, [repeating])[1];
success.empty_star_description := CALL_WITH_CATCH(
    StarOfMapBetweenProjectives, [identity, [], [1]])[1];

Print("{\"schema_version\":1,\"kind\":\"derivedcat_tau_errors\",");
Print("\"success\":{");
Print("\"dimension_vector_as_vertices\":", success.dimension_vector_as_vertices, ",");
Print("\"repeating_highest_known\":", success.repeating_highest_known, ",");
Print("\"empty_star_description\":", success.empty_star_description);
Print("}}\n");

QUIT_GAP(0);
