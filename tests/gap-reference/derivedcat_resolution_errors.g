if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

algebra := PathAlgebra(Rationals,
    Quiver(["1", "2"], [["1", "2", "a"]]));
cat := CatOfRightAlgebraModules(algebra);
simple := RightModuleOverPathAlgebra(algebra, [1, 0], []);
identityComplex := FiniteComplex(cat, 1, [IdentityMapping(simple)]);
zero := ZeroMapping(simple, simple);
periodic := Complex(cat, 0, [zero], ["repeat", [zero]], ["repeat", [zero]]);

success := rec();
success.identity_cancellation := CALL_WITH_CATCH(
    ProjectiveResolutionOfComplex, [identityComplex])[1];
success.nonfinite_lower_bound := CALL_WITH_CATCH(
    ProjectiveResolutionOfComplex, [periodic])[1];
success.zero_complex := CALL_WITH_CATCH(
    ProjectiveResolutionOfComplex, [ZeroComplex(cat)])[1];

Print("{\"schema_version\":1,\"kind\":\"derivedcat_resolution_errors\",");
Print("\"success\":{");
Print("\"identity_cancellation\":", success.identity_cancellation, ",");
Print("\"nonfinite_lower_bound\":", success.nonfinite_lower_bound, ",");
Print("\"zero_complex\":", success.zero_complex);
Print("}}\n");

QUIT_GAP(0);
