if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

algebra := PathAlgebra(Rationals,
    Quiver(["1", "2"], [["1", "2", "a"]]));
cat := CatOfRightAlgebraModules(algebra);
projectives := IndecProjectiveModules(algebra);
sum := DirectSumOfQPAModules([projectives[1], projectives[2]]);
zero := ZeroMapping(projectives[1], projectives[1]);
periodic := Complex(cat, 0, [zero], ["repeat", [zero]], ["repeat", [zero]]);

success := rec();
success.nonfinite_description := CALL_WITH_CATCH(
    DescriptionOfFiniteProjComplex, [periodic])[1];
success.non_direct_sum_complement := CALL_WITH_CATCH(
    DirectSumMinusSummands, [projectives[1], [1]])[1];
success.out_of_range_complement := CALL_WITH_CATCH(
    DirectSumMinusSummands, [sum, [3]])[1];
success.non_direct_sum_components := CALL_WITH_CATCH(
    FindAllMapComponents, [IdentityMapping(projectives[1])])[1];

Print("{\"schema_version\":1,\"kind\":\"derivedcat_foundation_errors\",");
Print("\"success\":{");
Print("\"nonfinite_description\":", success.nonfinite_description, ",");
Print("\"non_direct_sum_complement\":", success.non_direct_sum_complement, ",");
Print("\"out_of_range_complement\":", success.out_of_range_complement, ",");
Print("\"non_direct_sum_components\":", success.non_direct_sum_components);
Print("}}\n");

QUIT_GAP(0);
