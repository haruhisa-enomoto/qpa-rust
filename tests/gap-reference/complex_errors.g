if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

algebra := PathAlgebra(Rationals, Quiver(["v"], []));;
cat := CatOfRightAlgebraModules(algebra);;
zero := ZeroModule(algebra);;
one := RightModuleOverPathAlgebra(algebra, [1], []);;
two := RightModuleOverPathAlgebra(algebra, [2], []);;
f := RightModuleHomOverAlgebra(one, two, [[[1, 0]]]);;
g := RightModuleHomOverAlgebra(two, one, [[[0], [1]]]);;
ses := ShortExactSequence(cat, f, g);;
zeroEndo := ZeroMapping(one, one);;
periodic := Complex(cat, 0, [zeroEndo],
                    ["repeat", [zeroEndo]], ["repeat", [zeroEndo]]);;
good := GoodTruncation(ses, 0, 0);;

emptyFinite := CALL_WITH_CATCH(FiniteComplex, [cat, 0, []]);;
wrongObjects := CALL_WITH_CATCH(FiniteComplex,
    [cat, 0, [ZeroMapping(one, one), ZeroMapping(one, two)]]);;
nonzeroComposition := CALL_WITH_CATCH(FiniteComplex,
    [cat, 0, [IdentityMapping(one), IdentityMapping(one)]]);;
emptyRepeat := CALL_WITH_CATCH(Complex,
    [cat, 0, [], ["repeat", []], ["repeat", [zeroEndo]]]);;
notShortExact := CALL_WITH_CATCH(ShortExactSequence,
    [cat, zeroEndo, zeroEndo]);;
reversedGood := CALL_WITH_CATCH(GoodTruncation, [ses, -1, 1]);;
reversedBrutal := CALL_WITH_CATCH(BrutalTruncation, [ses, -1, 1]);;
reversedBoth := CALL_WITH_CATCH(SyzygyCosyzygyTruncation, [ses, -1, 1]);;
propertyFail := CALL_WITH_CATCH(IsShortExactSequence, [good]);;
yonedaFirstUnbounded := CALL_WITH_CATCH(YonedaProduct, [periodic, ses]);;
yonedaSecondUnbounded := CALL_WITH_CATCH(YonedaProduct, [ses, periodic]);;
twoStalk := StalkComplex(cat, two, 0);;
yonedaIncompatible := CALL_WITH_CATCH(YonedaProduct, [ses, twoStalk]);;

algebra3 := PathAlgebra(GF(3), Quiver(["v"], []));;
cat3 := CatOfRightAlgebraModules(algebra3);;
one3 := RightModuleOverPathAlgebra(algebra3, [1], []);;
stalk3 := StalkComplex(cat3, one3, 0);;
differentCategories := CALL_WITH_CATCH(ZeroChainMap, [ses, stalk3]);;
wrongChainSource := CALL_WITH_CATCH(FiniteChainMap,
    [ses, ses, 0, [IdentityMapping(one)]]);;
swap := RightModuleHomOverAlgebra(two, two, [[[0, 1], [1, 0]]]);;
noncommutingChain := CALL_WITH_CATCH(FiniteChainMap,
    [ses, ses, -1,
     [IdentityMapping(one), swap, IdentityMapping(one)]]);;
emptyCombined := CALL_WITH_CATCH(ComplexAndChainMaps,
    [[], [], 0, [], "zero", "zero"]);;

projective := ProjectiveResolution(one);;
wrongBoundComparison := CALL_WITH_CATCH(ComparisonLifting,
    [IdentityMapping(one), projective, StalkComplex(cat, one, 0)]);;
wrongObjectComparison := CALL_WITH_CATCH(ComparisonLifting,
    [IdentityMapping(two), projective, ProjectiveResolution(one)]);;

Print("{\"schema_version\":1,\"kind\":\"complex_errors\",\"success\":{");
Print("\"empty_finite\":", emptyFinite[1], ",");
Print("\"wrong_objects\":", wrongObjects[1], ",");
Print("\"nonzero_composition\":", nonzeroComposition[1], ",");
Print("\"empty_repeat\":", emptyRepeat[1], ",");
Print("\"not_short_exact\":", notShortExact[1], ",");
Print("\"reversed_good\":", reversedGood[1], ",");
Print("\"reversed_brutal\":", reversedBrutal[1], ",");
Print("\"reversed_both\":", reversedBoth[1], ",");
Print("\"property_fail\":", propertyFail[1], ",");
Print("\"yoneda_first_unbounded\":", yonedaFirstUnbounded[1], ",");
Print("\"yoneda_second_unbounded\":", yonedaSecondUnbounded[1], ",");
Print("\"yoneda_incompatible\":", yonedaIncompatible[1], ",");
Print("\"different_categories\":", differentCategories[1], ",");
Print("\"wrong_chain_source\":", wrongChainSource[1], ",");
Print("\"noncommuting_chain\":", noncommutingChain[1], ",");
Print("\"empty_combined\":", emptyCombined[1], ",");
Print("\"wrong_bound_comparison\":", wrongBoundComparison[1], ",");
Print("\"wrong_object_comparison\":", wrongObjectComparison[1], "}}\n");

QUIT_GAP(0);
