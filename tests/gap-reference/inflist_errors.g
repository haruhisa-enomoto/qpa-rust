if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

At := function(list, position)
    return list^position;
end;;

finite := FiniteInfList(0, [1, 2]);;
full := ConstantInfList(1);;
positive := MakeHalfInfList(1, 1, ["repeat", [1]], false, false);;
badPositive := MakeHalfInfList(2, 1, ["repeat", [1]], false, false);;
negative := MakeHalfInfList(-1, -1, ["repeat", [1]], false, false);;
emptyRepeat := MakeHalfInfList(0, 1, ["repeat", []], false, false);;

shiftFinite := CALL_WITH_CATCH(Shift, [finite, 1]);;
mapFinite := CALL_WITH_CATCH(InfList, [finite, value -> 2 * value]);;
emptyKnown := CALL_WITH_CATCH(LowestKnownPosition, [FiniteInfList(0, [])]);;
shiftedSplice := CALL_WITH_CATCH(ShiftedSplice, [full, 0, full, 0, 0]);;
illegalHalfPosition := CALL_WITH_CATCH(At, [positive, 0]);;
emptyRepeatValue := CALL_WITH_CATCH(At, [emptyRepeat, 0]);;
cutOutside := CALL_WITH_CATCH(Cut, [positive, 0]);;
finiteOutside := CALL_WITH_CATCH(At, [finite, 2]);;
badPositiveStart := CALL_WITH_CATCH(MakeInfListFromHalfInfLists,
                                    [0, [1], badPositive, fail]);;
badNegativeStart := CALL_WITH_CATCH(MakeInfListFromHalfInfLists,
                                    [1, [1], fail, negative]);;

Print("{\"schema_version\":1,\"kind\":\"inflist_errors\",\"success\":{");
Print("\"shift_finite\":", shiftFinite[1], ",");
Print("\"map_finite\":", mapFinite[1], ",");
Print("\"empty_known\":", emptyKnown[1], ",");
Print("\"shifted_splice\":", shiftedSplice[1], ",");
Print("\"illegal_half_position\":", illegalHalfPosition[1], ",");
Print("\"empty_repeat_value\":", emptyRepeatValue[1], ",");
Print("\"cut_outside\":", cutOutside[1], ",");
Print("\"finite_outside\":", finiteOutside[1], ",");
Print("\"bad_positive_start\":", badPositiveStart[1], ",");
Print("\"bad_negative_start\":", badNegativeStart[1], "}}\n");

QUIT_GAP(0);
