if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

SizeScreen([100000, 100000]);

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintPosition := function(value)
    if value = infinity then
        PrintJsonString("infinity");
    elif value = -infinity then
        PrintJsonString("-infinity");
    elif value = "none" then
        PrintJsonString("none");
    else
        Print(value);
    fi;
end;

ListAt := function(list, positions)
    return List(positions, position -> list^position);
end;

positiveRepeat := MakeHalfInfList(2, 1, ["repeat", [10, 20, 30]], false, false);;
negativeRepeat := MakeHalfInfList(-1, -1, ["repeat", [7, 8]], false, false);;

nextCallbacks := [];;
nextRepeatData := fail;;
nextCallback := function(position, direction, listType)
    Add(nextCallbacks, [position, direction, listType]);
end;;
nextRepeatCallback := function(tail, repeating)
    nextRepeatData := [tail, StartPosition(repeating), Direction(repeating),
                       RepeatingList(repeating)];
end;;
nextList := MakeHalfInfList(0, 1,
    ["next/repeat", value -> (value + 1) mod 4, -1],
    nextCallback, nextRepeatCallback);;
nextValue := nextList^10;;

positionCallbacks := [];;
positionCallback := function(position, direction, listType)
    Add(positionCallbacks, [position, direction, listType]);
end;;
positionList := MakeHalfInfList(-2, 1, ["pos", value -> value^2, true],
                                positionCallback, false);;
positionLowestBefore := LowestKnownPosition(positionList);;
positionValue := positionList^2;;
positionValueAgain := positionList^0;;

nonstoringCallbacks := [];;
nonstoringCallback := function(position, direction, listType)
    Add(nonstoringCallbacks, [position, direction, listType]);
end;;
nonstoring := MakeHalfInfList(5, -1, ["pos", value -> value + 100, false],
                              nonstoringCallback, false);;
nonstoringValue := nonstoring^1;;

full := MakeInfList(0, [100, 200], ["repeat", [1, 2, 3]],
                    ["repeat", [7, 8]], false);;
shifted := Shift(full, 3);;
mapped := InfList(full, value -> 2 * value);;
positiveCut := Cut(PositivePart(full), 3);;
negativeCut := Cut(NegativePart(full), -2);;
positiveFrom := PositivePartFrom(full, 3);;
negativeFrom := NegativePartFrom(full, -2);;

right := MakeInfList(0, [20, 21], ["repeat", [22, 23]],
                     ["repeat", [18, 19]], false);;
spliced := Splice(full, right, 0);;
concatenated := InfConcatenation(full, right);;

finite := FiniteInfList(-1, [4, 5, 6]);;
constant := ConstantInfList(9);;
halfMapped := HalfInfList(positiveRepeat, value -> value + 1);;

tailNext := function(value)
    if value = -1 then
        return 9;
    elif value = 9 then
        return 5;
    elif value = 5 then
        return 6;
    elif value = 6 then
        return 7;
    else
        return 5;
    fi;
end;;
positiveParent := MakeInfList(0, [99], ["next/repeat", tailNext, -1],
                              ["repeat", [8]], false);;
positiveParentValue := positiveParent^10;;
negativeParent := MakeInfList(0, [99], ["repeat", [8]],
                              ["next/repeat", tailNext, -1], false);;
negativeParentValue := negativeParent^(-10);;

Print("{\"schema_version\":1,\"kind\":\"inflist\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);

Print("},\"half_lists\":{");
Print("\"positive_repeat\":{\"start\":", StartPosition(positiveRepeat),
      ",\"direction\":", Direction(positiveRepeat), ",\"type\":");
PrintJsonString(InfListType(positiveRepeat));
Print(",\"repeating\":", RepeatingList(positiveRepeat),
      ",\"storing\":", IsStoringValues(positiveRepeat), ",\"initial_is_fail\":",
      InitialValue(positiveRepeat) = fail, ",\"values\":",
      ListAt(positiveRepeat, [2..8]), ",\"lowest\":");
PrintPosition(LowestKnownPosition(positiveRepeat)); Print(",\"highest\":");
PrintPosition(HighestKnownPosition(positiveRepeat)); Print("},");
Print("\"negative_repeat\":{\"values\":", ListAt(negativeRepeat, [-6..-1]),
      ",\"lowest\":"); PrintPosition(LowestKnownPosition(negativeRepeat));
Print(",\"highest\":"); PrintPosition(HighestKnownPosition(negativeRepeat)); Print("},");
Print("\"next_repeat\":{\"value_10\":", nextValue,
      ",\"stored\":", nextList!.values, ",\"callbacks\":", nextCallbacks,
      ",\"repeat_data\":", nextRepeatData, "},");
Print("\"position_storing\":{\"lowest_before\":");
PrintPosition(positionLowestBefore);
Print(",\"value_2\":", positionValue, ",\"value_0_again\":", positionValueAgain,
      ",\"stored\":", positionList!.values, ",\"callbacks\":", positionCallbacks,
      ",\"lowest_after\":"); PrintPosition(LowestKnownPosition(positionList));
Print(",\"highest_after\":"); PrintPosition(HighestKnownPosition(positionList)); Print("},");
Print("\"position_nonstoring\":{\"value_1\":", nonstoringValue,
      ",\"stored\":", nonstoring!.values, ",\"callbacks\":", nonstoringCallbacks,
      ",\"lowest\":"); PrintPosition(LowestKnownPosition(nonstoring)); Print("}},");

Print("\"full\":{\"middle_start\":", MiddleStart(full),
      ",\"middle_end\":", MiddleEnd(full), ",\"middle\":", MiddlePart(full),
      ",\"positive_start\":", StartPosition(PositivePart(full)),
      ",\"negative_start\":", StartPosition(NegativePart(full)),
      ",\"values\":", ListAt(full, [-5,-2,-1,0,1,2,5]),
      ",\"lowest_known\":"); PrintPosition(LowestKnownPosition(full));
Print(",\"highest_known\":"); PrintPosition(HighestKnownPosition(full));
Print(",\"lower_bound\":"); PrintPosition(LowerBound(full));
Print(",\"upper_bound\":"); PrintPosition(UpperBound(full)); Print("},");

Print("\"transformations\":{");
Print("\"shift\":{\"middle_start\":", MiddleStart(shifted),
      ",\"values\":", ListAt(shifted, [-5,-3,-2,-1,2]), "},");
Print("\"map\":", ListAt(mapped, [-2,-1,0,1,2]), ",");
Print("\"positive_cut\":{\"middle_start\":", MiddleStart(positiveCut),
      ",\"middle\":", MiddlePart(positiveCut), ",\"positive_start\":",
      StartPosition(PositivePart(positiveCut)), ",\"values\":",
      ListAt(positiveCut, [3..8]), "},");
Print("\"negative_cut\":{\"middle_start\":", MiddleStart(negativeCut),
      ",\"middle\":", MiddlePart(negativeCut), ",\"negative_start\":",
      StartPosition(NegativePart(negativeCut)), ",\"values\":",
      ListAt(negativeCut, [-6..-2]), "},");
Print("\"positive_from_values\":", ListAt(positiveFrom, [3..8]), ",");
Print("\"negative_from_values\":", ListAt(negativeFrom, [-6..-2]), ",");
Print("\"splice\":", ListAt(spliced, [-2..3]), ",");
Print("\"concatenation\":{\"middle_start\":", MiddleStart(concatenated),
      ",\"middle\":", MiddlePart(concatenated), ",\"positive_start\":",
      StartPosition(PositivePart(concatenated)), ",\"negative_start\":",
      StartPosition(NegativePart(concatenated)), ",\"values\":",
      ListAt(concatenated, [-4..3]), "},");
Print("\"half_map\":", ListAt(halfMapped, [2..7]), "},");

Print("\"constructors\":{");
Print("\"finite\":{\"middle_start\":", MiddleStart(finite),
      ",\"middle_end\":", MiddleEnd(finite), ",\"values\":",
      ListAt(finite, [-1..1]), ",\"lower_bound\":"); PrintPosition(LowerBound(finite));
Print(",\"upper_bound\":"); PrintPosition(UpperBound(finite)); Print("},");
Print("\"constant\":", ListAt(constant, [-3,0,4]),
      ",\"integers\":", ListAt(IntegersList, [-3,0,4]), "},");

Print("\"repeat_conversion\":{");
Print("\"positive_parent\":{\"value_10\":", positiveParentValue,
      ",\"middle_start\":", MiddleStart(positiveParent),
      ",\"middle\":", MiddlePart(positiveParent), ",\"positive_start\":",
      StartPosition(PositivePart(positiveParent)), ",\"period\":",
      RepeatingList(PositivePart(positiveParent)), "},");
Print("\"negative_parent\":{\"value_minus_10\":", negativeParentValue,
      ",\"middle_start\":", MiddleStart(negativeParent),
      ",\"middle\":", MiddlePart(negativeParent), ",\"negative_start\":",
      StartPosition(NegativePart(negativeParent)), ",\"period\":",
      RepeatingList(NegativePart(negativeParent)), "}}}", "\n");

QUIT_GAP(0);
