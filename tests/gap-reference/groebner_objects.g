if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintRationalElement := function(element)
    local data, i;
    data := CoefficientsAndMagmaElements(element);
    Print("[");
    for i in [1, 3 .. Length(data) - 1] do
        if i > 1 then Print(","); fi;
        Print("{\"path\":", ExtRepOfObj(data[i]));
        Print(",\"coefficient\":[", NumeratorRat(data[i + 1]), ",",
              DenominatorRat(data[i + 1]), "]}");
    od;
    Print("]");
end;

PrintRationalElementList := function(elements)
    local i;
    Print("[");
    for i in [1 .. Length(elements)] do
        if i > 1 then Print(","); fi;
        PrintRationalElement(elements[i]);
    od;
    Print("]");
end;

PrintPrimeElement := function(element)
    local data, i;
    data := CoefficientsAndMagmaElements(element);
    Print("[");
    for i in [1, 3 .. Length(data) - 1] do
        if i > 1 then Print(","); fi;
        Print("{\"path\":", ExtRepOfObj(data[i]));
        Print(",\"coefficient\":", IntFFE(data[i + 1]), "}");
    od;
    Print("]");
end;

PrintPrimeElementList := function(elements)
    local i;
    Print("[");
    for i in [1 .. Length(elements)] do
        if i > 1 then Print(","); fi;
        PrintPrimeElement(elements[i]);
    od;
    Print("]");
end;

PrintPathList := function(paths)
    local i;
    Print("[");
    for i in [1 .. Length(paths)] do
        if i > 1 then Print(","); fi;
        Print(ExtRepOfObj(paths[i]));
    od;
    Print("]");
end;

Q := Quiver(1, [[1, 1, "x"]]);;
KQ := PathAlgebra(Rationals, Q);;
v := ElementOfPathAlgebra(KQ, Q.v1);;
x := ElementOfPathAlgebra(KQ, Q.x);;
redundant := [x, x * x];;
I := Ideal(KQ, redundant);;
gb := GroebnerBasis(I, redundant);;
initial := ShallowCopy(AsList(gb));;
initial_view := StringView(gb);;
initial_print := StringPrint(gb);;
initial_flags := [HasIsTipReducedGroebnerBasis(gb),
                  HasIsCompletelyReducedGroebnerBasis(gb),
                  IsCompleteGroebnerBasis(gb)];;
gb := TipReduceGroebnerBasis(gb);;
tip_reduced := ShallowCopy(AsList(gb));;
tip_reduction := TipReduce(gb, x * x * x + v);;
complete_reduction := CompletelyReduce(gb, x * x * x + x + v);;
gb := CompletelyReduceGroebnerBasis(gb);;
completely_reduced := ShallowCopy(AsList(gb));;
final_flags := [IsTipReducedGroebnerBasis(gb),
                IsCompletelyReducedGroebnerBasis(gb),
                IsCompleteGroebnerBasis(gb)];;

cubic := x * x * x;;
J := Ideal(KQ, [cubic]);;
cubic_gb := GroebnerBasis(J, [cubic]);;
nontips := Nontips(cubic_gb);;
nontip_size := NontipSize(cubic_gb);;
finite_admits := AdmitsFinitelyManyNontips(cubic_gb);;
prefix_short := IsPrefixOfTipInTipIdeal(cubic_gb, x * x);;
prefix_equal := IsPrefixOfTipInTipIdeal(cubic_gb, cubic);;
prefix_long := IsPrefixOfTipInTipIdeal(cubic_gb, cubic * x);;
right_gb := RightGroebnerBasis(J);;
ideal_after_right_completely_reduced :=
    IsCompletelyReducedGroebnerBasis(GroebnerBasisOfIdeal(J));;
right_gb := TipReduceGroebnerBasis(right_gb);;
right_relations := ShallowCopy(AsList(right_gb));;
right_reduction := TipReduce(right_gb, cubic * x + x);;
right_view := StringView(right_gb);;
right_print := CALL_WITH_CATCH(StringPrint, [right_gb]);;
right_print_success := right_print[1];;

homogeneous_call := CALL_WITH_CATCH(IsHomogeneousGroebnerBasis, [gb]);;
homogeneous_success := homogeneous_call[1];;
SetIsHomogeneousGroebnerBasis(gb, true);
homogeneous_after_set :=
    [HasIsHomogeneousGroebnerBasis(gb), IsHomogeneousGroebnerBasis(gb)];;
zero_tip_reduction := CALL_WITH_CATCH(TipReduce, [gb, Zero(KQ)]);;
zero_tip_reduction_success := zero_tip_reduction[1];;
zero_prefix := CALL_WITH_CATCH(IsPrefixOfTipInTipIdeal,
    [cubic_gb, Zero(KQ)]);;
zero_prefix_success := zero_prefix[1];;
zero_prefix_value := zero_prefix[2];;
gb_iterator := Iterator(gb);;
iterator_values := [];;
while not IsDoneIterator(gb_iterator) do
    Add(iterator_values, NextIterator(gb_iterator));
od;

unreduced_setter := GroebnerBasis(I, [x * x]);;
SetEnumerator(unreduced_setter, [v]);
unreduced_enumerator := ShallowCopy(Enumerator(unreduced_setter));;
unreduced_as_list := ShallowCopy(AsList(unreduced_setter));;
unreduced_setter_flags :=
    [HasEnumerator(unreduced_setter), HasAsList(unreduced_setter)];;

enumerator_setter := CompletelyReduceGroebnerBasis(
    GroebnerBasis(I, [x * x]));;
SetEnumerator(enumerator_setter, [v]);
enumerator_setter_enumerator := ShallowCopy(Enumerator(enumerator_setter));;
enumerator_setter_as_list := ShallowCopy(AsList(enumerator_setter));;
enumerator_setter_flags :=
    [HasEnumerator(enumerator_setter), HasAsList(enumerator_setter)];;

as_list_setter := CompletelyReduceGroebnerBasis(
    GroebnerBasis(I, [x * x]));;
SetAsList(as_list_setter, [v]);
as_list_setter_enumerator := ShallowCopy(Enumerator(as_list_setter));;
as_list_setter_as_list := ShallowCopy(AsList(as_list_setter));;
as_list_setter_flags :=
    [HasEnumerator(as_list_setter), HasAsList(as_list_setter)];;

duplicate_ideal := Ideal(KQ, [x * x]);;
duplicate_gb := CompletelyReduceGroebnerBasis(
    GroebnerBasis(duplicate_ideal, [x * x]));;
SetAsList(duplicate_gb, [x * x, x * x]);
duplicate_right_relations :=
    ShallowCopy(AsList(RightGroebnerBasis(duplicate_ideal)));;

Qinf := Quiver(1, [[1, 1, "x"], [1, 1, "y"]]);;
KQinf := PathAlgebra(Rationals, Qinf);;
xinf := ElementOfPathAlgebra(KQinf, Qinf.x);;
yinf := ElementOfPathAlgebra(KQinf, Qinf.y);;
infinite_ideal := Ideal(KQinf, [xinf * xinf]);;
infinite_gb := GroebnerBasis(infinite_ideal, [xinf * xinf]);;
infinite_nontips := Nontips(infinite_gb);;
infinite_size := NontipSize(infinite_gb);;
infinite_right := CALL_WITH_CATCH(RightGroebnerBasis, [infinite_ideal]);;
infinite_right_success := infinite_right[1];;

MakeUnreducedRightBasis := function(relations)
    return Objectify(
        NewType(FamilyObj(relations),
            IsRightGroebnerBasis and IsGroebnerBasisDefaultRep),
        rec(relations := relations));
end;;
manual_right := MakeUnreducedRightBasis([xinf * xinf, yinf * yinf]);;
manual_right_result := CALL_WITH_CATCH(TipReduceGroebnerBasis, [manual_right]);;
manual_right_success := manual_right_result[1];;
manual_right_relations := ShallowCopy(AsList(manual_right_result[2]));;
manual_right_flags := [IsTipReducedGroebnerBasis(manual_right_result[2]),
                       HasIsCompleteGroebnerBasis(manual_right_result[2])];;
manual_prefix := MakeUnreducedRightBasis([xinf, xinf * xinf]);;
manual_prefix_result := CALL_WITH_CATCH(TipReduceGroebnerBasis, [manual_prefix]);;
manual_prefix_success := manual_prefix_result[1];;
manual_empty := Objectify(
    NewType(FamilyObj([xinf]),
        IsRightGroebnerBasis and IsGroebnerBasisDefaultRep),
    rec(relations := []));;
manual_empty_result := CALL_WITH_CATCH(TipReduceGroebnerBasis, [manual_empty]);;
manual_empty_success := manual_empty_result[1];;

F5 := GF(5);;
KQ5 := PathAlgebra(F5, Q);;
v5 := ElementOfPathAlgebra(KQ5, Q.v1);;
x5 := ElementOfPathAlgebra(KQ5, Q.x);;
redundant5 := [x5, x5 * x5];;
I5 := Ideal(KQ5, redundant5);;
gb5 := GroebnerBasis(I5, redundant5);;
gb5 := TipReduceGroebnerBasis(gb5);;
tip_reduction5 := TipReduce(gb5, 3 * x5 * x5 + v5);;

F2 := GF(2);;
KQ2 := PathAlgebra(F2, Q);;
v2 := ElementOfPathAlgebra(KQ2, Q.v1);;
x2 := ElementOfPathAlgebra(KQ2, Q.x);;
redundant2 := [x2, x2 * x2];;
I2 := Ideal(KQ2, redundant2);;
gb2 := GroebnerBasis(I2, redundant2);;
gb2 := TipReduceGroebnerBasis(gb2);;
tip_reduction2 := TipReduce(gb2, x2 * x2 + v2);;

Print("{\"schema_version\":1,\"kind\":\"groebner_objects\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},\"rational\":{\"initial\":");
PrintRationalElementList(initial);
Print(",\"initial_flags\":", initial_flags);
Print(",\"tip_reduced\":");
PrintRationalElementList(tip_reduced);
Print(",\"tip_reduction\":");
PrintRationalElement(tip_reduction);
Print(",\"complete_reduction\":");
PrintRationalElement(complete_reduction);
Print(",\"completely_reduced\":");
PrintRationalElementList(completely_reduced);
Print(",\"final_flags\":", final_flags);
Print(",\"nontips\":");
PrintPathList(nontips);
Print(",\"nontip_size\":", nontip_size);
Print(",\"finite_admits\":", finite_admits);
Print(",\"prefix\":[", prefix_short, ",", prefix_equal, ",", prefix_long, "]");
Print(",\"right_relations\":");
PrintRationalElementList(right_relations);
Print(",\"right_reduction\":");
PrintRationalElement(right_reduction);
Print(",\"surface\":{\"view\":");
PrintJsonString(initial_view);
Print(",\"print\":");
PrintJsonString(initial_print);
Print(",\"homogeneous_success\":", homogeneous_success);
Print(",\"homogeneous_after_set\":", homogeneous_after_set);
Print(",\"zero_tip_reduction_success\":", zero_tip_reduction_success);
Print(",\"zero_prefix\":[", zero_prefix_success, ",", zero_prefix_value, "]");
Print(",\"iterator\":");
PrintRationalElementList(iterator_values);
Print(",\"unreduced_setter\":{\"flags\":", unreduced_setter_flags,
      ",\"enumerator\":");
PrintRationalElementList(unreduced_enumerator);
Print(",\"as_list\":");
PrintRationalElementList(unreduced_as_list);
Print("},\"enumerator_setter\":{\"flags\":", enumerator_setter_flags,
      ",\"enumerator\":");
PrintRationalElementList(enumerator_setter_enumerator);
Print(",\"as_list\":");
PrintRationalElementList(enumerator_setter_as_list);
Print("},\"as_list_setter\":{\"flags\":", as_list_setter_flags,
      ",\"enumerator\":");
PrintRationalElementList(as_list_setter_enumerator);
Print(",\"as_list\":");
PrintRationalElementList(as_list_setter_as_list);
Print("},\"right_view\":");
PrintJsonString(right_view);
Print(",\"right_print_success\":", right_print_success);
Print(",\"ideal_after_right_completely_reduced\":",
      ideal_after_right_completely_reduced);
Print(",\"duplicate_right_relations\":");
PrintRationalElementList(duplicate_right_relations);
Print(",\"manual_right\":{\"success\":", manual_right_success,
      ",\"relations\":");
PrintRationalElementList(manual_right_relations);
Print(",\"flags\":", manual_right_flags,
      ",\"prefix_success\":", manual_prefix_success,
      ",\"empty_success\":", manual_empty_success, "}");
Print(",\"infinite\":{\"admits\":",
      AdmitsFinitelyManyNontips(infinite_gb),
      ",\"nontips_is_fail\":", infinite_nontips = fail,
      ",\"size_is_infinity\":", infinite_size = infinity,
      ",\"right_success\":", infinite_right_success, "}}");
Print("},\"f5\":{\"tip_reduced\":");
PrintPrimeElementList(AsList(gb5));
Print(",\"tip_reduction\":");
PrintPrimeElement(tip_reduction5);
Print("},\"f2\":{\"tip_reduced\":");
PrintPrimeElementList(AsList(gb2));
Print(",\"tip_reduction\":");
PrintPrimeElement(tip_reduction2);
Print("}}\n");

QUIT_GAP(0);
