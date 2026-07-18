if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

SizeScreen([100000, 100000]);

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    escaped := ReplacedString(escaped, "\n", "\\n");
    Print("\"", escaped, "\"");
end;

PrintElement := function(element)
    local data, i;
    data := CoefficientsAndMagmaElements(element);
    Print("[");
    for i in [1, 3 .. Length(data) - 1] do
        if i > 1 then
            Print(",");
        fi;
        Print("{\"path\":");
        PrintJsonString(StringPrint(data[i]));
        Print(",\"coefficient\":", IntFFE(data[i + 1]), "}");
    od;
    Print("]");
end;

Q := Quiver(2, [[1, 1, "a"], [1, 2, "b"]]);;
R := PathRing(GF(2), Q);;
A := PathAlgebra(GF(2), Q);;
a := A.a;;
b := A.b;;
mixed := a + b;;
zero := Zero(A);;

I := Ideal(A, [a^2]);;
gb := GBNPGroebnerBasis([a^2], A);;
SetGroebnerBasisOfIdeal(I, GroebnerBasis(I, gb));
B := A / I;;
ba := B.a;;
bb := B.b;;
bmixed := ba + bb;;

left_list := CALL_WITH_CATCH(IsLeftUniform, [[a, b], Q.v1]);;
right_list := CALL_WITH_CATCH(IsRightUniform, [[a, b], Q.v1]);;

Print("{\"schema_version\":1,\"kind\":\"pathalg_surface\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"path_ring\":{");
Print("\"categories\":[", IsPathRing(R), ",", IsPathAlgebra(R), ",",
    IsQuiverAlgebra(R), "]");
Print(",\"separate_identity\":", IsIdenticalObj(R, A));
Print(",\"quiver_identical\":", IsIdenticalObj(QuiverOfPathRing(R), Q));
Print("},\"path_algebra\":{");
Print("\"categories\":[", IsPathRing(A), ",", IsPathAlgebra(A), ",",
    IsQuiverAlgebra(A), "]");
Print(",\"quiver_identical\":", IsIdenticalObj(QuiverOfPathAlgebra(A), Q));
Print(",\"ordering_identical\":",
    IsIdenticalObj(OrderingOfAlgebra(A), OrderingOfQuiver(Q)));
Print(",\"view\":"); PrintJsonString(StringView(A));
Print(",\"print\":"); PrintJsonString(StringPrint(A));
Print(",\"cyclic_basis_is_fail\":", CanonicalBasis(A) = fail);
Print(",\"dot_matches_embedding\":[", a = ElementOfPathAlgebra(A, Q.a), ",",
    A.v1 = ElementOfPathAlgebra(A, Q.v1), "]");
Print(",\"zero_shortcut_equal\":[",
    zero * a = zero, ",", a * zero = zero, ",",
    zero + a = a, ",", a + zero = a, "]");
Print(",\"comparison\":[", zero < a, ",", a < a^2, ",", a^2 < a, ",",
    a < mixed, ",", mixed < a, "]");
Print(",\"leading\":{");
Print("\"term\":"); PrintElement(LeadingTerm(mixed));
Print(",\"tip_equal\":", Tip(mixed) = LeadingTerm(mixed));
Print(",\"coefficient\":", IntFFE(TipCoefficient(mixed)));
Print(",\"monomial\":"); PrintJsonString(StringPrint(TipMonomial(mixed)));
Print("},\"uniformity\":{");
Print("\"zero\":[", IsLeftUniform(zero), ",", IsRightUniform(zero), ",",
    IsUniform(zero), "]");
Print(",\"mixed\":[", IsLeftUniform(mixed), ",", IsRightUniform(mixed), ",",
    IsUniform(mixed), "]");
Print(",\"list_success\":[", left_list[1], ",", right_list[1], "]}");
Print("},\"quotient\":{");
Print("\"categories\":[", IsPathRing(B), ",", IsPathAlgebra(B), ",",
    IsQuiverAlgebra(B), ",", IsQuotientOfPathAlgebra(B), ",",
    IsFullFpPathAlgebra(B), "]");
Print(",\"quiver_identical\":", IsIdenticalObj(QuiverOfPathAlgebra(B), Q));
Print(",\"ordering_identical\":",
    IsIdenticalObj(OrderingOfAlgebra(B), OrderingOfAlgebra(A)));
Print(",\"view\":"); PrintJsonString(StringView(B));
Print(",\"print\":"); PrintJsonString(StringPrint(B));
Print(",\"element_print\":"); PrintJsonString(StringPrint(bmixed));
Print(",\"dot_matches_generator\":[", ba = GeneratorsOfAlgebraWithOne(B)[3], ",",
    bb = GeneratorsOfAlgebraWithOne(B)[4], "]");
Print(",\"equality\":[", ba = ba + ba^2, ",", ba = bb, "]");
Print(",\"comparison\":[", Zero(B) < ba, ",", ba < bb, ",", bb < ba, "]");
Print(",\"leading\":{");
Print("\"term\":"); PrintElement(LeadingTerm(bmixed)![1]);
Print(",\"coefficient\":", IntFFE(LeadingCoefficient(bmixed)));
Print(",\"monomial\":"); PrintJsonString(StringPrint(LeadingMonomial(bmixed)));
Print("},\"uniformity\":[", IsLeftUniform(bmixed), ",",
    IsRightUniform(bmixed), ",", IsUniform(bmixed), "]");
Print(",\"normal_form_function_bound\":",
    HasNormalFormFunction(ElementsFamily(FamilyObj(B))));
Print(",\"arithmetic\":{");
Print("\"sum\":"); PrintElement((ba + bb)![1]);
Print(",\"difference\":"); PrintElement((ba - bb)![1]);
Print(",\"product\":"); PrintElement((ba * bb)![1]);
Print(",\"left_scalar\":"); PrintElement((Z(2)^0 * bmixed)![1]);
Print(",\"right_scalar\":"); PrintElement((bmixed * Z(2)^0)![1]);
Print(",\"left_path\":"); PrintElement((Q.v1 * bmixed)![1]);
Print(",\"right_path\":"); PrintElement((bmixed * Q.v2)![1]);
Print(",\"negative\":"); PrintElement((-bmixed)![1]);
Print(",\"one\":"); PrintElement(One(bmixed)![1]);
Print(",\"zero\":"); PrintElement(Zero(bmixed)![1]);
Print(",\"zero_power\":"); PrintElement((bmixed^0)![1]);
Print("}}");
Print("}\n");

QUIT_GAP(0);
