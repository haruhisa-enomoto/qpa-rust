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
    local data, index;
    if IsElementOfQuotientOfPathAlgebra(element) then
        data := CoefficientsAndMagmaElements(element![1]);
    else
        data := CoefficientsAndMagmaElements(element);
    fi;
    Print("[");
    for index in [1, 3 .. Length(data) - 1] do
        if index > 1 then Print(","); fi;
        Print("{\"path\":", ExtRepOfObj(data[index]));
        Print(",\"coefficient\":[", NumeratorRat(data[index + 1]), ",",
              DenominatorRat(data[index + 1]), "]}");
    od;
    Print("]");
end;

PrintElementList := function(elements)
    local index;
    Print("[");
    for index in [1..Length(elements)] do
        if index > 1 then Print(","); fi;
        PrintRationalElement(elements[index]);
    od;
    Print("]");
end;

PrintPrimeElement := function(element)
    local data, index;
    if IsElementOfQuotientOfPathAlgebra(element) then
        data := CoefficientsAndMagmaElements(element![1]);
    else
        data := CoefficientsAndMagmaElements(element);
    fi;
    Print("[");
    for index in [1, 3 .. Length(data) - 1] do
        if index > 1 then Print(","); fi;
        Print("{\"path\":", ExtRepOfObj(data[index]));
        Print(",\"coefficient\":", IntFFE(data[index + 1]), "}");
    od;
    Print("]");
end;

PrintPaths := function(elements)
    local index, terms;
    Print("[");
    for index in [1..Length(elements)] do
        if index > 1 then Print(","); fi;
        terms := PathAlgebraElementTerms(elements[index]);
        if Length(terms) <> 1 or terms[1].coeff <> One(LeftActingDomain(PathAlgebraContainingElement(elements[index]))) then
            Error("expected a monic monomial element");
        fi;
        Print(ExtRepOfObj(terms[1].monom));
    od;
    Print("]");
end;

Q := Quiver(3, [[1, 1, "a"], [1, 2, "b"], [2, 3, "c"],
                 [3, 1, "d"], [2, 2, "e"]]);;
KQ := PathAlgebra(Rationals, Q, HighLevelGroebnerBasis);;
a := ElementOfPathAlgebra(KQ, Q.a);;
b := ElementOfPathAlgebra(KQ, Q.b);;
c := ElementOfPathAlgebra(KQ, Q.c);;
d := ElementOfPathAlgebra(KQ, Q.d);;
e := ElementOfPathAlgebra(KQ, Q.e);;

uniform := MakeUniformOnRight([a+b, a]);;
extended := GeneratorsTimesArrowsOnRight([a+b]);;
power0 := NthPowerOfArrowIdeal(KQ, 0);;
power1 := NthPowerOfArrowIdeal(KQ, 1);;
power2 := NthPowerOfArrowIdeal(KQ, 2);;
power3 := NthPowerOfArrowIdeal(KQ, 3);;
mutable_relations := [a+b];;
returned_relations := AddNthPowerToRelations(KQ, mutable_relations, 2);;

I := Ideal(KQ, [a^2]);;
projection := NaturalHomomorphismByIdeal(KQ, I);;
natural_range := Range(projection);;
image := ElementOfQuotientOfPathAlgebra(FamilyObj(One(natural_range)), a+a^2, false);;
lazy_sum := ElementOfQuotientOfPathAlgebra(FamilyObj(One(natural_range)), a+b, false);;
lazy_uniform := MakeUniformOnRight([lazy_sum]);;
lazy_extended := GeneratorsTimesArrowsOnRight([lazy_sum]);;

Iattached := Ideal(KQ, [a^2]);;
attached_basis := GroebnerBasis(Iattached, [a^2]);;
attached_projection := NaturalHomomorphismByIdeal(KQ, Iattached);;
attached_range := Range(attached_projection);;
attached_image := ElementOfQuotientOfPathAlgebra(
    FamilyObj(One(attached_range)), a+a^2, false);;

KQop := OppositePathAlgebra(KQ);;
opposite_element := OppositePathAlgebraElement(2*a*b+3*b*c);;
double_opposite_element := OppositePathAlgebraElement(opposite_element);;
opposite_relations := OppositeRelations([a^2+b*c, c*d]);;

truncated := TruncatedPathAlgebra(Rationals, Q, 3);;
truncated_opposite := OppositePathAlgebra(truncated);;
truncated_double_opposite := OppositePathAlgebra(truncated_opposite);;
qa := One(truncated) * Q.a;;
qb := One(truncated) * Q.b;;
qc := One(truncated) * Q.c;;
quotient_element := qa+qa^2+qb*qc;;
opposite_quotient_element := OppositePathAlgebraElement(quotient_element);;

KQ2 := PathAlgebra(GF(2), Q, HighLevelGroebnerBasis);;
a2 := ElementOfPathAlgebra(KQ2, Q.a);;
b2 := ElementOfPathAlgebra(KQ2, Q.b);;
c2 := ElementOfPathAlgebra(KQ2, Q.c);;
I2 := Ideal(KQ2, [a2^2]);;
projection2 := NaturalHomomorphismByIdeal(KQ2, I2);;
range2 := Range(projection2);;
image2 := ElementOfQuotientOfPathAlgebra(FamilyObj(One(range2)), a2+a2^2, false);;
I2attached := Ideal(KQ2, [a2^2]);;
attached_basis2 := GroebnerBasis(I2attached, [a2^2]);;
projection2attached := NaturalHomomorphismByIdeal(KQ2, I2attached);;
range2attached := Range(projection2attached);;
attached_image2 := ElementOfQuotientOfPathAlgebra(
    FamilyObj(One(range2attached)), a2+a2^2, false);;
opposite_element2 := OppositePathAlgebraElement(a2*b2+b2*c2);;

Print("{\"schema_version\":1,\"kind\":\"pathalg_foundational\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},\"right_uniform\":");
PrintElementList(uniform);
Print(",\"times_arrows\":");
PrintPaths(extended);
Print(",\"arrow_powers\":{\"n0\":");
PrintPaths(power0);
Print(",\"n1\":");
PrintPaths(power1);
Print(",\"n2\":");
PrintPaths(power2);
Print(",\"n3\":");
PrintPaths(power3);
Print("},\"add_power\":{\"same_object\":", IsIdenticalObj(mutable_relations, returned_relations));
Print(",\"relations\":");
PrintElementList(returned_relations);
Print("},\"natural_map\":{\"source_identical\":", IsIdenticalObj(Source(projection), KQ));
Print(",\"range_original_identical\":", IsIdenticalObj(OriginalPathAlgebra(natural_range), KQ));
Print(",\"surjective\":", IsSurjective(projection));
Print(",\"image\":");
PrintRationalElement(image);
Print(",\"image_is_normal_form\":", IsNormalForm(image));
Print(",\"right_uniform\":");
PrintElementList(lazy_uniform);
Print(",\"times_arrows\":");
PrintElementList(lazy_extended);
Print(",\"attached_image\":");
PrintRationalElement(attached_image);
Print(",\"attached_image_is_normal_form\":", IsNormalForm(attached_image));
Print("},\"opposite_parent\":{\"double_identical\":",
      IsIdenticalObj(OppositePathAlgebra(KQop), KQ));
Print(",\"quiver_double_identical\":",
      IsIdenticalObj(OppositeQuiver(QuiverOfPathAlgebra(KQop)), Q));
Print(",\"element\":");
PrintRationalElement(opposite_element);
Print(",\"double_element\":");
PrintRationalElement(double_opposite_element);
Print(",\"relations\":");
PrintElementList(opposite_relations);
Print("},\"truncated\":{\"dimension\":", Dimension(truncated));
Print(",\"opposite_dimension\":", Dimension(truncated_opposite));
Print(",\"double_identical\":", IsIdenticalObj(truncated_double_opposite, truncated));
Print(",\"opposite_element\":");
PrintRationalElement(opposite_quotient_element);
Print("},\"f2\":{\"lazy_image\":");
PrintPrimeElement(image2);
Print(",\"lazy_is_normal_form\":", IsNormalForm(image2));
Print(",\"attached_image\":");
PrintPrimeElement(attached_image2);
Print(",\"attached_is_normal_form\":", IsNormalForm(attached_image2));
Print(",\"opposite_element\":");
PrintPrimeElement(opposite_element2);
Print("}}");
Print("\n");

QUIT_GAP(0);
