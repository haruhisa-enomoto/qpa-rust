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
    data := CoefficientsAndMagmaElements(element);
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
    data := CoefficientsAndMagmaElements(element);
    Print("[");
    for index in [1, 3 .. Length(data) - 1] do
        if index > 1 then Print(","); fi;
        Print("{\"path\":", ExtRepOfObj(data[index]));
        Print(",\"coefficient\":", IntFFE(data[index + 1]), "}");
    od;
    Print("]");
end;

PrintPrimeElementList := function(elements)
    local index;
    Print("[");
    for index in [1..Length(elements)] do
        if index > 1 then Print(","); fi;
        PrintPrimeElement(elements[index]);
    od;
    Print("]");
end;

Q := Quiver(3, [[1, 2, "a"], [1, 2, "b"], [2, 2, "c"],
                 [2, 3, "d"], [3, 1, "e"]]);;
KQ := PathAlgebra(Rationals, Q, HighLevelGroebnerBasis);;
a := ElementOfPathAlgebra(KQ, Q.a);;
b := ElementOfPathAlgebra(KQ, Q.b);;
c := ElementOfPathAlgebra(KQ, Q.c);;
d := ElementOfPathAlgebra(KQ, Q.d);;
e := ElementOfPathAlgebra(KQ, Q.e);;
relations := [c^2, d*e, e*a, a*c*d-b*d, e*b*d];;
I := Ideal(KQ, relations);;
A := KQ / relations;;
regression_admissible := IsAdmissibleIdeal(I);;
regression_minimal := MinimalGeneratingSetOfIdeal(I);;

Qloop := Quiver(1, [[1, 1, "x"]]);;
KQloop := PathAlgebra(Rationals, Qloop, HighLevelGroebnerBasis);;
x := ElementOfPathAlgebra(KQloop, Qloop.x);;
zero := Zero(KQloop);;
trivial_loop := Ideal(KQloop, []);;
linear := Ideal(KQloop, [x]);;
monomial := Ideal(KQloop, [x^2]);;
derived_monomial := Ideal(KQloop, [x^2+x^3, x^3]);;
right_factor := Ideal(KQloop, [x^2]);;
right_gb := GroebnerBasis(right_factor, [x^2]);;
left_factor := Ideal(KQloop, [x]);;
product := ProductOfIdeals(left_factor, right_factor);;
redundant := Ideal(KQloop, [x^2, x^3]);;
redundant_admissible := IsAdmissibleIdeal(redundant);;
minimal := MinimalGeneratingSetOfIdeal(redundant);;

Qacyclic := Quiver(2, [[1, 2, "u"]]);;
KQacyclic := PathAlgebra(Rationals, Qacyclic, HighLevelGroebnerBasis);;
trivial_acyclic := Ideal(KQacyclic, []);;
homogeneous_empty := IsHomogeneousListOfElements([], 2);;

KQloop2 := PathAlgebra(GF(2), Qloop, HighLevelGroebnerBasis);;
x2 := ElementOfPathAlgebra(KQloop2, Qloop.x);;
monomial2 := Ideal(KQloop2, [x2^2]);;
redundant2 := Ideal(KQloop2, [x2^2, x2^3]);;
redundant_admissible2 := IsAdmissibleIdeal(redundant2);;
minimal2 := MinimalGeneratingSetOfIdeal(redundant2);;

Print("{\"schema_version\":1,\"kind\":\"pathalg_ideal\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},\"regression\":{");
Print("\"admissible\":", regression_admissible);
Print(",\"monomial\":", IsMonomialIdeal(I));
Print(",\"quotient_admissible\":", IsAdmissibleQuotientOfPathAlgebra(A));
Print(",\"minimal_generators\":");
PrintElementList(regression_minimal);
Print("},\"boundary\":{");
Print("\"trivial_loop_admissible\":", IsAdmissibleIdeal(trivial_loop));
Print(",\"trivial_acyclic_admissible\":", IsAdmissibleIdeal(trivial_acyclic));
Print(",\"linear_admissible\":", IsAdmissibleIdeal(linear));
Print(",\"monomial_is_monomial\":", IsMonomialIdeal(monomial));
Print(",\"derived_monomial\":", IsMonomialIdeal(derived_monomial));
Print("},\"homogeneous\":{");
Print("\"empty_degree_2\":", homogeneous_empty);
Print(",\"degree_2\":", IsHomogeneousListOfElements([x^2, 2*x^2, zero], 2));
Print(",\"mixed_degree_2\":", IsHomogeneousListOfElements([x^2+x^3], 2));
Print("},\"product_generators\":");
PrintElementList(GeneratorsOfTwoSidedIdeal(product));
Print(",\"minimal_generators\":");
PrintElementList(minimal);
Print(",\"f2\":{\"admissible\":", IsAdmissibleIdeal(monomial2));
Print(",\"monomial\":", IsMonomialIdeal(monomial2));
Print(",\"minimal_generators\":");
PrintPrimeElementList(minimal2);
Print("}");
Print("}\n");

QUIT_GAP(0);
