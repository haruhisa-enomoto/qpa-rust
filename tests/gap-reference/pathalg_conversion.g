if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

UnderlyingRepresentative := function(element)
    if IsElementOfQuotientOfPathAlgebra(element) then
        return element![1];
    fi;
    return element;
end;

PrintRationalElement := function(element)
    local data, index;
    data := CoefficientsAndMagmaElements(UnderlyingRepresentative(element));
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

PrintRationalVector := function(vector)
    local index;
    Print("[");
    for index in [1..Length(vector)] do
        if index > 1 then Print(","); fi;
        Print("[", NumeratorRat(vector[index]), ",", DenominatorRat(vector[index]), "]");
    od;
    Print("]");
end;

PrintPrimeElement := function(element)
    local data, index;
    data := CoefficientsAndMagmaElements(UnderlyingRepresentative(element));
    Print("[");
    for index in [1, 3 .. Length(data) - 1] do
        if index > 1 then Print(","); fi;
        Print("{\"path\":", ExtRepOfObj(data[index]));
        Print(",\"coefficient\":", IntFFE(data[index + 1]), "}");
    od;
    Print("]");
end;

PrintPrimeVector := function(vector)
    local index;
    Print("[");
    for index in [1..Length(vector)] do
        if index > 1 then Print(","); fi;
        Print(IntFFE(vector[index]));
    od;
    Print("]");
end;

Q := Quiver(3, [[1, 2, "a"], [2, 3, "b"], [1, 3, "c"]]);;
KQ := PathAlgebra(Rationals, Q, HighLevelGroebnerBasis);;
a := ElementOfPathAlgebra(KQ, Q.a);;
b := ElementOfPathAlgebra(KQ, Q.b);;
c := ElementOfPathAlgebra(KQ, Q.c);;
relation := a*b-c;;
A := KQ/[relation];;
B := CanonicalBasis(A);;
basis_vectors := BasisVectors(B);;
element := 2*A.v1+3*A.a+5*A.c;;
coefficients := Coefficients(B, element);;
rebuilt := ObjByExtRep(FamilyObj(Zero(A)), ExtRepOfObj(2*a+3*a*b));;
space := VectorSpace(Rationals, [A.v1, A.a, A.c], Zero(A));;
nice := NiceVector(space, element);;
ugly := UglyVector(space, nice);;

KQ2 := PathAlgebra(GF(2), Q, HighLevelGroebnerBasis);;
a2 := ElementOfPathAlgebra(KQ2, Q.a);;
b2 := ElementOfPathAlgebra(KQ2, Q.b);;
c2 := ElementOfPathAlgebra(KQ2, Q.c);;
A2 := KQ2/[a2*b2-c2];;
B2 := CanonicalBasis(A2);;
element2 := A2.v1+A2.a+A2.c;;
coefficients2 := Coefficients(B2, element2);;
rebuilt2 := ObjByExtRep(FamilyObj(Zero(A2)), ExtRepOfObj(a2+a2*b2));;

QL := Quiver(1, [[1, 1, "x"]]);;
KQL := PathAlgebra(Rationals, QL);;
xL := ElementOfPathAlgebra(KQL, QL.x);;
ILEarly := Ideal(KQL, [xL^2]);;
LEarly := KQL/ILEarly;;
lazy_cached_basis_before := HasBasisVectors(CanonicalBasis(LEarly));;
Dimension(LEarly);;
lazy_cached_basis_after := HasBasisVectors(CanonicalBasis(LEarly));;
IL := Ideal(KQL, [xL^2]);;
L := KQL/IL;;
famL := ElementsFamily(FamilyObj(L));;
lazy_before_basis := HasGroebnerBasisOfIdeal(IL);;
lazy_before_normal_function := HasNormalFormFunction(famL);;
lazy_before := ElementOfQuotientOfPathAlgebra(famL, xL+xL^2, false);;
lazy_dimension := Dimension(L);;
lazy_after_basis := HasGroebnerBasisOfIdeal(IL);;
lazy_tip_reduced := IsTipReducedGroebnerBasis(GroebnerBasisOfIdeal(IL));;
lazy_after_normal_function := HasNormalFormFunction(famL);;
lazy_after := ElementOfQuotientOfPathAlgebra(famL, xL+xL^2, false);;
lazy_basis := BasisVectors(CanonicalBasis(L));;

Print("{\"schema_version\":1,\"kind\":\"pathalg_conversion\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},\"rational\":{\"parent_relations\":");
PrintElementList(RelationsOfAlgebra(KQ));
Print(",\"quotient_relations\":");
PrintElementList(RelationsOfAlgebra(A));
Print(",\"ideal_identical\":", IsIdenticalObj(IdealOfQuotient(A),
      ElementsFamily(FamilyObj(A))!.ideal));
Print(",\"finite\":", IsFiniteDimensional(A));
Print(",\"dimension\":", Dimension(A));
Print(",\"basis\":");
PrintElementList(basis_vectors);
Print(",\"element\":");
PrintRationalElement(element);
Print(",\"coefficients\":");
PrintRationalVector(coefficients);
Print(",\"rebuilt\":");
PrintRationalElement(rebuilt);
Print(",\"rebuilt_is_normal_form\":", IsNormalForm(rebuilt));
Print(",\"nice\":");
PrintRationalVector(nice);
Print(",\"ugly\":");
PrintRationalElement(ugly);
Print(",\"vertex_positions\":[", VertexPosition(KQ.v2), ",", VertexPosition(A.v2), "]");
Print("},\"f2\":{\"basis\":");
Print("[");
for index in [1..Length(BasisVectors(B2))] do
    if index > 1 then Print(","); fi;
    PrintPrimeElement(BasisVectors(B2)[index]);
od;
Print("]");
Print(",\"coefficients\":");
PrintPrimeVector(coefficients2);
Print(",\"rebuilt\":");
PrintPrimeElement(rebuilt2);
Print(",\"rebuilt_is_normal_form\":", IsNormalForm(rebuilt2));
Print("},\"lazy_default\":{\"before_basis\":", lazy_before_basis);
Print(",\"before_normal_function\":", lazy_before_normal_function);
Print(",\"before_element\":");
PrintRationalElement(lazy_before);
Print(",\"before_element_is_normal\":", IsNormalForm(lazy_before));
Print(",\"cached_basis_before\":", lazy_cached_basis_before);
Print(",\"cached_basis_after\":", lazy_cached_basis_after);
Print(",\"dimension\":", lazy_dimension);
Print(",\"after_basis\":", lazy_after_basis);
Print(",\"tip_reduced\":", lazy_tip_reduced);
Print(",\"after_normal_function\":", lazy_after_normal_function);
Print(",\"after_element\":");
PrintRationalElement(lazy_after);
Print(",\"after_element_is_normal\":", IsNormalForm(lazy_after));
Print(",\"basis\":");
PrintElementList(lazy_basis);
Print(",\"basis_normal_flags\":[");
for index in [1..Length(lazy_basis)] do
    if index > 1 then Print(","); fi;
    Print(IsNormalForm(lazy_basis[index]));
od;
Print("]}}");
Print("\n");

QUIT_GAP(0);
