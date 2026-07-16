if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

ProductOfFields := function(field, count)
    local table, position;
    table := EmptySCTable(count, Zero(field));
    for position in [1..count] do
        SetEntrySCTable(table, position, position, [One(field), position]);
    od;
    return AlgebraByStructureConstants(field, table);
end;

DualNumbers := function(field)
    local table;
    table := EmptySCTable(2, Zero(field));
    SetEntrySCTable(table, 1, 1, [One(field), 1]);
    SetEntrySCTable(table, 1, 2, [One(field), 2]);
    SetEntrySCTable(table, 2, 1, [One(field), 2]);
    return AlgebraByStructureConstants(field, table);
end;

ZeroAlgebra := function(field)
    return AlgebraByStructureConstants(field, EmptySCTable(1, Zero(field)));
end;

IsZeroElement := function(algebra, element)
    return element = Zero(algebra);
end;

PolynomialCoefficients := function(polynomial)
    local data;
    data := CoefficientsOfLaurentPolynomial(polynomial);
    return Concatenation(ListWithIdenticalEntries(data[2],
        Zero(data[1][1])), data[1]);
end;

PairwiseOrthogonal := function(idempotents)
    local left, right;
    for left in [1..Length(idempotents)] do
        for right in [1..Length(idempotents)] do
            if left <> right and idempotents[left] * idempotents[right] <>
              Zero(Parent(idempotents[left])) then
                return false;
            fi;
        od;
    od;
    return true;
end;

PrintMatrixCase := function(field, size)
    local algebra, basis, e11_polynomial, irreducible_polynomial, zero_divisor,
          singular, has_before, idempotents, second, right_dimensions;
    algebra := FullMatrixAlgebra(field, size);
    basis := Basis(algebra);
    e11_polynomial := HOPF_MinPoly(algebra, BasisVectors(basis)[1]);
    irreducible_polynomial := fail;
    if Characteristic(field) = 2 and size = 2 then
        irreducible_polynomial := HOPF_MinPoly(algebra,
            BasisVectors(basis)[2] + BasisVectors(basis)[3] +
            BasisVectors(basis)[4]);
    fi;
    Reset(GlobalMersenneTwister, 1);
    zero_divisor := HOPF_ZeroDivisor(algebra);
    Reset(GlobalMersenneTwister, 1);
    singular := HOPF_SingularIdempotent(algebra);
    has_before := HasPrimitiveIdempotents(algebra);
    Reset(GlobalMersenneTwister, 1);
    idempotents := PrimitiveIdempotents(algebra);
    second := PrimitiveIdempotents(algebra);
    right_dimensions := List(idempotents, element -> Dimension(element * algebra));

    Print("{\"dimension\":", Dimension(algebra));
    Print(",\"e11_minpoly\":",
          List(PolynomialCoefficients(e11_polynomial), IntFFE));
    if irreducible_polynomial <> fail then
        Print(",\"irreducible_minpoly\":",
              List(PolynomialCoefficients(irreducible_polynomial), IntFFE));
    fi;
    Print(",\"left_identity\":",
          HOPF_LeftIdentity(algebra) = One(algebra));
    Print(",\"zero_divisor\":{\"nonzero\":",
          ForAll(zero_divisor, element -> not IsZeroElement(algebra, element)),
          ",\"product_zero\":", IsZeroElement(algebra,
              zero_divisor[1] * zero_divisor[2]), "}");
    Print(",\"singular\":{\"idempotent\":", singular * singular = singular,
          ",\"right_dimension\":", Dimension(singular * algebra), "}");
    Print(",\"primitive\":{\"has_before\":", has_before,
          ",\"has_after\":", HasPrimitiveIdempotents(algebra),
          ",\"cache_identical\":", IsIdenticalObj(idempotents, second),
          ",\"count\":", Length(idempotents),
          ",\"complete\":", Sum(idempotents) = One(algebra),
          ",\"idempotent\":",
          ForAll(idempotents, element -> element * element = element),
          ",\"orthogonal\":", PairwiseOrthogonal(idempotents),
          ",\"right_dimensions\":", right_dimensions, "}}");
end;

commutative_q := ProductOfFields(Rationals, 2);;
commutative_f2 := ProductOfFields(GF(2), 2);;
simple_f2 := ProductOfFields(GF(2), 1);;
zero_f3 := ZeroAlgebra(GF(3));;
nonsemisimple := CALL_WITH_CATCH(PrimitiveIdempotents, [DualNumbers(GF(2))]);;
nonsimple := CALL_WITH_CATCH(PrimitiveIdempotents, [commutative_f2]);;
nonfinite := CALL_WITH_CATCH(PrimitiveIdempotents,
    [FullMatrixAlgebra(Rationals, 2)]);;
nonfinite_zero_divisor := CALL_WITH_CATCH(HOPF_ZeroDivisor,
    [FullMatrixAlgebra(Rationals, 2)]);;

Print("{\"schema_version\":1,\"kind\":\"idempotent\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"f2_m2\":"); PrintMatrixCase(GF(2), 2);
Print(",\"f3_m2\":"); PrintMatrixCase(GF(3), 2);
Print(",\"f2_m3\":"); PrintMatrixCase(GF(2), 3);
Print(",\"commutative\":{");
Print("\"rational_zero_divisor_empty\":", HOPF_ZeroDivisor(commutative_q) = []);
Print(",\"binary_zero_divisor_empty\":", HOPF_ZeroDivisor(commutative_f2) = []);
Print(",\"rational_singular_is_one\":",
      HOPF_SingularIdempotent(commutative_q) = One(commutative_q));
Print(",\"simple_binary_primitive_count\":", Length(PrimitiveIdempotents(simple_f2)));
Print(",\"zero_left_identity_is_fail\":", HOPF_LeftIdentity(zero_f3) = fail, "}");
Print(",\"errors\":{");
Print("\"nonsemisimple_success\":", nonsemisimple[1]);
Print(",\"nonsimple_success\":", nonsimple[1]);
Print(",\"nonfinite_success\":", nonfinite[1]);
Print(",\"nonfinite_zero_divisor_success\":", nonfinite_zero_divisor[1], "}");
Print("}\n");

QUIT_GAP(0);
