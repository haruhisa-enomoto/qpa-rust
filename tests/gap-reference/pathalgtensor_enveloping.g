if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintJsonStringList := function(values)
    local first, value;
    Print("[");
    first := true;
    for value in values do
        if not first then
            Print(",");
        fi;
        PrintJsonString(value);
        first := false;
    od;
    Print("]");
end;

PathText := function(path)
    if IsQuiverVertex(path) then
        return String(path);
    elif IsZeroPath(path) then
        return "0";
    fi;
    return JoinStringsWithSeparator(List(WalkOfPathOrVertex(path), String), "*");
end;

PrintTerms := function(element)
    local first, term;
    Print("[");
    first := true;
    for term in PathAlgebraElementTerms(element) do
        if not first then
            Print(",");
        fi;
        Print("{\"monomial\":");
        PrintJsonString(PathText(term.monom));
        Print(",\"coefficient\":");
        PrintJsonString(String(term.coeff));
        Print("}");
        first := false;
    od;
    Print("]");
end;

PrintElementList := function(elements)
    local element, first;
    Print("[");
    first := true;
    for element in elements do
        if not first then
            Print(",");
        fi;
        PrintTerms(element);
        first := false;
    od;
    Print("]");
end;

AttachedBasis := function(algebra)
    local family;
    family := ElementsFamily(FamilyObj(algebra));
    return GroebnerBasisOfIdeal(family!.ideal);
end;

qloop := Quiver(["v"], [["v", "v", "x"]]);;
parent := PathAlgebra(Rationals, qloop);;
algebra := parent / [parent.x * parent.x];;
input_admissible_cached_before :=
    HasIsAdmissibleQuotientOfPathAlgebra(algebra);;
opposite_cached_before := HasOppositePathAlgebra(algebra);;
enveloping := EnvelopingAlgebra(algebra);;
enveloping_again := EnvelopingAlgebra(algebra);;
decomposition := TensorProductDecomposition(enveloping);;
opposite := OppositeAlgebra(algebra);;

nonadmissible_parent := PathAlgebra(Rationals, qloop);;
nonadmissible := nonadmissible_parent / [nonadmissible_parent.x];;
nonadmissible_cached_before :=
    HasIsAdmissibleQuotientOfPathAlgebra(nonadmissible);;
nonadmissible_enveloping := EnvelopingAlgebra(nonadmissible);;

qpath := Quiver(["u1", "u2"], [["u1", "u2", "a"]]);;
path_algebra := PathAlgebra(Rationals, qpath);;
path_enveloping := EnvelopingAlgebra(path_algebra);;
path_enveloping_again := EnvelopingAlgebra(path_algebra);;
path_decomposition := TensorProductDecomposition(path_enveloping);;

binary_parent := PathAlgebra(GF(2), qloop);;
binary_algebra := binary_parent / [binary_parent.x * binary_parent.x];;
binary_enveloping := EnvelopingAlgebra(binary_algebra);;

Print("{\"schema_version\":1,\"kind\":\"pathalgtensor_enveloping\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":");
PrintJsonString(PackageInfo("GBNP")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"rational_quotient\":{");
Print("\"property\":[", IsEnvelopingAlgebra(algebra), ",",
      IsEnvelopingAlgebra(enveloping), ",",
      IsIdenticalObj(enveloping, enveloping_again), "]");
Print(",\"input_admissibility\":[", input_admissible_cached_before, ",",
      HasIsAdmissibleQuotientOfPathAlgebra(algebra), ",",
      IsAdmissibleQuotientOfPathAlgebra(algebra), "]");
Print(",\"enveloping_admissibility\":[",
      HasIsAdmissibleQuotientOfPathAlgebra(enveloping), ",",
      IsAdmissibleQuotientOfPathAlgebra(enveloping), "]");
Print(",\"opposite_cache\":[", opposite_cached_before, ",",
      HasOppositePathAlgebra(algebra), ",",
      IsIdenticalObj(decomposition[1], opposite), ",",
      IsIdenticalObj(OppositeAlgebra(opposite), algebra), "]");
Print(",\"factor_identity\":[",
      IsIdenticalObj(decomposition[2], algebra), ",",
      IsIdenticalObj(decomposition[1], opposite), "]");
Print(",\"vertices\":");
PrintJsonStringList(List(VerticesOfQuiver(QuiverOfPathAlgebra(enveloping)), String));
Print(",\"arrows\":");
PrintJsonStringList(List(ArrowsOfQuiver(QuiverOfPathAlgebra(enveloping)), String));
Print(",\"relations\":");
PrintElementList(RelatorsOfFpAlgebra(enveloping));
Print(",\"groebner_basis\":");
PrintElementList(AttachedBasis(enveloping));
Print(",\"nonadmissible\":[", nonadmissible_cached_before, ",",
      HasIsAdmissibleQuotientOfPathAlgebra(nonadmissible), ",",
      IsAdmissibleQuotientOfPathAlgebra(nonadmissible), ",",
      HasIsAdmissibleQuotientOfPathAlgebra(nonadmissible_enveloping), ",",
      IsEnvelopingAlgebra(nonadmissible_enveloping), "]");
Print("},\"rational_path\":{");
Print("\"property\":[", IsEnvelopingAlgebra(path_algebra), ",",
      IsEnvelopingAlgebra(path_enveloping), ",",
      IsIdenticalObj(path_enveloping, path_enveloping_again), "]");
Print(",\"admissibility_cached\":",
      HasIsAdmissibleQuotientOfPathAlgebra(path_enveloping));
Print(",\"factor_identity\":[",
      IsIdenticalObj(path_decomposition[2], path_algebra), ",",
      IsIdenticalObj(path_decomposition[1], OppositeAlgebra(path_algebra)), "]");
Print(",\"relations\":");
PrintElementList(RelatorsOfFpAlgebra(path_enveloping));
Print("},\"binary\":{");
Print("\"property\":", IsEnvelopingAlgebra(binary_enveloping));
Print(",\"admissibility\":[",
      HasIsAdmissibleQuotientOfPathAlgebra(binary_enveloping), ",",
      IsAdmissibleQuotientOfPathAlgebra(binary_enveloping), "]");
Print(",\"relations\":");
PrintElementList(RelatorsOfFpAlgebra(binary_enveloping));
Print("}}\n");

QUIT_GAP(0);
