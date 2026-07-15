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

q1 := Quiver(["u1", "u2"], [["u1", "u2", "a"]]);;
q2 := Quiver(["v"], [["v", "v", "b"]]);;

left := PathAlgebra(Rationals, q1);;
right_parent := PathAlgebra(Rationals, q2);;
right := right_parent / [right_parent.b * right_parent.b];;
tensor := TensorProductOfAlgebras(left, right);;
tensor_quiver := QuiverOfPathAlgebra(tensor);;
decomposition := TensorProductDecomposition(tensor);;

binary_left := PathAlgebra(GF(2), q1);;
binary_right := PathAlgebra(GF(2), q2);;
binary_tensor := TensorProductOfAlgebras(binary_left, binary_right);;

Print("{\"schema_version\":1,\"kind\":\"pathalgtensor_full\",\"oracle\":{");
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
Print("},\"rational\":{");
Print("\"vertices\":");
PrintJsonStringList(List(VerticesOfQuiver(tensor_quiver), String));
Print(",\"arrows\":");
PrintJsonStringList(List(ArrowsOfQuiver(tensor_quiver), String));
Print(",\"relations\":");
PrintElementList(RelatorsOfFpAlgebra(tensor));
Print(",\"groebner_basis\":");
PrintElementList(AttachedBasis(tensor));
Print(",\"factor_identity\":[", IsIdenticalObj(decomposition[1], left), ",",
      IsIdenticalObj(decomposition[2], right), "]");
Print(",\"simple_tensor\":");
PrintTerms(SimpleTensor([left.a, right.b], tensor));
Print("},\"binary\":{");
Print("\"relations\":");
PrintElementList(RelatorsOfFpAlgebra(binary_tensor));
Print(",\"groebner_basis\":");
PrintElementList(AttachedBasis(binary_tensor));
Print("}}\n");

QUIT_GAP(0);
