if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
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

q1 := Quiver(["u1", "u2"], [["u1", "u2", "a"]]);;
q2 := Quiver(["v"], [["v", "v", "b"]]);;
left := PathAlgebra(Rationals, q1);;
right_parent := PathAlgebra(Rationals, q2);;
right := right_parent / [right_parent.b * right_parent.b];;
tensor := TensorProductOfAlgebras(left, right);;
left_inclusion := TensorAlgebraInclusion(tensor, 1);;
right_inclusion := TensorAlgebraInclusion(tensor, 2);;
zero_inclusion := TensorAlgebraInclusion(tensor, 0);;
left_element := left.u1 + 2 * left.a;;
right_element := right.v + 3 * right.b;;

binary_left := PathAlgebra(GF(2), q1);;
binary_right_parent := PathAlgebra(GF(2), q2);;
binary_right := binary_right_parent /
    [binary_right_parent.b * binary_right_parent.b];;
binary_tensor := TensorProductOfAlgebras(binary_left, binary_right);;
binary_left_inclusion := TensorAlgebraInclusion(binary_tensor, 1);;
binary_right_inclusion := TensorAlgebraInclusion(binary_tensor, 2);;

Print("{\"schema_version\":1,\"kind\":\"pathalgtensor_inclusion\",\"oracle\":{");
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
Print("\"left_generators\":");
PrintElementList(left_inclusion!.generators);
Print(",\"left_images\":");
PrintElementList(left_inclusion!.genimages);
Print(",\"right_generators\":");
PrintElementList(right_inclusion!.generators);
Print(",\"right_images\":");
PrintElementList(right_inclusion!.genimages);
Print(",\"zero_selects_second\":[",
      zero_inclusion!.generators = right_inclusion!.generators, ",",
      zero_inclusion!.genimages = right_inclusion!.genimages, "]");
Print(",\"source_range_identity\":[",
      IsIdenticalObj(Source(left_inclusion), left), ",",
      IsIdenticalObj(Source(right_inclusion), right), ",",
      IsIdenticalObj(Range(left_inclusion), tensor), ",",
      IsIdenticalObj(Range(right_inclusion), tensor), "]");
Print(",\"left_element_image\":");
PrintTerms(ImageElm(left_inclusion, left_element));
Print(",\"right_element_image\":");
PrintTerms(ImageElm(right_inclusion, right_element));
Print(",\"zero_index_image\":");
PrintTerms(ImageElm(zero_inclusion, right_element));
Print("},\"binary\":{");
Print("\"left_images\":");
PrintElementList(binary_left_inclusion!.genimages);
Print(",\"right_images\":");
PrintElementList(binary_right_inclusion!.genimages);
Print("}}\n");

QUIT_GAP(0);
