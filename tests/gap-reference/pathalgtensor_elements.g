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

q := Quiver(["u", "v"], [["u", "u", "a"], ["u", "v", "b"]]);;
parent := PathAlgebra(Rationals, q);;
quotient := parent / [parent.a * parent.a];;
parent_element := 2 * parent.u + 3 * parent.a + 5 * parent.a * parent.b;;
quotient_element := 2 * quotient.u + 3 * quotient.a
                    + 5 * quotient.a * quotient.b;;

q2 := Quiver(2, [[1, 2, "c"]]);;
binary := PathAlgebra(GF(2), q2);;
binary_element := binary.v1 + binary.c;;

Print("{\"schema_version\":1,\"kind\":\"pathalgtensor_elements\",\"oracle\":{");
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
Print("\"include_parent\":");
PrintTerms(IncludeInPathAlgebra(q.a * q.b, parent));
Print(",\"include_quotient_zero\":");
PrintTerms(IncludeInPathAlgebra(q.a * q.a, quotient));
Print(",\"parent_vertices\":");
PrintElementList(VerticesOfPathAlgebra(parent));
Print(",\"quotient_vertices\":");
PrintElementList(VerticesOfPathAlgebra(quotient));
Print(",\"parent_terms\":");
PrintTerms(parent_element);
Print(",\"quotient_terms\":");
PrintTerms(quotient_element);
Print("},\"binary\":{");
Print("\"include\":");
PrintTerms(IncludeInPathAlgebra(q2.c, binary));
Print(",\"vertices\":");
PrintElementList(VerticesOfPathAlgebra(binary));
Print(",\"terms\":");
PrintTerms(binary_element);
Print("}}\n");

QUIT_GAP(0);
