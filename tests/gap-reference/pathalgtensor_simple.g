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

q1 := Quiver(["u1", "u2"], [["u1", "u2", "a"]]);;
q2 := Quiver(["v"], [["v", "v", "b"]]);;
product := QuiverProduct(q1, q2);;

left := PathAlgebra(Rationals, q1);;
right_parent := PathAlgebra(Rationals, q2);;
right := right_parent / [right_parent.b * right_parent.b];;
target_parent := PathAlgebra(Rationals, product);;
target_quotient := target_parent / [target_parent.a_v * target_parent.u2_b];;
x := 2 * left.u1 + 3 * left.a;;
y := 5 * right.v + 7 * right.b;;

binary_left := PathAlgebra(GF(2), q1);;
binary_right := PathAlgebra(GF(2), q2);;
binary_target := PathAlgebra(GF(2), product);;
binary_x := binary_left.u1 + binary_left.a;;
binary_y := binary_right.v + binary_right.b;;

Print("{\"schema_version\":1,\"kind\":\"pathalgtensor_simple\",\"oracle\":{");
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
Print("\"quotient_factor\":");
PrintTerms(SimpleTensor([x, y], target_parent));
Print(",\"zero_short_circuit\":");
PrintTerms(SimpleTensor([Zero(left), y], target_parent));
Print(",\"quotient_target_reduction\":");
PrintTerms(SimpleTensor([left.a, right.b], target_quotient));
Print("},\"binary\":{");
Print("\"expansion\":");
PrintTerms(SimpleTensor([binary_x, binary_y], binary_target));
Print("}}\n");

QUIT_GAP(0);
