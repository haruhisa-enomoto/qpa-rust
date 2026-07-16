if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintRationalCohen := function(polynomial)
    local i;
    Print("[");
    for i in [1 .. Length(polynomial[1])] do
        if i > 1 then Print(","); fi;
        Print("{\"monomial\":", polynomial[1][i]);
        Print(",\"coefficient\":[", NumeratorRat(polynomial[2][i]), ",",
              DenominatorRat(polynomial[2][i]), "]}");
    od;
    Print("]");
end;

PrintRationalCohenList := function(polynomials)
    local i;
    Print("[");
    for i in [1 .. Length(polynomials)] do
        if i > 1 then Print(","); fi;
        PrintRationalCohen(polynomials[i]);
    od;
    Print("]");
end;

PrintPrimeCohen := function(polynomial)
    local i;
    Print("[");
    for i in [1 .. Length(polynomial[1])] do
        if i > 1 then Print(","); fi;
        Print("{\"monomial\":", polynomial[1][i]);
        Print(",\"coefficient\":", IntFFE(polynomial[2][i]), "}");
    od;
    Print("]");
end;

PrintPrimeCohenList := function(polynomials)
    local i;
    Print("[");
    for i in [1 .. Length(polynomials)] do
        if i > 1 then Print(","); fi;
        PrintPrimeCohen(polynomials[i]);
    od;
    Print("]");
end;

PrintRationalElement := function(element)
    local data, i;
    data := CoefficientsAndMagmaElements(element);
    Print("[");
    for i in [1, 3 .. Length(data) - 1] do
        if i > 1 then Print(","); fi;
        Print("{\"path\":", ExtRepOfObj(data[i]));
        Print(",\"coefficient\":[", NumeratorRat(data[i + 1]), ",",
              DenominatorRat(data[i + 1]), "]}");
    od;
    Print("]");
end;

PrintRationalElementList := function(elements)
    local i;
    Print("[");
    for i in [1 .. Length(elements)] do
        if i > 1 then Print(","); fi;
        PrintRationalElement(elements[i]);
    od;
    Print("]");
end;

PrintPrimeElement := function(element)
    local data, i;
    data := CoefficientsAndMagmaElements(element);
    Print("[");
    for i in [1, 3 .. Length(data) - 1] do
        if i > 1 then Print(","); fi;
        Print("{\"path\":", ExtRepOfObj(data[i]));
        Print(",\"coefficient\":", IntFFE(data[i + 1]), "}");
    od;
    Print("]");
end;

PrintPrimeElementList := function(elements)
    local i;
    Print("[");
    for i in [1 .. Length(elements)] do
        if i > 1 then Print(","); fi;
        PrintPrimeElement(elements[i]);
    od;
    Print("]");
end;

Q := Quiver(1, [[1, 1, "a"], [1, 1, "b"]]);;
KQ := PathAlgebra(Rationals, Q);;
v := ElementOfPathAlgebra(KQ, Q.v1);;
a := ElementOfPathAlgebra(KQ, Q.a);;
b := ElementOfPathAlgebra(KQ, Q.b);;
p1 := a * a * b - v;;
p2 := a * b * b - v;;
raw := QPA_Path2Cohen([p1, p2]);;
free := QPA_Path2CohenFree([p1, p2]);;
nc_basis := GBNPGroebnerBasisNC([p1, p2], KQ);;
checked_basis := GBNPGroebnerBasis([a * a, a * b], KQ);;
checked_ideal := TwoSidedIdealByGenerators(KQ, [a * a, a * b]);;
ideal_basis := GBNPGroebnerBasis(checked_ideal);;

F5 := GF(5);;
KQ5 := PathAlgebra(F5, Q);;
v5 := ElementOfPathAlgebra(KQ5, Q.v1);;
a5 := ElementOfPathAlgebra(KQ5, Q.a);;
b5 := ElementOfPathAlgebra(KQ5, Q.b);;
p15 := a5 * a5 * b5 - v5;;
p25 := a5 * b5 * b5 - v5;;
free5 := QPA_Path2CohenFree([p15, p25]);;
nc_basis5 := GBNPGroebnerBasisNC([p15, p25], KQ5);;

F2 := GF(2);;
KQ2 := PathAlgebra(F2, Q);;
v2 := ElementOfPathAlgebra(KQ2, Q.v1);;
a2 := ElementOfPathAlgebra(KQ2, Q.a);;
b2 := ElementOfPathAlgebra(KQ2, Q.b);;
p12 := a2 * a2 * b2 - v2;;
p22 := a2 * b2 * b2 - v2;;
free2 := QPA_Path2CohenFree([p12, p22]);;
nc_basis2 := GBNPGroebnerBasisNC([p12, p22], KQ2);;

Line := Quiver(2, [[1, 2, "c"]]);;
KLine := PathAlgebra(Rationals, Line);;
c := ElementOfPathAlgebra(KLine, Line.c);;
line_v1 := ElementOfPathAlgebra(KLine, Line.v1);;
line_v2 := ElementOfPathAlgebra(KLine, Line.v2);;
line_nonuniform := line_v1 + c + line_v2;;
uniform_element := MakeUniform(line_nonuniform);;
uniform_list := MakeUniform([line_nonuniform, c]);;
structure := QPA_RelationsForPathAlgebra(KLine);;
structure_as_paths := QPA_Cohen2Path(structure, KLine);;
line_nc_basis := GBNPGroebnerBasisNC([c], KLine);;

Print("{\"schema_version\":1,\"kind\":\"qpa_gbnp\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":");
PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"rational\":{\"raw\":");
PrintRationalCohenList(raw);
Print(",\"free\":");
PrintRationalCohenList(free);
Print(",\"nc_basis\":");
PrintRationalElementList(nc_basis);
Print(",\"checked_basis\":");
PrintRationalElementList(checked_basis);
Print(",\"ideal_basis\":");
PrintRationalElementList(ideal_basis);
Print(",\"vertex_is_in_arrow_ideal\":", QPA_InArrowIdeal([v], KQ));
Print("},\"f5\":{\"free\":");
PrintPrimeCohenList(free5);
Print(",\"nc_basis\":");
PrintPrimeElementList(nc_basis5);
Print("},\"f2\":{\"free\":");
PrintPrimeCohenList(free2);
Print(",\"nc_basis\":");
PrintPrimeElementList(nc_basis2);
Print("},\"multi_vertex\":{\"structure_relations\":");
PrintRationalCohenList(structure);
Print(",\"uniform_element\":");
PrintRationalElementList(uniform_element);
Print(",\"uniform_list\":");
PrintRationalElementList(uniform_list);
Print(",\"structure_as_paths\":");
PrintRationalElementList(structure_as_paths);
Print(",\"nc_basis\":");
PrintRationalElementList(line_nc_basis);
Print("}}\n");

QUIT_GAP(0);
