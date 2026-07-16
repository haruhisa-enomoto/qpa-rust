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
    local data, i;
    data := CoefficientsAndMagmaElements(element);
    Print("[");
    for i in [1, 3 .. Length(data) - 1] do
        if i > 1 then
            Print(",");
        fi;
        Print("{\"path\":");
        PrintJsonString(StringPrint(data[i]));
        Print(",\"coefficient\":[", NumeratorRat(data[i + 1]), ",",
              DenominatorRat(data[i + 1]), "]}");
    od;
    Print("]");
end;

PrintElementMatrix := function(matrix)
    local row, firstrow, entry, firstentry;
    Print("[");
    firstrow := true;
    for row in matrix do
        if not firstrow then
            Print(",");
        fi;
        Print("[");
        firstentry := true;
        for entry in row do
            if not firstentry then
                Print(",");
            fi;
            PrintRationalElement(entry);
            firstentry := false;
        od;
        Print("]");
        firstrow := false;
    od;
    Print("]");
end;

Q := Quiver(3, [[1, 2, "a"], [1, 3, "b"], [2, 3, "c"]]);;
A := PathAlgebra(Rationals, Q);;
matrix := [[A.a+A.b, Zero(A)], [A.c, A.v2]];;

presentation := rec();;
ObjectifyWithAttributes(
    presentation,
    NewType(ProjectivePresentationFamily, IsProjectivePresentationDefaultRep),
    ParentAlgebra, A,
    MatrixRepresentation, matrix,
    MatrixNumRows, 2,
    MatrixNumCols, 2,
    P0VertexList, [Q.v1, Q.v2]
);;

Print("{\"schema_version\":1,\"kind\":\"present2\",\"oracle\":{");
Print("\"gap_version\":");
PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":");
PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print("},\"category\":", IsProjectivePresentation(presentation));
Print(",\"default_rep\":", IsProjectivePresentationDefaultRep(presentation));
Print(",\"global_family\":", FamilyObj(presentation) = ProjectivePresentationFamily);
Print(",\"parent_identical\":", IsIdenticalObj(ParentAlgebra(presentation), A));
Print(",\"matrix_equal\":", MatrixRepresentation(presentation) = matrix);
Print(",\"matrix\":");
PrintElementMatrix(MatrixRepresentation(presentation));
Print(",\"shape\":[", MatrixNumRows(presentation), ",", MatrixNumCols(presentation), "]");
Print(",\"p0\":");
Print("[");
PrintJsonString(StringPrint(P0VertexList(presentation)[1]));
Print(",");
PrintJsonString(StringPrint(P0VertexList(presentation)[2]));
Print("]");
Print(",\"has_p1\":", HasP1VertexList(presentation));
Print("}\n");

QUIT_GAP(0);
