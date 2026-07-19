if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

SizeScreen([100000, 100000]);

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

PrintCoefficient := function(field, coefficient)
    if IsFinite(field) then
        if Size(field) = 2 then
            if coefficient = Zero(field) then Print("0"); else Print("1"); fi;
        else
            Print(IntFFE(coefficient));
        fi;
    else
        Print("[", NumeratorRat(coefficient), ",", DenominatorRat(coefficient), "]");
    fi;
end;

PrintRelation := function(field, relation)
    local data, i;
    data := CoefficientsAndMagmaElements(relation);
    Print("[");
    for i in [1, 3 .. Length(data) - 1] do
        if i > 1 then Print(","); fi;
        Print("{\"path\":"); PrintJsonString(StringPrint(data[i]));
        Print(",\"coefficient\":"); PrintCoefficient(field, data[i + 1]);
        Print("}");
    od;
    Print("]");
end;

PrintAlgebra := function(algebra)
    local quiver, adjacency, relations, i;
    quiver := QuiverOfPathAlgebra(OriginalPathAlgebra(algebra));
    adjacency := AdjacencyMatrixOfQuiver(quiver);
    relations := RelationsOfAlgebra(algebra);
    Print("{\"is_path\":", IsPathAlgebra(algebra));
    Print(",\"dimension\":", Dimension(algebra));
    Print(",\"adjacency\":", adjacency, ",\"relations\":[");
    for i in [1 .. Length(relations)] do
        if i > 1 then Print(","); fi;
        PrintRelation(LeftActingDomain(algebra), relations[i]);
    od;
    Print("]}");
end;

quiver := Quiver(3, [[1, 2, "a"], [2, 3, "b"]]);;
parent2 := PathAlgebra(GF(2), quiver);;
algebra2 := parent2 / [parent2.a * parent2.b];;
corner12_2 := QuiverAlgebraOfeAe(algebra2, algebra2.v1 + algebra2.v2);;
corner13_2 := QuiverAlgebraOfeAe(algebra2, algebra2.v1 + algebra2.v3);;
corner_conjugate_2 := QuiverAlgebraOfeAe(algebra2, algebra2.v1 + algebra2.a);;
corner_whole_2 := QuiverAlgebraOfeAe(algebra2, One(algebra2));;

minimal_quiver := Quiver(4,
    [[1, 2, "a"], [2, 3, "b"], [3, 4, "c"]]);;
minimal_parent := PathAlgebra(GF(2), minimal_quiver);;
minimal_algebra := minimal_parent / [minimal_parent.a * minimal_parent.b];;
minimal_corner := QuiverAlgebraOfeAe(minimal_algebra, One(minimal_algebra));;

parentq := PathAlgebra(Rationals, quiver);;
algebraq := parentq / [parentq.a * parentq.b];;
corner12_q := CALL_WITH_CATCH(QuiverAlgebraOfeAe,
    [algebraq, algebraq.v1 + algebraq.v2]);;

path_result := CALL_WITH_CATCH(QuiverAlgebraOfeAe, [parent2, parent2.v1]);;
other := PathAlgebra(GF(2), Quiver(1, []));;
foreign_result := CALL_WITH_CATCH(QuiverAlgebraOfeAe, [algebra2, other.v1]);;
nonidempotent_result := CALL_WITH_CATCH(QuiverAlgebraOfeAe, [algebra2, algebra2.a]);;
zero_result := CALL_WITH_CATCH(QuiverAlgebraOfeAe, [algebra2, Zero(algebra2)]);;

Print("{\"schema_version\":1,\"kind\":\"pathalg_eae\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION); Print("}");
Print(",\"binary\":{\"corner12\":"); PrintAlgebra(corner12_2);
Print(",\"corner13\":"); PrintAlgebra(corner13_2);
Print(",\"conjugate\":"); PrintAlgebra(corner_conjugate_2);
Print(",\"whole\":"); PrintAlgebra(corner_whole_2);
Print(",\"minimal\":"); PrintAlgebra(minimal_corner); Print("}");
Print(",\"rational\":{\"success\":", corner12_q[1]);
if corner12_q[1] then Print(",\"corner12\":"); PrintAlgebra(corner12_q[2]); fi;
Print("}");
Print(",\"errors\":{\"path_success\":", path_result[1]);
Print(",\"foreign_success\":", foreign_result[1]);
Print(",\"nonidempotent_success\":", nonidempotent_result[1]);
Print(",\"zero_success\":", zero_result[1], "}}\n");

QUIT_GAP(0);
