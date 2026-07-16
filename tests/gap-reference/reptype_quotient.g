if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    escaped := ReplacedString(escaped, "\n", "\\n");
    Print("\"", escaped, "\"");
end;

PrintCaught := function(caught)
    Print("{\"success\":", caught[1], ",\"value\":");
    if caught[1] and caught[2] <> fail then
        Print(caught[2]);
    else
        Print("null");
    fi;
    Print("}");
end;

semi_quiver := Quiver(["s"], []);;
semi_parent := PathAlgebra(Rationals, semi_quiver);;
semi_ideal := TwoSidedIdealByGenerators(semi_parent, [Zero(semi_parent)]);;
semi := FactorPathAlgebraByRelators(
    semi_parent, semi_ideal, OrderingOfAlgebra(semi_parent));;

dynkin_quiver := Quiver(["1", "2", "3"],
    [["1", "2", "a"], ["2", "3", "b"]]);;
dynkin_parent := PathAlgebra(GF(2), dynkin_quiver);;
dynkin := dynkin_parent / [dynkin_parent.a * dynkin_parent.b];;

loop_quiver := Quiver(["v"], [["v", "v", "x"]]);;
loop_parent := PathAlgebra(Rationals, loop_quiver);;
dual_numbers := loop_parent / [loop_parent.x^2];;
truncated_cube := loop_parent / [loop_parent.x^3];;
nonadmissible := loop_parent / [loop_parent.x];;

wild_quiver := Quiver(["w"],
    [["w", "w", "x"], ["w", "w", "y"]]);;
wild_parent := PathAlgebra(GF(2), wild_quiver);;
wild := wild_parent /
    [wild_parent.x^2, wild_parent.x * wild_parent.y,
     wild_parent.y * wild_parent.x, wild_parent.y^2];;

diamond_quiver := Quiver(["1", "2", "3", "4", "5", "6", "7"],
    [["1", "2", "a"], ["1", "3", "b"],
     ["2", "4", "c"], ["3", "4", "d"],
     ["5", "6", "u"], ["6", "7", "z"]]);;
diamond_parent := PathAlgebra(Rationals, diamond_quiver);;
diamond := diamond_parent / [diamond_parent.u * diamond_parent.z];;

star_quiver := Quiver(["1", "2", "c", "4", "5"],
    [["1", "c", "a"], ["2", "c", "b"],
     ["c", "4", "d"], ["c", "5", "e"]]);;
star_parent := PathAlgebra(Rationals, star_quiver);;
star := star_parent /
    [star_parent.a * star_parent.e, star_parent.b * star_parent.d];;

Print("CASE semisimple\n");
semi_value := IsFiniteTypeAlgebra(semi);;
Print("CASE dynkin\n");
dynkin_value := IsFiniteTypeAlgebra(dynkin);;
Print("CASE radical-square-zero\n");
dual_value := IsFiniteTypeAlgebra(dual_numbers);;
Print("CASE separated-nondynkin\n");
wild_value := IsFiniteTypeAlgebra(wild);;
Print("CASE nondistributive\n");
diamond_value := IsFiniteTypeAlgebra(diamond);;
Print("CASE special-biserial\n");
star_value := IsFiniteTypeAlgebra(star);;
Print("CASE inconclusive\n");
cube_value := IsFiniteTypeAlgebra(truncated_cube);;
Print("CASE nonadmissible\n");
nonadmissible_value := CALL_WITH_CATCH(IsFiniteTypeAlgebra, [nonadmissible]);;

Print("{\"schema_version\":1,\"kind\":\"reptype_quotient\",\"oracle\":{");
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
Print("},\"values\":{");
Print("\"semisimple\":", semi_value);
Print(",\"dynkin\":", dynkin_value);
Print(",\"radical_square_zero\":", dual_value);
Print(",\"separated_nondynkin\":", wild_value);
Print(",\"nondistributive\":", diamond_value);
Print(",\"special_biserial\":", star_value);
Print(",\"inconclusive\":");
if cube_value = fail then Print("null"); else Print(cube_value); fi;
Print(",\"nonadmissible\":"); PrintCaught(nonadmissible_value);
Print("},\"distributive\":{");
Print("\"diamond\":", IsDistributiveAlgebra(diamond));
Print(",\"star\":", IsDistributiveAlgebra(star));
Print(",\"truncated_cube\":", IsDistributiveAlgebra(truncated_cube));
Print(",\"wild_gf2\":", IsDistributiveAlgebra(wild));
Print(",\"path_tree\":", IsDistributiveAlgebra(star_parent));
Print(",\"path_disconnected\":", IsDistributiveAlgebra(diamond_parent));
Print("},\"is_quotient\":[", IsQuotientOfPathAlgebra(semi), ",",
    IsQuotientOfPathAlgebra(diamond), "]}\n");

QUIT_GAP(0);
