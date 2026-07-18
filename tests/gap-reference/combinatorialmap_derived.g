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

PrintGenus := function(genus)
    Print("[", NumeratorRat(genus), ",", DenominatorRat(genus), "]");
end;

DerivedData := function(algebra)
    local map, boundary, basis_winding, gcd_list, entry, pair;
    map := CombinatorialMapOfGentleAlgebra(algebra);
    boundary := List(
        BoundaryCurvesOfCombinatorialMap(map),
        entry -> [-WindingNumber(map, entry[1]), entry[2]]);
    basis_winding := List(
        HomologyBasisOfCombinatorialMap(map),
        pair -> [WindingNumber(map, pair[1]), WindingNumber(map, pair[2])]);
    gcd_list := [];
    if GenusOfCombinatorialMap(map) = 1 then
        for entry in boundary do
            Add(gcd_list, entry[1] + 2);
        od;
        for pair in basis_winding do
            Append(gcd_list, pair);
        od;
    fi;
    Print("{\"genus\":"); PrintGenus(GenusOfCombinatorialMap(map));
    Print(",\"face_count\":", Length(FacesOfCombinatorialMap(map)));
    Print(",\"vertex_count\":", Length(OrbitsPerms(
        [OrderingOfCombinatorialMap(map)], [1..Size(map)])));
    Print(",\"boundary\":", boundary);
    Print(",\"basis_winding\":", basis_winding);
    if Length(gcd_list) = 0 then
        Print(",\"gcd\":null");
    else
        Print(",\"gcd\":", AbsInt(Gcd(Integers, gcd_list)));
    fi;
    Print("}");
end;

GenusOneA := function(field)
    local quiver, algebra;
    quiver := Quiver(6, [
        [2,1,"aa"], [3,2,"bb"], [1,4,"cc"],
        [4,2,"dd"], [3,5,"ee"], [5,4,"ff"],
        [4,6,"gg"], [6,5,"hh"], [5,6,"ii"]]);
    algebra := PathAlgebra(field, quiver);
    return algebra/[
        algebra.aa*algebra.cc, algebra.cc*algebra.dd,
        algebra.dd*algebra.aa, algebra.ee*algebra.ff,
        algebra.ff*algebra.gg, algebra.gg*algebra.hh,
        algebra.hh*algebra.ii];
end;

GenusOneB := function(field)
    local quiver, algebra;
    quiver := Quiver(6, [
        [2,1,"aaa"], [1,3,"bbb"], [3,2,"ccc"],
        [1,5,"ddd"], [2,6,"eee"], [3,4,"fff"],
        [4,5,"ggg"], [5,6,"hhh"], [6,4,"iii"]]);
    algebra := PathAlgebra(field, quiver);
    return algebra/[
        algebra.aaa*algebra.bbb, algebra.bbb*algebra.ccc,
        algebra.ccc*algebra.aaa, algebra.fff*algebra.ggg,
        algebra.ggg*algebra.hhh, algebra.hhh*algebra.iii];
end;

Chain := function(field)
    return PathAlgebra(field, Quiver(3, [
        [1,2,"chain_a"], [2,3,"chain_b"]]));
end;

Branch := function(field)
    local algebra;
    algebra := PathAlgebra(field, Quiver(4, [
        [1,2,"branch_a"], [2,3,"branch_b"], [2,4,"branch_c"]]));
    return algebra/[algebra.branch_a*algebra.branch_b];
end;

genus_one_a := GenusOneA(Rationals);
genus_one_b := GenusOneB(Rationals);
chain_q := Chain(Rationals);
chain_q_copy := Chain(Rationals);
chain_f2 := Chain(GF(2));
branch_q := Branch(Rationals);

Print("{\"schema_version\":1,\"kind\":\"combinatorialmap_derived\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"genus_one_a\":"); DerivedData(genus_one_a);
Print(",\"genus_one_b\":"); DerivedData(genus_one_b);
Print(",\"results\":{");
Print("\"source_genus_one\":", AreDerivedEquivalent(genus_one_a, genus_one_b));
Print(",\"reflexive_genus_one\":", AreDerivedEquivalent(genus_one_a, genus_one_a));
Print(",\"genus_zero_disc\":", AreDerivedEquivalent(chain_q, chain_q_copy));
Print(",\"different_field\":", AreDerivedEquivalent(chain_q, chain_f2));
Print(",\"different_vertex_count\":", AreDerivedEquivalent(chain_q, branch_q));
Print(",\"different_genus\":", AreDerivedEquivalent(chain_q, genus_one_a));
Print("}}\n");

QUIT_GAP(0);
