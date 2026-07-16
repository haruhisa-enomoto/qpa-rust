if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

PrintJsonString := function(value)
    local escaped;
    escaped := ReplacedString(value, "\\", "\\\\");
    escaped := ReplacedString(escaped, "\"", "\\\"");
    Print("\"", escaped, "\"");
end;

space := FullRowSpace(Rationals, 2);;
basis := Basis(space);;
identity := LeftModuleHomomorphismByMatrix(
    basis, [[1, 0], [0, 1]], basis);;
nilpotent := LeftModuleHomomorphismByMatrix(
    basis, [[0, 1], [0, 0]], basis);;
algebra := AlgebraWithOneByGenerators(
    Rationals, [identity, nilpotent]);;

had_basis := HasBasis(algebra);;
radical := RadicalOfAlgebra(algebra);;
radical_basis := BasisVectors(Basis(radical));;

Print("{\"schema_version\":1,\"kind\":\"algebra_mapping_radical\",\"oracle\":{");
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
Print("},\"had_basis_before\":", had_basis);
Print(",\"dimension\":", Dimension(radical));
Print(",\"basis_matrices\":", List(radical_basis, x -> x!.matrix));
Print(",\"all_square_zero\":", ForAll(radical_basis, x -> IsZero(x * x)));
Print("}\n");

QUIT_GAP(0);
