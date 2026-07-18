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

q := Quiver(["v", "w", "x"],
    [["v", "w", "a"], ["w", "x", "b"], ["v", "x", "c"]]);
algebra := PathAlgebra(GF(2), q);
map := [[algebra.a, algebra.w], [algebra.c, algebra.b]];
module := VertexProjectivePresentation(algebra, map);
other_q := Quiver(["z"], []);
other_algebra := PathAlgebra(GF(2), other_q);

hom_result := CALL_WITH_CATCH(Hom, [algebra, module, module]);
end_result := CALL_WITH_CATCH(End, [algebra, module]);
zero_module_result := CALL_WITH_CATCH(VertexProjectivePresentation,
    [algebra, [[algebra.v, Zero(algebra)], [Zero(algebra), algebra.w]]]);
foreign_hom_result := CALL_WITH_CATCH(Hom, [other_algebra, module, module]);
foreign_end_result := CALL_WITH_CATCH(End, [other_algebra, module]);
module_generators := GeneratorsOfAlgebraModule(module);
by_images := LeftModuleHomomorphismByImagesNC(module, module,
    module_generators, module_generators);
SetFilterObj(by_images, IsAlgebraModuleHomomorphism);
images_by_images_result := CALL_WITH_CATCH(ImagesSource, [by_images]);
by_matrix := LeftModuleHomomorphismByMatrix(Basis(module),
    IdentityMat(Dimension(module), GF(2)), Basis(module));
SetFilterObj(by_matrix, IsAlgebraModuleHomomorphism);
images_by_matrix_result := CALL_WITH_CATCH(ImagesSource, [by_matrix]);

Print("{\"schema_version\":1,\"kind\":\"projpathalgmodule_fp_hom\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"dimension\":", Dimension(module));
Print(",\"success\":[", hom_result[1], ",", end_result[1], ",",
    zero_module_result[1], ",",
    foreign_hom_result[1], ",", foreign_end_result[1], ",",
    images_by_images_result[1], ",", images_by_matrix_result[1], "]");
Print(",\"basisimage_bound\":[", IsBound(by_images!.basisimage), ",",
    IsBound(by_matrix!.basisimage), "]");
if images_by_images_result[1] then
    Print(",\"images_by_images_dimension\":", Dimension(images_by_images_result[2]));
fi;
if images_by_matrix_result[1] then
    Print(",\"images_by_matrix_dimension\":", Dimension(images_by_matrix_result[2]));
fi;
if hom_result[1] then
    Print(",\"hom_dimension\":", Dimension(hom_result));
fi;
if end_result[1] then
    Print(",\"end_dimension\":", Dimension(end_result));
fi;
Print("}\n");

QUIT_GAP(0);
