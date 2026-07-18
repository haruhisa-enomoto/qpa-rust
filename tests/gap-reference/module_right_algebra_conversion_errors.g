if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

Convert := function(module)
    return RightAlgebraModuleToPathAlgebraMatModule(module);
end;;
BuildAndConvertRegular := function(algebra)
    local module;
    module := RightAlgebraModule(algebra, \*, algebra);
    return RightAlgebraModuleToPathAlgebraMatModule(module);
end;;

infinite_algebra := PathAlgebra(Rationals,
    Quiver(["v"], [["v", "v", "x"]]));;
nonfinite := CALL_WITH_CATCH(BuildAndConvertRegular,
    [infinite_algebra]);;

matrix_algebra := FullMatrixAlgebra(Rationals, 2);;
matrix_module := RightAlgebraModule(matrix_algebra, \*, matrix_algebra);;
wrong_algebra := CALL_WITH_CATCH(Convert, [matrix_module]);;

Print("{\"schema_version\":1,\"kind\":\"module_right_algebra_conversion_errors\",\"success\":{");
Print("\"nonfinite\":", nonfinite[1], ",");
Print("\"wrong_algebra\":", wrong_algebra[1], "}}\n");

QUIT_GAP(0);
