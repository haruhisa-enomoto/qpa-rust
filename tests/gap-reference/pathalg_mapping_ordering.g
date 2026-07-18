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

PrintElement := function(element)
    local data, i;
    data := CoefficientsAndMagmaElements(element);
    Print("[");
    for i in [1, 3 .. Length(data) - 1] do
        if i > 1 then
            Print(",");
        fi;
        Print("{\"path\":"); PrintJsonString(StringPrint(data[i]));
        Print(",\"coefficient\":", IntFFE(data[i + 1]), "}");
    od;
    Print("]");
end;

source_quiver := Quiver(1, [[1, 1, "a"], [1, 1, "b"]]);;
range_quiver := Quiver(1, [[1, 1, "x"], [1, 1, "y"]]);;
source := PathAlgebra(GF(5), source_quiver);;
range := PathAlgebra(GF(5), range_quiver);;
source_generators := GeneratorsOfAlgebra(source);;
range_images := [range.v1, range.y, range.x];;
element := 2 * source.v1 + 3 * source.a^2
    + 4 * source.a * source.b + source.b;;
mapped := MappedExpression(element, source_generators, range_images);;
mapped_zero := MappedExpression(Zero(source), source_generators, range_images);;

path_hom := AlgebraWithOneHomomorphismByImagesNC(source, range,
    source_generators, range_images);;
unretained := CALL_WITH_CATCH(ImagesRepresentative, [path_hom, element]);;
path_hom!.generators := source_generators;;
path_hom!.genimages := range_images;;
retained := CALL_WITH_CATCH(ImagesRepresentative, [path_hom, element]);;

assigned := CALL_WITH_CATCH(AssignGeneratorVariables, [source]);;
invalid_quiver := Quiver(1, [[1, 1, "not-valid"]]);;
invalid_algebra := PathAlgebra(GF(5), invalid_quiver);;
invalid_assignment := CALL_WITH_CATCH(AssignGeneratorVariables,
    [invalid_algebra]);;
read_only_quiver := Quiver(1, [[1, 1, "Z"]]);;
read_only_algebra := PathAlgebra(GF(5), read_only_quiver);;
read_only_assignment := CALL_WITH_CATCH(AssignGeneratorVariables,
    [read_only_algebra]);;

source_quotient := source / [source.a^2];;
range_quotient := range / [range.y^2];;
quotient_generators := GeneratorsOfAlgebraWithOne(source_quotient);;
quotient_images := [range_quotient.v1, range_quotient.y, range_quotient.x];;
quotient_element := 2 * source_quotient.v1
    + 3 * source_quotient.a * source_quotient.b
    + 4 * source_quotient.b^2;;
quotient_mapped := MappedExpression(quotient_element,
    quotient_generators, quotient_images);;
quotient_hom := AlgebraWithOneHomomorphismByImagesNC(source_quotient,
    range_quotient, quotient_generators, quotient_images);;
quotient_hom!.generators := quotient_generators;;
quotient_hom!.genimages := quotient_images;;
quotient_retained := CALL_WITH_CATCH(ImagesRepresentative,
    [quotient_hom, quotient_element]);;
ordered := CALL_WITH_CATCH(OrderedBy,
    [source_quotient, OrderingOfAlgebra(source_quotient)]);;

Print("{\"schema_version\":1,\"kind\":\"pathalg_mapping_ordering\",\"oracle\":{");
Print("\"gap_version\":"); PrintJsonString(GAPInfo.Version);
Print(",\"gap_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION);
Print(",\"qpa_version\":"); PrintJsonString(PackageInfo("qpa")[1].Version);
Print(",\"qpa_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION);
Print(",\"gbnp_version\":"); PrintJsonString(PackageInfo("gbnp")[1].Version);
Print(",\"gbnp_revision\":");
PrintJsonString(GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION);
Print("},\"path\":{");
Print("\"mapped\":"); PrintElement(mapped);
Print(",\"mapped_zero\":"); PrintElement(mapped_zero);
Print(",\"unretained_success\":", unretained[1]);
Print(",\"retained_success\":", retained[1]);
if retained[1] then
    Print(",\"retained_image\":"); PrintElement(retained[2]);
fi;
Print("},\"quotient\":{");
Print("\"mapped\":"); PrintElement(quotient_mapped![1]);
Print(",\"retained_success\":", quotient_retained[1]);
if quotient_retained[1] then
    Print(",\"retained_image\":"); PrintElement(quotient_retained[2]![1]);
fi;
Print(",\"ordered_success\":", ordered[1]);
Print("},\"assignment\":{");
Print("\"success\":", assigned[1]);
Print(",\"invalid_success\":", invalid_assignment[1]);
Print(",\"read_only_success\":", read_only_assignment[1]);
Print(",\"names\":[\"v1\",\"a\",\"b\"]");
Print(",\"values\":[");
PrintElement(v1); Print(","); PrintElement(a); Print(","); PrintElement(b);
Print("]}}\n");

QUIT_GAP(0);
